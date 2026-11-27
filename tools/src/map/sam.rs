//! Turns alignments into SAM text.
//! That covers the flag and column constants, the MD, cs and NM tags, quality strings, read groups and the short-read record lines.

use super::*;

/// SAM flag bits (SAM specification, section 1.4).
pub(crate) const SAM_FLAG_PAIRED: u16 = 0x1;

pub(crate) const SAM_FLAG_PROPER_PAIR: u16 = 0x2;

pub(crate) const SAM_FLAG_UNMAPPED: u16 = 0x4;

pub(crate) const SAM_FLAG_MATE_UNMAPPED: u16 = 0x8;

pub(crate) const SAM_FLAG_REVERSE: u16 = 0x10;

pub(crate) const SAM_FLAG_MATE_REVERSE: u16 = 0x20;

pub(crate) const SAM_FLAG_READ1: u16 = 0x40;

pub(crate) const SAM_FLAG_READ2: u16 = 0x80;

pub(crate) const SAM_FLAG_SUPPLEMENTARY: u16 = 0x800;

pub(super) const SAM_MAX_MAPQ: i32 = 60;

/// SAM alignment columns, 0-based (SAM specification, section 1.4).
pub(crate) const SAM_COL_QNAME: usize = 0;
pub(crate) const SAM_COL_FLAG: usize = 1;

pub(crate) const SAM_COL_RNAME: usize = 2;

pub(crate) const SAM_COL_POS: usize = 3;

pub(crate) const SAM_COL_MAPQ: usize = 4;

pub(crate) const SAM_COL_CIGAR: usize = 5;

pub(crate) const SAM_COL_RNEXT: usize = 6;

pub(crate) const SAM_COL_PNEXT: usize = 7;

pub(crate) const SAM_COL_TLEN: usize = 8;

pub(crate) const SAM_COL_SEQ: usize = 9;

pub(crate) const SAM_COL_QUAL: usize = 10;

/// Mandatory columns of a SAM alignment line, QNAME to QUAL, after which optional tags follow.
pub(crate) const SAM_MANDATORY_COLUMNS: usize = 11;

/// An optional SAM field is TAG:TYPE:VALUE, so it has at least this many bytes ("XY:" plus a value) and a ':' after the two-letter tag.
pub(super) const SAM_TAG_MIN_LEN: usize = 4;

pub(super) const SAM_TAG_SEP_POS: usize = 2;

/// Phred quality encoding, Q = -10 log10(p), written as ASCII Q + 33 (Sanger/Illumina 1.8+).
pub(super) const PHRED_SCALE: f64 = 10.0;

pub(super) const PHRED_ASCII_OFFSET: u8 = 33;

/// The SAM record of a short read at a verified placement.
#[allow(clippy::too_many_arguments)]
pub(super) fn short_format(
    read_name: &str,
    seq: &[u8],
    rev_seq: &[u8],
    hit: &ShortHit,
    aln: &ShortAln,
    mapq: i32,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
) -> String {
    let align_seq: &[u8] = if hit.rev { rev_seq } else { seq };
    let flag = if hit.rev { SAM_FLAG_REVERSE } else { 0 };
    // SAM SEQ is the read as aligned to the forward reference strand.
    let out_seq = std::str::from_utf8(align_seq).unwrap_or("*");
    let tags = diff_tag_fields(
        &aln.cig,
        align_seq,
        &ctx.ref_seqs[hit.rid],
        hit.win_start + aln.offset,
        ctx.cs_mode,
    );

    format!(
        "{read_name}\t{flag}\t{}\t{}\t{mapq}\t{}\t*\t0\t0\t{out_seq}\t*\tNM:i:{}{tags}",
        ref_names[hit.rid],
        hit.win_start + aln.offset + 1,
        aln.cig,
        aln.edit_dist
    )
}

/// The SAM record of a short read at its verified placement `hit`, with a MAPQ from its own evidence that `mapq_cap` caps.
pub(super) fn short_record(
    read_name: &str,
    seq: &[u8],
    hit: &ShortHit,
    ctx: &AlignCtx<'_>,
    ref_names: &[String],
    mapq_cap: Option<i32>,
) -> String {
    let rev_seq = revcomp(seq);
    let aln = short_align(seq, &rev_seq, hit, ctx);
    let mut mapq = final_mapq(hit, aln.events, seq.len(), ctx.rho);

    if let Some(cap) = mapq_cap {
        mapq = mapq.min(cap);
    }

    let mut rec = short_format(read_name, seq, &rev_seq, hit, &aln, mapq, ctx, ref_names);

    // OPTIK_MAPQ_DEBUG appends every MAPQ input as an XD tag for calibration scripts.
    if mapq_debug() {
        let comps: Vec<String> = hit.comp.iter().map(|comp| comp.to_string()).collect();
        let [n_loci_all, n_min, n_rep, n_best, n_other] = hit.seed_ev;

        rec.push_str(&format!(
            "\tXD:Z:{},{},{},{},{},{},{},{}:{}",
            hit.edit_dist,
            aln.events,
            n_loci_all,
            n_min,
            n_rep,
            n_best,
            n_other,
            comps.len(),
            comps.join(";")
        ));
    }

    rec
}

/// Writes the read's base qualities into the QUAL column of its records.
/// They are reversed on reverse-strand records, whose SEQ is reverse-complemented, and FASTA reads keep '*'.
/// Supplementary records are hard-clipped, holding only the aligned bases in SEQ as minimap2 and BWA-MEM write them, so their QUAL is cut by the CIGAR's leading and trailing H the same way.
pub(super) fn fill_qual(lines: &mut [String], qual: &[u8]) {
    if qual.is_empty() {
        return;
    }

    for line in lines.iter_mut() {
        let mut fields: Vec<&str> = line.splitn(SAM_MANDATORY_COLUMNS + 1, '\t').collect();

        if fields.len() < SAM_MANDATORY_COLUMNS || fields[SAM_COL_SEQ] == "*" {
            continue;
        }

        let flag: u16 = fields[SAM_COL_FLAG].parse().unwrap_or(0);
        // SEQ of a reverse record is reverse-complemented, so its qualities run backwards.
        let oriented: Vec<u8> = if flag & SAM_FLAG_REVERSE != 0 {
            qual.iter().rev().copied().collect()
        } else {
            qual.to_vec()
        };
        let (lead, trail) = hard_clips(fields[SAM_COL_CIGAR]);

        if lead + trail > oriented.len() {
            continue;
        }

        let qual_str: String = oriented[lead..oriented.len() - trail]
            .iter()
            .map(|&qual_byte| qual_byte as char)
            .collect();

        // A length mismatch would make invalid SAM. Such a record keeps QUAL as *.
        if qual_str.len() != fields[SAM_COL_SEQ].len() {
            continue;
        }

        fields[SAM_COL_QUAL] = &qual_str;
        let joined = fields.join("\t");

        *line = joined;
    }
}

/// minimap2's cs tag, off by default, short with `:n` runs, or long with `=SEQ`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum CsMode {
    Off,
    Short,
    Long,
}

/// MD, per the SAM spec, and optionally minimap2's cs of one record.
/// The CIGAR's clips are skipped, and aligned columns are compared as NM counts them, by equal bases, where a non-ACGT reference base never matches.
/// Deletions are spelled from the reference and, in cs only, insertions from the read.
/// MD keeps a number between any two elements, as in `10A0^TT5`, the way samtools calmd writes it.
pub(super) fn diff_tags(
    cig: &str,
    seq: &[u8],
    ref_seq: &[u8],
    pos: usize,
    cs_mode: CsMode,
) -> (String, Option<String>) {
    let lower = |base: u8| (base as char).to_ascii_lowercase();
    let mut md_tag = String::new();
    let mut run = 0usize; // MD matching run
    let mut csv = String::new();
    let mut ident = 0usize; // cs short identical run
    let mut ident_long = String::new(); // cs long identical bases
    // Identical bases are buffered and written only when a difference or the end interrupts them, because cs gives a run as one :N (short form) or =BASES (long form).
    let flush_cs = |csv: &mut String, ident: &mut usize, ident_long_buf: &mut String| {
        if *ident > 0 {
            csv.push_str(&format!(":{ident}"));
            *ident = 0;
        }

        if !ident_long_buf.is_empty() {
            csv.push('=');
            csv.push_str(ident_long_buf);
            ident_long_buf.clear();
        }
    };
    let (mut read_i, mut ref_i, mut run_len) = (0usize, pos, 0usize);

    for byte in cig.bytes() {
        if byte.is_ascii_digit() {
            run_len = push_digit(run_len, byte);
            continue;
        }

        match byte {
            // Soft-clipped bases sit in SEQ and are skipped. Hard-clipped ones are absent.
            b'S' => read_i += run_len,
            b'H' => {}
            b'M' | b'=' | b'X' => {
                for _ in 0..run_len {
                    let (read_base, ref_base) = (
                        seq.get(read_i).copied().unwrap_or(b'N'),
                        ref_seq.get(ref_i).copied().unwrap_or(b'N'),
                    );

                    // An N on either side is written as a mismatch in MD and cs, as minimap2 does, even when both bases are N.
                    if read_base == ref_base && matches!(ref_base, b'A' | b'C' | b'G' | b'T') {
                        run += 1;

                        match cs_mode {
                            CsMode::Short => ident += 1,
                            CsMode::Long => ident_long.push(read_base as char),
                            CsMode::Off => {}
                        }
                    } else {
                        md_tag.push_str(&run.to_string());
                        md_tag.push(ref_base as char);
                        run = 0;

                        if cs_mode != CsMode::Off {
                            flush_cs(&mut csv, &mut ident, &mut ident_long);
                            csv.push('*');
                            csv.push(lower(ref_base));
                            csv.push(lower(read_base));
                        }
                    }

                    read_i += 1;
                    ref_i += 1;
                }
            }
            // MD has no notation for insertions, so only cs records the inserted bases.
            b'I' => {
                if cs_mode != CsMode::Off {
                    flush_cs(&mut csv, &mut ident, &mut ident_long);
                    csv.push('+');
                    let end = (read_i + run_len).min(seq.len());

                    if read_i < end {
                        for &base in &seq[read_i..end] {
                            csv.push(lower(base));
                        }
                    }
                }

                read_i += run_len;
            }
            // A deletion breaks the MD match run and lists the deleted reference bases after ^.
            b'D' | b'N' => {
                md_tag.push_str(&run.to_string());
                md_tag.push('^');
                run = 0;

                if cs_mode != CsMode::Off {
                    flush_cs(&mut csv, &mut ident, &mut ident_long);
                    csv.push('-');
                }

                for j in ref_i..ref_i + run_len {
                    let ref_base = ref_seq.get(j).copied().unwrap_or(b'N');

                    md_tag.push(ref_base as char);

                    if cs_mode != CsMode::Off {
                        csv.push(lower(ref_base));
                    }
                }

                ref_i += run_len;
            }
            _ => {}
        }

        run_len = 0;
    }

    // MD always ends with a match count, which may be 0.
    md_tag.push_str(&run.to_string());

    if cs_mode == CsMode::Off {
        return (md_tag, None);
    }

    flush_cs(&mut csv, &mut ident, &mut ident_long);
    (md_tag, Some(csv))
}

/// "\tMD:Z:..[\tcs:Z:..]" for a record.
pub(super) fn diff_tag_fields(cig: &str, seq: &[u8], ref_seq: &[u8], pos: usize, cs_mode: CsMode) -> String {
    let (md_tag, cs_tag) = diff_tags(cig, seq, ref_seq, pos, cs_mode);

    match cs_tag {
        Some(cs_tag) => format!("\tMD:Z:{md_tag}\tcs:Z:{cs_tag}"),
        None => format!("\tMD:Z:{md_tag}"),
    }
}

/// SAM QNAME length limit.
pub(crate) const MAX_QNAME: usize = 254;

/// Parses a `-R` read-group header line (the bwa and minimap2 form, with `\t` written as the two characters or as a tab), returning the normalized line and its ID.
pub(crate) fn parse_read_group(read_group: &str) -> Result<(String, String), String> {
    // Shell users type \t literally, as with bwa -R.
    let line = read_group.replace("\\t", "\t");

    if !line.starts_with("@RG\t") {
        return Err(format!(
            "-R must be a read-group line starting with '@RG\\t', got '{read_group}'"
        ));
    }

    if line.contains('\n') || line.contains('\r') {
        return Err("-R may not contain a line break".into());
    }

    let rg_id = line
        .split('\t')
        .skip(1)
        .find_map(|field| field.strip_prefix("ID:"))
        .filter(|rg_id| !rg_id.is_empty())
        .ok_or_else(|| format!("-R read group has no ID field: '{read_group}'"))?;

    for field in line.split('\t').skip(1) {
        if field.len() < SAM_TAG_MIN_LEN || field.as_bytes()[SAM_TAG_SEP_POS] != b':' {
            return Err(format!("-R field '{field}' is not TAG:VALUE"));
        }
    }

    Ok((line.clone(), rg_id.to_string()))
}

/// Mate names, where a trailing /1 or /2 is not part of the read's name.
pub(super) fn mate_name(name: &str) -> &str {
    name.strip_suffix("/1")
        .or_else(|| name.strip_suffix("/2"))
        .unwrap_or(name)
}
