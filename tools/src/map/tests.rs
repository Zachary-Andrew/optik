#[test]
fn local_pieces_clips_and_splits() {
    use super::{MM_ASM20, MM_READS, asm_score, local_pieces};
    let offs = |cols: &[u8]| {
        let (mut read_off, mut ref_off) = (vec![0usize], vec![0usize]);

        for &col_op in cols {
            read_off.push(read_off.last().unwrap() + (col_op != b'D') as usize);
            ref_off.push(ref_off.last().unwrap() + (col_op != b'I') as usize);
        }

        (read_off, ref_off)
    };
    let kmer_len = 15;
    // An error-free path is one piece over every column.
    let read = rand_dna(400, 1);
    let cols = vec![b'M'; 400];
    let (read_off, ref_off) = offs(&cols);
    let anchor_offs = [0, 100, 200, 300];

    assert_eq!(
        local_pieces(
            &cols,
            &read,
            &read,
            &read_off,
            &ref_off,
            &anchor_offs,
            kmer_len,
            MM_READS
        ),
        vec![(0, 400, 800)]
    );
    // A junk end past the last anchor is clipped back to the last match before it.
    let mut reference = read.clone();

    for i in 386..400 {
        if i % 7 != 0 {
            reference[i] = if read[i] == b'A' { b'C' } else { b'A' };
        }
    }

    let pieces = local_pieces(
        &cols,
        &read,
        &reference,
        &read_off,
        &ref_off,
        &anchor_offs,
        kmer_len,
        MM_READS,
    );

    assert!(
        pieces.len() == 1 && pieces[0].0 == 0 && pieces[0].1 <= 386,
        "{pieces:?}"
    );
    // One long deletion between anchors is neither a Z-drop nor clipped, so it stays one piece with all columns.
    let mut cols2 = vec![b'M'; 200];

    cols2.extend(std::iter::repeat_n(b'D', 300));
    cols2.extend(std::iter::repeat_n(b'M', 200));
    let ref_del: Vec<u8> = read[..200]
        .iter()
        .chain(rand_dna(300, 2).iter())
        .chain(read[200..].iter())
        .copied()
        .collect();
    let (qp2, rp2) = offs(&cols2);
    let pieces_del = local_pieces(
        &cols2,
        &read,
        &ref_del,
        &qp2,
        &rp2,
        &[0, 100, 250, 350],
        kmer_len,
        MM_READS,
    );

    assert!(
        pieces_del.len() == 1 && pieces_del[0].0 == 0 && pieces_del[0].1 == 700,
        "{pieces_del:?}"
    );
    // 300 columns of unrelated sequence between two anchored blocks split into two pieces.
    let read_junk: Vec<u8> = read[..150]
        .iter()
        .chain(rand_dna(300, 3).iter())
        .chain(read[150..300].iter())
        .copied()
        .collect();
    let ref_junk: Vec<u8> = read[..150]
        .iter()
        .chain(rand_dna(300, 4).iter())
        .chain(read[150..300].iter())
        .copied()
        .collect();
    let cols3 = vec![b'M'; 600];
    let (qp3, rp3) = offs(&cols3);
    let pieces_junk = local_pieces(
        &cols3,
        &read_junk,
        &ref_junk,
        &qp3,
        &rp3,
        &[0, 60, 120, 460, 520, 580],
        kmer_len,
        MM_READS,
    );

    assert!(
        pieces_junk.len() == 2
            && pieces_junk[0].0 == 0
            && pieces_junk[0].1 <= 154
            && pieces_junk[1].0 >= 446
            && pieces_junk[1].1 == 600,
        "{pieces_junk:?}"
    );
    // Under asm20 scoring the long deletion is still one piece, thanks to the two-piece gap's e2 allowance.
    let pieces_asm = local_pieces(
        &cols2,
        &read,
        &ref_del,
        &qp2,
        &rp2,
        &[0, 100, 250, 350],
        kmer_len,
        MM_ASM20,
    );

    assert!(
        pieces_asm.len() == 1 && pieces_asm[0].0 == 0 && pieces_asm[0].1 == 700,
        "{pieces_asm:?}"
    );
    // An error-free 400 bp block scores 400, and the same block at 19% substitutions scores under 200.
    assert_eq!(asm_score("400M", &read, &read, 0), 400);
    let mut diverged = read.clone();

    for i in (0..400).filter(|i| i % 16 < 3) {
        diverged[i] = if read[i] == b'A' { b'C' } else { b'A' };
    }

    assert!(asm_score("400M", &read, &diverged, 0) < super::ASM_MIN_DP);
    assert_eq!(
        asm_score(
            "10S100M30D100M",
            read[..210].to_vec().as_slice(),
            &[&read[10..110], &diverged[..30], &read[110..210]].concat(),
            0
        ),
        200 - 26 - 30
    );
}

#[test]
fn hard_clip_supplementary() {
    use super::hard_clip;
    assert_eq!(
        hard_clip("3S4M2S", "AAACCCCGG"),
        ("3H4M2H".to_string(), "CCCC")
    );
    assert_eq!(
        hard_clip("5M2I3M", "ACGTACGTAC"),
        ("5M2I3M".to_string(), "ACGTACGTAC")
    );
    assert_eq!(
        hard_clip("10S5M", "AAAAAAAAAACCCCC"),
        ("10H5M".to_string(), "CCCCC")
    );
    assert_eq!(
        hard_clip("5M1D2M12S", "CCCCCGGAAAAAAAAAAAA"),
        ("5M1D2M12H".to_string(), "CCCCCGG")
    );
}

/// QUAL follows SEQ on hard-clipped supplementary records, cut by the CIGAR's H in the record's orientation, with qualities reversed on reverse records.
#[test]
fn qual_of_hard_clipped_records() {
    use super::{fill_qual, hard_clips};
    assert_eq!(hard_clips("3H4M2H"), (3, 2));
    assert_eq!(hard_clips("10H5M"), (10, 0));
    assert_eq!(hard_clips("5M12H"), (0, 12));
    assert_eq!(hard_clips("3S4M2S"), (0, 0));
    // Read qualities are "ABCDEFGHI", and forward record 3H4M2H holds bases 3..7.
    let mut lines = vec!["r\t2048\tc\t1\t60\t3H4M2H\t*\t0\t0\tCCCC\t*".to_string()];

    fill_qual(&mut lines, b"ABCDEFGHI");
    assert!(lines[0].ends_with("\tCCCC\tDEFG"), "{}", lines[0]);
    // A reverse record has its qualities reversed ("IHGFEDCBA"), then cut 2H..1H.
    let mut lines = vec!["r\t2064\tc\t1\t60\t2H6M1H\t*\t0\t0\tAAAAAA\t*".to_string()];

    fill_qual(&mut lines, b"ABCDEFGHI");
    assert!(lines[0].ends_with("\tAAAAAA\tGFEDCB"), "{}", lines[0]);
    // An unclipped primary keeps the whole quality string.
    let mut lines = vec!["r\t0\tc\t1\t60\t9M\t*\t0\t0\tAAAAAAAAA\t*".to_string()];

    fill_qual(&mut lines, b"ABCDEFGHI");
    assert!(lines[0].ends_with("\tABCDEFGHI"), "{}", lines[0]);
}

use super::*;

/// Insert model from quartiles +- 3 IQR, with none below 30 pairs.
#[test]
fn insert_model_from_quartiles() {
    assert_eq!(insert_model(vec![400; 29]), None);
    let tlens: Vec<i64> = (0..100).map(|i| 300 + 2 * i).collect(); // 300..498
    let model = insert_model(tlens).unwrap();

    // Q1 = 350 and Q3 = 448 by rounded ranks, so IQR = 98.
    assert_eq!((model.min_tlen, model.max_tlen), (350 - 294, 448 + 294));
}

/// FR template length, from forward start to reverse end on the same reference.
#[test]
fn fr_template_length() {
    let fwd_place = Place {
        rid: 0,
        rev: false,
        start: 100,
        end: 250,
    };
    let rev_place = Place {
        rid: 0,
        rev: true,
        start: 350,
        end: 500,
    };

    assert_eq!(fr_tlen(fwd_place, rev_place), Some(400));
    assert_eq!(fr_tlen(rev_place, fwd_place), Some(400));
    assert_eq!(
        fr_tlen(
            fwd_place,
            Place {
                rev: false,
                ..rev_place
            }
        ),
        None
    ); // FF
    assert_eq!(
        fr_tlen(
            fwd_place,
            Place {
                rid: 1,
                ..rev_place
            }
        ),
        None
    ); // other contig
    assert_eq!(
        fr_tlen(
            Place {
                start: 600,
                end: 750,
                ..fwd_place
            },
            rev_place
        ),
        None
    ); // RF (outward)
    let model = InsertModel {
        min_tlen: 200,
        max_tlen: 600,
    };

    assert!(proper(fwd_place, rev_place, &model));
    assert!(!proper(
        fwd_place,
        Place {
            start: 950,
            end: 1100,
            ..rev_place
        },
        &model
    ));
}

#[test]
fn mate_names() {
    assert_eq!(mate_name("r7/1"), "r7");
    assert_eq!(mate_name("r7/2"), "r7");
    assert_eq!(mate_name("r7"), "r7");
    assert_eq!(mate_name("r7/3"), "r7/3");
}

/// SAM pair fields, i.e. flags, RNEXT/PNEXT, signed TLEN, MC/MQ, and an unmapped mate placed at its partner.
#[test]
fn pair_fields() {
    let mut lines_a = vec!["p\t0\tc\t101\t60\t150M\t*\t0\t0\tACGT\t*\tNM:i:0".to_string()];
    let mut lines_b = vec!["p\t16\tc\t351\t40\t140M10S\t*\t0\t0\tACGT\t*\tNM:i:1".to_string()];

    set_pair_fields(&mut lines_a, &mut lines_b, true);
    let fields_a: Vec<&str> = lines_a[0].split('\t').collect();
    let fields_b: Vec<&str> = lines_b[0].split('\t').collect();

    assert_eq!(fields_a[1], (1 | 2 | 0x20 | 0x40).to_string());
    assert_eq!(fields_b[1], (1 | 2 | 16 | 0x80).to_string());
    assert_eq!((fields_a[6], fields_a[7], fields_a[8]), ("=", "351", "390")); // 101..490
    assert_eq!(
        (fields_b[6], fields_b[7], fields_b[8]),
        ("=", "101", "-390")
    );
    assert!(
        lines_a[0].ends_with("\tMC:Z:140M10S\tMQ:i:40")
            && lines_b[0].ends_with("\tMC:Z:150M\tMQ:i:60")
    );
    // With one mate unmapped it takes the partner's position, with flags 0x8 / 0x4.
    let mut lines_a = vec!["p\t0\tc\t101\t60\t150M\t*\t0\t0\tACGT\t*\tNM:i:0".to_string()];
    let mut lines_unmapped = vec!["p\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*".to_string()];

    set_pair_fields(&mut lines_a, &mut lines_unmapped, false);
    let fields_a: Vec<&str> = lines_a[0].split('\t').collect();
    let fields_unmapped: Vec<&str> = lines_unmapped[0].split('\t').collect();

    assert_eq!(fields_a[1], (1 | 8 | 0x40).to_string());
    assert_eq!((fields_a[6], fields_a[7], fields_a[8]), ("=", "101", "0"));
    assert_eq!(fields_unmapped[1], (1 | 4 | 0x80).to_string());
    assert_eq!(
        (
            fields_unmapped[2],
            fields_unmapped[3],
            fields_unmapped[6],
            fields_unmapped[7],
            fields_unmapped[8]
        ),
        ("c", "101", "=", "101", "0")
    );
}

/// Deterministic xorshift64 PRNG, avoiding an external dev-dependency.
fn rng_state(seed: u64) -> u64 {
    let mut state = seed;

    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    state
}

fn next_rng(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

const BASES: &[u8] = b"ACGT";

/// Random DNA window of `len` bases from a deterministic PRNG.
fn rand_dna(len: usize, seed: u64) -> Vec<u8> {
    let mut state = rng_state(seed);

    (0..len)
        .map(|_| BASES[(next_rng(&mut state) % 4) as usize])
        .collect()
}

/// Reconstructs the semi-global edit distance a CIGAR implies by walking it from (query=0, text=text_start).
/// Returns None if the CIGAR overshoots.
fn cigar_edit_distance(
    cigar: &str,
    query: &[u8],
    text: &[u8],
    text_start: usize,
) -> Option<usize> {
    let mut query_pos = 0usize; // query offset
    let mut text_pos = text_start; // text offset
    let mut edits = 0usize;
    let mut num = 0usize;

    for byte in cigar.bytes() {
        match byte {
            b'0'..=b'9' => num = num * 10 + (byte - b'0') as usize,
            _ => {
                if num == 0 {
                    return None;
                }

                match byte {
                    b'M' => {
                        for _ in 0..num {
                            if query_pos >= query.len() || text_pos >= text.len() {
                                return None;
                            }

                            if query[query_pos] != text[text_pos] {
                                edits += 1;
                            }

                            query_pos += 1;
                            text_pos += 1;
                        }
                    }
                    b'I' => {
                        for _ in 0..num {
                            if query_pos >= query.len() {
                                return None;
                            }

                            edits += 1;
                            query_pos += 1;
                        }
                    }
                    b'D' => {
                        for _ in 0..num {
                            if text_pos >= text.len() {
                                return None;
                            }

                            edits += 1;
                            text_pos += 1;
                        }
                    }
                    _ => return None,
                }

                num = 0;
            }
        }
    }

    if query_pos != query.len() {
        return None;
    }

    Some(edits)
}

/// Independent O(mn) semi-global edit distance (query end-to-end, text ends free).
fn naive_semi_global_ed(query: &[u8], text: &[u8]) -> usize {
    let mut prev = vec![0usize; text.len() + 1]; // free text start

    for i in 1..=query.len() {
        let mut cur = vec![i; text.len() + 1];

        for j in 1..=text.len() {
            let sub = prev[j - 1] + (query[i - 1] != text[j - 1]) as usize;

            cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
        }

        prev = cur;
    }

    *prev.iter().min().unwrap() // free text end
}

#[test]
fn left_align_normalizes_and_preserves() {
    // Random reads with planted indels in repeat contexts.
    // left_align keeps the same columns, i.e. the ops multiset, the mismatches and both strings fully consumed, and leaves no gap that could still move left.
    let expand = |cigar_str: &str| -> Vec<u8> {
        let mut cols = Vec::new();
        let mut run_len = 0usize;

        for byte in cigar_str.bytes() {
            if byte.is_ascii_digit() {
                run_len = run_len * 10 + (byte - b'0') as usize;
            } else {
                cols.extend(std::iter::repeat_n(byte, run_len));
                run_len = 0;
            }
        }

        cols
    };
    let mut rng = rng_state(77);

    for case in 0..400 {
        // Text with homopolymers and short tandem repeats.
        let mut text = Vec::new();

        while text.len() < 200 {
            let unit: Vec<u8> = (0..(1 + next_rng(&mut rng) % 3))
                .map(|_| b"ACGT"[(next_rng(&mut rng) % 4) as usize])
                .collect();

            for _ in 0..(1 + next_rng(&mut rng) % 4) {
                text.extend_from_slice(&unit);
            }
        }

        let read_start = 20 + (next_rng(&mut rng) % 40) as usize;
        let mut read: Vec<u8> = text[read_start..read_start + 120].to_vec();
        let indel_at = 30 + (next_rng(&mut rng) % 60) as usize;
        let indel_len = 1 + (next_rng(&mut rng) % 4) as usize;
        let del = case % 2 == 0;

        if del {
            read.drain(indel_at..indel_at + indel_len);
        } else {
            let ins: Vec<u8> = read[indel_at - indel_len..indel_at].to_vec();

            for (j, base) in ins.into_iter().enumerate() {
                read.insert(indel_at + j, base);
            }
        }

        let Some((cig, off)) = crate::wfa_affine_cigar(&read, &text) else {
            continue;
        };
        let out = left_align(&cig, &read, &text, off);
        let (cols_before, cols_after) = (expand(&cig), expand(&out));
        let count = |cols: &[u8], op: u8| cols.iter().filter(|&&col| col == op).count();

        for op in *b"MID" {
            assert_eq!(
                count(&cols_before, op),
                count(&cols_after, op),
                "case {case}: {cig} -> {out}"
            );
        }

        assert_eq!(
            cigar_ed_implied(&cig, &read, &text, off),
            cigar_ed_implied(&out, &read, &text, off),
            "case {case}"
        );
        // No gap can still shift, i.e. the column before each gap run is not a match equal to the gap's last base.
        let (mut query_i, mut text_i) = (0usize, off);
        let mut col = 0;

        while col < cols_after.len() {
            let op = cols_after[col];

            if op == b'M' {
                query_i += 1;
                text_i += 1;
                col += 1;
                continue;
            }

            let mut run_end = col;

            while run_end < cols_after.len() && cols_after[run_end] == op {
                run_end += 1;
            }

            let len = run_end - col;

            if col > 0 && cols_after[col - 1] == b'M' {
                let (was, last) = if op == b'D' {
                    (text[text_i - 1], text[text_i + len - 1])
                } else {
                    (read[query_i - 1], read[query_i + len - 1])
                };

                assert_ne!(last, was, "case {case}: gap not left-aligned in {out}");
            }

            if op == b'D' {
                text_i += len;
            } else {
                query_i += len;
            }

            col = run_end;
        }
    }
}

#[test]
fn left_align_moves_gap_past_mismatch_in_run() {
    // HG002 chr20:59112821-59112830 holds a 1-base deletion and an A>T SNP in an A run.
    // The reference reads GATT A12 TGT and the read GATT A7 T A3 TGT.
    // The aligner may write the T as a mismatch beside the deletion, and normalization puts the deletion at the start of the run and the T at the SNP, as GIAB and BWA-MEM write it.
    let text = b"GATTAAAAAAAAAAAATGT";
    let read = b"GATTAAAAAAATAAATGT";

    assert_eq!(left_align("12M1D6M", read, text, 0), "4M1D14M");
    assert_eq!(
        cigar_ed_implied("12M1D6M", read, text, 0),
        cigar_ed_implied("4M1D14M", read, text, 0)
    );
}

#[test]
fn polish_gap_clusters_replaces_cancelling_indels() {
    // The read is the reference with 10 bases inserted at 100.
    // A path forced through 5D ... 15I is replaced by a cheaper free alignment of the window that consumes the same bases with fewer edits.
    let mut rng = rng_state(5);
    let reference: Vec<u8> = (0..200)
        .map(|_| b"ACGT"[(next_rng(&mut rng) % 4) as usize])
        .collect();
    let mut read = reference[..100].to_vec();

    read.extend((0..10).map(|_| b"ACGT"[(next_rng(&mut rng) % 4) as usize]));
    read.extend_from_slice(&reference[100..]);
    let ops = vec![(95u32, b'M'), (5, b'D'), (10, b'M'), (15, b'I'), (90, b'M')];
    let (out, saved) = polish_gap_clusters(&ops, &read, &reference, 0, 0);

    assert!(saved > 0, "{out:?}");
    let cig = |run_ops: &[(u32, u8)]| {
        run_ops
            .iter()
            .map(|(run_len, op)| format!("{run_len}{}", *op as char))
            .collect::<String>()
    };
    let read_bases: u32 = out
        .iter()
        .filter(|run| run.1 != b'D')
        .map(|run| run.0)
        .sum();
    let ref_bases: u32 = out
        .iter()
        .filter(|run| run.1 != b'I')
        .map(|run| run.0)
        .sum();

    assert_eq!(
        (read_bases as usize, ref_bases as usize),
        (read.len(), reference.len())
    );
    let (before, after) = (
        cigar_ed_implied(&cig(&ops), &read, &reference, 0).unwrap(),
        cigar_ed_implied(&cig(&out), &read, &reference, 0).unwrap(),
    );

    assert!(after < before, "{before} -> {after}: {}", cig(&out));
    // An alignment with a single gap is left alone.
    let one = vec![(100u32, b'M'), (10, b'I'), (100, b'M')];

    assert_eq!(polish_gap_clusters(&one, &read, &reference, 0, 0).1, 0);
    // 1-2 bp gaps from sequencing noise do not form clusters by themselves.
    let noise = vec![
        (95u32, b'M'),
        (1, b'D'),
        (20, b'M'),
        (2, b'I'),
        (20, b'M'),
        (1, b'D'),
        (63, b'M'),
    ];

    assert_eq!(polish_gap_clusters(&noise, &read, &reference, 0, 0).1, 0);
}

/// The scalar CIGAR DP (the fallback when the wavefront aligner gives up) on hand-checked cases.
#[test]
fn scalar_cigar_trivials() {
    // Exact window placement -> single M run, zero edits, start = 5.
    let window = b"AAAAACGGTACTTCCG".to_vec();
    let (cig, start) = cigar(b"CGGTAC", &window);

    assert_eq!((cig.as_str(), start), ("6M", 5));

    // For a query with no exact placement the CIGAR must consume the query and imply the optimal (minimum) semi-global edit distance.
    let read2 = b"ACGTACGT".to_vec();
    let (sc2, ss2) = cigar(&read2, &window);

    assert_eq!(
        cigar_edit_distance(&sc2, &read2, &window, ss2),
        Some(naive_semi_global_ed(&read2, &window)),
        "cigar={sc2} start={ss2}"
    );

    // Empty query is a special-cased "*".
    let (sc3, ss3) = cigar(b"", &window);

    assert_eq!((sc3.as_str(), ss3), ("*", 0));
}

/// MAPQ vectors from the Lean model of posterior_mapq_w (lean/Vectors/Main.lean), regenerated with `pixi run vectors`.
/// The Lean side proves the posterior identity, its monotonicity and the MAPQ 0 rule for that model, and this test ties the model to the Rust function.
/// Cases within 1e-6 of a rounding boundary are left out by the generator, because Lean's and Rust's powi may differ in the last bit.
#[test]
fn posterior_mapq_matches_lean_vectors() {
    use super::posterior_mapq_w;
    let cases: Vec<serde_json::Value> =
        serde_json::from_str(include_str!("../../../lean/vectors/mapq.json")).expect("valid vector JSON");

    assert!(!cases.is_empty(), "no MAPQ vectors");

    for case in &cases {
        let best = case["best_ed"].as_i64().expect("best_ed") as i32;
        let competitors: Vec<(i32, f64)> = case["competitors"]
            .as_array()
            .expect("competitors")
            .iter()
            .map(|pair| (pair[0].as_i64().expect("ed") as i32, pair[1].as_f64().expect("weight")))
            .collect();
        let len = case["len"].as_u64().expect("len") as usize;
        let hidden = case["hidden"].as_f64().expect("hidden");
        let expect = case["mapq"].as_i64().expect("mapq") as i32;

        assert_eq!(posterior_mapq_w(best, &competitors, len, hidden), expect, "vector {case}");
    }
}

/// `lean/vectors/cigar.json` lists `Cigar.edImplied`, `Cigar.refSpan` and `Cigar.colCount` for random M/I/D CIGARs, some of which the code must reject.
/// The CIGARs are strings, so the test covers the code's digit parsing as well.
#[test]
fn cigar_helpers_match_lean_vectors() {
    use super::paired::ref_span;
    let cases: Vec<serde_json::Value> =
        serde_json::from_str(include_str!("../../../lean/vectors/cigar.json")).expect("valid vector JSON");

    assert!(!cases.is_empty(), "no CIGAR vectors");

    for case in &cases {
        let cig = case["cigar"].as_str().expect("cigar");
        let query = case["query"].as_str().expect("query").as_bytes();
        let text = case["text"].as_str().expect("text").as_bytes();
        let toff = case["toff"].as_u64().expect("toff") as usize;
        // A null NM is a CIGAR the code must reject.
        let expect_nm = case["nm"].as_i64().map(|nm| nm as i32);

        assert_eq!(cigar_ed_implied(cig, query, text, toff), expect_nm, "NM of vector {case}");
        assert_eq!(ref_span(cig), case["ref_span"].as_i64().expect("ref_span"), "span of vector {case}");
        assert_eq!(cigar_cols(cig) as u64, case["cols"].as_u64().expect("cols"), "columns of vector {case}");
    }
}

/// `lean/vectors/bpm.json` gives each text's least semi-global edit distance by the recurrence of `Bpm.D`, which `Bpm.bpm_eq_dp` proves `bpm_semiglobal_ed` returns.
/// Ends are compared only within `max_ed`, since past the cap the code reports no placement.
#[test]
fn placements_match_lean_vectors() {
    use crate::native::placements_x4;
    let cases: Vec<serde_json::Value> =
        serde_json::from_str(include_str!("../../../lean/vectors/bpm.json")).expect("valid vector JSON");

    assert!(!cases.is_empty(), "no placement vectors");

    for case in &cases {
        let query = case["query"].as_str().expect("query").as_bytes();
        let max_ed = case["max_ed"].as_i64().expect("max_ed") as i32;
        let lanes = case["texts"].as_array().expect("texts");
        let texts: Vec<&[u8]> = lanes.iter().map(|lane| lane["text"].as_str().expect("text").as_bytes()).collect();
        // One text takes the scalar path, and two or more the AVX2 lanes.
        let got = placements_x4(query, &texts, max_ed, 0);

        for (i, lane) in lanes.iter().enumerate() {
            let best = lane["best"].as_i64().expect("best") as i32;

            assert_eq!(got[i].best, best, "distance of text {i} in vector {case}");

            if best <= max_ed {
                assert_eq!(got[i].end as i64, lane["end"].as_i64().expect("end"), "end of text {i} in vector {case}");
            }
        }
    }
}
