#!/usr/bin/env python3
"""test_harness.py - known-answer tests of the validation and benchmark harness.

The validations and tool_bench.py's parity checks are code themselves.
Each helper they rely on is tested here against hand-computed answers, so that a harness bug cannot silently turn a real difference into a PASS.

Run with python3 -m unittest bench/test_harness.py, or python3 bench/test_harness.py.
"""

import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import dbg_oracles
import edge_cases
import gen_reads
import hw
import validate_dbg
import validate_mapping as vm
import validate_official as vo
import validate_query as vq
import validate_unitigs as vu


def write(path, text, gz=False):
    """Write text to path, gzip-compressed when gz is set."""
    op = gzip.open if gz else open
    with op(path, "wt") as fh:
        fh.write(text)


class Tmp(unittest.TestCase):
    """Test case with a fresh temporary directory per test."""

    def setUp(self):
        self.td = tempfile.mkdtemp(prefix="harness_")

    def tearDown(self):
        shutil.rmtree(self.td, ignore_errors=True)

    def p(self, name):
        """Path of `name` inside the test's temporary directory."""
        return os.path.join(self.td, name)


class Readers(Tmp):
    """FASTA and FASTQ readers shared by the validators."""

    def test_validate_dbg_reader_fasta_fastq_gz_crlf(self):
        write(
            self.p("a.fq.gz"), "@r1\nACGT\n+\n@@@@\n@r2\nGG\nTT\n+\n@@\n>>\n", gz=True
        )
        self.assertEqual(
            validate_dbg.read_fasta(self.p("a.fq.gz")), [("r1", "ACGT"), ("r2", "GGTT")]
        )
        with open(self.p("b.fa"), "w", newline="") as fh:
            fh.write(">x desc\r\nAC\r\nGT\r\n>y\r\nTT\r\n")
        self.assertEqual(
            validate_unitigs_reader(self.p("b.fa")), [("x", "ACGT"), ("y", "TT")]
        )


def validate_unitigs_reader(p):
    """Records of p as validate_unitigs reads them."""
    return list(vu.read_fasta(p))


class KmerOracles(Tmp):
    """k-mer enumeration and the k-mer and edge filter oracles."""

    def test_kmers_of_counts_invalid_and_canonical(self):
        total, ks = edge_cases.kmers_of("ACGTNac", 2)
        self.assertEqual(total, 6)  # 7 - 2 + 1 positions
        self.assertEqual(
            ks, ["AC", "CG", "AC", "AC"]
        )  # GT->AC, TN and NA invalid, ac->AC

    def test_filter_oracles_kmer_vs_edge_models(self):
        # At k = 3 and c = 2, AAA occurs once in each read but every 4-mer only once overall, so the k-mer filter keeps AAA while the edge filter keeps nothing.
        write(self.p("r.fa"), ">a\nAAAC\n>b\nGAAA\n")
        kf, ef = dbg_oracles.filter_oracles(self.p("r.fa"), 3, 2)
        self.assertEqual(kf, {"AAA"})
        self.assertEqual(ef, set())
        # The shared edge AAAC (twice) keeps both endpoints through the edge filter.
        write(self.p("s.fa"), ">a\nAAAC\n>b\nAAAC\n>c\nAAAG\n")
        kf, ef = dbg_oracles.filter_oracles(self.p("s.fa"), 3, 2)
        self.assertEqual(kf, {"AAA", "AAC"})
        self.assertEqual(ef, {"AAA", "AAC"})
        # Wrapped FASTA records are joined, not split at line ends.
        write(self.p("w.fa"), ">a\nAA\nAC\n>b\nAAAC\n")
        kf, _ = dbg_oracles.filter_oracles(self.p("w.fa"), 3, 2)
        self.assertEqual(kf, {"AAA", "AAC"})


class UnitigValidator(Tmp):
    """validate_unitigs on correct GFAs and on GFAs with one planted defect each."""

    def test_correct_and_broken_gfas(self):
        k = 7
        rng = random.Random(3)
        while (
            True
        ):  # a read whose canonical k-mers are all distinct forms one simple path
            seq = "".join(rng.choice("ACGT") for _ in range(40))
            ks = [edge_cases.canon(seq[i : i + k]) for i in range(len(seq) - k + 1)]
            if len(set(ks)) == len(ks):
                break
        write(self.p("in.fa"), f">a\n{seq}\n")
        res = vu.validate(self.p("in.fa"), self._gfa([seq]), k, 1)
        self.assertTrue(
            all(c["pass"] for c in res.values()),
            {n: c for n, c in res.items() if not c["pass"]},
        )
        # Split into two overlapping segments joined by a link, which is not maximal
        cut = 20
        split = vu.validate(
            self.p("in.fa"),
            self._gfa([seq[: cut + k - 1], seq[cut:]], [(0, "+", 1, "+")]),
            k,
            1,
        )
        self.assertFalse(split["maximal"]["pass"])
        # A missing k-mer, a duplicated k-mer and a spurious link are each caught.
        self.assertFalse(
            vu.validate(self.p("in.fa"), self._gfa([seq[:-1]]), k, 1)["kmer_set_equal"][
                "pass"
            ]
        )
        self.assertFalse(
            vu.validate(self.p("in.fa"), self._gfa([seq, seq[3 : 3 + k]]), k, 1)[
                "no_duplicates"
            ]["pass"]
        )
        self.assertFalse(
            vu.validate(self.p("in.fa"), self._gfa([seq], [(0, "+", 0, "+")]), k, 1)[
                "links_exact"
            ]["pass"]
        )
        # A wrong KC, a missing tag and a wrong LN are each caught.
        self.assertFalse(
            vu.validate(self.p("in.fa"), self._gfa([seq], kc_delta=1), k, 1)[
                "tags_exact"
            ]["pass"]
        )
        self.assertFalse(
            vu.validate(self.p("in.fa"), self._gfa([seq], tags=False), k, 1)[
                "tags_exact"
            ]["pass"]
        )
        self.assertFalse(
            vu.validate(self.p("in.fa"), self._gfa([seq], ln_delta=1), k, 1)[
                "tags_exact"
            ]["pass"]
        )
        # A repeated segment name, which a dict reader would drop, and a segment with a non-ACGT base fail well-formedness instead of crashing.
        dup = self._gfa([seq[: cut + k - 1], seq[cut:]], names=["0", "0"])
        self.assertFalse(
            vu.validate(self.p("in.fa"), dup, k, 1)["segments_wellformed"]["pass"]
        )
        self.assertFalse(
            vu.validate(self.p("in.fa"), self._gfa([seq[:10] + "R" + seq[11:]]), k, 1)[
                "segments_wellformed"
            ]["pass"]
        )

    def test_kdbg_node_table(self):
        import struct

        k = 5
        seq = "ACGTTGCATGCCATAGGCTTAACG"
        write(self.p("in.fa"), f">a\n{seq}\n>b\n{seq[3:15]}\n")
        solid, kp1, _linked, counts = vu.oracle(self.p("in.fa"), k)
        kmask = (1 << (2 * k)) - 1

        def has_edge(cur, nxt):
            """True when edge cur -> nxt is in the oracle's edge-centric graph."""
            m = (cur << 2) | (nxt & 3)
            return (m in kp1 or vu.rc(m, k + 1) in kp1) and vu.canon(nxt, k) in solid

        keys = sorted(solid)
        masks = []
        for c in keys:
            out = sum(1 << b for b in range(4) if has_edge(c, ((c << 2) | b) & kmask))
            inn = sum(
                1 << b
                for b in range(4)
                if has_edge(vu.rc(c, k), vu.rc((b << (2 * (k - 1))) | (c >> 2), k))
            )
            masks.append((inn << 4) | out)

        def table(ab=None, mk=None, ks=None):
            """Bytes of an OKKDBG node table with the given keys, abundances and masks."""
            ks = ks or keys
            body = (
                b"OKKDBG"
                + bytes([k])
                + struct.pack("<Q", len(ks))
                + struct.pack(f"<{len(ks)}Q", *ks)
                + struct.pack(f"<{len(ks)}I", *(ab or [counts[c] for c in ks]))
                + bytes(mk or masks)
            )
            write_bin = self.p("t.kdbg")
            with open(write_bin, "wb") as fh:
                fh.write(body)
            return vu.validate_kdbg(write_bin, k, solid, counts, has_edge)

        self.assertTrue(table()["pass"], table())
        bad_ab = [counts[c] for c in keys]
        bad_ab[3] += 1
        self.assertFalse(table(ab=bad_ab)["pass"])
        bad_mk = list(masks)
        bad_mk[0] ^= 1
        self.assertFalse(table(mk=bad_mk)["pass"])
        self.assertFalse(
            table(ks=keys[:-1], ab=[counts[c] for c in keys[:-1]], mk=masks[:-1])[
                "pass"
            ]
        )

    def test_kdbg_masks_with_cutoff(self):
        """With -c 2 a predecessor below the cutoff is not a node, so no mask bit may point to it, as an in-mask built from raw windows would."""
        k = 5
        core = "ACGTTGCATGCCATAG"
        write(self.p("in.fa"), f">a\nT{core}\n>b\n{core}\n")  # T+core[:4] seen once
        solid, kp1, _, counts = vu.oracle(self.p("in.fa"), k, 2)
        kmask = (1 << (2 * k)) - 1

        def has_edge(cur, nxt):
            """True when edge cur -> nxt is in the oracle's edge-centric graph."""
            m = (cur << 2) | (nxt & 3)
            return (m in kp1 or vu.rc(m, k + 1) in kp1) and vu.canon(nxt, k) in solid

        # The oracle's own masks must point only at solid neighbours.
        import struct

        keys = sorted(solid)
        orig = vu.validate_kdbg
        for c in keys:
            for b in range(4):
                nb_out = ((c << 2) | b) & kmask
                nb_in = (b << (2 * (k - 1))) | (c >> 2)
                if has_edge(c, nb_out):
                    self.assertIn(vu.canon(nb_out, k), solid)
                if has_edge(vu.rc(c, k), vu.rc(nb_in, k)):
                    self.assertIn(vu.canon(nb_in, k), solid)
        # The validator must reject a table whose in-mask includes the predecessor T+core[:4] -> core[:5], which is below the cutoff.
        first = vu.canon(vu.pack(core[:5]), k)
        self.assertIn(first, solid)
        self.assertNotIn(vu.canon(vu.pack("T" + core[:4]), k), solid)
        masks = []
        for c in keys:
            out = sum(1 << b for b in range(4) if has_edge(c, ((c << 2) | b) & kmask))
            inn = sum(
                1 << b
                for b in range(4)
                if has_edge(vu.rc(c, k), vu.rc((b << (2 * (k - 1))) | (c >> 2), k))
            )
            masks.append((inn << 4) | out)

        def run(mk):
            body = (
                b"OKKDBG"
                + bytes([k])
                + struct.pack("<Q", len(keys))
                + struct.pack(f"<{len(keys)}Q", *keys)
                + struct.pack(f"<{len(keys)}I", *[counts[c] for c in keys])
                + bytes(mk)
            )
            with open(self.p("c.kdbg"), "wb") as fh:
                fh.write(body)
            return orig(self.p("c.kdbg"), k, solid, counts, has_edge)["pass"]

        self.assertTrue(run(masks))
        bad = list(masks)
        i = keys.index(first)
        bad[i] |= (
            0xF0 if first == vu.pack(core[:5]) else 0x0F
        )  # edge to the dropped k-mer
        self.assertFalse(run(bad))

    def _gfa(self, segs, links=(), tags=True, kc_delta=0, ln_delta=0, names=None, k=7):
        """Write a GFA of segs and links, with tags optionally skewed by kc_delta and ln_delta."""
        counts = vu.oracle(self.p("in.fa"), k)[3]
        lines = ["H\tVN:Z:1.0"]
        for i, s in enumerate(segs):
            name = names[i] if names else str(i)
            line = f"S\t{name}\t{s}"
            if tags and set(s) <= set("ACGT") and len(s) >= k:
                kc = (
                    sum(
                        counts.get(vu.canon(vu.pack(s[j : j + k]), k), 0)
                        for j in range(len(s) - k + 1)
                    )
                    + kc_delta
                )
                line += f"\tLN:i:{len(s) + ln_delta}\tKC:i:{kc}\tkm:f:{kc / (len(s) - k + 1):.2f}"
            lines.append(line)
        lines += [f"L\t{a}\t{oa}\t{b}\t{ob}\t{6}M" for a, oa, b, ob in links]
        write(self.p("g.gfa"), "\n".join(lines) + "\n")
        return self.p("g.gfa")


class SamHelpers(Tmp):
    """SAM helpers of the validators, such as NM recomputation, PAF conversion and the record audit."""

    REF = "ACGTACGTTTGCAAGCTTAGGCATCGATCGGATCCA"

    def test_sam_nm_from_cigar(self):
        ref = self.REF
        read = ref[4:14]
        self.assertEqual(vo.sam_nm_from_cigar(read, "10M", ref, 5, "+"), 0)
        mut = read[:3] + ("A" if read[3] != "A" else "C") + read[4:]
        self.assertEqual(vo.sam_nm_from_cigar(mut, "10M", ref, 5, "+"), 1)
        self.assertEqual(
            vo.sam_nm_from_cigar(read[:5] + "GG" + read[5:], "5M2I5M", ref, 5, "+"), 2
        )
        self.assertEqual(vo.sam_nm_from_cigar("N" + read[1:], "10M", ref, 5, "+"), 1)
        self.assertIsNone(vo.sam_nm_from_cigar(read, "10M", ref, len(ref) - 3, "+"))

    def test_one_primary_per_read_keeps_longest(self):
        sam = (
            "@HD\tVN:1.6\n"
            "r\t0\tc\t1\t60\t50M100S\t*\t0\t0\t*\t*\n"
            "r\t0\tc\t900\t60\t60S90M\t*\t0\t0\t*\t*\n"
            "s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
        )
        write(self.p("in.sam"), sam)
        split = vo.one_primary_per_read(self.p("in.sam"), self.p("out.sam"))
        self.assertEqual(split, 1)
        with open(self.p("out.sam")) as fh:
            rows = [line.split("\t") for line in fh if not line.startswith("@")]
        self.assertEqual(
            [(t[0], int(t[1]), t[3]) for t in rows],
            [("r", 0x800, "1"), ("r", 0, "900"), ("s", 4, "0")],
        )

    def test_faidx_matches_samtools(self):
        write(self.p("f.fa"), ">a x\nACGTACGT\nACG\n>b\nTTTT\n")
        vo.faidx(self.p("f.fa"))
        mine = Path(self.p("f.fa.fai")).read_text()
        samtools = shutil.which("samtools")
        if samtools:
            os.remove(self.p("f.fa.fai"))
            subprocess.run([samtools, "faidx", self.p("f.fa")], check=True)
            self.assertEqual(mine, Path(self.p("f.fa.fai")).read_text())
        else:
            self.assertEqual(mine, "a\t11\t5\t8\t9\nb\t4\t25\t4\t5\n")

    def test_sam_to_paf_minus_strand_clip_coordinates(self):
        # A read of length 10 on the reverse strand whose SAM SEQ, in reference orientation, has 3 leading soft-clipped bases.
        # On the read as sequenced those are the last 3 bases, so the aligned query interval is [0, 7).
        sam = "q\t16\tc\t5\t60\t3S7M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        paf = vm.sam_to_paf(sam, {"c": 100})
        t = paf[0].split("\t")
        self.assertEqual((t[2], t[3], t[4], t[7], t[8]), ("0", "7", "-", "4", "11"))
        sam_f = "q\t0\tc\t5\t60\t3S7M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        t = vm.sam_to_paf(sam_f, {"c": 100})[0].split("\t")
        self.assertEqual((t[2], t[3], t[4]), ("3", "10", "+"))

    def test_audit_records_catches_each_defect(self):
        ref = {"c": "ACGTACGTTAGCCATGCATGACGTTTGACCATG"}
        read = ref["c"][4:24]
        good = f"r\t0\tc\t5\t60\t20M\t*\t0\t0\t{read}\t*\tNM:i:0"
        rcread = read.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        a = vm.audit_records(good, {"r": read}, ref)
        self.assertEqual(
            (a["records"], a["seq_bad"], a["cigar_bad"], a["nm_bad"]), (1, 0, 0, 0)
        )
        # On the reverse strand SEQ is the reverse complement of the read.
        self.assertEqual(
            vm.audit_records(good.replace("\t0\tc", "\t16\tc"), {"r": rcread}, ref)[
                "seq_bad"
            ],
            0,
        )
        bad = {
            "seq": good.replace(
                read, read[:5] + ("A" if read[5] != "A" else "C") + read[6:]
            ),
            "seq*": good.replace(read, "*"),
            "cigar": good.replace("20M", "19M"),
            "nm_missing": good.replace("\tNM:i:0", ""),
            "nm_wrong": good.replace("NM:i:0", "NM:i:1"),
            "span": good.replace("\tc\t5\t", "\tc\t30\t"),
        }
        a = {n: vm.audit_records(line, {"r": read}, ref) for n, line in bad.items()}
        self.assertTrue(a["seq"]["seq_bad"] and a["seq"]["nm_bad"])  # mismatching base
        self.assertTrue(a["seq*"]["seq_bad"])
        self.assertTrue(a["cigar"]["cigar_bad"])
        self.assertTrue(a["nm_missing"]["nm_bad"])
        self.assertTrue(a["nm_wrong"]["nm_bad"])
        self.assertTrue(a["span"]["span_bad"])

    def test_sam_to_paf_reports_unconvertible_records(self):
        sam = (
            "r1\t0\tc\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
            "r2\t0\tc\t1\t60\t8M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        dropped = []
        lines = vm.sam_to_paf(sam, {"c": 100}, dropped)
        self.assertEqual((len(lines), dropped), (1, ["r2"]))

    def test_parse_sam_skips_secondary_and_supplementary(self):
        sam = (
            "r\t2048\tc\t1\t60\t5M\t*\t0\t0\tACGTA\t*\n"
            "r\t256\tc\t9\t60\t5M\t*\t0\t0\tACGTA\t*\n"
            "r\t16\tc\t20\t60\t5M\t*\t0\t0\tACGTA\t*\n"
        )
        recs = vm.parse_sam(sam)
        self.assertEqual((recs["r"]["pos"], recs["r"]["is_reverse"]), (20, True))


class TruthScoring(Tmp):
    """Truth scoring of mapped reads in the align category and in scaling_bench."""

    def test_align_truth_metrics(self):
        sys.path.insert(0, os.path.join(HERE, "categories"))
        from categories.align import _truth_metrics

        reads = self.p("reads.fa")
        write(reads, ">a\nX\n")
        write(
            reads + ".truth.tsv",
            "read_name\tref_name\tpos\tstrand\tread_length\n"
            "a\tc\t100\t+\t150\n"  # correct, soft-clipped start
            "b\tc\t500\t-\t150\n"  # wrong strand
            "d\tc\t900\t+\t150\n"  # off by 40 > tolerance 15
            "e\tc\t50\t+\t150\n",
        )  # unmapped
        sam = (
            "a\t0\tc\t106\t60\t5S145M\t*\t0\t0\t*\t*\n"
            "b\t0\tc\t501\t60\t150M\t*\t0\t0\t*\t*\n"
            "d\t16\tc\t941\t35\t150M\t*\t0\t0\t*\t*\n"
            "e\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
        )
        m = _truth_metrics(sam, reads)
        self.assertEqual(
            (m["correct_reads"], m["wrong_reads"], m["unmapped_reads"]), (1, 2, 1)
        )
        self.assertEqual(m["mapq30_wrong"], 2)
        self.assertEqual(m["mapq30_reads"], 3)

    def test_scaling_map_accuracy(self):
        import scaling_bench

        tsv = self.p("r.tsv")
        write(
            tsv,
            "read_name\tref\tpos\tstrand\na\tc\t100\t+\nb\tc\t500\t-\nd\tc\t900\t+\ne\tc\t50\t+\n",
        )
        sam = self.p("o.sam")
        write(
            sam,
            "@HD\tVN:1.6\n"
            "a\t0\tc\t106\t60\t5S145M\t*\t0\t0\t*\t*\n"  # correct, MAPQ 60
            "a\t2048\tc\t1\t60\t5M145S\t*\t0\t0\t*\t*\n"  # supplementary, so ignored
            "b\t0\tc\t501\t60\t150M\t*\t0\t0\t*\t*\n"  # wrong strand, confident
            "d\t16\tc\t941\t12\t150M\t*\t0\t0\t*\t*\n"  # wrong, not confident
            "e\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n",
        )
        # (correct, total, MAPQ>=30 wrong, MAPQ>=30)
        self.assertEqual(scaling_bench.map_accuracy(sam, tsv), (1, 4, 1, 2))


class ParityChecks(Tmp):
    """tool_bench's query and graph parity checks."""

    def test_query_parity(self):
        import tool_bench

        same = {
            "a": {"query_kmers": 10, "query_positive": 4},
            "b": {"query_kmers": 10, "query_positive": 4},
        }
        self.assertTrue(tool_bench.query_parity(same)["pass"])
        diff = {
            "a": {"query_kmers": 10, "query_positive": 4},
            "b": {"query_kmers": 10, "query_positive": 5},
        }
        self.assertFalse(tool_bench.query_parity(diff)["pass"])
        self.assertFalse(tool_bench.query_parity({"a": {"query_kmers": 1}})["pass"])

    def test_graph_parity_per_model(self):
        import tool_bench

        write(self.p("r.fa"), ">a\nAAAC\n>b\nGAAA\n")
        ds = {
            "k": 3,
            "dbg_input": self.p("r.fa"),
            "queries": self.p("r.fa"),
            "dbg_min_count": 2,
        }

        def tool(name, kmers):
            """A graph tool's result dict whose dumped k-mer set is `kmers`."""
            f = self.p(name + ".kmers")
            write(f, "\n".join(sorted(kmers)))
            return {
                "graph_kmers_file": f,
                "graph_kmers": len(kmers),
                "graph_kmer_occurrences": len(kmers),
                "min_count": 2,
            }

        ok = {
            "optik": tool("optik", {"AAA"}),
            "bifrost": tool("bifrost", {"AAA"}),
            "cuttlefish": tool("cuttlefish", set()),
        }
        self.assertTrue(tool_bench.graph_parity(ok, ds)["pass"])
        bad = dict(ok, optik=tool("optik2", set()))
        self.assertFalse(tool_bench.graph_parity(bad, ds)["pass"])
        dup = dict(ok, optik=dict(tool("optik3", {"AAA"}), graph_kmer_occurrences=2))
        self.assertFalse(tool_bench.graph_parity(dup, ds)["pass"])


class ReportParsing(unittest.TestCase):
    """The query-report parsers of the different harnesses agree."""

    TEXT = (
        "==== query report:\nnum_kmers = 120\n"
        "num_positive_kmers = 100 (83.3333%)\nnum_negative_kmers = 20 (16.6%)\n"
    )

    def test_parsers_agree(self):
        sys.path.insert(0, os.path.join(HERE, "categories"))
        from categories.ssh import _parse_report_field

        self.assertEqual(_parse_report_field(self.TEXT, "num_positive_kmers"), 100)
        self.assertEqual(vq.parse_counts(self.TEXT)["num_positive_kmers"], 100)
        self.assertEqual(vo.parse_optik_report(self.TEXT)["num_kmers"], 120)


class ReadSimulator(Tmp):
    """gen_reads.py's truth file."""

    def test_gen_reads_truth_is_exact_for_error_free_reads(self):
        rng = random.Random(1)
        ref = "".join(rng.choice("ACGT") for _ in range(3000))
        ref = ref[:1000] + "N" * 20 + ref[1020:]
        reads, truth = gen_reads.generate_reads([("c", ref)], 300, 100, 0.0, 7)
        self.assertEqual(len(reads), len(truth))
        strands = set()
        for (_, read), (_, _rname, pos, strand, span) in zip(
            reads, truth, strict=False
        ):
            window = ref[pos : pos + span]
            self.assertNotIn("N", window)
            self.assertEqual(
                read, window if strand == "+" else gen_reads.revcomp(window)
            )
            strands.add(strand)
        self.assertEqual(strands, {"+", "-"})


def undefined_globals(src, name="<src>"):
    """Global names a module references but never defines, imports or takes from builtins.

    This is a small pyflakes built on the standard library's symtable.
    The bench runners execute only inside a full benchmark run, so a name lost in a refactor, as fmsi's _ensure_uncompressed once was, would otherwise show up only as a crashed benchmark row."""
    import builtins
    import symtable

    top = symtable.symtable(src, name, "exec")
    defined = {
        sym.get_name()
        for sym in top.get_symbols()
        if sym.is_assigned() or sym.is_imported() or sym.is_namespace()
    }
    if "*" in defined or "import *" in src:
        return set()
    known = (
        defined
        | set(dir(builtins))
        | {
            "__file__",
            "__name__",
            "__doc__",
            "__spec__",
            "__builtins__",
            "__loader__",
            "__package__",
        }
    )
    missing = set()

    def walk(table):
        """Collect referenced globals that are not known, recursing into child scopes."""
        for sym in table.get_symbols():
            if sym.is_global() and sym.is_referenced() and sym.get_name() not in known:
                missing.add(sym.get_name())
        for child in table.get_children():
            walk(child)

    for child in top.get_children():
        walk(child)
    for sym in top.get_symbols():  # uses at module level
        if (
            sym.is_referenced()
            and not sym.is_assigned()
            and not sym.is_imported()
            and sym.get_name() not in known
        ):
            missing.add(sym.get_name())
    return missing


class StaticNames(unittest.TestCase):
    """The undefined-name checker, and the bench modules checked with it."""

    def test_checker_catches_a_lost_helper(self):
        bad = "from ._base import run_cmd\ndef run(p):\n    return _ensure_uncompressed(p, run_cmd)\n"
        self.assertEqual(undefined_globals(bad), {"_ensure_uncompressed"})
        self.assertEqual(
            undefined_globals("import os\ndef f():\n    return os.sep\n"), set()
        )

    def test_no_bench_module_uses_an_undefined_name(self):
        import glob

        bad = {}
        for path in sorted(
            glob.glob(os.path.join(HERE, "*.py"))
            + glob.glob(os.path.join(HERE, "categories", "*.py"))
        ):
            miss = undefined_globals(Path(path).read_text(), path)
            if miss:
                bad[os.path.relpath(path, HERE)] = sorted(miss)
        self.assertEqual(bad, {}, "undefined global names (crash when that code runs)")


class OracleKernels(unittest.TestCase):
    """The fast reverse-complement and pack kernels of the dBG oracles against base-by-base versions."""

    def test_table_reverse_complement_matches_base_by_base(self):
        # The dBG oracles' rc() is a byte-table shortcut, which must equal the plain per-base definition at every width the validators use.
        import random

        import validate_unitigs as vu

        rnd = random.Random(7)
        for w in range(1, 33):
            vals = [0, (1 << (2 * w)) - 1] + [
                rnd.getrandbits(2 * w) for _ in range(2000)
            ]
            for v in vals:
                self.assertEqual(vu.rc(v, w), vu._rc_slow(v, w), (w, v))

    def test_base4_pack_matches_base_by_base(self):
        import random

        import validate_unitigs as vu

        rnd = random.Random(11)
        for n in range(40):
            for _ in range(200):
                s = "".join(rnd.choice("ACGT") for _ in range(n))
                self.assertEqual(vu.pack(s), vu._pack_slow(s), s)


class MicroBenchGuard(unittest.TestCase):
    """micro_bench must refuse to report kernel numbers when the counters are unusable or the calibration kernels fall out of range."""

    GOOD = (
        {"kernel": "nop", "instructions": 1.03, "cycles": 0.17, "l1d_misses": 0},
        {
            "kernel": "add_dependent",
            "instructions": 1.03,
            "cycles": 1.0,
            "l1d_misses": 0,
        },
        {
            "kernel": "add_independent",
            "instructions": 4.03,
            "cycles": 1.0,
            "l1d_misses": 0,
        },
        {"kernel": "chase_l1", "instructions": 1.03, "cycles": 5.0, "l1d_misses": 0.0},
        {
            "kernel": "chase_dram",
            "instructions": 1.03,
            "cycles": 450.0,
            "l1d_misses": 1.0,
        },
    )

    def test_good_calibration_is_usable(self):
        import micro_bench

        self.assertTrue(micro_bench.evaluate([dict(r) for r in self.GOOD])["usable"])

    def test_counter_error_is_unusable(self):
        import micro_bench

        rows = [dict(r) for r in self.GOOD] + [
            {"kernel": "x", "error": "cycles: Permission denied"}
        ]
        r = micro_bench.evaluate(rows)
        self.assertFalse(r["usable"])
        self.assertIn("counters unavailable", r["reason"])

    def test_off_calibration_is_unusable(self):
        import micro_bench

        rows = [dict(r) for r in self.GOOD]
        rows[1]["cycles"] = 2.5  # as with multiplexed or scaled counts
        r = micro_bench.evaluate(rows)
        self.assertFalse(r["usable"])
        self.assertIn("add_dependent", r["reason"])

    def test_cached_dram_chase_is_unusable(self):
        import micro_bench

        rows = [dict(r) for r in self.GOOD]
        rows[4]["cycles"] = 62.0  # the L3-resident chase, restarting at the head
        r = micro_bench.evaluate(rows)
        self.assertFalse(r["usable"])
        self.assertIn("chase_dram", r["reason"])


def sampled_kmers(path, k, mod=97):
    """Canonical k-mers of a FASTA(.gz) whose hash % mod == 0, per record with wrapped lines joined.

    This is a fixed 1/mod sample, the same k-mers in every file."""
    import gzip as _gz
    import zlib

    comp = str.maketrans("ACGT", "TGCA")
    out, seqs, cur = set(), [], []
    op = _gz.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    for s in seqs:
        r = s.translate(comp)[::-1]
        n = len(s)
        for i in range(n - k + 1):
            a, b = s[i : i + k], r[n - k - i : n - i]
            c = min(b, a)
            if "N" not in c and zlib.crc32(c.encode()) % mod == 0:
                out.add(c)
    return out


class DatasetConsistency(unittest.TestCase):
    """Every dataset's query, graph and mapping genome must come from the organism its unitigs index.

    The SE datasets once pointed at E. coli (1.8% of its k-mers are in the S. enterica unitigs), so an 'SE' run repeated the E. coli rows under another name."""

    def test_query_genome_matches_unitigs(self):
        from config import DATASETS

        for name, ds in DATASETS.items():
            if not (Path(ds["unitigs"]).exists() and Path(ds["queries"]).exists()):
                continue
            k = ds["k"]
            u = sampled_kmers(ds["unitigs"], k)
            q = sampled_kmers(ds["queries"], k)
            frac = len(u & q) / max(1, len(q))
            self.assertGreater(
                frac,
                0.9,
                f"{name}: only {frac:.1%} of the query genome's k-mers are in the unitigs",
            )
            self.assertEqual(ds.get("dbg_input"), ds["queries"], name)


class Measure(Tmp):
    """bench/measure.py drops warmups, keeps exactly --runs runs, returns the stdout and stderr of the final run, and reports failures untimed."""

    def test_warmups_dropped_runs_kept(self):
        import measure

        cnt = self.p("n")
        # Each run appends one line, so the last run's stdout shows the run count.
        cmd = ["sh", "-c", f"echo x >> {cnt}; wc -l < {cnt}; echo err >&2"]
        r = measure.measure(cmd, 2, 3, scratch=self.p("m"))
        self.assertTrue(r.ok, r.err)
        self.assertEqual((r.warmup, r.runs, len(r.per_run)), (2, 3, 3))
        with open(r.out) as f:
            self.assertEqual(
                f.read().strip(), "5"
            )  # 2 warmups + 3 runs, the last one's output
        self.assertEqual(r.err.strip(), "err")
        self.assertIsNotNone(r.wall["sd"])
        self.assertLessEqual(r.wall["min"], r.wall["mean"])
        self.assertLessEqual(r.wall["mean"], r.wall["max"])

    def test_single_run_has_no_sd(self):
        import measure

        r = measure.measure(["true"], 0, 1, scratch=self.p("m"))
        self.assertTrue(r.ok)
        self.assertIsNone(r.wall["sd"])

    def test_failure_not_timed(self):
        import measure

        r = measure.measure(
            ["sh", "-c", "echo boom >&2; exit 3"], 1, 2, scratch=self.p("m")
        )
        self.assertFalse(r.ok)
        self.assertEqual(r.code, 3)
        self.assertIn("boom", r.err)
        self.assertEqual(r.wall, {})


class MachineLimits(unittest.TestCase):
    """hw.py, the CPU and memory limits every bench sizes itself by."""

    def test_thread_points(self):
        self.assertEqual(hw.thread_points(1, 1), [1])
        self.assertEqual(hw.thread_points(1, 12), [1, 2, 4, 8, 12])
        self.assertEqual(hw.thread_points(1, 16), [1, 2, 4, 8, 16])
        self.assertEqual(hw.thread_points(4, 32), [4, 8, 16, 32])
        self.assertEqual(hw.thread_points(3, 32), [3, 4, 8, 16, 32])
        self.assertEqual(hw.thread_points(8, 4), [4])

    def test_cpus_follow_the_affinity_mask(self):
        # A child pinned to two CPUs sees two, as optik's available_parallelism() would.
        if len(os.sched_getaffinity(0)) < 2:
            self.skipTest("needs two CPUs")
        two = sorted(os.sched_getaffinity(0))[:2]
        code = f"import os, sys; sys.path.insert(0, {HERE!r}); import hw; "
        code += f"os.sched_setaffinity(0, {two!r}); print(hw.cpus())"
        out = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=False
        )
        self.assertEqual(out.stdout.strip(), "2", out.stderr)

    def test_busy_windows_count_cpu_time_not_runnable_threads(self):
        import thread_check

        # Two threads each on a CPU for the whole 0.1 s, sampled every 10 ms: 2 CPUs busy.
        two = [(i / 100, {1: i * 10**7, 2: i * 10**7}) for i in range(11)]
        self.assertAlmostEqual(thread_check.peak_busy(two, 0.05), 2.0)
        # A third thread that is runnable but never gets a CPU adds nothing.
        idle = [(w, {**c, 3: 0}) for w, c in two]
        self.assertAlmostEqual(thread_check.peak_busy(idle, 0.05), 2.0)
        # A thread exiting mid-run loses its last share, never inflating the result.
        gone = [(w, c if w < 0.05 else {1: c[1]}) for w, c in two]
        self.assertLessEqual(thread_check.peak_busy(gone, 0.05), 2.0)

    def test_default_cap_fits_the_machine(self):
        mem = hw.mem_bytes()
        self.assertGreater(mem, 0)
        self.assertLessEqual(hw.default_cap_bytes(), mem)
        self.assertGreaterEqual(hw.default_cap_bytes(), mem // 2)


if __name__ == "__main__":
    unittest.main(verbosity=1)
