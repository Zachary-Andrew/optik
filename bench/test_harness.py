#!/usr/bin/env python3
"""test_harness.py - known-answer tests of the validation and benchmark harness.

The validations and tool_bench.py's parity checks are code themselves.
Each helper they rely on is tested here against hand-computed answers, so that a harness bug cannot silently turn a real difference into a PASS.

Run with python3 -m unittest bench/test_harness.py, or python3 bench/test_harness.py.
"""

import contextlib
import gzip
import hashlib
import http.server
import io
import os
import random
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
import unittest
import unittest.mock
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import dbg_oracles
import edge_cases
import gen_reads
import hw
import validate_dbg
import validate_mapping
import validate_official
import validate_query
import validate_unitigs


def write(path, text, compressed=False):
    """Write text to path, gzip-compressed when gz is set."""
    opener = gzip.open if compressed else open
    with opener(path, "wt") as handle:
        handle.write(text)


class Tmp(unittest.TestCase):
    """Test case with a fresh temporary directory per test."""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp(prefix="harness_")

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def tmp_path(self, name):
        """Path of `name` inside the test's temporary directory."""
        return os.path.join(self.tmp_dir, name)


class Readers(Tmp):
    """FASTA and FASTQ readers shared by the validators."""

    def test_validate_dbg_reader_fasta_fastq_gz_crlf(self):
        write(
            self.tmp_path("a.fq.gz"), "@r1\nACGT\n+\n@@@@\n@r2\nGG\nTT\n+\n@@\n>>\n", compressed=True
        )
        self.assertEqual(
            list(validate_dbg.read_fasta(self.tmp_path("a.fq.gz"))), [("r1", "ACGT"), ("r2", "GGTT")]
        )
        with open(self.tmp_path("b.fa"), "w", newline="") as handle:
            handle.write(">x desc\r\nAC\r\nGT\r\n>y\r\nTT\r\n")
        self.assertEqual(
            validate_unitigs_reader(self.tmp_path("b.fa")), [("x", "ACGT"), ("y", "TT")]
        )


def validate_unitigs_reader(path):
    """Records of p as validate_unitigs reads them."""
    return list(validate_unitigs.read_fasta(path))


class KmerOracles(Tmp):
    """k-mer enumeration and the k-mer and edge filter oracles."""

    def test_kmers_of_counts_invalid_and_canonical(self):
        total, kmers = edge_cases.kmers_of("ACGTNac", 2)
        self.assertEqual(total, 6)  # 7 - 2 + 1 positions
        self.assertEqual(
            kmers, ["AC", "CG", "AC", "AC"]
        )  # GT->AC, TN and NA invalid, ac->AC

    def test_filter_oracles_kmer_vs_edge_models(self):
        # At k = 3 and c = 2, AAA occurs once in each read but every 4-mer only once overall, so the k-mer filter keeps AAA while the edge filter keeps nothing.
        write(self.tmp_path("r.fa"), ">a\nAAAC\n>b\nGAAA\n")
        kmer_filter, edge_filter = dbg_oracles.filter_oracles(self.tmp_path("r.fa"), 3, 2)
        self.assertEqual(kmer_filter, {"AAA"})
        self.assertEqual(edge_filter, set())
        # The shared edge AAAC (twice) keeps both endpoints through the edge filter.
        write(self.tmp_path("s.fa"), ">a\nAAAC\n>b\nAAAC\n>c\nAAAG\n")
        kmer_filter, edge_filter = dbg_oracles.filter_oracles(self.tmp_path("s.fa"), 3, 2)
        self.assertEqual(kmer_filter, {"AAA", "AAC"})
        self.assertEqual(edge_filter, {"AAA", "AAC"})
        # Wrapped FASTA records are joined, not split at line ends.
        write(self.tmp_path("w.fa"), ">a\nAA\nAC\n>b\nAAAC\n")
        kmer_filter, _ = dbg_oracles.filter_oracles(self.tmp_path("w.fa"), 3, 2)
        self.assertEqual(kmer_filter, {"AAA", "AAC"})


class UnitigValidator(Tmp):
    """validate_unitigs on correct GFAs and on GFAs with one planted defect each."""

    def test_correct_and_broken_gfas(self):
        kmer_len = 7
        rng = random.Random(3)
        while (
            True
        ):  # a read whose canonical k-mers are all distinct forms one simple path
            seq = "".join(rng.choice("ACGT") for _ in range(40))
            kmers = [edge_cases.canon(seq[i : i + kmer_len]) for i in range(len(seq) - kmer_len + 1)]
            if len(set(kmers)) == len(kmers):
                break
        write(self.tmp_path("in.fa"), f">a\n{seq}\n")
        res = validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq]), kmer_len, 1)
        self.assertTrue(
            all(check["pass"] for check in res.values()),
            {name: check for name, check in res.items() if not check["pass"]},
        )
        # Split into two overlapping segments joined by a link, which is not maximal
        cut = 20
        split = validate_unitigs.validate(
            self.tmp_path("in.fa"),
            self._gfa([seq[: cut + kmer_len - 1], seq[cut:]], [(0, "+", 1, "+")]),
            kmer_len,
            1,
        )
        self.assertFalse(split["maximal"]["pass"])
        # A missing k-mer, a duplicated k-mer and a spurious link are each caught.
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq[:-1]]), kmer_len, 1)["kmer_set_equal"][
                "pass"
            ]
        )
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq, seq[3 : 3 + kmer_len]]), kmer_len, 1)[
                "no_duplicates"
            ]["pass"]
        )
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq], [(0, "+", 0, "+")]), kmer_len, 1)[
                "links_exact"
            ]["pass"]
        )
        # A wrong KC, a missing tag and a wrong LN are each caught.
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq], kc_delta=1), kmer_len, 1)[
                "tags_exact"
            ]["pass"]
        )
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq], tags=False), kmer_len, 1)[
                "tags_exact"
            ]["pass"]
        )
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq], ln_delta=1), kmer_len, 1)[
                "tags_exact"
            ]["pass"]
        )
        # A repeated segment name, which a dict reader would drop, and a segment with a non-ACGT base fail well-formedness instead of crashing.
        dup = self._gfa([seq[: cut + kmer_len - 1], seq[cut:]], names=["0", "0"])
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), dup, kmer_len, 1)["segments_wellformed"]["pass"]
        )
        self.assertFalse(
            validate_unitigs.validate(self.tmp_path("in.fa"), self._gfa([seq[:10] + "R" + seq[11:]]), kmer_len, 1)[
                "segments_wellformed"
            ]["pass"]
        )

    def test_kdbg_node_table(self):
        import struct

        kmer_len = 5
        seq = "ACGTTGCATGCCATAGGCTTAACG"
        write(self.tmp_path("in.fa"), f">a\n{seq}\n>b\n{seq[3:15]}\n")
        solid, kp1, _linked, counts = validate_unitigs.oracle(self.tmp_path("in.fa"), kmer_len)
        kmask = (1 << (2 * kmer_len)) - 1

        def has_edge(cur, nxt):
            """True when edge cur -> nxt is in the oracle's edge-centric graph."""
            edge = (cur << 2) | (nxt & 3)
            return (edge in kp1 or validate_unitigs.revcomp_packed(edge, kmer_len + 1) in kp1) and validate_unitigs.canon(nxt, kmer_len) in solid

        keys = sorted(solid)
        masks = []
        for key in keys:
            out = sum(1 << b for b in range(4) if has_edge(key, ((key << 2) | b) & kmask))
            inn = sum(
                1 << b
                for b in range(4)
                if has_edge(validate_unitigs.revcomp_packed(key, kmer_len), validate_unitigs.revcomp_packed((b << (2 * (kmer_len - 1))) | (key >> 2), kmer_len))
            )
            masks.append((inn << 4) | out)

        def table(abundances=None, masks_arg=None, keys_arg=None):
            """Bytes of an OKKDBG node table with the given keys, abundances and masks."""
            keys_arg = keys_arg or keys
            body = (
                b"OKKDBG"
                + bytes([kmer_len])
                + struct.pack("<Q", len(keys_arg))
                + struct.pack(f"<{len(keys_arg)}Q", *keys_arg)
                + struct.pack(f"<{len(keys_arg)}I", *(abundances or [counts[key] for key in keys_arg]))
                + bytes(masks_arg or masks)
            )
            write_bin = self.tmp_path("t.kdbg")
            with open(write_bin, "wb") as handle:
                handle.write(body)
            return validate_unitigs.validate_kdbg(write_bin, kmer_len, solid, counts, has_edge)

        self.assertTrue(table()["pass"], table())
        bad_ab = [counts[key] for key in keys]
        bad_ab[3] += 1
        self.assertFalse(table(abundances=bad_ab)["pass"])
        bad_mk = list(masks)
        bad_mk[0] ^= 1
        self.assertFalse(table(masks_arg=bad_mk)["pass"])
        self.assertFalse(
            table(keys_arg=keys[:-1], abundances=[counts[key] for key in keys[:-1]], masks_arg=masks[:-1])[
                "pass"
            ]
        )

    def test_kdbg_masks_with_cutoff(self):
        """With -c 2 a predecessor below the cutoff is not a node, so no mask bit may point to it, as an in-mask built from raw windows would."""
        kmer_len = 5
        core = "ACGTTGCATGCCATAG"
        write(self.tmp_path("in.fa"), f">a\nT{core}\n>b\n{core}\n")  # T+core[:4] seen once
        solid, kp1, _, counts = validate_unitigs.oracle(self.tmp_path("in.fa"), kmer_len, 2)
        kmask = (1 << (2 * kmer_len)) - 1

        def has_edge(cur, nxt):
            """True when edge cur -> nxt is in the oracle's edge-centric graph."""
            edge = (cur << 2) | (nxt & 3)
            return (edge in kp1 or validate_unitigs.revcomp_packed(edge, kmer_len + 1) in kp1) and validate_unitigs.canon(nxt, kmer_len) in solid

        # The oracle's own masks must point only at solid neighbours.
        import struct

        keys = sorted(solid)
        orig = validate_unitigs.validate_kdbg
        for key in keys:
            for b in range(4):
                nb_out = ((key << 2) | b) & kmask
                nb_in = (b << (2 * (kmer_len - 1))) | (key >> 2)
                if has_edge(key, nb_out):
                    self.assertIn(validate_unitigs.canon(nb_out, kmer_len), solid)
                if has_edge(validate_unitigs.revcomp_packed(key, kmer_len), validate_unitigs.revcomp_packed(nb_in, kmer_len)):
                    self.assertIn(validate_unitigs.canon(nb_in, kmer_len), solid)
        # The validator must reject a table whose in-mask includes the predecessor T+core[:4] -> core[:5], which is below the cutoff.
        first = validate_unitigs.canon(validate_unitigs.pack(core[:5]), kmer_len)
        self.assertIn(first, solid)
        self.assertNotIn(validate_unitigs.canon(validate_unitigs.pack("T" + core[:4]), kmer_len), solid)
        masks = []
        for key in keys:
            out = sum(1 << b for b in range(4) if has_edge(key, ((key << 2) | b) & kmask))
            inn = sum(
                1 << b
                for b in range(4)
                if has_edge(validate_unitigs.revcomp_packed(key, kmer_len), validate_unitigs.revcomp_packed((b << (2 * (kmer_len - 1))) | (key >> 2), kmer_len))
            )
            masks.append((inn << 4) | out)

        def run(masks_arg):
            body = (
                b"OKKDBG"
                + bytes([kmer_len])
                + struct.pack("<Q", len(keys))
                + struct.pack(f"<{len(keys)}Q", *keys)
                + struct.pack(f"<{len(keys)}I", *[counts[key] for key in keys])
                + bytes(masks_arg)
            )
            with open(self.tmp_path("c.kdbg"), "wb") as handle:
                handle.write(body)
            return orig(self.tmp_path("c.kdbg"), kmer_len, solid, counts, has_edge)["pass"]

        self.assertTrue(run(masks))
        bad = list(masks)
        first_idx = keys.index(first)
        bad[first_idx] |= (
            0xF0 if first == validate_unitigs.pack(core[:5]) else 0x0F
        )  # edge to the dropped k-mer
        self.assertFalse(run(bad))

    def _gfa(self, segs, links=(), tags=True, kc_delta=0, ln_delta=0, names=None, kmer_len=7):
        """Write a GFA of segs and links, with tags optionally skewed by kc_delta and ln_delta."""
        counts = validate_unitigs.oracle(self.tmp_path("in.fa"), kmer_len)[3]
        lines = ["H\tVN:Z:1.0"]
        for i, seq in enumerate(segs):
            name = names[i] if names else str(i)
            line = f"S\t{name}\t{seq}"
            if tags and set(seq) <= set("ACGT") and len(seq) >= kmer_len:
                kc_value = (
                    sum(
                        counts.get(validate_unitigs.canon(validate_unitigs.pack(seq[j : j + kmer_len]), kmer_len), 0)
                        for j in range(len(seq) - kmer_len + 1)
                    )
                    + kc_delta
                )
                line += f"\tLN:i:{len(seq) + ln_delta}\tKC:i:{kc_value}\tkm:f:{kc_value / (len(seq) - kmer_len + 1):.2f}"
            lines.append(line)
        lines += [f"L\t{from_seg}\t{from_orient}\t{to_seg}\t{to_orient}\t{6}M" for from_seg, from_orient, to_seg, to_orient in links]
        write(self.tmp_path("g.gfa"), "\n".join(lines) + "\n")
        return self.tmp_path("g.gfa")


class SamHelpers(Tmp):
    """SAM helpers of the validators, such as NM recomputation, PAF conversion and the record audit."""

    REF = "ACGTACGTTTGCAAGCTTAGGCATCGATCGGATCCA"

    def test_sam_nm_from_cigar(self):
        ref = self.REF
        read = ref[4:14]
        self.assertEqual(validate_official.sam_nm_from_cigar(read, "10M", ref, 5, "+"), 0)
        mut = read[:3] + ("A" if read[3] != "A" else "C") + read[4:]
        self.assertEqual(validate_official.sam_nm_from_cigar(mut, "10M", ref, 5, "+"), 1)
        self.assertEqual(
            validate_official.sam_nm_from_cigar(read[:5] + "GG" + read[5:], "5M2I5M", ref, 5, "+"), 2
        )
        self.assertEqual(validate_official.sam_nm_from_cigar("N" + read[1:], "10M", ref, 5, "+"), 1)
        self.assertIsNone(validate_official.sam_nm_from_cigar(read, "10M", ref, len(ref) - 3, "+"))

    def test_one_primary_per_read_keeps_longest(self):
        sam = (
            "@HD\tVN:1.6\n"
            "r\t0\tc\t1\t60\t50M100S\t*\t0\t0\t*\t*\n"
            "r\t0\tc\t900\t60\t60S90M\t*\t0\t0\t*\t*\n"
            "s\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
        )
        write(self.tmp_path("in.sam"), sam)
        split = validate_official.one_primary_per_read(self.tmp_path("in.sam"), self.tmp_path("out.sam"))
        self.assertEqual(split, 1)
        with open(self.tmp_path("out.sam")) as handle:
            rows = [line.split("\t") for line in handle if not line.startswith("@")]
        self.assertEqual(
            [(row[0], int(row[1]), row[3]) for row in rows],
            [("r", 0x800, "1"), ("r", 0, "900"), ("s", 4, "0")],
        )

    def test_faidx_matches_samtools(self):
        write(self.tmp_path("f.fa"), ">a x\nACGTACGT\nACG\n>b\nTTTT\n")
        validate_official.faidx(self.tmp_path("f.fa"))
        mine = Path(self.tmp_path("f.fa.fai")).read_text()
        samtools = shutil.which("samtools")
        if samtools:
            os.remove(self.tmp_path("f.fa.fai"))
            subprocess.run([samtools, "faidx", self.tmp_path("f.fa")], check=True)
            self.assertEqual(mine, Path(self.tmp_path("f.fa.fai")).read_text())
        else:
            self.assertEqual(mine, "a\t11\t5\t8\t9\nb\t4\t25\t4\t5\n")

    def test_sam_to_paf_minus_strand_clip_coordinates(self):
        # A read of length 10 on the reverse strand whose SAM SEQ, in reference orientation, has 3 leading soft-clipped bases.
        # On the read as sequenced those are the last 3 bases, so the aligned query interval is [0, 7).
        sam = "q\t16\tc\t5\t60\t3S7M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        paf = validate_mapping.sam_to_paf(sam, {"c": 100})
        fields = paf[0].split("\t")
        self.assertEqual((fields[2], fields[3], fields[4], fields[7], fields[8]), ("0", "7", "-", "4", "11"))
        sam_f = "q\t0\tc\t5\t60\t3S7M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        fields = validate_mapping.sam_to_paf(sam_f, {"c": 100})[0].split("\t")
        self.assertEqual((fields[2], fields[3], fields[4]), ("3", "10", "+"))

    def test_audit_records_catches_each_defect(self):
        ref = {"c": "ACGTACGTTAGCCATGCATGACGTTTGACCATG"}
        read = ref["c"][4:24]
        good = f"r\t0\tc\t5\t60\t20M\t*\t0\t0\t{read}\t*\tNM:i:0"
        rcread = read.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        audit = validate_mapping.audit_records(good, {"r": read}, ref)
        self.assertEqual(
            (audit["records"], audit["seq_bad"], audit["cigar_bad"], audit["nm_bad"]), (1, 0, 0, 0)
        )
        # On the reverse strand SEQ is the reverse complement of the read.
        self.assertEqual(
            validate_mapping.audit_records(good.replace("\t0\tc", "\t16\tc"), {"r": rcread}, ref)[
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
        audit = {defect: validate_mapping.audit_records(line, {"r": read}, ref) for defect, line in bad.items()}
        self.assertTrue(audit["seq"]["seq_bad"] and audit["seq"]["nm_bad"])  # mismatching base
        self.assertTrue(audit["seq*"]["seq_bad"])
        self.assertTrue(audit["cigar"]["cigar_bad"])
        self.assertTrue(audit["nm_missing"]["nm_bad"])
        self.assertTrue(audit["nm_wrong"]["nm_bad"])
        self.assertTrue(audit["span"]["span_bad"])

    def test_sam_to_paf_reports_unconvertible_records(self):
        sam = (
            "r1\t0\tc\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
            "r2\t0\tc\t1\t60\t8M\t*\t0\t0\tACGTACGTAC\t*\tNM:i:0\n"
        )
        dropped = []
        lines = validate_mapping.sam_to_paf(sam, {"c": 100}, dropped)
        self.assertEqual((len(lines), dropped), (1, ["r2"]))

    def test_parse_sam_skips_secondary_and_supplementary(self):
        sam = (
            "r\t2048\tc\t1\t60\t5M\t*\t0\t0\tACGTA\t*\n"
            "r\t256\tc\t9\t60\t5M\t*\t0\t0\tACGTA\t*\n"
            "r\t16\tc\t20\t60\t5M\t*\t0\t0\tACGTA\t*\n"
        )
        recs = validate_mapping.parse_sam(sam)
        self.assertEqual((recs["r"]["pos"], recs["r"]["is_reverse"]), (20, True))


class TruthScoring(Tmp):
    """Truth scoring of mapped reads in the align category and in scaling_bench."""

    def test_align_truth_metrics(self):
        sys.path.insert(0, os.path.join(HERE, "categories"))
        from categories.align import _truth_metrics

        reads = self.tmp_path("reads.fa")
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
        metrics = _truth_metrics(sam, reads)
        self.assertEqual(
            (metrics["correct_reads"], metrics["wrong_reads"], metrics["unmapped_reads"]), (1, 2, 1)
        )
        self.assertEqual(metrics["mapq30_wrong"], 2)
        self.assertEqual(metrics["mapq30_reads"], 3)

    def test_scaling_map_accuracy(self):
        import scaling_bench

        tsv = self.tmp_path("r.tsv")
        write(
            tsv,
            "read_name\tref\tpos\tstrand\na\tc\t100\t+\nb\tc\t500\t-\nd\tc\t900\t+\ne\tc\t50\t+\n",
        )
        sam = self.tmp_path("o.sam")
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

        write(self.tmp_path("r.fa"), ">a\nAAAC\n>b\nGAAA\n")
        dataset = {
            "k": 3,
            "dbg_input": self.tmp_path("r.fa"),
            "queries": self.tmp_path("r.fa"),
            "dbg_min_count": 2,
        }

        def tool(name, kmers):
            """A graph tool's result dict whose dumped k-mer set is `kmers`."""
            kmers_path = self.tmp_path(name + ".kmers")
            write(kmers_path, "\n".join(sorted(kmers)))
            return {
                "graph_kmers_file": kmers_path,
                "graph_kmers": len(kmers),
                "graph_kmer_occurrences": len(kmers),
                "min_count": 2,
            }

        passing = {
            "optik": tool("optik", {"AAA"}),
            "bifrost": tool("bifrost", {"AAA"}),
            "cuttlefish": tool("cuttlefish", set()),
        }
        self.assertTrue(tool_bench.graph_parity(passing, dataset)["pass"])
        bad = dict(passing, optik=tool("optik2", set()))
        self.assertFalse(tool_bench.graph_parity(bad, dataset)["pass"])
        dup = dict(passing, optik=dict(tool("optik3", {"AAA"}), graph_kmer_occurrences=2))
        self.assertFalse(tool_bench.graph_parity(dup, dataset)["pass"])


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
        self.assertEqual(validate_query.parse_counts(self.TEXT)["num_positive_kmers"], 100)
        self.assertEqual(validate_official.parse_optik_report(self.TEXT)["num_kmers"], 120)


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
        # The dBG oracles' revcomp_packed() is a byte-table shortcut, which must equal the plain per-base definition at every width the validators use.
        import random

        import validate_unitigs

        rnd = random.Random(7)
        for w in range(1, 33):
            vals = [0, (1 << (2 * w)) - 1] + [
                rnd.getrandbits(2 * w) for _ in range(2000)
            ]
            for value in vals:
                self.assertEqual(validate_unitigs.revcomp_packed(value, w), validate_unitigs._rc_slow(value, w), (w, value))

    def test_base4_pack_matches_base_by_base(self):
        import random

        import validate_unitigs

        rnd = random.Random(11)
        for n in range(40):
            for _ in range(200):
                seq = "".join(rnd.choice("ACGT") for _ in range(n))
                self.assertEqual(validate_unitigs.pack(seq), validate_unitigs._pack_slow(seq), seq)


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

        self.assertTrue(micro_bench.evaluate([dict(row) for row in self.GOOD])["usable"])

    def test_counter_error_is_unusable(self):
        import micro_bench

        rows = [dict(row) for row in self.GOOD] + [
            {"kernel": "x", "error": "cycles: Permission denied"}
        ]
        result = micro_bench.evaluate(rows)
        self.assertFalse(result["usable"])
        self.assertIn("counters unavailable", result["reason"])

    def test_off_calibration_is_unusable(self):
        import micro_bench

        rows = [dict(row) for row in self.GOOD]
        rows[1]["cycles"] = 2.5  # as with multiplexed or scaled counts
        result = micro_bench.evaluate(rows)
        self.assertFalse(result["usable"])
        self.assertIn("add_dependent", result["reason"])

    def test_cached_dram_chase_is_unusable(self):
        import micro_bench

        rows = [dict(row) for row in self.GOOD]
        rows[4]["cycles"] = 62.0  # the L3-resident chase, restarting at the head
        result = micro_bench.evaluate(rows)
        self.assertFalse(result["usable"])
        self.assertIn("chase_dram", result["reason"])


def sampled_kmers(path, kmer_len, mod=97):
    """Canonical k-mers of a FASTA(.gz) whose hash % mod == 0, per record with wrapped lines joined.

    This is a fixed 1/mod sample, the same k-mers in every file."""
    import gzip as _gz
    import zlib

    comp = str.maketrans("ACGT", "TGCA")
    out, seqs, cur = set(), [], []
    opener = _gz.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    for seq in seqs:
        seq_rc = seq.translate(comp)[::-1]
        seq_len = len(seq)
        for i in range(seq_len - kmer_len + 1):
            fwd_kmer, rev_kmer = seq[i : i + kmer_len], seq_rc[seq_len - kmer_len - i : seq_len - i]
            canon_kmer = min(rev_kmer, fwd_kmer)
            if "N" not in canon_kmer and zlib.crc32(canon_kmer.encode()) % mod == 0:
                out.add(canon_kmer)
    return out


class DatasetConsistency(unittest.TestCase):
    """Every dataset's query, graph and mapping genome must come from the organism its unitigs index.

    The SE datasets once pointed at E. coli (1.8% of its k-mers are in the S. enterica unitigs), so an 'SE' run repeated the E. coli rows under another name."""

    def test_query_genome_matches_unitigs(self):
        from config import DATASETS

        for name, dataset in DATASETS.items():
            if not (Path(dataset["unitigs"]).exists() and Path(dataset["queries"]).exists()):
                continue
            kmer_len = dataset["k"]
            unitig_kmers = sampled_kmers(dataset["unitigs"], kmer_len)
            query_kmers = sampled_kmers(dataset["queries"], kmer_len)
            frac = len(unitig_kmers & query_kmers) / max(1, len(query_kmers))
            self.assertGreater(
                frac,
                0.9,
                f"{name}: only {frac:.1%} of the query genome's k-mers are in the unitigs",
            )
            self.assertEqual(dataset.get("dbg_input"), dataset["queries"], name)


class Measure(Tmp):
    """bench/measure.py drops warmups, keeps exactly --runs runs, returns the stdout and stderr of the final run, and reports failures untimed."""

    def test_warmups_dropped_runs_kept(self):
        import measure

        cnt = self.tmp_path("n")
        # Each run appends one line, so the last run's stdout shows the run count.
        cmd = ["sh", "-c", f"echo x >> {cnt}; wc -l < {cnt}; echo err >&2"]
        result = measure.measure(cmd, 2, 3, scratch=self.tmp_path("m"))
        self.assertTrue(result.succeeded, result.err)
        self.assertEqual((result.warmup, result.runs, len(result.per_run)), (2, 3, 3))
        with open(result.out) as handle:
            self.assertEqual(
                handle.read().strip(), "5"
            )  # 2 warmups + 3 runs, the last one's output
        self.assertEqual(result.err.strip(), "err")
        self.assertIsNotNone(result.wall["sd"])
        self.assertLessEqual(result.wall["min"], result.wall["mean"])
        self.assertLessEqual(result.wall["mean"], result.wall["max"])

    def test_single_run_has_no_sd(self):
        import measure

        result = measure.measure(["true"], 0, 1, scratch=self.tmp_path("m"))
        self.assertTrue(result.succeeded)
        self.assertIsNone(result.wall["sd"])

    def test_failure_not_timed(self):
        import measure

        result = measure.measure(
            ["sh", "-c", "echo boom >&2; exit 3"], 1, 2, scratch=self.tmp_path("m")
        )
        self.assertFalse(result.succeeded)
        self.assertEqual(result.code, 3)
        self.assertIn("boom", result.err)
        self.assertEqual(result.wall, {})


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
        # A child pinned to two CPUs reports two, as optik's available_parallelism() would.
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
        idle = [(when, {**on_cpu, 3: 0}) for when, on_cpu in two]
        self.assertAlmostEqual(thread_check.peak_busy(idle, 0.05), 2.0)
        # A thread exiting mid-run loses its last share, never inflating the result.
        gone = [(when, on_cpu if when < 0.05 else {1: on_cpu[1]}) for when, on_cpu in two]
        self.assertLessEqual(thread_check.peak_busy(gone, 0.05), 2.0)

    def test_default_cap_fits_the_machine(self):
        mem = hw.mem_bytes()
        self.assertGreater(mem, 0)
        self.assertLessEqual(hw.default_cap_bytes(), mem)
        self.assertGreaterEqual(hw.default_cap_bytes(), mem // 2)


class RecordFormats(unittest.TestCase):
    """The record helpers that replaced bare column indexes."""

    def test_sam_record_reads_fields_by_name_and_still_indexes(self):
        from sam_spec import SamCol, SamFlag, SamRecord

        rec = SamRecord.parse("r1\t16\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tNM:i:0")
        self.assertEqual((rec.qname, rec.flag, rec.pos, rec.mapq, rec.cigar), ("r1", 16, 100, 60, "5M"))
        self.assertEqual(rec[SamCol.RNAME], "chr1")
        self.assertEqual(rec.tags, ("NM:i:0",))
        self.assertTrue(rec.is_reverse and rec.is_primary and not rec.is_unmapped)
        self.assertEqual(str(SamFlag.REVERSE | SamFlag.PAIRED), "17")
        self.assertFalse(SamRecord.parse("r\t2048\tchr1\t1\t0\t*\t*\t0\t0\t*\t*").is_primary)

    def test_fastq_helpers_count_and_extract(self):
        import io

        from fastq import count_fastq_records, fastq_sequences

        text = "@a\nACGT\n+\nIIII\n@b\nGG\n+\nII\n"
        self.assertEqual(count_fastq_records(io.StringIO(text)), 2)
        self.assertEqual(fastq_sequences(text.splitlines(keepends=True)), ["ACGT", "GG"])

    def test_overconfident_rule(self):
        from mapq_calibration import overconfident

        self.assertFalse(overconfident(7, 2.0))
        self.assertTrue(overconfident(8, 2.0))

    def test_rusage_line_with_and_without_io_fields(self):
        import measure

        full = measure.RusageLine.parse("0 1.5 2048 4096 512")
        self.assertEqual((full.exit_code, full.wall_s, full.peak_rss_kib, full.bytes_written, full.bytes_read), (0, 1.5, 2048, 4096, 512))
        short = measure.RusageLine.parse("137 0.25 100")
        self.assertIsNone(short.bytes_written)

    def test_wfa_bench_summary_and_fai_rows(self):
        import compare_tools_small_sets as small

        summary = small.parse_wfa_bench_summary("align_seconds 1.5 invalid 0 pairs 1000")
        self.assertEqual(summary, {"align_seconds": 1.5, "invalid": 0, "pairs": 1000})
        with tempfile.TemporaryDirectory() as tmp:
            fasta = os.path.join(tmp, "t.fa")
            Path(fasta).write_text(">a desc\nACGTACGT\nACG\n>b\nTTTT\n")
            small.write_fai(fasta)
            # The same rows samtools faidx writes for this file.
            self.assertEqual(Path(fasta + ".fai").read_text(), "a\t11\t8\t8\t9\nb\t4\t24\t4\t5\n")


class _RangeServer(http.server.BaseHTTPRequestHandler):
    """Serves DATA with Range support, plus the one fault named by FAULT, so rangefetch can be tested offline.

    The faults copy what proxies, VPNs and queueing servers did to real fetches.
    full200 answers every fourth range with the whole file, shifted every fourth with a 206 for a later offset, and drop cuts every response short.
    slow_first holds each response back for SLOW_FIRST_DELAY_S, and bigger reports a total one byte larger than DATA.
    """

    DATA = b""
    FAULT = None
    requests = 0
    lock = threading.Lock()
    EVERY = 4  # full200 and shifted hit every fourth request
    SHIFT = 7  # bytes by which a shifted response starts late
    DROP_SHARE = 0.4  # drop sends this share of the requested range, at least one byte
    SLOW_FIRST_DELAY_S = 1.5

    def log_message(self, *args):
        pass  # keep the test output clean

    def do_HEAD(self):
        self.send_response(200)
        self.send_header("Content-Length", str(len(self.DATA)))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()

    def do_GET(self):
        with self.lock:
            type(self).requests += 1
            nth = self.requests
        lo, hi = 0, len(self.DATA) - 1
        match = re.match(r"bytes=(\d+)-(\d+)", self.headers.get("Range", ""))
        if match:
            lo, hi = int(match.group(1)), min(int(match.group(2)), len(self.DATA) - 1)
        if self.FAULT == "slow_first":
            time.sleep(self.SLOW_FIRST_DELAY_S)
        if not match or (self.FAULT == "full200" and nth % self.EVERY == 0):
            self.send_response(200)
            self.send_header("Content-Length", str(len(self.DATA)))
            self.end_headers()
            self.wfile.write(self.DATA)
            return
        if self.FAULT == "shifted" and nth % self.EVERY == 0:
            lo = min(lo + self.SHIFT, hi)
        total = len(self.DATA) + (1 if self.FAULT == "bigger" else 0)
        body = self.DATA[lo : hi + 1]
        self.send_response(206)
        self.send_header("Content-Range", f"bytes {lo}-{hi}/{total}")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if self.FAULT == "drop":
            body = body[: max(1, int(len(body) * self.DROP_SHARE))]
        self.wfile.write(body)


class _QuietServer(http.server.ThreadingHTTPServer):
    """Threading server that ignores clients hanging up mid-response, which killed curls do on purpose."""

    daemon_threads = True

    def handle_error(self, request, client_address):
        pass


class RangeFetch(Tmp):
    """rangefetch.Scheduler against a local server that misbehaves the way real servers and proxies did.

    Before the fix, full200 and shifted produced a file of the right size with the wrong md5, which is the cf2-gut failure.
    drop exhausted RANGE_MAX_ATTEMPTS even though every attempt made progress, which is the cf2-rnaseq-1kg failure.
    """

    SIZE = 1 << 20
    FETCH_TIMEOUT_S = 90

    def setUp(self):
        super().setUp()
        import rangefetch

        self.rangefetch = rangefetch
        _RangeServer.DATA = random.Random(7).randbytes(self.SIZE)
        _RangeServer.FAULT, _RangeServer.requests = None, 0
        self.server = _QuietServer(("127.0.0.1", 0), _RangeServer)
        threading.Thread(target=self.server.serve_forever, daemon=True).start()
        self.url = f"http://127.0.0.1:{self.server.server_port}/data.bin"
        # Small ranges and short timers keep every case to a few seconds.
        for name, value in {
            "CHUNK": 64 << 10,
            "CHUNK_MAX": 64 << 10,
            "MIN_SPLIT": 8 << 10,
            "TAKEOVER_S": 0.5,
            "RETRY_SLEEP_PER_ATTEMPT_S": 0,
            "SYNC_INTERVAL_S": 0.2,
        }.items():
            patcher = unittest.mock.patch.object(rangefetch, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        super().tearDown()

    def fetch(self, fault):
        """Run one fetch under `fault`, returning (scheduler, error text or None, md5 matches)."""
        _RangeServer.FAULT = fault
        sched = self.rangefetch.Scheduler(max_conns=4, log=lambda *args, **kwargs: None)
        part = self.tmp_path("data.bin.part")
        outcome = {}

        def run():
            try:
                sched.fetch(self.url, part, self.SIZE)
            except SystemExit as err:
                outcome["err"] = str(err)

        worker = threading.Thread(target=run, daemon=True)
        worker.start()
        worker.join(self.FETCH_TIMEOUT_S)
        self.assertFalse(worker.is_alive(), f"fetch under {fault} did not finish")
        same = os.path.exists(part) and Path(part).read_bytes() == _RangeServer.DATA
        return sched, outcome.get("err"), same

    def test_honest_server(self):
        _, err, same = self.fetch(None)
        self.assertIsNone(err)
        self.assertTrue(same)

    def test_full_body_answers_are_discarded(self):
        _, err, same = self.fetch("full200")
        self.assertIsNone(err)
        self.assertTrue(same, "a 200 with the whole file was written into a range")

    def test_shifted_ranges_are_discarded(self):
        _, err, same = self.fetch("shifted")
        self.assertIsNone(err)
        self.assertTrue(same, "a 206 for another offset was written into a range")

    def test_progressing_attempts_do_not_fail_the_file(self):
        _, err, same = self.fetch("drop")
        self.assertIsNone(err, "attempts that each made progress still failed the file")
        self.assertTrue(same)

    def test_waiting_for_a_first_byte_is_not_taken_over(self):
        sched, err, same = self.fetch("slow_first")
        self.assertIsNone(err)
        self.assertTrue(same)
        self.assertEqual(sched.takeovers, 0, "ranges still waiting for a first byte were taken over")

    def test_changed_size_fails_at_once(self):
        _, err, _ = self.fetch("bigger")
        self.assertIsNotNone(err)
        self.assertIn("serves", err)

    def test_part_lock_refuses_a_second_writer(self):
        part = self.tmp_path("x.part")
        with self.rangefetch.part_lock(part), self.assertRaises(SystemExit), self.rangefetch.part_lock(part):
            pass

    def test_preallocated_part_always_has_done(self):
        # A crash before the first sync must not leave a full-size .part that _adopt would take as complete.
        part = self.tmp_path("data.bin.part")
        sched = self.rangefetch.Scheduler(max_conns=1, log=lambda *args, **kwargs: None)
        _RangeServer.FAULT = "slow_first"
        worker = threading.Thread(target=sched.fetch, args=(self.url, part, self.SIZE), daemon=True)
        worker.start()
        deadline = time.time() + self.FETCH_TIMEOUT_S
        while not os.path.exists(part) and time.time() < deadline:
            time.sleep(0.01)
        self.assertTrue(os.path.exists(part + ".done"), ".part exists without .part.done")
        # The fetch must end before tearDown removes its directory.
        worker.join(self.FETCH_TIMEOUT_S)

    def test_plain_path_drops_stale_ranged_state(self):
        # A .part preallocated by the scheduler, with its .done, must not be resumed or trusted by the single-connection path.
        from published.fetch import download

        dest = self.tmp_path("small.bin")
        Path(dest + ".part").write_bytes(bytes(self.SIZE))
        Path(dest + ".part.done").write_text("[[0, 1048576]]")
        digest = hashlib.md5(_RangeServer.DATA).hexdigest()
        with contextlib.redirect_stdout(io.StringIO()):
            download(self.url, dest, "md5", digest, self.SIZE)
        self.assertEqual(Path(dest).read_bytes(), _RangeServer.DATA)
        self.assertFalse(os.path.exists(dest + ".part.done"))
        self.assertFalse(os.path.exists(dest + ".part.lock"))

    def test_slow_answers_keep_a_minimum_of_connections(self):
        # On a high-latency path every first byte is slow, which once cut ENA to one connection at 0.1 MB/s.
        sched = self.rangefetch.Scheduler(max_conns=16, log=lambda *args, **kwargs: None)
        host = "slow.example"
        for _ in range(50):
            sched._host(host).cut = 0.0  # let every call count, as in a long fetch
            sched._answered(host, self.rangefetch.SLOW_FIRST_S * 3, busy=1)
        self.assertEqual(sched._host(host).cap, self.rangefetch.SLOW_ANSWER_MIN_CAP)
        sched._rejected(host, retry_after=0)
        self.assertLess(sched._host(host).cap, self.rangefetch.SLOW_ANSWER_MIN_CAP, "a 429/503 may still cut below the floor")

    def test_unanswered_probe_is_unknown(self):
        # A probe that gets no answer is not a "no", so it must not send the file down the path that discards a partial .part.
        closed = socket.socket()
        closed.bind(("127.0.0.1", 0))
        port = closed.getsockname()[1]
        closed.close()
        self.assertIsNone(self.rangefetch.supports_ranges(f"http://127.0.0.1:{port}/x"))
        self.assertTrue(self.rangefetch.supports_ranges(self.url))


if __name__ == "__main__":
    unittest.main(verbosity=1)
