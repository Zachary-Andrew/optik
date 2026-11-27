"""Runs the de Bruijn graph section, where every graph builder runs and their unitigs are checked against each other."""

import gzip
import os
import subprocess
import sys

from diagnostics import DETAIL_CHARS
from dna import n50, read_fasta, write_fasta
from external_oracles import judge_unitigs, kmc_count, kmc_db, same_unitigs
from sam_spec import GfaSegmentCol
from small_common import BIN, TEMP_DIR

from small_sets.common import HERE, KMC_MEM_GB, SET_K

# The read set runs at cutoff 2, which drops most error k-mers.
READS_CUTOFF = 2


def gfa_or_fasta_seqs(path):
    """Segment sequences of a GFA file, or the records of a FASTA file."""
    opener = gzip.open if path.endswith(".gz") else open
    if ".gfa" in path:
        with opener(path, "rt") as handle:
            return [line.split("\t")[GfaSegmentCol.SEQ] for line in handle if line.startswith("S\t")]
    return [seq for _, seq in read_fasta(path)]


def section_dbg(bench, mode, work, threads, quick):
    """Build compacted de Bruijn graphs with every graph tool and judge each against KMC."""
    print(f"\n=== dbg ({mode})", flush=True)
    dbg_dir = os.path.join(work, "dbg")
    os.makedirs(dbg_dir, exist_ok=True)
    genome = os.path.join(TEMP_DIR, "salmonella_LT2.fa")
    reads = os.path.join(dbg_dir, "reads30x.fa")
    if not os.path.exists(reads):
        subprocess.run(
            [
                sys.executable,
                os.path.join(HERE, "gen_reads.py"),
                "-r",
                genome,
                "-c",
                "30",
                "-l",
                "150",
                "-e",
                "0.005",
                "--ins",
                "0.0005",
                "--del",
                "0.0005",
                "--rc",
                "0.5",
                "-s",
                "31",
                "-o",
                reads,
            ],
            check=True,
        )
    cases = [("genome", genome, 1)] + ([] if quick else [("reads_30x", reads, READS_CUTOFF)])
    for name, inp, min_count in cases:
        reads_in = name.startswith("reads")
        case_dir = os.path.join(dbg_dir, f"{name}_{mode}")
        os.makedirs(case_dir, exist_ok=True)
        truth = kmc_db([inp], SET_K, min_count, os.path.join(case_dir, "truth"), threads, True, KMC_MEM_GB)
        print(f"  {name}: KMC solid k-mers (c >= {min_count}) = {kmc_count(truth)}", flush=True)
        cf_mode = "--read" if reads_in else "--ref"
        tools = {
            # All tools write the same artifact, uncompressed FASTA unitigs without GFA links or gzip, so the timings compare graph construction rather than output formats.
            "optik": (
                [
                    BIN["optik"],
                    "dbruijn",
                    "-i",
                    inp,
                    "-k",
                    "31",
                    "-c",
                    str(min_count),
                    "-t",
                    str(threads),
                    "-o",
                    f"{case_dir}/optik.fa",
                ],
                f"{case_dir}/optik.fa",
                f"k=31 -c {min_count} -t {threads}; FASTA unitigs",
            ),
            "cuttlefish": (
                [
                    BIN["cuttlefish"],
                    "build",
                    "-s",
                    inp,
                    "-k",
                    "31",
                    "-t",
                    str(threads),
                    "-o",
                    f"{case_dir}/cf",
                    cf_mode,
                ]
                + (["-c", str(min_count)] if cf_mode == "--read" else []),
                f"{case_dir}/cf.fa",
                f"k=31 {cf_mode}"
                + (f" -c {min_count} (edge cutoff)" if cf_mode == "--read" else "")
                + f" -t {threads}",
            ),
            # For Bifrost, -r keeps every k-mer and -s drops k-mers seen once through its Bloom filter.
            "bifrost": (
                [
                    BIN["bifrost"],
                    "build",
                    "-r" if min_count == 1 else "-s",
                    inp,
                    "-k",
                    "31",
                    "-t",
                    str(threads),
                    "-o",
                    f"{case_dir}/bf",
                    "-f",
                    "-n",
                    "-N",
                ],
                f"{case_dir}/bf.fasta",
                (
                    f"k=31 {'-r (all k-mers)' if min_count == 1 else '-s (drop singletons)'} -t {threads}; "
                    "-f -n -N (FASTA, uncompressed, no index file)"
                ),
            ),
        }
        written = {}
        for tool, (cmd, outp, params) in tools.items():
            # The KMC judging evicts the input from the cache, so it is read in again before every run.
            result = bench.measure(f"{mode}/dbg/{name}/{tool}", cmd, cwd=case_dir, inputs=[inp])
            rc_, err = (0 if result.succeeded else result.code), result.err
            if rc_ != 0 or not os.path.exists(outp):
                alt = outp.replace(".gfa.gz", ".gfa")
                if rc_ == 0 and os.path.exists(alt):
                    outp = alt
                else:
                    bench.add(
                        section="dbg",
                        mode=mode,
                        dataset=name,
                        tool=tool,
                        gate=False,
                        gate_detail=f"exit {rc_}, output {outp} missing: {err[-DETAIL_CHARS:]}",
                        params=params,
                    )
                    continue
            seqs = gfa_or_fasta_seqs(outp)
            unitigs_fa = os.path.join(case_dir, f"{tool}.unitigs.fa")
            write_fasta(unitigs_fa, ((str(i), seq) for i, seq in enumerate(seqs)))
            written[tool] = unitigs_fa
            # KMC judges each tool under that tool's own documented graph model.
            passed, det, num_distinct, n_model = judge_unitigs(
                tool, unitigs_fa, [inp], SET_K, min_count, reads_in, case_dir, threads, KMC_MEM_GB, truth
            )
            lens = [len(seq) for seq in seqs]
            unitig_n50 = n50(lens)
            bench.add(
                section="dbg",
                mode=mode,
                dataset=name,
                tool=tool,
                wall_s=result.wall_row(),
                rss_mb=round(result.rss(), 1),
                unitigs=len(lens),
                total_bp=sum(lens),
                n50=unitig_n50,
                kmers=num_distinct,
                truth_kmers=n_model,
                graph="edge-centric" if tool == "cuttlefish" else "node-centric",
                gate=passed,
                gate_detail=det,
                params=params,
            )
        # Maximal unitigs are unique, so on the genome optik's must equal Cuttlefish's string for string.
        if (
            not reads_in
            and "optik" in written
            and "cuttlefish" in written
            and not same_unitigs(written["optik"], written["cuttlefish"])
        ):
            for result in bench.rows:
                if (
                    result.get("section"),
                    result.get("dataset"),
                    result.get("mode"),
                    result.get("tool"),
                ) == (
                    "dbg",
                    name,
                    mode,
                    "optik",
                ):
                    result["gate"] = False
                    result["gate_detail"] = (
                        result["gate_detail"] + "; " if result["gate_detail"] else ""
                    ) + "unitigs differ from Cuttlefish's (maximal unitigs are unique)"
