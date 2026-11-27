"""External correctness oracles of the tool benches.

Every correctness number in both benches comes from an external tool, not from the bench scripts.
KMC 3 provides k-mer sets, counts and set operations, Picard ValidateSamFile checks the SAM spec, samtools calmd recomputes NM/MD against the reference, and minimap2's paftools.js (mapeval, pafcmp, call) and samtools' wgsim_eval.pl score mappings.
compare_tools_small_sets.py re-exports these helpers, so the published-sets bench reaches them as pb.<name>.
"""

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import hw
import memcap
from diagnostics import DETAIL_CHARS
from dna import read_fasta, revcomp
from small_common import BIN, TEMP_DIR
from units import BYTES_PER_GIB

KMC_TOOLS = os.path.join(TEMP_DIR, "KMC/bin/kmc_tools")


# memory budget of tools that take one: the SSHash 2026 authors' 16 GB, or this share of the usable memory, at least 2 GB
MAX_TOOL_BUDGET_GB = 16
TOOL_MEMORY_SHARE = 0.8
MIN_TOOL_BUDGET_GB = 2
# KMC -m defaults of a database build and of a hit count, and the cap headroom above -m
KMC_DEFAULT_MEM_GB = 8
KMC_HITS_MEM_GB = 4
KMC_CAP_HEADROOM_GB = 4
# Picard error lines quoted in a failed check
PICARD_ERRORS_SHOWN = 4

# A k-mer that occurs twice in one tool's unitigs is a duplicate.
DUPLICATE_MIN_COUNT = 2
# Bifrost on reads keeps k-mers seen at least this often, its only cutoff.
BIFROST_CUTOFF = 2
# Failure details quote this many example records.
SHOWN_EXAMPLES = 2


def mem_budget_gb():
    """Memory budget in GB for tools that take one (SSHash -g, SBWT -m).

    This is the SSHash 2026 authors' 16 GB, or 80% of this machine's usable memory, hw.mem_bytes() with cgroup limits applied, when that is less."""
    return max(MIN_TOOL_BUDGET_GB, min(MAX_TOOL_BUDGET_GB, int(hw.mem_bytes() / BYTES_PER_GIB * TOOL_MEMORY_SHARE)))


def find_tool(name, *local):
    """Path of a helper tool, searched first on PATH for any Linux install, then where bench/build_tools.sh puts it under temp/.
    None when neither has it."""
    return shutil.which(name) or next(
        (candidate for candidate in local if candidate and os.path.exists(candidate)), None
    )


SAMTOOLS = find_tool("samtools", os.path.join(TEMP_DIR, "samtools/samtools"))

BWA_BIN = find_tool(
    "bwa", os.path.join(TEMP_DIR, "bwa/bwa"), os.path.join(TEMP_DIR, "gcroots/bwa/bin/bwa")
)


def picard_cmd():
    """Picard as a command list, from `picard` on PATH or else build_tools.sh's picard.jar run under java.
    None when neither exists."""
    picard_path = find_tool("picard", os.path.join(TEMP_DIR, "gcroots/picard/bin/picard"))
    if picard_path:
        return [picard_path]
    jar = os.path.join(TEMP_DIR, "picard/picard.jar")
    return (
        ["java", "-jar", jar] if os.path.exists(jar) and shutil.which("java") else None
    )


def kmc_db(inputs, kmer_len, cutoff, out, threads, fasta_multi=True, mem_gb=KMC_DEFAULT_MEM_GB, canonical=True):
    """Build a KMC 3 canonical k-mer database of `inputs` keeping counts >= cutoff.

    Counters go up to 2^32 - 1, so multiplicities are exact.
    An existing database is reused only when its stamp (the inputs' paths, sizes and mtimes, k and cutoff) matches, because a database built from a tool's earlier output must never judge its new output."""
    stamp = json.dumps(
        {
            "inputs": [
                (os.path.abspath(input_path), os.path.getsize(input_path), os.path.getmtime(input_path))
                for input_path in inputs
            ],
            "k": kmer_len,
            "cutoff": cutoff,
            "fmt": fasta_multi,
            "canonical": canonical,
        }
    )
    stamp_file = out + ".stamp"
    if not (
        os.path.exists(out + ".kmc_suf")
        and os.path.exists(stamp_file)
        and Path(stamp_file).read_text() == stamp
    ):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        for ext in (".kmc_pre", ".kmc_suf", ".stamp"):
            if os.path.exists(out + ext):
                os.remove(out + ext)
        lst = out + ".lst"
        with open(lst, "w") as handle:
            handle.write("\n".join(inputs) + "\n")
        fmt = "-fm" if fasta_multi else "-fq"
        # -sm keeps KMC within -m.
        # Without it, KMC's second stage on large read sets went past -m15 to 28 GB (cf2-gut, 2026-10-18).
        # memcap backs this up with a cap of -m plus headroom.
        cmd = (
            [
                BIN["kmc"],
                f"-k{kmer_len}",
                f"-ci{cutoff}",
                # Counts saturate at 2^32 - 1, matching optik's own counter.
                "-cs4294967295",
                fmt,
                f"-t{threads}",
                f"-m{mem_gb}",
                "-sm",
            ]
            # -b turns off canonical counting for inputs where strands must stay apart.
            + ([] if canonical else ["-b"])
            + [f"@{lst}", out, os.path.dirname(out)]
        )
        result = memcap.run(cmd, (mem_gb + KMC_CAP_HEADROOM_GB) * BYTES_PER_GIB, capture_output=True)
        if result.returncode:
            raise subprocess.CalledProcessError(
                result.returncode,
                cmd,
                (
                    f"KMC passed the {mem_gb + KMC_CAP_HEADROOM_GB} GiB cap (-m{mem_gb} -sm)"
                    if result.capped
                    else result.stderr.decode(errors="replace")[-DETAIL_CHARS:]
                ),
                result.stderr,
            )
        with open(stamp_file, "w") as handle:
            handle.write(stamp)
    return out


def kmc_count(db_path):
    """Number of distinct k-mers in a KMC database, as kmc_tools info reports."""
    out = subprocess.run(
        [KMC_TOOLS, "info", db_path], capture_output=True, text=True, check=True
    ).stdout
    return int(re.search(r"total k-mers\s*:\s*(\d+)", out).group(1))


def kmc_diff(db_a, db_b, out):
    """Size of |A - B| computed with KMC's own set operations (kmc_tools), listing nothing."""
    subprocess.run(
        [KMC_TOOLS, "simple", db_a, db_b, "kmers_subtract", out],
        check=True,
        capture_output=True,
    )
    return kmc_count(out)


def kmc_weight(db_path, work_dir):
    """Sum of a KMC database's counters, i.e. its k-mer occurrences.

    kmc_tools' histogram covers counts up to 65,535, and a dump adds the rarer higher counts."""
    hist, big = os.path.join(work_dir, "w.hist"), os.path.join(work_dir, "w.big")
    subprocess.run(
        [KMC_TOOLS, "transform", db_path, "histogram", hist, "-ci1", "-cx65535"],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [KMC_TOOLS, "transform", db_path, "-ci65536", "dump", big],
        check=True,
        capture_output=True,
    )
    tot = 0
    with open(hist) as handle:
        for line in handle:
            count, num_kmers = line.split()
            tot += int(count) * int(num_kmers)
    with open(big) as handle:
        for line in handle:
            tot += int(line.split()[1])
    return tot


def kmc_hits(set_db, query_fa, kmer_len, work_dir, threads, mem_gb=KMC_HITS_MEM_GB):
    """Count query k-mer windows, with multiplicity, whose canonical k-mer is in the set.

    KMC intersects the queries' k-mer counts with the set, keeping the query counters, and the counters are summed.
    A window holding a non-ACGT base is never in a set."""
    os.makedirs(work_dir, exist_ok=True)
    query_db = kmc_db([query_fa], kmer_len, 1, os.path.join(work_dir, "q"), threads, True, mem_gb)
    hit = os.path.join(work_dir, "hit")
    if not os.path.exists(hit + ".kmc_suf"):
        subprocess.run(
            [KMC_TOOLS, "simple", query_db, set_db, "intersect", hit, "-ocleft"],
            check=True,
            capture_output=True,
        )
    return kmc_weight(hit, work_dir)


def is_fasta(paths):
    """True when the inputs are FASTA rather than FASTQ, judged by their names (.fq, .fastq, optionally .gz, are FASTQ).
    Read sets come both ways, simulated ones as FASTA."""
    return not all(re.search(r"\.f(ast)?q(\.gz)?$", path) for path in paths)


def edge_model_db(inputs, kmer_len, cutoff, work_dir, threads, fasta_multi, mem_gb):
    """Cuttlefish's read model, counted by KMC.

    The cutoff applies to edges, so the model is the k-mers of the (k+1)-mers that KMC counts c or more times."""
    out = os.path.join(work_dir, "edge_model")
    edges_db = kmc_db(
        inputs, kmer_len + 1, cutoff, os.path.join(work_dir, "edges"), threads, fasta_multi, mem_gb
    )  # reused only when its stamp matches
    if not os.path.exists(out + ".kmc_suf") or os.path.getmtime(
        out + ".kmc_suf"
    ) < os.path.getmtime(edges_db + ".kmc_suf"):
        edges_fa = os.path.join(work_dir, "edges.fa")
        subprocess.run(
            f"{KMC_TOOLS} transform {edges_db} dump {edges_fa}.txt && "
            f'awk \'{{print ">"NR"\\n"$1}}\' {edges_fa}.txt > {edges_fa} && rm {edges_fa}.txt',
            shell=True,
            check=True,
        )
        kmc_db([edges_fa], kmer_len, 1, out, threads, True, mem_gb)
        os.remove(edges_fa)
    return out


def judge_unitigs(
    tool, unitigs_fa, inputs, kmer_len, cutoff, reads, work_dir, threads, mem_gb, truth_db, shared=None
):
    """Judge a graph tool's unitig FASTA `unitigs_fa` with KMC under the tool's own documented model.

    The unitigs' k-mer set must equal KMC's k-mers with count >= cutoff.
    For Cuttlefish on reads the model is the k-mers of KMC's (k+1)-mers with count >= cutoff.
    Bifrost on reads uses its fixed cutoff, KMC count >= 2, and any extra k-mer must be an input k-mer its Bloom filter let through.
    No k-mer may lie in two unitigs, i.e. KMC count >= 2 over the unitigs must be empty.
    Returns (ok, detail, distinct k-mers, model k-mers).
    `work_dir` holds this run's databases, and `shared`, default work_dir, the input-only model databases, reusable across runs of the same input, k and cutoff."""
    shared = shared or work_dir
    unitigs_db = kmc_db([unitigs_fa], kmer_len, 1, os.path.join(work_dir, f"{tool}.kmc"), threads, True, mem_gb)
    # A k-mer that KMC counts twice in the unitigs appears in two places, which a set of maximal unitigs never does.
    dups = kmc_count(
        kmc_db([unitigs_fa], kmer_len, DUPLICATE_MIN_COUNT, os.path.join(work_dir, f"{tool}.dup"), threads, True, mem_gb)
    )
    # Each tool is judged under its own documented model, so Cuttlefish's cutoff applies to (k+1)-mer edges and Bifrost's to its fixed count of two.
    ref_db, model_txt = truth_db, f"KMC k-mers with count >= {cutoff}"
    if tool == "cuttlefish" and reads and cutoff > 1:
        ref_db = edge_model_db(inputs, kmer_len, cutoff, shared, threads, is_fasta(inputs), mem_gb)
        model_txt = f"KMC (k+1)-mers with count >= {cutoff} (Cuttlefish's edge cutoff)"
    if tool == "bifrost" and reads:
        ref_db = kmc_db(
            inputs,
            kmer_len,
            BIFROST_CUTOFF,
            os.path.join(shared, "truth_c2"),
            threads,
            is_fasta(inputs),
            mem_gb,
        )
        model_txt = "KMC k-mers with count >= 2 (Bifrost's fixed cutoff)"
    # The two set differences give the k-mers a tool lost and the ones it invented.
    missing = kmc_diff(ref_db, unitigs_db, os.path.join(work_dir, f"{tool}.missing"))
    extra = kmc_diff(unitigs_db, ref_db, os.path.join(work_dir, f"{tool}.extra"))
    det = []
    # Bifrost's Bloom filter lets some single-occurrence k-mers through by design, so only extras absent from the input count against it.
    # Bifrost's Bloom filter may admit singletons. Its extra k-mers pass if they are at least real input k-mers.
    if tool == "bifrost" and reads and extra:
        every = kmc_db(
            inputs,
            kmer_len,
            1,
            os.path.join(shared, "input_all"),
            threads,
            is_fasta(inputs),
            mem_gb,
        )
        foreign = kmc_diff(
            os.path.join(work_dir, f"{tool}.extra"), every, os.path.join(work_dir, f"{tool}.foreign")
        )
        det.append(
            f"{extra} input k-mers seen once let through by the Bloom filter (Bifrost's model)"
        )
        same = missing == 0 and foreign == 0
        if foreign:
            det.append(f"{foreign} k-mers not in the input")
    else:
        same = missing == 0 and extra == 0
    if not same:
        det.append(f"k-mer set != {model_txt}: {missing} missing, {extra} extra")
    if dups:
        det.append(f"{dups} k-mers in two unitigs")
    return same and dups == 0, "; ".join(det), kmc_count(unitigs_db), kmc_count(ref_db)


def same_unitigs(fa_a, fa_b):
    """Check that two edge-centric tools wrote the same unitigs for one genome.

    Maximal unitigs are unique, so the two sets must match.
    They are compared as sorted canonical strings with GNU sort and cmp."""
    keys = []
    for unitigs_fa in (fa_a, fa_b):
        out = unitigs_fa + ".canon"
        with open(out + ".tmp", "w") as out_handle:
            for _, seq in read_fasta(unitigs_fa):
                seq_rc = revcomp(seq)
                out_handle.write((min(seq, seq_rc)) + "\n")
        subprocess.run(
            ["sort", "-S", "1G", "-o", out, out + ".tmp"],
            check=True,
            env=dict(os.environ, LC_ALL="C"),
        )
        os.remove(out + ".tmp")
        keys.append(out)
    return subprocess.run(["cmp", "-s", keys[0], keys[1]], check=False).returncode == 0


def external_sam_checks(sam, ref):
    """Judge a SAM with the tools that pipelines run on it.

    Picard ValidateSamFile checks the SAM spec (without -R it does not require read groups, nor qualities from FASTA input), and samtools calmd recomputes NM and MD from the reference.
    Returns (ok, detail)."""
    det, passed = [], True
    picard = picard_cmd()
    if picard:
        # No -R, since with a reference Picard needs coordinate-sorted input and SEQ on secondary records.
        # calmd's check below covers NM/MD against the reference.
        proc = subprocess.run(
            picard
            + [
                "ValidateSamFile",
                "-I",
                sam,
                "-MODE",
                "SUMMARY",
                "-IGNORE",
                "MISSING_READ_GROUP",
                "-IGNORE",
                "RECORD_MISSING_READ_GROUP",
                "-IGNORE",
                "QUALITY_NOT_STORED",
                "-MAX_OUTPUT",
                "10",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        errs = [
            line
            for line in (proc.stdout + proc.stderr).splitlines()
            if line.startswith(("ERROR", "WARNING"))
            # ValidateSamFile is run without -R, so its NM warning is about the missing reference rather than the SAM.
            and "NM validation cannot be performed" not in line
        ]
        if errs:
            passed = False
            det.append("Picard: " + "; ".join(errs[:PICARD_ERRORS_SHOWN]))
        elif proc.returncode != 0:
            # Picard itself failed, e.g. no java, JVM heap or unreadable file, so nothing was validated.
            passed = False
            det.append(
                f"Picard exit {proc.returncode} without a report: {(proc.stderr or proc.stdout)[-DETAIL_CHARS:]}"
            )
    else:
        passed = False
        det.append(
            "Picard not found (PATH, or bench/build_tools.sh picard): SAM spec not checked"
        )
    # calmd recomputes NM and MD from the reference and names each record whose own tag differs.
    # NM is checked for every mapper, since all write it, and MD where the mapper writes it.
    # The BAM goes to /dev/null and the messages to a file, one per differing record.
    errf = sam + ".calmd.err"
    with open(errf, "w") as err_handle:
        calmd_exit = subprocess.run(
            [SAMTOOLS or "samtools", "calmd", "-b", sam, ref],
            stdout=subprocess.DEVNULL,
            stderr=err_handle,
            check=False,
        ).returncode
    nm_diff, md_diff, md_n = {}, [], set()
    with open(errf, errors="replace") as err_handle:
        for line in err_handle:
            match = re.search(r"different (?P<tag>NM|MD) for read '(?P<read>[^']*)'", line)
            if not match:
                continue
            if match.group("tag") == "NM":
                nm_diff[match.group("read")] = line.strip()
            elif md_differs_only_at_n(line):
                md_n.add(match.group("read"))
            else:
                md_diff.append(line.strip())
    os.remove(errf)
    # BWA's index replaces reference ambiguity codes by random bases.
    # Its MD therefore names a base where calmd, reading the FASTA, writes N, and its NM may count one edit fewer there.
    # That is a convention rather than an error, for exactly the records whose MD differs only at an N.
    # Any other NM/MD difference fails.
    real = md_diff + [message for read_name, message in nm_diff.items() if read_name not in md_n]
    if calmd_exit != 0 or real:
        passed = False
        det.append(
            f"samtools calmd: {len(real)} records with NM/MD differing from the reference"
            + (f" (exit {calmd_exit})" if calmd_exit else "")
            + f": {real[:SHOWN_EXAMPLES]}"
        )
    if md_n:
        det.append(
            f"{len(md_n)} MD strings name a base at a reference N (the mapper's index randomized it)"
        )
    return passed, "; ".join(det)


def md_differs_only_at_n(msg):
    """True for a calmd 'different MD' message whose old and new MD agree at every reference position except where the new one has N.
    At such a position the old MD names a base or counts a match against the randomized base."""
    match = re.search(r"different MD for read '[^']*': '(?P<old>[^']*)' -> '(?P<new>[^']*)'", msg)
    if not match:
        return False

    def expand(md_str):
        """MD string expanded to one entry per reference position."""
        # '=' marks a match, a base a mismatch and '^'+base a deletion.
        out = []
        for token in re.findall(r"\d+|\^[A-Z]+|[A-Z]", md_str):
            out += (
                ["="] * int(token)
                if token.isdigit()
                else ["^" + base for base in token[1:]]
                if token[0] == "^"
                else [token]
            )
        return out

    old_md, new_md = expand(match.group("old")), expand(match.group("new"))
    return (
        len(old_md) == len(new_md)
        and old_md != new_md
        and all(
            old_pos == new_pos or (new_pos == "N" and old_pos != "N" and not old_pos.startswith("^"))
            for old_pos, new_pos in zip(old_md, new_md, strict=False)
        )
    )
