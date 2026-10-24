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
from small_common import BIN, T, rc, read_fasta

KMC_TOOLS = os.path.join(T, "KMC/bin/kmc_tools")


def mem_budget_gb():
    """Memory budget in GB for tools that take one (SSHash -g, SBWT -m).

    This is the SSHash 2026 authors' 16 GB, or 80% of this machine's usable memory, hw.mem_bytes() with cgroup limits applied, when that is less."""
    return max(2, min(16, int(hw.mem_bytes() / 2**30 * 0.8)))


def find_tool(name, *local):
    """Path of a helper tool, searched first on PATH for any Linux install, then where bench/build_tools.sh puts it under temp/.
    None when neither has it."""
    return shutil.which(name) or next(
        (x for x in local if x and os.path.exists(x)), None
    )


SAMTOOLS = find_tool("samtools", os.path.join(T, "samtools/samtools"))

BWA_BIN = find_tool(
    "bwa", os.path.join(T, "bwa/bwa"), os.path.join(T, "gcroots/bwa/bin/bwa")
)


def picard_cmd():
    """Picard as a command list, from `picard` on PATH or else build_tools.sh's picard.jar run under java.
    None when neither exists."""
    p = find_tool("picard", os.path.join(T, "gcroots/picard/bin/picard"))
    if p:
        return [p]
    jar = os.path.join(T, "picard/picard.jar")
    return (
        ["java", "-jar", jar] if os.path.exists(jar) and shutil.which("java") else None
    )


def kmc_db(inputs, k, cutoff, out, threads, fasta_multi=True, mem_gb=8, canonical=True):
    """Build a KMC 3 canonical k-mer database of `inputs` keeping counts >= cutoff.

    Counters go up to 2^32 - 1, so multiplicities are exact.
    An existing database is reused only when its stamp (the inputs' paths, sizes and mtimes, k and cutoff) matches, because a database built from a tool's earlier output must never judge its new output."""
    stamp = json.dumps(
        {
            "inputs": [
                (os.path.abspath(x), os.path.getsize(x), os.path.getmtime(x))
                for x in inputs
            ],
            "k": k,
            "cutoff": cutoff,
            "fmt": fasta_multi,
            "canonical": canonical,
        }
    )
    sf = out + ".stamp"
    if not (
        os.path.exists(out + ".kmc_suf")
        and os.path.exists(sf)
        and Path(sf).read_text() == stamp
    ):
        os.makedirs(os.path.dirname(out), exist_ok=True)
        for ext in (".kmc_pre", ".kmc_suf", ".stamp"):
            if os.path.exists(out + ext):
                os.remove(out + ext)
        lst = out + ".lst"
        with open(lst, "w") as f:
            f.write("\n".join(inputs) + "\n")
        fmt = "-fm" if fasta_multi else "-fq"
        # -sm keeps KMC within -m.
        # Without it, KMC's second stage on large read sets went past -m15 to 28 GB (cf2-gut, 2026-10-18).
        # memcap backs this up with a cap of -m plus headroom.
        cmd = (
            [
                BIN["kmc"],
                f"-k{k}",
                f"-ci{cutoff}",
                "-cs4294967295",
                fmt,
                f"-t{threads}",
                f"-m{mem_gb}",
                "-sm",
            ]
            + ([] if canonical else ["-b"])
            + [f"@{lst}", out, os.path.dirname(out)]
        )
        r = memcap.run(cmd, (mem_gb + 4) << 30, capture_output=True)
        if r.returncode:
            raise subprocess.CalledProcessError(
                r.returncode,
                cmd,
                (
                    f"KMC passed the {mem_gb + 4} GiB cap (-m{mem_gb} -sm)"
                    if r.capped
                    else r.stderr.decode(errors="replace")[-300:]
                ),
                r.stderr,
            )
        with open(sf, "w") as f:
            f.write(stamp)
    return out


def kmc_count(db):
    """Number of distinct k-mers in a KMC database, as kmc_tools info reports."""
    out = subprocess.run(
        [KMC_TOOLS, "info", db], capture_output=True, text=True, check=True
    ).stdout
    return int(re.search(r"total k-mers\s*:\s*(\d+)", out).group(1))


def kmc_diff(a, b, out):
    """Size of |A - B| computed with KMC's own set operations (kmc_tools), listing nothing."""
    subprocess.run(
        [KMC_TOOLS, "simple", a, b, "kmers_subtract", out],
        check=True,
        capture_output=True,
    )
    return kmc_count(out)


def kmc_weight(db, d):
    """Sum of a KMC database's counters, i.e. its k-mer occurrences.

    kmc_tools' histogram covers counts up to 65,535, and a dump adds the rarer higher counts."""
    hist, big = os.path.join(d, "w.hist"), os.path.join(d, "w.big")
    subprocess.run(
        [KMC_TOOLS, "transform", db, "histogram", hist, "-ci1", "-cx65535"],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [KMC_TOOLS, "transform", db, "-ci65536", "dump", big],
        check=True,
        capture_output=True,
    )
    tot = 0
    with open(hist) as fh:
        for line in fh:
            c, n = line.split()
            tot += int(c) * int(n)
    with open(big) as fh:
        for line in fh:
            tot += int(line.split()[1])
    return tot


def kmc_hits(set_db, query_fa, k, d, threads, mem_gb=4):
    """Count query k-mer windows, with multiplicity, whose canonical k-mer is in the set.

    KMC intersects the queries' k-mer counts with the set, keeping the query counters, and the counters are summed.
    A window holding a non-ACGT base is never in a set."""
    os.makedirs(d, exist_ok=True)
    q = kmc_db([query_fa], k, 1, os.path.join(d, "q"), threads, True, mem_gb)
    hit = os.path.join(d, "hit")
    if not os.path.exists(hit + ".kmc_suf"):
        subprocess.run(
            [KMC_TOOLS, "simple", q, set_db, "intersect", hit, "-ocleft"],
            check=True,
            capture_output=True,
        )
    return kmc_weight(hit, d)


def is_fasta(paths):
    """True when the inputs are FASTA rather than FASTQ, judged by their names (.fq, .fastq, optionally .gz, are FASTQ).
    Read sets come both ways, simulated ones as FASTA."""
    return not all(re.search(r"\.f(ast)?q(\.gz)?$", p) for p in paths)


def edge_model_db(inputs, k, c, d, threads, fasta_multi, mem_gb):
    """Cuttlefish's read model, counted by KMC.

    The cutoff applies to edges, so the model is the k-mers of the (k+1)-mers that KMC counts c or more times."""
    out = os.path.join(d, "edge_model")
    e = kmc_db(
        inputs, k + 1, c, os.path.join(d, "edges"), threads, fasta_multi, mem_gb
    )  # reused only when its stamp matches
    if not os.path.exists(out + ".kmc_suf") or os.path.getmtime(
        out + ".kmc_suf"
    ) < os.path.getmtime(e + ".kmc_suf"):
        fa = os.path.join(d, "edges.fa")
        subprocess.run(
            f"{KMC_TOOLS} transform {e} dump {fa}.txt && "
            f'awk \'{{print ">"NR"\\n"$1}}\' {fa}.txt > {fa} && rm {fa}.txt',
            shell=True,
            check=True,
        )
        kmc_db([fa], k, 1, out, threads, True, mem_gb)
        os.remove(fa)
    return out


def judge_unitigs(
    tool, fa, inputs, k, c, reads, d, threads, mem_gb, truth_db, shared=None
):
    """Judge a graph tool's unitig FASTA `fa` with KMC under the tool's own documented model.

    The unitigs' k-mer set must equal KMC's k-mers with count >= c.
    For Cuttlefish on reads the model is the k-mers of KMC's (k+1)-mers with count >= c.
    Bifrost on reads uses its fixed cutoff, KMC count >= 2, and any extra k-mer must be an input k-mer its Bloom filter let through.
    No k-mer may lie in two unitigs, i.e. KMC count >= 2 over the unitigs must be empty.
    Returns (ok, detail, distinct k-mers, model k-mers).
    `d` holds this run's databases, and `shared`, default d, the input-only model databases, reusable across runs of the same input, k and cutoff."""
    shared = shared or d
    db = kmc_db([fa], k, 1, os.path.join(d, f"{tool}.kmc"), threads, True, mem_gb)
    dups = kmc_count(
        kmc_db([fa], k, 2, os.path.join(d, f"{tool}.dup"), threads, True, mem_gb)
    )
    ref_db, model_txt = truth_db, f"KMC k-mers with count >= {c}"
    if tool == "cuttlefish" and reads and c > 1:
        ref_db = edge_model_db(inputs, k, c, shared, threads, is_fasta(inputs), mem_gb)
        model_txt = f"KMC (k+1)-mers with count >= {c} (Cuttlefish's edge cutoff)"
    if tool == "bifrost" and reads:
        ref_db = kmc_db(
            inputs,
            k,
            2,
            os.path.join(shared, "truth_c2"),
            threads,
            is_fasta(inputs),
            mem_gb,
        )
        model_txt = "KMC k-mers with count >= 2 (Bifrost's fixed cutoff)"
    missing = kmc_diff(ref_db, db, os.path.join(d, f"{tool}.missing"))
    extra = kmc_diff(db, ref_db, os.path.join(d, f"{tool}.extra"))
    det = []
    if tool == "bifrost" and reads and extra:
        every = kmc_db(
            inputs,
            k,
            1,
            os.path.join(shared, "input_all"),
            threads,
            is_fasta(inputs),
            mem_gb,
        )
        foreign = kmc_diff(
            os.path.join(d, f"{tool}.extra"), every, os.path.join(d, f"{tool}.foreign")
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
    return same and dups == 0, "; ".join(det), kmc_count(db), kmc_count(ref_db)


def same_unitigs(fa_a, fa_b):
    """Check that two edge-centric tools wrote the same unitigs for one genome.

    Maximal unitigs are unique, so the two sets must match.
    They are compared as sorted canonical strings with GNU sort and cmp."""
    keys = []
    for fa in (fa_a, fa_b):
        out = fa + ".canon"
        with open(out + ".tmp", "w") as o:
            for _, s in read_fasta(fa):
                r = rc(s)
                o.write((min(s, r)) + "\n")
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
    det, ok = [], True
    picard = picard_cmd()
    if picard:
        # No -R, since with a reference Picard needs coordinate-sorted input and SEQ on secondary records.
        # calmd's check below covers NM/MD against the reference.
        p = subprocess.run(
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
            x
            for x in (p.stdout + p.stderr).splitlines()
            if x.startswith(("ERROR", "WARNING"))
            and "NM validation cannot be performed" not in x
        ]
        if errs:
            ok = False
            det.append("Picard: " + "; ".join(errs[:4]))
        elif p.returncode != 0:
            # Picard itself failed, e.g. no java, JVM heap or unreadable file, so nothing was validated.
            ok = False
            det.append(
                f"Picard exit {p.returncode} without a report: {(p.stderr or p.stdout)[-200:]}"
            )
    else:
        ok = False
        det.append(
            "Picard not found (PATH, or bench/build_tools.sh picard): SAM spec not checked"
        )
    # calmd recomputes NM and MD from the reference and names each record whose own tag differs.
    # NM is checked for every mapper, since all write it, and MD where the mapper writes it.
    # The BAM goes to /dev/null and the messages to a file, one per differing record.
    errf = sam + ".calmd.err"
    with open(errf, "w") as ef:
        rc = subprocess.run(
            [SAMTOOLS or "samtools", "calmd", "-b", sam, ref],
            stdout=subprocess.DEVNULL,
            stderr=ef,
            check=False,
        ).returncode
    nm_diff, md_diff, md_n = {}, [], set()
    with open(errf, errors="replace") as ef:
        for x in ef:
            m = re.search(r"different (NM|MD) for read '([^']*)'", x)
            if not m:
                continue
            if m.group(1) == "NM":
                nm_diff[m.group(2)] = x.strip()
            elif md_differs_only_at_n(x):
                md_n.add(m.group(2))
            else:
                md_diff.append(x.strip())
    os.remove(errf)
    # BWA's index replaces reference ambiguity codes by random bases.
    # Its MD therefore names a base where calmd, reading the FASTA, writes N, and its NM may count one edit fewer there.
    # That is a convention rather than an error, for exactly the records whose MD differs only at an N.
    # Any other NM/MD difference fails.
    real = md_diff + [v for r_, v in nm_diff.items() if r_ not in md_n]
    if rc != 0 or real:
        ok = False
        det.append(
            f"samtools calmd: {len(real)} records with NM/MD differing from the reference"
            + (f" (exit {rc})" if rc else "")
            + f": {real[:2]}"
        )
    if md_n:
        det.append(
            f"{len(md_n)} MD strings name a base at a reference N (the mapper's index randomized it)"
        )
    return ok, "; ".join(det)


def md_differs_only_at_n(msg):
    """True for a calmd 'different MD' message whose old and new MD agree at every reference position except where the new one has N.
    At such a position the old MD names a base or counts a match against the randomized base."""
    m = re.search(r"different MD for read '[^']*': '([^']*)' -> '([^']*)'", msg)
    if not m:
        return False

    def expand(md):
        """MD string expanded to one entry per reference position."""
        # '=' marks a match, a base a mismatch and '^'+base a deletion.
        out = []
        for x in re.findall(r"\d+|\^[A-Z]+|[A-Z]", md):
            out += (
                ["="] * int(x)
                if x.isdigit()
                else ["^" + c for c in x[1:]]
                if x[0] == "^"
                else [x]
            )
        return out

    a, b = expand(m.group(1)), expand(m.group(2))
    return (
        len(a) == len(b)
        and a != b
        and all(
            x == y or (y == "N" and x != "N" and not x.startswith("^"))
            for x, y in zip(a, b, strict=False)
        )
    )
