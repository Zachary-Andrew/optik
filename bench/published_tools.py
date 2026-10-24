"""Compared tools of compare_tools_published_sets.py.

This module locates each tool, sets up the environment the bench builds its helper tools in, and holds the build and query command lines of every k-mer index that the query section compares.
"""

import os
import shutil
import subprocess

import compare_tools_small_sets as pb
from published_catalog import TOOLS, T
from small_common import BIN


def tool_path(name):
    """Path of a helper tool, taken from PATH or else from what `tools` built under TOOLS."""
    return shutil.which(name) or os.path.join(TOOLS, name)


def have_tool(name):
    """True when helper tool `name` is on PATH or was built under TOOLS."""
    return shutil.which(name) or os.path.exists(os.path.join(TOOLS, name))


def sshash_budget(mach):
    """SSHash's -g in GiB, the authors' 16 or less on a smaller machine."""
    return max(1, min(16, int(mach["mem_gib"] * 0.8)))


def ggcat_bin():
    """Path of the GGCAT binary that `tools` installs."""
    return os.path.join(TOOLS, "ggcat/bin/ggcat")


def ggcat_version():
    """Version of the installed GGCAT, reported with each SPSS."""
    r = subprocess.run(
        [ggcat_bin(), "--version"], capture_output=True, text=True, check=False
    )
    return (
        r.stdout.strip().split()[-1] if r.returncode == 0 and r.stdout.strip() else "?"
    )


def dev_env(pkgs):
    """Build environment that exposes these libraries' headers, libraries and pkg-config files.

    The system's copies are used if the compiler finds them.
    Otherwise the nixpkgs ones are used, which needs no root, kept alive by GC roots under temp/gcroots.
    Winnowmap's bundled meryl builds htslib, which links curl and openssl and asks pkg-config for them as well as xz, zlib and bzip2."""
    env = dict(os.environ)
    probe = {"curl": "curl/curl.h", "openssl": "openssl/hmac.h"}
    missing = [
        p
        for p in pkgs
        if subprocess.run(
            [
                "bash",
                "-c",
                f"echo '#include <{probe[p]}>' | g++ -E -x c++ - >/dev/null 2>&1",
            ],
            check=False,
        ).returncode
        != 0
    ]
    if not missing:
        return env
    if not shutil.which("nix"):
        raise SystemExit(
            f"headers of {', '.join(missing)} not found: install them (Debian/Ubuntu: libcurl4-openssl-dev "
            "libssl-dev zlib1g-dev liblzma-dev libbz2-dev; Fedora: libcurl-devel openssl-devel "
            "zlib-devel xz-devel bzip2-devel), or have nix available to fetch them"
        )
    inc, lib, pc = [], [], []
    for p in dict.fromkeys(missing + ["xz", "zlib", "bzip2"]):
        dev = os.path.join(T, "gcroots", f"{p}-dev")
        subprocess.run(
            ["nix", "build", f"nixpkgs#{p}.dev", "--out-link", dev], check=True
        )
        dev += "-dev"  # nix names a dev output's link LINK-dev
        if p in missing:
            out = os.path.join(T, "gcroots", f"{p}-out")
            subprocess.run(
                ["nix", "build", f"nixpkgs#{p}.out", "--out-link", out], check=True
            )
            inc.append(os.path.join(dev, "include"))
            lib.append(os.path.join(out, "lib"))
        pc += [
            x
            for x in (
                os.path.join(dev, "lib/pkgconfig"),
                os.path.join(dev, "share/pkgconfig"),
            )
            if os.path.isdir(x)
        ]

    def join(xs, var):
        """xs followed by the current value of search path `var`, colon-separated."""
        return ":".join(xs + ([env[var]] if env.get(var) else []))

    env["CPATH"] = join(inc, "CPATH")
    env["LIBRARY_PATH"] = join(lib, "LIBRARY_PATH")
    env["PKG_CONFIG_PATH"] = join(pc, "PKG_CONFIG_PATH")
    env["LDFLAGS"] = (
        " ".join(f"-Wl,-rpath,{x}" for x in lib) + " " + env.get("LDFLAGS", "")
    )
    return env


# Tools that write per-k-mer answers to their own file (-o, --ids).
# Their stdout must never be captured into that file.
OWN_OUTPUT = {"sbwt", "optik-spss --ids"}

# Tools with no thread option for the step are measured once, at t = 1.
# These are FMSI's build, i.e. KmerCamel compute and fmsi index, and the query of every compared tool, SSHash, SBWT and FMSI.
NO_BUILD_THREADS = {"fmsi"}


def query_tools(mode, d, unitigs_fa, k, m, threads, mach, cap_gib, qthreads=1):
    """Build and query commands of every compared k-mer index, keyed by tool name.

    Each entry gives the build command (None for a query view of another tool's index, named by `reuse`), the index files, a query function that maps a query file to (command, answer file), the kind of answers the tool writes and a parameter string for the report.
    `mode` is "parity" or "optimal", and indexes go to directory `d`."""
    par = mode == "parity"
    bt = threads
    qt = qthreads  # every paper ran single-threaded queries (1), and a sweep scales optik's
    ms = str(m)
    budget = sshash_budget(mach)  # the authors' 16 GiB, or less on a machine with less
    sbwt_gb = max(2, budget)  # SBWT's -m must be at least 2, and the authors used 16
    return {
        "optik-spss": {
            "build": [
                BIN["optik"],
                "build-spss",
                "-i",
                unitigs_fa,
                "-k",
                str(k),
                "-m",
                ms if par else "0",
                "-t",
                str(bt),
                "-o",
                f"{d}/optik.spss",
            ],
            "files": [f"{d}/optik.spss"],
            "query": lambda q: (
                [
                    BIN["optik"],
                    "query",
                    "-i",
                    f"{d}/optik.spss",
                    "-q",
                    q,
                    "-t",
                    str(qt),
                ],
                None,
            ),
            "answers": "counts",
            "params": f"k={k} m={ms if par else 'auto'} -t {bt}; query -t {qt}",
        },
        # The same index answering per k-mer, with a dense id per window and -1 when absent.
        # This is the like-for-like comparison with the per-k-mer output of SBWT and FMSI.
        # Timed runs write to /dev/null.
        "optik-spss --ids": {
            "build": None,
            "reuse": "optik-spss",
            "query": lambda q: (
                [
                    BIN["optik"],
                    "query",
                    "-i",
                    f"{d}/optik.spss",
                    "-q",
                    q,
                    "-t",
                    str(qt),
                    "--ids",
                    f"{d}/ids.out",
                ],
                f"{d}/ids.out",
            ),
            "timed": lambda q: [
                BIN["optik"],
                "query",
                "-i",
                f"{d}/optik.spss",
                "-q",
                q,
                "-t",
                str(qt),
                "--ids",
                "/dev/null",
            ],
            "answers": "per-k-mer",
            "params": f"the optik-spss index; query --ids (per-k-mer ids) -t {qt}",
        },
        "optik-fm": {
            "build": [
                BIN["optik"],
                "build-fm",
                "-i",
                unitigs_fa,
                "-k",
                str(k),
                "-c",
                "1",
                "-t",
                str(bt),
                "-o",
                f"{d}/optik.fm",
            ],
            "files": [f"{d}/optik.fm"],
            "query": lambda q: (
                [
                    BIN["optik"],
                    "query",
                    "-i",
                    f"{d}/optik.fm",
                    "-q",
                    q,
                    "-t",
                    str(qt),
                ],
                None,
            ),
            "answers": "counts",
            "params": f"k={k} -t {bt}; query -t {qt}",
        },
        "sshash": {
            # The authors' own command, with -m per collection and a 16 GiB RAM budget, -g 16
            "build": [
                BIN["sshash"],
                "build",
                "-i",
                unitigs_fa,
                "-k",
                str(k),
                "-m",
                ms,
                "-g",
                str(budget),
                "-t",
                str(bt),
                "-d",
                d,
                "-o",
                f"{d}/idx.sshash",
            ],
            "files": [f"{d}/idx.sshash"],
            "query": lambda q: (
                [BIN["sshash"], "query", "-i", f"{d}/idx.sshash", "-q", q],
                None,
            ),
            "answers": "counts",
            "params": f"k={k} m={ms} -g {budget} -t {bt}; query single-threaded",
        },
        "sbwt": {
            # The authors' command in kmer_sets_benchmark's build_sbwt_indexes.sh uses -m 16 and the default plain-matrix variant with --add-reverse-complements, the index their streaming queries use.
            # Their index without reverse complements is also built, for its size (extra).
            "build": [
                BIN["sbwt"],
                "build",
                "-i",
                unitigs_fa,
                "-o",
                f"{d}/idx.sbwt",
                "-k",
                str(k),
                "--variant",
                "plain-matrix",
                "--add-reverse-complements",
                "-m",
                str(sbwt_gb),
                "-t",
                str(bt),
                "--temp-dir",
                d,
            ],
            "files": [f"{d}/idx.sbwt"],
            "extra": (
                [
                    BIN["sbwt"],
                    "build",
                    "-i",
                    unitigs_fa,
                    "-o",
                    f"{d}/idx_norc.sbwt",
                    "-k",
                    str(k),
                    "--variant",
                    "plain-matrix",
                    "-m",
                    str(sbwt_gb),
                    "-t",
                    str(bt),
                    "--temp-dir",
                    d,
                ],
                f"{d}/idx_norc.sbwt",
            ),
            "query": lambda q: (
                [
                    BIN["sbwt"],
                    "search",
                    "-i",
                    f"{d}/idx.sbwt",
                    "-q",
                    q,
                    "-o",
                    f"{d}/sbwt.out",
                ],
                f"{d}/sbwt.out",
            ),
            "timed": lambda q: [
                BIN["sbwt"],
                "search",
                "-i",
                f"{d}/idx.sbwt",
                "-q",
                q,
                "-o",
                "/dev/null",
            ],
            "answers": "per-k-mer",
            "params": f"k={k} plain-matrix +rc -m {sbwt_gb} -t {bt}; search single-threaded",
        },
        "fmsi": {
            # The authors run kmercamel compute -S on SPSS input, then fmsi index.
            "build": [
                "bash",
                "-c",
                (
                    f"{BIN['kmercamel']} compute -S -k {k} -M {d}/ms.fa -o /dev/null {unitigs_fa} && "
                    f"{BIN['fmsi']} index -k {k} {d}/ms.fa"
                ),
            ],
            "files_glob": f"{d}/ms.fa.fmsi.",
            # Both modes query with -S -O, the paper's settings.
            # kmer_sets_benchmark calls fmsi_lookup_*_query<minimal_hash=false, maximized_ones=true> over the kLCP array, which the CLI exposes as -O and -S.
            # These are FMSI's own optimizations rather than a parameter to match.
            "query": lambda q: (
                [
                    BIN["fmsi"],
                    "query",
                    "-S",
                    "-O",
                    "-q",
                    q,
                    "-k",
                    str(k),
                    f"{d}/ms.fa",
                ],
                f"{d}/fmsi.out",
            ),
            "answers": "per-k-mer",  # in timed runs hyperfine sends stdout to /dev/null
            "params": f"k={k} kmercamel compute -S + fmsi index; query -S -O (the paper's)",
        },
    }


def idx_files(t, d):
    """Index files of one query tool's build.

    These are the tool's explicit file list, or else every file in d whose name starts with its glob prefix."""
    if "files" in t:
        return [f for f in t["files"] if os.path.exists(f)]
    return [
        os.path.join(d, f)
        for f in os.listdir(d)
        if os.path.join(d, f).startswith(t["files_glob"])
    ]


def count_of(tool, outf, capt, err):
    """Positive k-mer count a query run reported, parsed from its answer file or captured output."""
    if tool in ("sbwt", "optik-spss --ids"):  # one id per window, -1 when absent
        return pb.sbwt_counts(outf)
    if tool == "fmsi":
        return pb.fmsi_counts(capt)
    return pb.parse_counts(open(capt).read() + err)
