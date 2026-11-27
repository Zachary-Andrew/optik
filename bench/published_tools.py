"""Compared tools of compare_tools_published_sets.py.

This module locates each tool, sets up the environment the bench builds its helper tools in, and holds the build and query command lines of every k-mer index that the query section compares.
"""

import os
import shutil
import subprocess

import compare_tools_small_sets
from published_catalog import (
    TEMP_DIR,
    TOOLS,
)
from small_common import BIN

# SSHash's -g follows the authors' 16 GiB, or takes this share of memory on a smaller machine.
SSHASH_MAX_BUDGET_GIB = 16
SSHASH_MEMORY_SHARE = 0.8


# SBWT's -m must be at least 2 GB.
SBWT_MIN_MEM_GB = 2


def tool_path(name):
    """Path of a helper tool, taken from PATH or else from what `tools` built under TOOLS."""
    return shutil.which(name) or os.path.join(TOOLS, name)


def have_tool(name):
    """True when helper tool `name` is on PATH or was built under TOOLS."""
    return shutil.which(name) or os.path.exists(os.path.join(TOOLS, name))


def sshash_budget(mach):
    """SSHash's -g in GiB, the authors' 16 or less on a smaller machine."""
    return max(1, min(SSHASH_MAX_BUDGET_GIB, int(mach["mem_gib"] * SSHASH_MEMORY_SHARE)))


def ggcat_bin():
    """Path of the GGCAT binary that `tools` installs."""
    return os.path.join(TOOLS, "ggcat/bin/ggcat")


def ggcat_version():
    """Version of the installed GGCAT, reported with each SPSS."""
    result = subprocess.run(
        [ggcat_bin(), "--version"], capture_output=True, text=True, check=False
    )
    return (
        result.stdout.strip().split()[-1] if result.returncode == 0 and result.stdout.strip() else "?"
    )


def dev_env(pkgs):
    """Build environment that exposes these libraries' headers, libraries and pkg-config files.

    The system's copies are used if the compiler finds them.
    Otherwise the nixpkgs ones are used, which needs no root, kept alive by GC roots under temp/gcroots.
    Winnowmap's bundled meryl builds htslib, which links curl and openssl and asks pkg-config for them as well as xz, zlib and bzip2."""
    env = dict(os.environ)
    probe = {"curl": "curl/curl.h", "openssl": "openssl/hmac.h"}
    missing = [
        pkg
        for pkg in pkgs
        if subprocess.run(
            [
                "bash",
                "-c",
                f"echo '#include <{probe[pkg]}>' | g++ -E -x c++ - >/dev/null 2>&1",
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
    inc, lib, pkgconfig_dirs = [], [], []
    for pkg in dict.fromkeys(missing + ["xz", "zlib", "bzip2"]):
        dev = os.path.join(TEMP_DIR, "gcroots", f"{pkg}-dev")
        subprocess.run(
            ["nix", "build", f"nixpkgs#{pkg}.dev", "--out-link", dev], check=True
        )
        dev += "-dev"  # nix names a dev output's link LINK-dev
        if pkg in missing:
            out = os.path.join(TEMP_DIR, "gcroots", f"{pkg}-out")
            subprocess.run(
                ["nix", "build", f"nixpkgs#{pkg}.out", "--out-link", out], check=True
            )
            inc.append(os.path.join(dev, "include"))
            lib.append(os.path.join(out, "lib"))
        pkgconfig_dirs += [
            pc_dir
            for pc_dir in (
                os.path.join(dev, "lib/pkgconfig"),
                os.path.join(dev, "share/pkgconfig"),
            )
            if os.path.isdir(pc_dir)
        ]

    def join(paths, var):
        """xs followed by the current value of search path `var`, colon-separated."""
        return ":".join(paths + ([env[var]] if env.get(var) else []))

    env["CPATH"] = join(inc, "CPATH")
    env["LIBRARY_PATH"] = join(lib, "LIBRARY_PATH")
    env["PKG_CONFIG_PATH"] = join(pkgconfig_dirs, "PKG_CONFIG_PATH")
    env["LDFLAGS"] = (
        " ".join(f"-Wl,-rpath,{lib_dir}" for lib_dir in lib) + " " + env.get("LDFLAGS", "")
    )
    return env


# Tools that write per-k-mer answers to their own file (-o, --ids).
# Their stdout must never be captured into that file.
OWN_OUTPUT = {"sbwt", "optik-spss --ids"}

# Tools with no thread option for the step are measured once, at t = 1.
# These are FMSI's build, i.e. KmerCamel compute and fmsi index, and the query of every compared tool, SSHash, SBWT and FMSI.
NO_BUILD_THREADS = {"fmsi"}


def query_tools(mode, work_dir, unitigs_fa, kmer_len, minimizer_len, threads, mach, cap_gib, qthreads=1):
    """Build and query commands of every compared k-mer index, keyed by tool name.

    Each entry gives the build command (None for a query view of another tool's index, named by `reuse`), the index files, a query function that maps a query file to (command, answer file), the kind of answers the tool writes and a parameter string for the report.
    `mode` is "parity" or "optimal", and indexes go to directory `d`."""
    par = mode == "parity"
    build_threads = threads
    query_threads = qthreads  # every paper ran single-threaded queries (1), and a sweep scales optik's
    m_str = str(minimizer_len)
    budget = sshash_budget(mach)  # the authors' 16 GiB, or less on a machine with less
    sbwt_gb = max(SBWT_MIN_MEM_GB, budget)  # the authors used 16
    return {
        "optik-spss": {
            "build": [
                BIN["optik"],
                "build-spss",
                "-i",
                unitigs_fa,
                "-k",
                str(kmer_len),
                "-m",
                # Parity runs pass the competitor's m, and -m 0 elsewhere lets choose_m() pick one from the input.
                m_str if par else "0",
                "-t",
                str(build_threads),
                "-o",
                f"{work_dir}/optik.spss",
            ],
            "files": [f"{work_dir}/optik.spss"],
            "query": lambda query_file: (
                [
                    BIN["optik"],
                    "query",
                    "-i",
                    f"{work_dir}/optik.spss",
                    "-q",
                    query_file,
                    "-t",
                    str(query_threads),
                ],
                None,
            ),
            "answers": "counts",
            "params": f"k={kmer_len} m={m_str if par else 'auto'} -t {build_threads}; query -t {query_threads}",
        },
        # The same index answering per k-mer, with a dense id per window and -1 when absent.
        # This is the like-for-like comparison with the per-k-mer output of SBWT and FMSI.
        # Timed runs write to /dev/null.
        "optik-spss --ids": {
            "build": None,
            "reuse": "optik-spss",
            "query": lambda query_file: (
                [
                    BIN["optik"],
                    "query",
                    "-i",
                    f"{work_dir}/optik.spss",
                    "-q",
                    query_file,
                    "-t",
                    str(query_threads),
                    "--ids",
                    f"{work_dir}/ids.out",
                ],
                f"{work_dir}/ids.out",
            ),
            "timed": lambda query_file: [
                BIN["optik"],
                "query",
                "-i",
                f"{work_dir}/optik.spss",
                "-q",
                query_file,
                "-t",
                str(query_threads),
                "--ids",
                "/dev/null",
            ],
            "answers": "per-k-mer",
            "params": f"the optik-spss index; query --ids (per-k-mer ids) -t {query_threads}",
        },
        "optik-fm": {
            "build": [
                BIN["optik"],
                "build-fm",
                "-i",
                unitigs_fa,
                "-k",
                str(kmer_len),
                "-c",
                "1",
                "-t",
                str(build_threads),
                "-o",
                f"{work_dir}/optik.fm",
            ],
            "files": [f"{work_dir}/optik.fm"],
            "query": lambda query_file: (
                [
                    BIN["optik"],
                    "query",
                    "-i",
                    f"{work_dir}/optik.fm",
                    "-q",
                    query_file,
                    "-t",
                    str(query_threads),
                ],
                None,
            ),
            "answers": "counts",
            "params": f"k={kmer_len} -t {build_threads}; query -t {query_threads}",
        },
        "sshash": {
            # The authors' own command, with -m per collection and a 16 GiB RAM budget, -g 16
            "build": [
                BIN["sshash"],
                "build",
                "-i",
                unitigs_fa,
                "-k",
                str(kmer_len),
                "-m",
                m_str,
                "-g",
                str(budget),
                "-t",
                str(build_threads),
                "-d",
                work_dir,
                "-o",
                f"{work_dir}/idx.sshash",
            ],
            "files": [f"{work_dir}/idx.sshash"],
            "query": lambda query_file: (
                [BIN["sshash"], "query", "-i", f"{work_dir}/idx.sshash", "-q", query_file],
                None,
            ),
            "answers": "counts",
            "params": f"k={kmer_len} m={m_str} -g {budget} -t {build_threads}; query single-threaded",
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
                f"{work_dir}/idx.sbwt",
                "-k",
                str(kmer_len),
                "--variant",
                "plain-matrix",
                "--add-reverse-complements",
                "-m",
                str(sbwt_gb),
                "-t",
                str(build_threads),
                "--temp-dir",
                work_dir,
            ],
            "files": [f"{work_dir}/idx.sbwt"],
            "extra": (
                [
                    BIN["sbwt"],
                    "build",
                    "-i",
                    unitigs_fa,
                    "-o",
                    f"{work_dir}/idx_norc.sbwt",
                    "-k",
                    str(kmer_len),
                    "--variant",
                    "plain-matrix",
                    "-m",
                    str(sbwt_gb),
                    "-t",
                    str(build_threads),
                    "--temp-dir",
                    work_dir,
                ],
                f"{work_dir}/idx_norc.sbwt",
            ),
            "query": lambda query_file: (
                [
                    BIN["sbwt"],
                    "search",
                    "-i",
                    f"{work_dir}/idx.sbwt",
                    "-q",
                    query_file,
                    "-o",
                    f"{work_dir}/sbwt.out",
                ],
                f"{work_dir}/sbwt.out",
            ),
            "timed": lambda query_file: [
                BIN["sbwt"],
                "search",
                "-i",
                f"{work_dir}/idx.sbwt",
                "-q",
                query_file,
                "-o",
                "/dev/null",
            ],
            "answers": "per-k-mer",
            "params": f"k={kmer_len} plain-matrix +rc -m {sbwt_gb} -t {build_threads}; search single-threaded",
        },
        "fmsi": {
            # The authors run kmercamel compute -S on SPSS input, then fmsi index.
            "build": [
                "bash",
                "-c",
                (
                    f"{BIN['kmercamel']} compute -S -k {kmer_len} -M {work_dir}/ms.fa -o /dev/null {unitigs_fa} && "
                    f"{BIN['fmsi']} index -k {kmer_len} {work_dir}/ms.fa"
                ),
            ],
            "files_glob": f"{work_dir}/ms.fa.fmsi.",
            # Both modes query with -S -O, the paper's settings.
            # kmer_sets_benchmark calls fmsi_lookup_*_query<minimal_hash=false, maximized_ones=true> over the kLCP array, which the CLI exposes as -O and -S.
            # These are FMSI's own optimizations rather than a parameter to match.
            "query": lambda query_file: (
                [
                    BIN["fmsi"],
                    "query",
                    "-S",
                    "-O",
                    "-q",
                    query_file,
                    "-k",
                    str(kmer_len),
                    f"{work_dir}/ms.fa",
                ],
                f"{work_dir}/fmsi.out",
            ),
            "answers": "per-k-mer",  # in timed runs hyperfine sends stdout to /dev/null
            "params": f"k={kmer_len} kmercamel compute -S + fmsi index; query -S -O (the paper's)",
        },
    }


def idx_files(tool, work_dir):
    """Index files of one query tool's build.

    These are the tool's explicit file list, or else every file in d whose name starts with its glob prefix."""
    if "files" in tool:
        return [index_file for index_file in tool["files"] if os.path.exists(index_file)]
    return [
        os.path.join(work_dir, index_file)
        for index_file in os.listdir(work_dir)
        if os.path.join(work_dir, index_file).startswith(tool["files_glob"])
    ]


def count_of(tool, outf, capt, err):
    """Positive k-mer count a query run reported, parsed from its answer file or captured output."""
    if tool in ("sbwt", "optik-spss --ids"):  # one id per window, -1 when absent
        return compare_tools_small_sets.sbwt_counts(outf)
    if tool == "fmsi":
        return compare_tools_small_sets.fmsi_counts(capt)
    return compare_tools_small_sets.parse_counts(open(capt).read() + err)
