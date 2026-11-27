"""Works out what this machine can run.

It estimates each tool's memory, adds up download and disk needs, and decides which datasets fit.
"""

import functools
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import hw
from external_oracles import SAMTOOLS, picard_cmd
from published_catalog import CATALOG, DATA_ROOT, TIERS
from published_fetch import ena_files, ena_prefetch, zenodo_file
from published_tools import have_tool
from small_common import BIN
from units import BITS_PER_BYTE, BYTES_PER_GB, BYTES_PER_GIB, BYTES_PER_MIB

from published.common import (
    _FALLBACK_BPK,
    _PHASE_TOOL,
    _SIZES,
    BUDGETED,
    DECOMPRESSED_PER_GB,
    EST_MEM_OVERHEAD_GB,
    FETCH_ERRORS,
    GENOME_WORK_FACTOR,
    INDEX_BITS,
    KMERS_PER_COMPRESSED_GB,
    MAP_DOWNLOAD_FACTOR,
    MAP_LOCAL_FACTOR,
    MATES_PER_PAIR,
    NCBI_ASSEMBLY_BYTES,
    READS_WORK_FACTOR,
    REFERENCE_K,
    SAM_BYTES_PER_READ,
    SCALING_LOGS,
)

# The estimate is printed in millions of k-mers.
KMERS_PER_MILLION = 1e6


def machine():
    """Usable memory in GiB, hw.mem_bytes(), i.e. MemTotal or the tightest cgroup memory.max on this process's path such as a systemd-run -p MemoryMax= scope, and the free disk in GB where the datasets go."""
    mem = hw.mem_bytes() / BYTES_PER_GIB
    os.makedirs(DATA_ROOT.path, exist_ok=True)
    return {
        "mem_gib": round(mem, 1),
        "disk_free_gb": round(shutil.disk_usage(DATA_ROOT.path).free / BYTES_PER_GB, 1),
    }


def bytes_per_kmer():
    """{k: {tool: bytes per distinct k-mer}} from the chr1 rows of the scaling logs, or from the fallback run above."""
    peaks = {}
    for log in SCALING_LOGS:
        if not os.path.exists(log):
            continue
        with open(log) as handle:
            for line in handle:
                # Rows are written by scaling_bench.py's row(): label (with k=), tool, phase, wall s, peak MB, status, key=value notes.
                fields = line.rstrip("\n").split("\t")
                try:
                    label, tool_name, phase, _wall_s, peak_mb, status, *notes = fields
                except ValueError:  # a header or truncated line
                    continue
                if not label.startswith("chr1") or status != "ok":
                    continue
                kmer_len = int(label.split("k=")[1]) if "k=" in label else REFERENCE_K
                row = peaks.setdefault(kmer_len, {})
                if tool_name == "optik-spss" and phase == "build (default m)":
                    match = re.search(r"kmers=(\d+)", " ".join(notes))
                    if match:
                        row["_kmers"] = int(match.group(1))
                tool = _PHASE_TOOL.get((tool_name, phase))
                if tool and peak_mb:
                    row[tool] = max(row.get(tool, 0), float(peak_mb))
    src = peaks if REFERENCE_K in peaks and "_kmers" in peaks[REFERENCE_K] else _FALLBACK_BPK
    out = {}
    for kmer_len, row in src.items():
        if "_kmers" in row:
            out[kmer_len] = {
                tool: peak_mb * BYTES_PER_MIB / row["_kmers"]
                for tool, peak_mb in row.items()
                if tool != "_kmers"
            }
    return out, (
        "scaling logs "
        + ", ".join(os.path.basename(log) for log in SCALING_LOGS if os.path.exists(log))
        if src is peaks
        else "chr1 scaling run of 2026-10-12 (fallback table)"
    )


@functools.lru_cache(maxsize=1)
def cached_bytes_per_kmer():
    """bytes_per_kmer(), read from the scaling logs only once per run."""
    return bytes_per_kmer()


def est_mem_gb(tool, kmer_len, kmers):
    """Estimated peak in GB of `tool` building over `kmers` distinct k-mers.

    The tool's measured bytes per k-mer are taken at the nearest measured k at or below k, since shorter k costs more per k-mer, and at k = 31 when nothing smaller was measured.
    None if the tool was never measured."""
    table, _ = cached_bytes_per_kmer()
    measured_ks = sorted(table)
    for table_k in [measured_k for measured_k in measured_ks if measured_k <= kmer_len][::-1] + [REFERENCE_K]:
        if tool in table.get(table_k, {}):
            return table[table_k][tool] * kmers / BYTES_PER_GB + EST_MEM_OVERHEAD_GB
    return None


def tool_skip(dataset, tool, mach, kmer_len=None, kmers=None):
    """Why `tool` cannot run on this dataset here, or None.

    A memory need the paper reports for this dataset takes precedence over the measured estimate.
    Tools that keep to a memory budget (BUDGETED) are given this machine's and are never skipped for it.
    SBWT's -m is "not strictly enforced", so SBWT is estimated like the others."""
    if tool in dataset.get("skip_tools", {}):
        return dataset["skip_tools"][tool]
    allow_gb = mach["mem_gib"] * BYTES_PER_GIB / BYTES_PER_GB
    pub = dataset.get("mem", {}).get(tool)
    if pub:
        if pub[0] > allow_gb:
            return (
                f"needs {pub[0]} GB ({pub[1]}); this machine allows {allow_gb:.1f} GB"
            )
        return None
    if tool in BUDGETED or not kmers:
        return None
    est = est_mem_gb(tool, kmer_len or REFERENCE_K, kmers)
    if est and est > allow_gb:
        return (
            f"estimated {est:.1f} GB ({kmers / KMERS_PER_MILLION:.0f} M k-mers x measured bytes per k-mer, "
            f"{cached_bytes_per_kmer()[1]}); this machine allows {allow_gb:.1f} GB"
        )
    return None


def download_gb(ds_id, dataset):
    """Bytes the sources publish for the dataset's files (Zenodo, ENA or the HTTP Content-Length), cached.
    None where a source gives no size."""
    cache = json.loads(Path(_SIZES).read_text()) if os.path.exists(_SIZES) else {}
    if ds_id in cache:
        return cache[ds_id]
    tot = 0
    try:
        ena_prefetch([file_spec["ena"] for file_spec in dataset["files"] if "ena" in file_spec])
        for file_spec in dataset["files"]:
            if "zenodo" in file_spec:
                tot += zenodo_file(file_spec["zenodo"], file_spec["key"]).size
            elif "ena" in file_spec:
                tot += sum(ena_file.size for ena_file in ena_files(file_spec["ena"]))
            elif "url" in file_spec:
                hdr = subprocess.run(
                    ["curl", "-fsSLI", "--retry", "4", file_spec["url"]],
                    capture_output=True,
                    text=True,
                    check=False,
                ).stdout
                sizes = re.findall(r"(?im)^content-length:\s*(\d+)", hdr)
                tot += int(sizes[-1]) if sizes else 0
            elif "ncbi" in file_spec:
                tot += NCBI_ASSEMBLY_BYTES  # a ~3 Gbp assembly, zipped
    except FETCH_ERRORS as err:
        print(f"  {ds_id}: size lookup failed ({err})")
        return None
    cache[ds_id] = round(tot / BYTES_PER_GB, 2)
    Path(_SIZES).write_text(json.dumps(cache, indent=1))
    return cache[ds_id]


def present_gb(ds_id):
    """GB already on disk under dataset `i`'s directory."""
    ds_dir = os.path.join(DATA_ROOT.path, ds_id)
    if not os.path.isdir(ds_dir):
        return 0.0
    return (
        sum(os.path.getsize(os.path.join(root, fname)) for root, _, files in os.walk(ds_dir) for fname in files)
        / BYTES_PER_GB
    )


def disk_parts(ds_id, dataset):
    """(download GB still to fetch, working GB while the dataset runs, how).

    Compressed FASTA expands about 4.5x, because gzip and xz both store DNA near 2 bits per base, i.e. 8 / 1.8 bytes.
    A k-mer set's SPSS holds about one base per k-mer, and SSHash and KMC write temporary files about the size of their input.
    A mapping run holds its reads, one SAM per mapper at about 2 kB per simulated 150 bp read pair end over all files, and a copy of the reference."""
    download_size = download_gb(ds_id, dataset) or 0.0
    left = max(0.0, download_size - present_gb(ds_id))
    if dataset.get("as_") and not dataset["files"]:
        # a scaling set reads its source's files, or a share of them, so it downloads nothing
        download_size = (download_gb(dataset["as_"], CATALOG[dataset["as_"]]) or 0.0) * dataset.get("runs_frac", 1.0)
        left = 0.0
    if dataset["section"] == "query":
        kmers = dataset.get("kmers31") or download_size * KMERS_PER_COMPRESSED_GB
        spss = kmers / BYTES_PER_GB
        idx = sum(INDEX_BITS.values()) * kmers / (BITS_PER_BYTE * BYTES_PER_GB)
        work = (spss + idx) * len(dataset.get("k", [REFERENCE_K])) + spss
        how = f"download {download_size:.1f} + SPSS {spss:.1f} + indexes {idx:.1f} + temp {spss:.1f} GB per k"
    elif dataset["section"] == "dbg":
        # Read sets, as measured on cf2-gut with 48 GB of FASTQ.gz.
        # The temporary files of optik and Cuttlefish 3 peaked at 103-106 GB, 2.2x the input, and the KMC databases the rows are judged by stay on disk beside them, about another 1x.
        work = download_size * READS_WORK_FACTOR if dataset.get("kind") == "reads" else download_size * DECOMPRESSED_PER_GB * GENOME_WORK_FACTOR
        how = (
            f"input {download_size:.1f} GB + graph temporaries and KMC databases (3x)"
            if dataset.get("kind") == "reads"
            else f"download {download_size:.1f} GB + decompressed genome and unitigs"
        )
    else:
        local = (
            os.path.getsize(dataset["local"]) / BYTES_PER_GB
            if dataset.get("local") and os.path.exists(dataset["local"])
            else 0.0
        )
        work = (
            download_size * MAP_DOWNLOAD_FACTOR
            + local * MAP_LOCAL_FACTOR
            + dataset.get("pairs", 0) * MATES_PER_PAIR * SAM_BYTES_PER_READ / BYTES_PER_GB
        )
        how = f"download {download_size:.1f} GB + reads, SAMs and reference copies"
    return left, work, how


def disk_need(ds_id, dataset):
    """(GB, how) one dataset needs free, counting its remaining download and its working files."""
    left, work, how = disk_parts(ds_id, dataset)
    return left + work, how


def preflight(ds_id, dataset, mach):
    """(runnable, reasons) of dataset `i` on this machine."""
    if "excluded" in dataset:
        return False, ["excluded: " + dataset["excluded"]]
    why = []
    miss = [tool for tool in dataset.get("needs", []) if not have_tool(tool)]
    if dataset.get("spss") and not (
        have_tool("ggcat/bin/ggcat") or os.path.exists(BIN["cuttlefish"])
    ):
        miss.append("GGCAT or Cuttlefish (the SPSS step)")
    if dataset["section"] == "map":
        # the external checks every mapper's SAM goes through
        if not SAMTOOLS:
            miss.append("samtools (PATH, or bench/build_tools.sh samtools)")
        if not picard_cmd():
            miss.append("Picard (PATH, or bench/build_tools.sh picard; needs java 17+)")
    if dataset.get("local") and not os.path.exists(dataset["local"]):
        why.append(f"needs {dataset['local']} (bench/build_tools.sh hg38 / scaling)")
    if (
        dataset.get("as_")
        and CATALOG[dataset["as_"]]["files"]
        and not os.path.exists(os.path.join(DATA_ROOT.path, dataset["as_"], "PROVENANCE.json"))
    ):
        why.append(f"needs {dataset['as_']} fetched (fetch --only {dataset['as_']})")
    if dataset.get("stream_from") and not os.path.exists(
        os.path.join(DATA_ROOT.path, dataset["stream_from"], dataset["stream"])
    ):
        why.append(f"needs {dataset['stream']} (fetch --only {dataset['stream_from']})")
    if miss:
        why.append(
            "not installed: "
            + ", ".join(miss)
            + " (on PATH, or `compare_tools_published_sets.py tools` "
            "for GGCAT, pbsim2, k8, Winnowmap/meryl and mason2)"
        )
    if dataset["files"] and download_gb(ds_id, dataset) is None:
        why.append(
            "download size could not be determined (source unreachable?); --force to try anyway"
        )
    need, how = disk_need(ds_id, dataset)
    if need > mach["disk_free_gb"]:
        why.append(
            f"needs ~{need:.0f} GB free disk ({how}), {mach['disk_free_gb']:.0f} GB free"
        )
    return not why, why


def selected(args, mach):
    """The (id, dataset) pairs this invocation runs.

    Applies --only, --skip, --section and --tier, then each dataset's preflight, and finally checks that the selection fits the free disk as a whole (unless --force)."""
    out, download_total, peak = [], 0.0, 0.0
    only = set(args.only.split(",")) if args.only else None
    skip = set(args.skip.split(",")) if getattr(args, "skip", None) else set()
    unknown = (skip | (only or set())) - set(CATALOG)
    if unknown:
        raise SystemExit("unknown dataset ids: " + ", ".join(sorted(unknown)))
    for ds_id, dataset in CATALOG.items():
        if (only is not None and ds_id not in only) or ds_id in skip:
            continue
        if args.section and dataset["section"] not in args.section.split(","):
            continue
        if "excluded" in dataset:
            continue
        # An opt-in dataset runs only when --only names it, never through --tier or --force.
        if "opt_in" in dataset and only is None:
            continue
        if only is None and args.tier != "auto" and TIERS[dataset["tier"]] > TIERS[args.tier]:
            continue
        runnable, why = preflight(ds_id, dataset, mach)
        if not runnable and not args.force:
            print(f"  {ds_id}: skipped on this machine -- {'; '.join(why)}")
            continue
        # The selection must fit as a whole.
        # Downloads stay and add up, while working files go after each dataset, so only the largest of those is needed at once.
        left, work, how = disk_parts(ds_id, dataset)
        if download_total + left + max(peak, work) > mach["disk_free_gb"] and not args.force:
            print(
                f"  {ds_id}: skipped -- needs ~{left + work:.0f} GB ({how}); "
                f"{mach['disk_free_gb'] - download_total - peak:.0f} GB left after the datasets selected before it"
            )
            continue
        download_total += left
        peak = max(peak, work)
        out.append((ds_id, dataset))
    return out
