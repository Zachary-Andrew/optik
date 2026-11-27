"""Download helpers of compare_tools_published_sets.py.

The module resolves Zenodo and ENA files, runs ranged parallel downloads (rangefetch.py), computes checksums and verifies fetched data against each source.
Functions that need the data root stay in compare_tools_published_sets.py, which can move that root at run time.
"""

import concurrent.futures
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
import urllib.parse
import zipfile
from typing import NamedTuple

import rangefetch
from diagnostics import DETAIL_CHARS
from published_inputs import _open_text

# download attempts, the backoff unit in seconds, and parallel downloads
MAX_ATTEMPTS = 4
RETRY_BACKOFF_S = 5
FETCH_THREADS = 8
# bytes hashed per read
HASH_CHUNK_BYTES = 1 << 22
# hex characters of an abbreviated commit
SHORT_SHA_LEN = 7
ZIP_SUFFIX = ".zip"


class SourceFile(NamedTuple):
    """A file a dataset is downloaded from, with the checksum its host publishes."""

    url: str
    checksum_algorithm: str
    digest: str
    size: int


def _http_text(url):
    """Fetch a metadata URL through curl with retries.

    ENA's portal answers urllib with truncated bodies, and with transient 500 errors under load."""
    result = subprocess.run(
        [
            "curl",
            "-fsSL",
            "--retry",
            "6",
            "--retry-all-errors",
            "--retry-delay",
            "5",
            "--max-time",
            "300",
            url,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise OSError(f"curl exit {result.returncode} for {url}: {result.stderr.strip()[-DETAIL_CHARS:]}")
    return result.stdout


def _get_json(url):
    """Parsed JSON body of url, fetched by _http_text."""
    return json.loads(_http_text(url))


_ZCACHE = {}


def zenodo_file(rec, key):
    """The SourceFile of file `key` in Zenodo record `rec`, cached per record."""
    if rec not in _ZCACHE:
        _ZCACHE[rec] = {
            file_meta["key"]: file_meta
            for file_meta in _get_json(f"https://zenodo.org/api/records/{rec}")["files"]
        }
    file_meta = _ZCACHE[rec][key]
    algo, digest = file_meta["checksum"].split(":", 1)
    return SourceFile(file_meta["links"]["self"], algo, digest, file_meta["size"])


_ENA = {}


def ena_files(acc):
    """The SourceFile of every FASTQ file of a run, experiment, study or project, with the MD5 and size ENA publishes.

    ENA's portal sometimes answers with only the header line.
    That case is retried instead of being read as "no files"."""
    if acc in _ENA:
        return _ENA[acc]
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession="
        + acc
        + "&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes&format=tsv&limit=0"
    )
    out = []
    # The ENA portal occasionally returns an empty report for a valid accession, so the query is retried with a growing delay.
    for attempt in range(MAX_ATTEMPTS):
        lines = _http_text(url).strip().split("\n")
        head = lines[0].split("\t")
        for line in lines[1:]:
            row = dict(zip(head, line.split("\t"), strict=False))
            if not row.get("fastq_ftp"):
                continue
            out += [
                SourceFile("https://" + url_path, "md5", md5, int(size_str))
                for url_path, md5, size_str in zip(
                    # Paired runs list their files, checksums and sizes as parallel semicolon lists.
                    row["fastq_ftp"].split(";"),
                    row["fastq_md5"].split(";"),
                    row["fastq_bytes"].split(";"),
                    strict=False,
                )
            ]
        if out:
            break
        time.sleep(RETRY_BACKOFF_S * (attempt + 1))
    if not out:
        raise ValueError(f"ENA lists no FASTQ files for {acc} (4 attempts)")
    _ENA[acc] = out
    return out


def ena_prefetch(accs):
    """Look up many ENA accessions concurrently, since each lookup takes about 25 s."""
    todo = [acc for acc in dict.fromkeys(accs) if acc not in _ENA]
    if todo:
        with concurrent.futures.ThreadPoolExecutor(FETCH_THREADS) as pool:
            list(pool.map(ena_files, todo))


def checksum(path, algo):
    """Hex digest of the file at path under hashlib algorithm `algo`, read in 4 MiB blocks."""
    digest = hashlib.new(algo)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(HASH_CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


# A connection that stalls without an error would otherwise hang curl indefinitely, so a transfer slower than 10 kB/s for a minute is dropped and resumed.
CURL = [
    "curl",
    "-fsSL",
    "--retry",
    "5",
    "--retry-all-errors",
    "--speed-limit",
    "10000",
    "--speed-time",
    "60",
]

RANGED_MIN = 8 << 20  # files below this size use a single connection

_RANGES = {}


def _ranged(url):
    """True when url's host serves byte ranges, False when it does not, and None when the probe got no answer.

    Only a definite answer is cached, so one dropped probe does not decide every later file of the host."""
    host = urllib.parse.urlsplit(url).netloc
    if host not in _RANGES:
        answer = rangefetch.supports_ranges(url)
        if answer is None:
            return None
        _RANGES[host] = answer
    return _RANGES[host]


def fasta_lengths(path):
    """Map each record's first header token to its sequence length in a plain or .gz FASTA file."""
    out, name = {}, None
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            elif name is not None:
                out[name] += len(line.strip())
    return out


def verify_by_source(url, path):
    """Check a file whose source publishes no checksum against what the source does publish.

    GitHub raw files at a commit are compared by the git blob hash from GitHub's contents API, NCBI nuccore records by esummary's slen length, NCBI assemblies by the length of every sequence in Datasets sequence_reports and ENA WGS sets by base count.
    A mismatch raises."""
    # A GitHub raw file at a pinned commit must hash to the blob that GitHub's contents API reports.
    match = re.match(r"https://github.com/([^/]+)/([^/]+)/raw/([0-9a-f]{40})/(.+)", url)
    if match:
        owner, repo, rev, file_path = match.groups()
        meta = _get_json(
            f"https://api.github.com/repos/{owner}/{repo}/contents/{file_path}?ref={rev}"
        )
        got = subprocess.run(
            ["git", "hash-object", path], capture_output=True, text=True, check=True
        ).stdout.strip()
        if got != meta["sha"]:
            raise SystemExit(
                f"{path}: git blob {got} != GitHub's {meta['sha']} at {rev[:SHORT_SHA_LEN]}"
            )
        return f"git blob {got} = GitHub contents API at {rev[:SHORT_SHA_LEN]}"
    # An NCBI nucleotide download must have the total length that esummary gives.
    match = re.search(r"efetch\.fcgi\?db=nuccore.*&id=([^&]+)", url)
    if match:
        acc = match.group(1)
        summary = _get_json(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={acc}&retmode=json"
        )
        want = int(summary["result"][summary["result"]["uids"][0]]["slen"])
        got = sum(fasta_lengths(path).values())
        if got != want:
            raise SystemExit(f"{path}: {got} bases != NCBI's {want} for {acc}")
        return f"length {got} = NCBI esummary slen"
    # An NCBI Datasets genome must match its sequence report, sequence by sequence.
    match = re.search(r"/genome/accession/([^/]+)/download", url)
    if match:
        acc = match.group(1)
        rep = _get_json(
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/sequence_reports"
        )
        want = {report["genbank_accession"]: int(report["length"]) for report in rep["reports"]}
        got = fasta_lengths(path[: -len(ZIP_SUFFIX)] + ".fa")
        if got != want:
            raise SystemExit(
                f"{path}: sequence lengths differ from NCBI Datasets' report for {acc}"
            )
        return f"all {len(want)} sequence lengths = NCBI Datasets sequence_reports"
    # An ENA WGS set must have the base count that the ENA portal reports.
    match = re.search(r"/ena/wgs/public/[^/]+/([A-Z]{4}\d{2})\.fasta\.gz$", url)
    if match:
        wgs_set = match.group(1)
        rows = (
            _http_text(
                "https://www.ebi.ac.uk/ena/portal/api/search?result=wgs_set&query=wgs_set=%22"
                f"{wgs_set}%22&fields=wgs_set,base_count&format=tsv"
            )
            .strip()
            .split("\n")
        )
        want = int(
            dict(zip(rows[0].split("\t"), rows[1].split("\t"), strict=False))["base_count"]
        )
        got = sum(fasta_lengths(path).values())
        if got != want:
            raise SystemExit(f"{path}: {got} bases != ENA's base_count {want} for {wgs_set}")
        return f"{got} bases = ENA wgs_set base_count"
    raise SystemExit(
        f"{path}: its source publishes neither a checksum nor a record to verify it against"
    )


def file_urls(dataset):
    """List a dataset's files as (url, local name, algo, digest, size) tuples."""
    out = []
    # The ENA accessions are looked up concurrently first, since each lookup takes about 25 s.
    ena_prefetch([file_spec["ena"] for file_spec in dataset["files"] if "ena" in file_spec])
    for file_spec in dataset["files"]:
        # Zenodo and ENA publish a checksum and size for each file, so those entries include them.
        if "zenodo" in file_spec:
            url, algo, dig, size = zenodo_file(file_spec["zenodo"], file_spec["key"])
            out.append((url, file_spec["key"], algo, dig, size))
        elif "ena" in file_spec:
            out += [
                (url, os.path.basename(url), algo, digest, size) for url, algo, digest, size in ena_files(file_spec["ena"])
            ]
        # Plain URLs and the NCBI services publish no checksum, and verify_by_source() checks those files instead.
        elif "url" in file_spec:
            out.append((file_spec["url"], os.path.basename(file_spec["url"]), None, None, None))
        # Sources without a checksum get None fields and are checked only by parsing.
        elif "efetch" in file_spec:
            acc = file_spec["efetch"]
            out.append(
                (
                    (
                        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta"
                        f"&retmode=text&id={acc}"
                    ),
                    acc + ".fa",
                    None,
                    None,
                    None,
                )
            )
        elif "ncbi" in file_spec:
            acc = file_spec["ncbi"]
            out.append(
                (
                    (
                        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
                        "?include_annotation_type=GENOME_FASTA"
                    ),
                    acc + ".zip",
                    None,
                    None,
                    None,
                )
            )
    return out


def unzip_ncbi(zip_path, fasta_path):
    """Concatenate the FASTA members of NCBI Datasets zip z into fa, once."""
    if not os.path.exists(fasta_path):
        with zipfile.ZipFile(zip_path) as zip_file, open(fasta_path + ".part", "wb") as out_file:
            for member in sorted(zip_file.namelist()):
                if member.endswith((".fna", ".fa", ".fasta")):
                    with zip_file.open(member) as src:
                        shutil.copyfileobj(src, out_file)
        os.replace(fasta_path + ".part", fasta_path)
