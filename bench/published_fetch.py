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
import zipfile

import rangefetch
from published_inputs import _open_text


def _http_text(url):
    """Fetch a metadata URL through curl with retries.

    ENA's portal answers urllib with truncated bodies, and with transient 500 errors under load."""
    r = subprocess.run(
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
    if r.returncode != 0:
        raise OSError(f"curl exit {r.returncode} for {url}: {r.stderr.strip()[-200:]}")
    return r.stdout


def _get_json(url):
    """Parsed JSON body of url, fetched by _http_text."""
    return json.loads(_http_text(url))


_ZCACHE = {}


def zenodo_file(rec, key):
    """(URL, checksum algorithm, digest, size) of file `key` in Zenodo record `rec`, cached per record."""
    if rec not in _ZCACHE:
        _ZCACHE[rec] = {
            f["key"]: f
            for f in _get_json(f"https://zenodo.org/api/records/{rec}")["files"]
        }
    f = _ZCACHE[rec][key]
    algo, digest = f["checksum"].split(":", 1)
    return f["links"]["self"], algo, digest, f["size"]


_ENA = {}


def ena_files(acc):
    """List every FASTQ file of a run, experiment, study or project, with the MD5 and size ENA publishes.

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
    for attempt in range(4):
        lines = _http_text(url).strip().split("\n")
        head = lines[0].split("\t")
        for line in lines[1:]:
            row = dict(zip(head, line.split("\t"), strict=False))
            if not row.get("fastq_ftp"):
                continue
            out += [
                ("https://" + u, "md5", m, int(b))
                for u, m, b in zip(
                    row["fastq_ftp"].split(";"),
                    row["fastq_md5"].split(";"),
                    row["fastq_bytes"].split(";"),
                    strict=False,
                )
            ]
        if out:
            break
        time.sleep(5 * (attempt + 1))
    if not out:
        raise ValueError(f"ENA lists no FASTQ files for {acc} (4 attempts)")
    _ENA[acc] = out
    return out


def ena_prefetch(accs):
    """Look up many ENA accessions concurrently, since each lookup takes about 25 s."""
    todo = [a for a in dict.fromkeys(accs) if a not in _ENA]
    if todo:
        with concurrent.futures.ThreadPoolExecutor(8) as ex:
            list(ex.map(ena_files, todo))


def checksum(path, algo):
    """Hex digest of the file at path under hashlib algorithm `algo`, read in 4 MiB blocks."""
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


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
    """True when url's host serves byte ranges, probed once per host."""
    host = url.split("/")[2]
    if host not in _RANGES:
        _RANGES[host] = rangefetch.supports_ranges(url)
    return _RANGES[host]


def fasta_lengths(path):
    """Map each record's first header token to its sequence length in a plain or .gz FASTA file."""
    out, name = {}, None
    with _open_text(path) as f:
        for line in f:
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
    m = re.match(r"https://github.com/([^/]+)/([^/]+)/raw/([0-9a-f]{40})/(.+)", url)
    if m:
        owner, repo, rev, fp = m.groups()
        meta = _get_json(
            f"https://api.github.com/repos/{owner}/{repo}/contents/{fp}?ref={rev}"
        )
        got = subprocess.run(
            ["git", "hash-object", path], capture_output=True, text=True, check=True
        ).stdout.strip()
        if got != meta["sha"]:
            raise SystemExit(
                f"{path}: git blob {got} != GitHub's {meta['sha']} at {rev[:7]}"
            )
        return f"git blob {got} = GitHub contents API at {rev[:7]}"
    m = re.search(r"efetch\.fcgi\?db=nuccore.*&id=([^&]+)", url)
    if m:
        acc = m.group(1)
        r = _get_json(
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={acc}&retmode=json"
        )
        want = int(r["result"][r["result"]["uids"][0]]["slen"])
        got = sum(fasta_lengths(path).values())
        if got != want:
            raise SystemExit(f"{path}: {got} bases != NCBI's {want} for {acc}")
        return f"length {got} = NCBI esummary slen"
    m = re.search(r"/genome/accession/([^/]+)/download", url)
    if m:
        acc = m.group(1)
        rep = _get_json(
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/sequence_reports"
        )
        want = {x["genbank_accession"]: int(x["length"]) for x in rep["reports"]}
        got = fasta_lengths(path[:-4] + ".fa")
        if got != want:
            raise SystemExit(
                f"{path}: sequence lengths differ from NCBI Datasets' report for {acc}"
            )
        return f"all {len(want)} sequence lengths = NCBI Datasets sequence_reports"
    m = re.search(r"/ena/wgs/public/[^/]+/([A-Z]{4}\d{2})\.fasta\.gz$", url)
    if m:
        ws = m.group(1)
        t = (
            _http_text(
                "https://www.ebi.ac.uk/ena/portal/api/search?result=wgs_set&query=wgs_set=%22"
                f"{ws}%22&fields=wgs_set,base_count&format=tsv"
            )
            .strip()
            .split("\n")
        )
        want = int(
            dict(zip(t[0].split("\t"), t[1].split("\t"), strict=False))["base_count"]
        )
        got = sum(fasta_lengths(path).values())
        if got != want:
            raise SystemExit(f"{path}: {got} bases != ENA's base_count {want} for {ws}")
        return f"{got} bases = ENA wgs_set base_count"
    raise SystemExit(
        f"{path}: its source publishes neither a checksum nor a record to verify it against"
    )


def file_urls(d):
    """List a dataset's files as (url, local name, algo, digest, size) tuples."""
    out = []
    ena_prefetch([f["ena"] for f in d["files"] if "ena" in f])
    for f in d["files"]:
        if "zenodo" in f:
            url, algo, dig, size = zenodo_file(f["zenodo"], f["key"])
            out.append((url, f["key"], algo, dig, size))
        elif "ena" in f:
            out += [
                (u, os.path.basename(u), a, g, s) for u, a, g, s in ena_files(f["ena"])
            ]
        elif "url" in f:
            out.append((f["url"], os.path.basename(f["url"]), None, None, None))
        elif "efetch" in f:
            acc = f["efetch"]
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
        elif "ncbi" in f:
            acc = f["ncbi"]
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


def unzip_ncbi(z, fa):
    """Concatenate the FASTA members of NCBI Datasets zip z into fa, once."""
    if not os.path.exists(fa):
        with zipfile.ZipFile(z) as zf, open(fa + ".part", "wb") as o:
            for n in sorted(zf.namelist()):
                if n.endswith((".fna", ".fa", ".fasta")):
                    with zf.open(n) as src:
                        shutil.copyfileobj(src, o)
        os.replace(fa + ".part", fa)
