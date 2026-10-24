#!/usr/bin/env python3
"""genomes.py - the small reference genomes the walkthrough and benches use, fetched on first use.

No genome is kept in the repository.
Each one here is downloaded from its archive into temp/genomes/, checked against a fixed SHA-256 of the uncompressed file, and reused afterwards.

  python3 bench/genomes.py ecoli        prints the path, downloading the genome first when needed
  python3 bench/genomes.py --list       lists the genomes and where they come from

  import genomes; genomes.path("ecoli") does the same from a bench.
"""

import argparse
import gzip
import hashlib
import os
import shutil
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEST = os.path.join(ROOT, "temp", "genomes")
NCBI_GCA = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"

# name -> (file name, gzipped source URL, SHA-256 of the uncompressed file, description)
GENOMES = {
    "ecoli": (
        "GCA_000005845.2_ASM584v2_genomic.fna",
        f"{NCBI_GCA}/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz",
        "4baa96ad53f849f485d7125913453b7ec34ebba03ea6489c16e839371f3867e9",
        "E. coli K-12 MG1655, GenBank GCA_000005845.2 (4.6 Mbp), from NCBI",
    ),
}


def _sha256(p):
    """SHA-256 of the file at `p`."""
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def path(name):
    """Path of genome `name` under temp/genomes/, downloaded and verified first when it is missing.

    Raises SystemExit when the download fails or the file does not match its checksum."""
    fname, url, sha, _ = GENOMES[name]
    dest = os.path.join(DEST, fname)
    if os.path.exists(dest):
        return dest
    os.makedirs(DEST, exist_ok=True)
    # Temporary names are per process, so two benches fetching at once never share a file.
    gz, part = f"{dest}.{os.getpid()}.gz.part", f"{dest}.{os.getpid()}.part"
    print(f"genomes: fetching {name} from {url}", file=sys.stderr)
    r = subprocess.run(
        ["curl", "-fsSL", "--retry", "5", "--retry-all-errors", "-o", gz, url],
        check=False,
    )
    if r.returncode != 0:
        if os.path.exists(gz):
            os.remove(gz)
        raise SystemExit(
            f"genomes: download of {name} failed (curl exit {r.returncode})"
        )
    with gzip.open(gz, "rb") as src, open(part, "wb") as out:
        shutil.copyfileobj(src, out)
    os.remove(gz)
    got = _sha256(part)
    if got != sha:
        os.remove(part)
        raise SystemExit(f"genomes: {name} has SHA-256 {got}, expected {sha}")
    os.replace(part, dest)
    return dest


def main():
    """Prints a genome's path after fetching it, or lists the genomes."""
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("name", nargs="?", choices=sorted(GENOMES))
    ap.add_argument(
        "--list", action="store_true", help="list the genomes and their sources"
    )
    a = ap.parse_args()
    if a.list or not a.name:
        for n, (fname, url, _, what) in GENOMES.items():
            have = (
                "present"
                if os.path.exists(os.path.join(DEST, fname))
                else "not fetched"
            )
            print(
                f"{n:8} {what}\n         {url}\n         temp/genomes/{fname} ({have})"
            )
        return
    print(path(a.name))


if __name__ == "__main__":
    main()
