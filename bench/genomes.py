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

from units import BYTES_PER_MIB

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEST = os.path.join(ROOT, "temp", "genomes")
NCBI_GCA = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"
# bytes read per block while hashing a genome
HASH_BLOCK_BYTES = 4 * BYTES_PER_MIB

# name -> (file name, gzipped source URL, SHA-256 of the uncompressed file, description)
GENOMES = {
    "ecoli": (
        "GCA_000005845.2_ASM584v2_genomic.fna",
        f"{NCBI_GCA}/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz",
        "4baa96ad53f849f485d7125913453b7ec34ebba03ea6489c16e839371f3867e9",
        "E. coli K-12 MG1655, GenBank GCA_000005845.2 (4.6 Mbp), from NCBI",
    ),
}


def _sha256(file_path):
    """SHA-256 of the file at `p`."""
    digest = hashlib.sha256()
    with open(file_path, "rb") as handle:
        for block in iter(lambda: handle.read(HASH_BLOCK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def path(name):
    """Path of genome `name` under temp/genomes/, downloaded and verified first when it is missing.

    Raises SystemExit when the download fails or the file does not match its checksum."""
    fname, url, sha, _ = GENOMES[name]
    dest = os.path.join(DEST, fname)
    if os.path.exists(dest):
        return dest
    os.makedirs(DEST, exist_ok=True)
    # Temporary names are per process, so two benches fetching at once never share a file.
    gz_part, part = f"{dest}.{os.getpid()}.gz.part", f"{dest}.{os.getpid()}.part"
    print(f"genomes: fetching {name} from {url}", file=sys.stderr)
    result = subprocess.run(
        ["curl", "-fsSL", "--retry", "5", "--retry-all-errors", "-o", gz_part, url],
        check=False,
    )
    if result.returncode != 0:
        if os.path.exists(gz_part):
            os.remove(gz_part)
        raise SystemExit(
            f"genomes: download of {name} failed (curl exit {result.returncode})"
        )
    with gzip.open(gz_part, "rb") as src, open(part, "wb") as out:
        shutil.copyfileobj(src, out)
    os.remove(gz_part)
    got = _sha256(part)
    if got != sha:
        os.remove(part)
        raise SystemExit(f"genomes: {name} has SHA-256 {got}, expected {sha}")
    os.replace(part, dest)
    return dest


def main():
    """Prints a genome's path after fetching it, or lists the genomes."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("name", nargs="?", choices=sorted(GENOMES))
    parser.add_argument(
        "--list", action="store_true", help="list the genomes and their sources"
    )
    args = parser.parse_args()
    if args.list or not args.name:
        for genome, (fname, url, _, what) in GENOMES.items():
            have = (
                "present"
                if os.path.exists(os.path.join(DEST, fname))
                else "not fetched"
            )
            print(
                f"{genome:8} {what}\n         {url}\n         temp/genomes/{fname} ({have})"
            )
        return
    print(path(args.name))


if __name__ == "__main__":
    main()
