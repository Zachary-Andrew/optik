"""Downloads, verifies and builds the compared tools.

The short prefix downloads that the smoke run relies on come from here as well.
"""

import concurrent.futures
import datetime
import glob
import json
import os
import subprocess
import tarfile
import time
import urllib.parse
from pathlib import Path

import rangefetch
from published_catalog import DATA_ROOT, K8_URL, SEQAN_TAG, TOOLS
from published_fetch import (
    CURL,
    RANGED_MIN,
    ZIP_SUFFIX,
    _ranged,
    checksum,
    file_urls,
    unzip_ncbi,
    verify_by_source,
)
from published_tools import dev_env, ggcat_bin
from units import BYTES_PER_GB

from published.common import (
    DOWNLOAD_ATTEMPTS,
    DOWNLOAD_RETRY_STEP_S,
    EXEC_MODE,
    PREFIX_EXT,
    SMOKE_PREFIX_BYTES,
)

# Two read files are enough to show a read set's format in a smoke slice.
SMOKE_READ_FILES = 2


def _verified_copy(url, name, dest):
    """A complete, verified copy of url's file under another dataset's directory whose PROVENANCE.json records the same URL as verified, or None."""
    for prov in glob.glob(os.path.join(DATA_ROOT.path, "*", "PROVENANCE.json")):
        try:
            with open(prov) as handle:
                recs = json.load(handle).get("files", [])
        except (OSError, ValueError):
            continue
        for rec in recs:
            path = rec.get("path", "")
            if (
                rec.get("url") == url
                and rec.get("verified")
                and os.path.basename(path) == name
                and os.path.isfile(path)
                and not os.path.samefile(os.path.dirname(path), os.path.dirname(dest))
            ):
                return path
    return None


def download(url, dest, algo=None, digest=None, size=None):
    """Fetch url to dest with curl, resuming partial files.

    The file is verified against the source's checksum when the source publishes one (Zenodo, ENA).
    Otherwise its MD5 is recorded."""
    if (
        os.path.exists(dest)
        and (size is None or os.path.getsize(dest) == size)
        and (digest is None or checksum(dest, algo) == digest)
    ):
        return {
            "url": url,
            "path": dest,
            "verified": digest is not None,
            "cached": True,
        }
    same = _verified_copy(url, os.path.basename(dest), dest)
    if same:
        # The same file of the same source was already verified for another dataset, as when CHM13's NCBI zip serves two, so it is linked instead of fetched again.
        for stale in (dest, dest + ".part"):
            if os.path.exists(stale):
                os.remove(stale)
        os.link(same, dest)
        print(
            f"    {os.path.basename(dest)}: linked from {os.path.relpath(same, DATA_ROOT.path)}",
            flush=True,
        )
        return {
            "url": url,
            "path": dest,
            "verified": True,
            "cached": True,
            "linked_from": same,
        }
    part = dest + ".part"
    # -sS because several downloads run at once (fetch --jobs) and progress bars would interleave
    print(
        f"    downloading {os.path.basename(dest)}"
        + (f" ({size / BYTES_PER_GB:.2f} GB)" if size else ""),
        flush=True,
    )
    total = size if size is not None else rangefetch.content_length(url)
    ranged = _ranged(url) if total is not None else False
    if ranged is None:
        # Guessing "no ranges" here would restart the file from zero and throw away a partial .part.
        raise SystemExit(f"{dest}: could not probe {urllib.parse.urlsplit(url).netloc} for range support; retry the fetch")
    with rangefetch.part_lock(part):
        if total is not None and total >= RANGED_MIN and ranged:
            # Many paths limit throughput per connection (2026-10-16 through a VPN, about 1.5 MB/s each and 10 MB/s over 8 or more).
            # Large files are therefore split into byte ranges spread over connections that every download shares (bench/rangefetch.py).
            rangefetch.scheduler(
                log=print, progress=os.path.join(DATA_ROOT.path, "_fetch_progress.json")
            ).fetch(url, part, total)
        else:
            # A .part with a .done beside it was preallocated by the range scheduler, so it is full-size and sparse, not a prefix curl could resume.
            # Its intervals would also outlive the rewrite and be trusted later, so both go.
            if os.path.exists(part + ".done"):
                for stale in (part, part + ".done"):
                    os.remove(stale)
            # A stalled transfer resumes with -C - where the server serves ranges.
            # Where it does not, as NCBI Datasets builds its zips per request, the transfer restarts from zero, because curl refuses to resume there and each retry would fail the same way.
            for attempt in range(DOWNLOAD_ATTEMPTS):
                if not ranged and os.path.exists(part):
                    os.remove(part)
                if (
                    subprocess.run(
                        CURL + (["-C", "-"] if ranged else []) + ["-o", part, url],
                        check=False,
                    ).returncode
                    == 0
                ):
                    break
                time.sleep(DOWNLOAD_RETRY_STEP_S * (attempt + 1))
            else:
                raise SystemExit(f"{dest}: download failed after {DOWNLOAD_ATTEMPTS} attempts")
        if size is not None and os.path.getsize(part) != size:
            raise SystemExit(f"{dest}: size {os.path.getsize(part)} != published {size}")
        got = checksum(part, algo or "md5")
        if digest is not None and got != digest:
            os.remove(part)
            raise SystemExit(f"{dest}: {algo} {got} != published {digest}")
        os.replace(part, dest)
    return {
        "url": url,
        "path": dest,
        algo or "md5": got,
        "verified": digest is not None,
        "date": datetime.datetime.now().astimezone().date().isoformat(),
    }


def fetch_one(ds_id, dataset):
    """Download, verify and unpack every file of dataset `i`, then write its PROVENANCE.json."""
    if "excluded" in dataset:
        print(f"  {ds_id}: excluded -- {dataset['excluded']}")
        return
    out = os.path.join(DATA_ROOT.path, ds_id)
    os.makedirs(out, exist_ok=True)
    prov = {
        "id": ds_id,
        "title": dataset.get("title"),
        "paper": dataset["paper"],
        "source": dataset["source"],
        "use": dataset.get("use"),
        "files": [],
    }
    files = file_urls(dataset)
    # All files go to the range scheduler together, so a connection moves on to the next file's ranges while the last ranges of a file finish.
    def fetch_one(file_tuple):
        """Download one (url, name, *options) entry into the dataset directory."""
        url, name, *options = file_tuple
        return download(url, os.path.join(out, name), *options)

    with concurrent.futures.ThreadPoolExecutor(max(1, len(files))) as pool:
        recs = list(pool.map(fetch_one, files))
    prov["files"] += recs
    for (_url, name, *_), rec in zip(files, recs, strict=False):
        if name.endswith(".zip"):
            fa_path = os.path.join(out, name[: -len(ZIP_SUFFIX)] + ".fa")
            src_fa = rec.get("linked_from", "")[: -len(ZIP_SUFFIX)] + ".fa"
            if (
                not os.path.exists(fa_path)
                and rec.get("linked_from")
                and os.path.isfile(src_fa)
            ):
                os.link(src_fa, fa_path)  # the linked zip's own extraction
            unzip_ncbi(os.path.join(out, name), fa_path)
    if dataset.get("local"):
        if not os.path.exists(dataset["local"]):
            raise SystemExit(
                f"{ds_id}: local input {dataset['local']} missing (bench/build_tools.sh hg38)"
            )
        prov["files"].append({"path": dataset["local"], "local": True})
    for rec in prov["files"]:
        if not rec.get("local") and not rec.get("verified"):
            rec["verified_by"] = verify_by_source(rec["url"], rec["path"])
            rec["verified"] = True
    Path(os.path.join(out, "PROVENANCE.json")).write_text(json.dumps(prov, indent=1))
    print(f"  {ds_id}: ok ({len(prov['files'])} files)", flush=True)


def build_tools():
    """Build the tools this bench needs beyond bench/build_tools.sh.

    These are GGCAT, PBSIM2, k8 to run paftools.js, Winnowmap with meryl, and Mason's simulator, all under TOOLS and linked there by name."""
    os.makedirs(TOOLS, exist_ok=True)
    if not os.path.exists(ggcat_bin()):
        # Latest default branch.
        # The deposit's eulertigs came from 14b2853, and the k-mer set is the same for every version, which the KMC gate checks.
        subprocess.run(
            [
                "cargo",
                "install",
                "--locked",
                "--git",
                "https://github.com/algbio/ggcat",
                "--root",
                os.path.join(TOOLS, "ggcat"),
                "ggcat_cmdline",
            ],
            check=True,
        )
    pb2 = os.path.join(TOOLS, "pbsim2")
    if not os.path.exists(os.path.join(pb2, "src/pbsim")):
        if not os.path.exists(pb2):
            subprocess.run(
                ["git", "clone", "-q", "https://github.com/yukiteruono/pbsim2", pb2],
                check=True,
            )
        # A fresh clone's timestamps make `make` rerun autoconf, which then lacks config.guess.
        subprocess.run(
            "autoreconf -fi && ./configure && make -j8", shell=True, cwd=pb2, check=True
        )
    k8_path = os.path.join(TOOLS, "k8")
    if not os.path.exists(k8_path):
        tarball = os.path.join(TOOLS, "k8.tar.bz2")
        download(K8_URL, tarball)
        with tarfile.open(tarball) as tar:
            member = next(member for member in tar.getmembers() if member.name.endswith("k8-x86_64-Linux"))
            member.name = "k8"
            tar.extract(member, TOOLS, filter="data")
        os.chmod(k8_path, EXEC_MODE)
    winnowmap_dir = os.path.join(TOOLS, "winnowmap_src")
    if not os.path.exists(os.path.join(winnowmap_dir, "bin/winnowmap")):
        if not os.path.exists(winnowmap_dir):
            subprocess.run(
                [
                    "git",
                    "clone",
                    "-q",
                    "--depth",
                    "1",
                    "https://github.com/marbl/Winnowmap",
                    winnowmap_dir,
                ],
                check=True,
            )
        # Older meryl copies include <sys/sysctl.h>, which glibc 2.32 and later no longer ship.
        # meryl uses it only on BSD and macOS, so each such include gets a guard.
        for root, _, files in os.walk(os.path.join(winnowmap_dir, "ext")):
            for fname in files:
                file_path = os.path.join(root, fname)
                if not fname.endswith((".C", ".H", ".c", ".h", ".cc", ".hpp", ".cpp")):
                    continue
                txt = Path(file_path).read_text(errors="replace")
                if (
                    "#include <sys/sysctl.h>\n" in txt
                    and "__has_include(<sys/sysctl.h>)" not in txt
                ):
                    Path(file_path).write_text(
                        txt.replace(
                            "#include <sys/sysctl.h>\n",
                            "#if __has_include(<sys/sysctl.h>)\n#include <sys/sysctl.h>\n#endif\n",
                        )
                    )
        subprocess.run(
            ["make", "-j8"], cwd=winnowmap_dir, check=True, env=dev_env(["curl", "openssl"])
        )
    seqan_dir = os.path.join(TOOLS, "seqan_src")
    mason_dir = os.path.join(TOOLS, "mason_build")
    if not os.path.exists(os.path.join(mason_dir, "bin/mason_simulator")):
        if not os.path.exists(seqan_dir):
            subprocess.run(
                [
                    "git",
                    "clone",
                    "-q",
                    "--depth",
                    "1",
                    "--branch",
                    SEQAN_TAG,
                    "https://github.com/seqan/seqan",
                    seqan_dir,
                ],
                check=True,
            )
        os.makedirs(mason_dir, exist_ok=True)
        subprocess.run(
            [
                "cmake",
                seqan_dir,
                "-DCMAKE_BUILD_TYPE=Release",
                "-DSEQAN_BUILD_SYSTEM=APP:mason2",
            ],
            cwd=mason_dir,
            check=True,
        )
        subprocess.run(["make", "-j8", "mason_simulator"], cwd=mason_dir, check=True)
    for name, target in (
        ("winnowmap", os.path.join(winnowmap_dir, "bin/winnowmap")),
        ("meryl", os.path.join(winnowmap_dir, "bin/meryl")),
        ("mason_simulator", os.path.join(mason_dir, "bin/mason_simulator")),
    ):
        lnk = os.path.join(TOOLS, name)
        if not os.path.lexists(lnk):
            os.symlink(target, lnk)
    print(
        f"ggcat {ggcat_bin()}\npbsim2 {pb2}/src/pbsim\nk8 {k8_path}\nwinnowmap/meryl {winnowmap_dir}/bin\n"
        f"mason_simulator {mason_dir}/bin/mason_simulator"
    )


def fetch_prefix(ds_id, dataset, nbytes=SMOKE_PREFIX_BYTES):
    """Fetch the first nbytes of each file not yet fetched, through HTTP Range.

    That is enough for a smoke slice, and a compressed prefix decompresses up to its cut."""
    src = os.path.join(DATA_ROOT.path, ds_id)
    os.makedirs(src, exist_ok=True)
    if dataset.get("kind") == "reads" or dataset.get("reads"):
        urls = file_urls(dict(dataset, files=dataset["files"][:1]))[:SMOKE_READ_FILES]
    else:
        urls = file_urls(dataset)
    for url, name, _, _, _ in urls:
        dest = os.path.join(src, name)
        if os.path.exists(dest) or os.path.exists(dest + ".prefix"):
            continue
        if name.endswith((".zip", ".tar.gz")):
            continue  # an archive is only readable whole
        subprocess.run(
            ["curl", "-fsL", "-r", f"0-{nbytes - 1}", "-o", dest + ".prefix", url],
            check=True,
        )
    for fname in os.listdir(src):  # prefixes keep the real names, in a side directory
        if fname.endswith(".prefix"):
            os.makedirs(os.path.join(src, "_prefix"), exist_ok=True)
            os.replace(os.path.join(src, fname), os.path.join(src, "_prefix", fname[: -len(PREFIX_EXT)]))
