"""Paths and tool binaries shared by compare_tools_small_sets.py, external_oracles.py and compare_tools_published_sets.py.

Sequence helpers live in dna.py.
"""

import os

from config import PATHS

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


TEMP_DIR = os.path.join(ROOT, "temp")

BIN = {
    "optik": os.path.join(ROOT, "target/release/optik"),
    "sshash": str(PATHS["sshash"]),
    "sbwt": str(PATHS["sbwt"]),
    "fmsi": str(PATHS["fmsi"]),
    "kmercamel": str(PATHS["kmercamel"]),
    "cuttlefish": str(PATHS["cuttlefish"]),
    "bifrost": str(PATHS["bifrost"]),
    "minimap2": str(PATHS["minimap2"]),
    "wfmash": str(PATHS["wfmash"]),
    "kmc": os.path.join(TEMP_DIR, "KMC/bin/kmc"),
    "kmc_dump": os.path.join(TEMP_DIR, "KMC/bin/kmc_dump"),
    "kmc_tools": os.path.join(TEMP_DIR, "KMC/bin/kmc_tools"),
}

DATA = os.path.join(TEMP_DIR, "sshash/data")
