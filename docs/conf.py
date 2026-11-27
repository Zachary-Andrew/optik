"""Sphinx configuration for optik's documentation site.

Breathe renders the C++ core from the Doxygen XML that docs/build.sh writes first.
autodoc renders the bench's shared Python modules straight from their docstrings.
The Rust CLI's API comes from cargo doc, which build.sh copies under rust/ in the built site.
lean_docs.py reads the Lean modules directly, and myst parses their Markdown docstrings along with lean/APPENDIX_B.md.
"""

import os
import sys

DOCS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(DOCS)

# autodoc imports the bench modules, which import their siblings by bare name.
sys.path.insert(0, os.path.join(ROOT, "bench"))
sys.path.insert(0, DOCS)

project = "optik"
author = "Zachary-Andrew"

extensions = [
    "breathe",
    "sphinx.ext.autodoc",
    "myst_parser",
    "lean_docs",
]

breathe_projects = {"optik": os.path.join(DOCS, "_build", "doxygen", "xml")}
breathe_default_project = "optik"
breathe_default_members = ("members", "private-members", "undoc-members")

autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": True}

# Every file page declares the namespace tpoptoa again, because the namespace spans all of them.
# Members are never duplicated, since doxygen_fixups.py drops the copies Doxygen makes for using-declarations.
suppress_warnings = ["duplicate_declaration.cpp"]

primary_domain = "cpp"
highlight_language = "cpp"

html_theme = "furo"
html_title = "optik"
exclude_patterns = ["_build", ".venv"]
