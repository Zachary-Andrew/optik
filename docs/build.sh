#!/usr/bin/env bash
# Builds the documentation site into docs/_build/html.
# Doxygen and Graphviz come from the system, and the pinned Sphinx side from docs/.venv, created on first use.
# Breathe reads the XML Doxygen writes, so Doxygen runs first, followed by the fix-ups that let Sphinx parse it.
# The Rust reference is copied in last, under rust/, where docs/rust.rst links to it.
set -euo pipefail

cd "$(dirname "$0")/.."

if [ ! -x docs/.venv/bin/sphinx-build ]; then
    python3 -m venv docs/.venv
fi

# A venv made before requirements.txt last changed, such as one without myst-parser, is brought up to date.
if [ docs/requirements.txt -nt docs/.venv/installed ]; then
    docs/.venv/bin/pip install --quiet -r docs/requirements.txt
    touch docs/.venv/installed
fi

rm -rf docs/_build
mkdir -p docs/_build/doxygen
doxygen docs/Doxyfile

# Doxygen exits 0 even with warnings, so a non-empty log fails the build.
if [ -s docs/_build/doxygen/warnings.log ]; then
    cat docs/_build/doxygen/warnings.log >&2
    exit 1
fi

python3 docs/doxygen_fixups.py docs/_build/doxygen/xml
docs/.venv/bin/sphinx-build -W --keep-going -q -b html docs docs/_build/html

cargo doc --no-deps --document-private-items
cp -r target/doc docs/_build/html/rust

echo "site: docs/_build/html/index.html"
