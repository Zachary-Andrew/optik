"""Sphinx directive that documents a Lean module of lean/OptikProofs from its source.

``.. lean-module:: Kmer`` turns the module's ``/-! -/`` header into a section and lists every declaration that has a ``/-- -/`` docstring.
Each entry shows the docstring and the signature up to ``:=``, while proof bodies are left out.
Docstrings are Markdown, as Lean's doc-gen4 reads them, so myst parses them instead of reStructuredText.
"""

import os
import re

from docutils import nodes
from docutils.utils import new_document
from myst_parser.parsers.sphinx_ import MystParser
from sphinx.util import logging
from sphinx.util.docutils import SphinxDirective

LOG = logging.getLogger(__name__)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEAN_DIR = os.path.join(ROOT, "lean")
SOURCE_DIR = os.path.join(LEAN_DIR, "OptikProofs")

# The library root imports every module, which lets the build list the ones the page leaves out.
LIBRARY_ROOT = os.path.join(LEAN_DIR, "OptikProofs.lean")

# Attributes such as @[simp] and modifiers such as noncomputable may precede the keyword.
PREFIX = re.compile(r"\s*(@\[[^\]]*\]\s*)*(noncomputable\s+|private\s+|protected\s+)*")
KEYWORD = re.compile(r"(theorem|lemma|def|abbrev|structure|inductive|instance)\s+(\S+)")

# A pattern-matching definition opens its equations with "|" on the line after the signature.
EQUATION = re.compile(r"\s*\|")

# A signature this long means the parser missed its end, which the build should report.
MAX_SIGNATURE_LINES = 12


def comment_text(lines):
    """Strips the comment markers and the common indent from a docstring's lines."""
    text = "\n".join(lines).strip()
    text = re.sub(r"^/-[-!]\s?", "", text)
    text = re.sub(r"\s*-/$", "", text)
    body = text.split("\n")
    indents = [len(s) - len(s.lstrip()) for s in body[1:] if s.strip()]
    cut = min(indents, default=0)
    return "\n".join([body[0]] + [s[cut:] for s in body[1:]])


def block_end(lines, start):
    """Index of the line that closes the comment opened on line ``start``."""
    end = start

    while "-/" not in lines[end]:
        end += 1

    return end


def signature(lines, start):
    """The declaration from line ``start`` up to its ``:=``, ``where`` or first equation, or None when no end is found."""
    sig = []

    for line in lines[start:start + MAX_SIGNATURE_LINES]:
        if sig and EQUATION.match(line):
            return "\n".join(sig)

        head, assign, _ = line.partition(":=")

        if assign:
            sig.append(head.rstrip())
            return "\n".join(sig)

        if line.rstrip().endswith(" where"):
            sig.append(line.rstrip()[: -len(" where")])
            return "\n".join(sig)

        sig.append(line.rstrip())

    return None


def parse_module(path):
    """The header and the documented declarations of one Lean file, plus the lines the parser could not read."""
    with open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")

    header = None
    decls = []
    unread = []
    i = 0

    while i < len(lines):
        stripped = lines[i].lstrip()

        if stripped.startswith("/-!") and header is None:
            end = block_end(lines, i)
            header = comment_text(lines[i:end + 1])
            i = end + 1
            continue

        if not stripped.startswith("/--"):
            i += 1
            continue

        end = block_end(lines, i)
        doc = comment_text(lines[i:end + 1])
        at = end + 1
        m = KEYWORD.match(lines[at][PREFIX.match(lines[at]).end():]) if at < len(lines) else None
        sig = signature(lines, at) if m else None

        if sig is None:
            unread.append(end + 2)
        else:
            decls.append((m.group(1), m.group(2), at + 1, doc, sig))

        i = at

    return header, decls, unread


def markdown(module, decls):
    """The Markdown for one module's declaration list, with a download link to its source."""
    source = os.path.relpath(os.path.join(SOURCE_DIR, module + ".lean"), os.path.join(ROOT, "docs"))
    parts = [f"Source: {{download}}`{module}.lean <{source}>`"]

    for kind, name, line, doc, sig in decls:
        parts.append(f"**`{module}.{name}`** ({kind}, line {line})")
        parts.append(doc)
        parts.append("```lean\n" + sig + "\n```")

    return "\n\n".join(parts)


class LeanModule(SphinxDirective):
    """``.. lean-module:: Name`` renders lean/OptikProofs/Name.lean."""

    required_arguments = 1

    def parse_md(self, text):
        doc = new_document(str(self.env.doc2path(self.env.docname)), self.state.document.settings)
        MystParser().parse(text, doc)
        return doc.children

    def run(self):
        module = self.arguments[0]
        path = os.path.join(SOURCE_DIR, module + ".lean")
        self.env.note_dependency(path)
        self.env.lean_documented = getattr(self.env, "lean_documented", set()) | {module}

        if not os.path.exists(path):
            return [self.state.document.reporter.warning(f"no Lean module {path}", line=self.lineno)]

        header, decls, unread = parse_module(path)

        for line in unread:
            LOG.warning("docstring not followed by a declaration this parser reads", location=(path, line))

        if header is None or not header.startswith("# "):
            return [self.state.document.reporter.warning(f"{path} has no /-! # Title -/ header", line=self.lineno)]

        title_text, _, body = header.partition("\n")
        section = nodes.section(ids=[nodes.make_id("lean-" + module)])

        # myst wraps the title in a paragraph, and the section title takes only its inline children.
        title_para = self.parse_md(title_text[len("# "):])[0]
        section += nodes.title("", "", *title_para.children)
        section += self.parse_md(body)
        section += self.parse_md(markdown(module, decls))
        return [section]


def check_coverage(app, env):
    """Warns about modules OptikProofs.lean imports that no lean-module directive documents."""
    with open(LIBRARY_ROOT, encoding="utf-8") as f:
        imported = set(re.findall(r"^import OptikProofs\.(\w+)", f.read(), re.M))

    for module in sorted(imported - getattr(env, "lean_documented", set())):
        LOG.warning(f"Lean module {module} is imported by OptikProofs.lean but not documented")


def merge_documented(app, env, docnames, other):
    env.lean_documented = getattr(env, "lean_documented", set()) | getattr(other, "lean_documented", set())


def setup(app):
    app.add_directive("lean-module", LeanModule)
    app.connect("env-check-consistency", check_coverage)
    app.connect("env-merge-info", merge_documented)
    return {"parallel_read_safe": True}
