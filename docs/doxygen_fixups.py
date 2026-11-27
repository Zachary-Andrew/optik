"""Rewrite Doxygen's XML in place so Breathe and Sphinx's C++ domain can read all of it.

Doxygen names an anonymous namespace anonymous_namespace{file.cpp}, which Sphinx cannot parse as a C++ name, so Breathe stops with an assertion.
Each one becomes anonymous_<file stem>, which is unique per file just like the original.
A `using ns::name;` line is also recorded as a member of the file that holds it.
Breathe would then declare that name a second time, so members whose declaration line is a using-declaration are dropped from file pages.

Usage: python3 docs/doxygen_fixups.py docs/_build/doxygen/xml
"""

import pathlib
import re
import sys

ANON = re.compile(r"anonymous_namespace\{(\w+)\.\w+\}")
FILE_MEMBER = re.compile(r'<memberdef kind="\w+" id="(?P<id>[^"]+)".*?<location file="(?P<file>[^"]+)" line="(?P<line>\d+)".*?</memberdef>\s*', re.DOTALL)


REPO = pathlib.Path(__file__).resolve().parent.parent
USING = re.compile(r"\s*using\s+[\w:]+::\w+\s*;")


def using_lines():
    """Set of (repository-relative file, line) pairs whose source line is a using-declaration."""
    found = set()
    for source_dir in ("include", "src"):
        for path in (REPO / source_dir).iterdir():
            for number, text in enumerate(path.read_text(errors="replace").splitlines(), 1):
                if USING.match(text):
                    found.add((f"{source_dir}/{path.name}", number))
    return found


def drop_using_members(text, usings):
    """`text` without the members whose declaration line is in `usings`, along with their member-list entries."""
    dropped = set()

    def keep_or_drop(match):
        if (match.group("file"), int(match.group("line"))) in usings:
            dropped.add(match.group("id"))
            return ""
        return match.group(0)

    text = FILE_MEMBER.sub(keep_or_drop, text)
    for member_id in dropped:
        text = re.sub(rf'\s*<member refid="{re.escape(member_id)}"[^>]*>.*?</member>', "", text, flags=re.DOTALL)
    # A section emptied of members would render as a bare heading.
    return re.sub(r'\s*<sectiondef kind="\w+">\s*</sectiondef>', "", text)


def main(xml_dir):
    usings = using_lines()
    for path in pathlib.Path(xml_dir).glob("*.xml"):
        text = path.read_text()
        fixed = ANON.sub(r"anonymous_\1", text)
        if 'kind="file"' in fixed[: fixed.find("<compoundname>")]:
            fixed = drop_using_members(fixed, usings)
        if fixed != text:
            path.write_text(fixed)


if __name__ == "__main__":
    main(sys.argv[1])
