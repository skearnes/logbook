# Copyright 2026 Steven Kearnes
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Check the metadata block every Markdown file under ``entries/`` opens with.

An entry states its date, author, status, tags, and license in a bullet list
under the title. Supporting documents beside it — plans, specs, generated tables
— state date, author, and license; status and tags describe the investigation,
which is the entry's job, so they are optional there and checked only when
present. This verifies the required fields are filled in, that the stated date
matches the entry directory, and that the license is the one the repository
actually grants. An entry's ``README.md`` must additionally open with a
``# YYYY-MM-DD — Topic`` title carrying that same date; supporting documents
keep whatever title they came with.

The repository's own ``README.md`` and ``TEMPLATE.md`` are not entries and are
not checked.

Run from the repository root:

    python3 .github/scripts/check_metadata.py
"""

import argparse
import re
import sys
from pathlib import Path

FIELD = re.compile(r"^- \*\*(?P<name>[A-Za-z]+):\*\*\s*(?P<value>.*)$")
PLACEHOLDER = re.compile(r"^`?<.*>`?$|^YYYY-MM-DD$|^draft \| final$")
SLUG = re.compile(r"^(?P<date>\d{4}-\d{2}-\d{2})-[a-z0-9-]+$")
TITLE = re.compile(r"^# (?P<date>\d{4}-\d{2}-\d{2}) — .+$")
LICENSE = "[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)"
REQUIRED = ("Date", "Author", "License")
# Only the entry itself has to be classified and tagged; supporting documents
# may carry these, and are checked when they do.
ENTRY_REQUIRED = ("Status", "Tags")


def parse(path: Path) -> tuple[str, dict[str, str]]:
    """Read a document's title line and the metadata bullets under it.

    Args:
        path: Markdown file to read.

    Returns:
        The H1 title line (empty if the file does not open with one) and the
        metadata fields, field name mapped to value. Scanning stops at the first
        section heading, so bolded bullets in the body are not mistaken for
        metadata.
    """
    title = ""
    found = {}
    for line in path.read_text().splitlines():
        if line.startswith("## "):
            break
        if line.startswith("# ") and not title:
            title = line
        match = FIELD.match(line)
        if match:
            found[match.group("name")] = match.group("value").strip()
    return title, found


def check(path: Path, root: Path) -> list[str]:
    """Check one document's metadata block.

    Args:
        path: Markdown file under ``entries/`` to check.
        root: Repository root, used to render paths in messages.

    Returns:
        One message per problem found, empty if the document is clean.
    """
    problems = []
    rel = path.relative_to(root)
    entry = (root / "entries" / rel.parts[1]) if len(rel.parts) > 1 else None
    is_entry = path.name == "README.md" and path.parent == entry
    title, found = parse(path)

    required = REQUIRED + ENTRY_REQUIRED if is_entry else REQUIRED
    for name in required:
        value = found.get(name)
        if value is None:
            problems.append(f"{rel}: missing metadata field: {name}")
        elif not value or PLACEHOLDER.match(value):
            problems.append(f"{rel}: {name} is still the template placeholder")

    # Optional fields still have to be filled in if a document states them.
    for name in () if is_entry else ENTRY_REQUIRED:
        value = found.get(name)
        if value is not None and (not value or PLACEHOLDER.match(value)):
            problems.append(f"{rel}: {name} is still the template placeholder")

    slug = SLUG.match(entry.name) if entry else None
    if not slug:
        problems.append(f"{rel}: entry directory is not named YYYY-MM-DD-topic")
    elif "Date" in found and found["Date"] != slug.group("date"):
        problems.append(
            f"{rel}: Date {found['Date']} does not match the entry date "
            f"{slug.group('date')}"
        )

    # Only the entry itself carries the date in its title; supporting documents
    # beside it keep whatever title they came with.
    if path.name == "README.md" and path.parent == entry:
        heading = TITLE.match(title)
        if not heading:
            problems.append(f"{rel}: title is not '# YYYY-MM-DD — Topic'")
        elif slug and heading.group("date") != slug.group("date"):
            problems.append(
                f"{rel}: title date {heading.group('date')} does not match the "
                f"entry date {slug.group('date')}"
            )
    elif not title:
        problems.append(f"{rel}: no H1 title")

    # Status wording is deliberately open — documents qualify it ("final (spec
    # pending)"), so only an unedited template value is a problem.
    license_ = found.get("License", "")
    if license_ and license_ != LICENSE:
        problems.append(f"{rel}: License should read {LICENSE}")

    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root", nargs="?", default=".", help="repository root to check"
    )
    args = parser.parse_args()
    root = Path(args.root).resolve()

    documents = sorted(root.glob("entries/*/**/*.md"))
    problems = []
    for path in documents:
        problems.extend(check(path, root))

    for problem in problems:
        print(problem)
    print(f"checked {len(documents)} documents, {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
