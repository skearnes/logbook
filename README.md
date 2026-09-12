# Logbook

A running work logbook. Each entry is a dated directory under
[`entries/`](entries/), named `YYYY-MM-DD-topic/`, holding a `README.md` with
the writeup plus any assets (images, scripts, data) for that entry.

## Conventions

- **All changes go through a pull request into `main` — never commit or push
  directly to `main`.** Work on a branch, open a PR, and merge it (CI runs the
  markdown lint, the link check, the metadata check, and the license-header
  check on every PR).
- Start a new entry by creating `entries/YYYY-MM-DD-topic/` and copying
  [`TEMPLATE.md`](TEMPLATE.md) into it as `README.md`.
- Name each directory `entries/YYYY-MM-DD-topic/` (e.g.
  `entries/2026-03-14-pi-day/`).
- Keep entry assets (images, scripts, data) in the same directory and link to
  them with relative paths from the `README.md`.
- Begin each entry's `README.md` with an H1 title covering the date and topic.
- Use H2 sections within an entry as needed.

Every Markdown file under `entries/` — an entry's `README.md` and any
supporting document beside it — opens with a metadata block under its title:

```markdown
- **Date:** 2026-03-14
- **Author:** Ada Lovelace
- **Status:** draft
- **Tags:** topic, another-topic
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
```

`Date` is the entry's date, matching the directory name (and, for a `README.md`,
the date in its title). `Status` wording is open — qualify it ("final (archive)")
when that says more than one word can. `License` must read exactly as above.
Supporting documents need only `Date`, `Author`, and `License`: `Status` and
`Tags` describe the investigation, which is the entry's job, and are optional
beside it. A document prepared with Claude adds an `Acknowledgments` line after
`Author` that names Claude Code and the model, as in
`- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 5)`;
delete the line from the template otherwise.
[`.github/scripts/check_metadata.py`](.github/scripts/check_metadata.py)
enforces all of this in CI and as a pre-commit hook; run it directly with:

```bash
python3 .github/scripts/check_metadata.py
```

Relative links are checked the same way.
[`.github/scripts/check_links.py`](.github/scripts/check_links.py) verifies that
every relative link across the tracked Markdown resolves, that any heading
anchor it names exists, and that link text spelling a path agrees with the path
it points at. External URLs are not fetched:

```bash
python3 .github/scripts/check_links.py
```

## Large assets

Assets too large or too regenerable for git (caches, embeddings, model
artifacts) live in the Google Cloud Storage bucket
[`gs://skearnes-logbook`](https://console.cloud.google.com/storage/browser/skearnes-logbook)
(project `skearnes`), mirroring the repo path. For example an entry's
`entries/<date-topic>/assets/cache/` maps to
`gs://skearnes-logbook/entries/<date-topic>/assets/cache/`. Keep these paths out
of git via a local `.gitignore` and sync them with:

```bash
gsutil rsync -r <entry>/assets/cache \
  gs://skearnes-logbook/<entry>/assets/cache
```

Note (macOS): do **not** use `gsutil -m` here — multiprocessing forks workers
that segfault when a child does a DNS lookup via Apple's non-fork-safe
Network.framework ("Python quit unexpectedly"). Set `parallel_process_count = 1`
(and e.g. `parallel_thread_count = 24`) under `[GSUtil]` in `~/.boto` so syncs
use threads, not forked processes. Also install the compiled crcmod C extension
(`python3 -m pip install --no-binary :all: crcmod`) so checksumming is fast.

## Linting

Markdown is linted with [markdownlint-cli2][cli2] in CI on every push and pull
request. To run the same check locally:

```bash
npx markdownlint-cli2 "**/*.md"
```

The same lint runs as a [pre-commit][pre-commit] hook. Enable it once with:

```bash
pre-commit install
```

Lint rules live in [`.markdownlint-cli2.jsonc`](.markdownlint-cli2.jsonc).

## License

Dual-licensed, following the [Open Reaction Database][ord] split:

- **Prose** — every `README.md` and other Markdown writeup, plus figures — is
  licensed [CC BY-SA 4.0](LICENSE-CC-BY-SA). Quote and adapt it with
  attribution; derivative writeups carry the same license.
- **Code** — the `.py` scripts under each entry's `assets/`, the workflows and
  scripts under [`.github/`](.github/), and the repository config files — is
  licensed [Apache-2.0](LICENSE-APACHE). Reuse it in your own projects under
  those terms.

Where the two could both plausibly apply (a fenced code block inside a
writeup), take the Apache-2.0 terms.

Every code file carries an Apache-2.0 header. [addlicense][addlicense] adds and
checks them, as a pre-commit hook and in CI; Markdown has no comment syntax it
can use, so each Markdown file states its license in the metadata block instead
(see [Conventions](#conventions)). To stamp new code files locally:

```bash
go install github.com/google/addlicense@v1.2.0
addlicense -c "Steven Kearnes" -l apache <files>
```

Third-party benchmark datasets and leaderboard exports (Polaris/ASAP,
OpenADMET) are **not** redistributed here — they stay with their upstream
sources under their own licenses, and entries link to where to fetch them.
Results computed from those datasets are covered by this repo's license.

[ord]: https://github.com/open-reaction-database

[pre-commit]: https://pre-commit.com

[cli2]: https://github.com/DavidAnson/markdownlint-cli2

[addlicense]: https://github.com/google/addlicense
