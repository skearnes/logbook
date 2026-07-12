# Logbook

A running work logbook. Each entry is a dated directory under
[`entries/`](entries/), named `YYYY-MM-DD-topic/`, holding a `README.md` with
the writeup plus any assets (images, scripts, data) for that entry.

## Conventions

- Start a new entry by creating `entries/YYYY-MM-DD-topic/` and copying
  [`TEMPLATE.md`](TEMPLATE.md) into it as `README.md`.
- Name each directory `entries/YYYY-MM-DD-topic/` (e.g.
  `entries/2026-03-14-pi-day/`).
- Keep entry assets (images, scripts, data) in the same directory and link to
  them with relative paths from the `README.md`.
- Begin each entry's `README.md` with an H1 title covering the date and topic.
- Use H2 sections within an entry as needed.

## Large assets

Assets too large or too regenerable for git (caches, embeddings, model
artifacts) live in the Google Cloud Storage bucket
[`gs://skearnes-logbook`](https://console.cloud.google.com/storage/browser/skearnes-logbook)
(project `skearnes`), mirroring the repo path. For example an entry's
`entries/<date-topic>/assets/cache/` maps to
`gs://skearnes-logbook/entries/<date-topic>/assets/cache/`. Keep these paths out
of git via a local `.gitignore` and sync them with:

```bash
gsutil -m rsync -r <entry>/assets/cache \
  gs://skearnes-logbook/<entry>/assets/cache
```

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

[pre-commit]: https://pre-commit.com

[cli2]: https://github.com/DavidAnson/markdownlint-cli2
