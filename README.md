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
