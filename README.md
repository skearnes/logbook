# Logbook

A running work logbook. Each entry is a dated Markdown file under
[`entries/`](entries/), named `YYYY-MM-DD-topic.md` so the topic is visible at
a glance.

## Conventions

- Start a new entry by copying [`TEMPLATE.md`](TEMPLATE.md) into `entries/`.
- Name each file `entries/YYYY-MM-DD-topic.md` (e.g.
  `entries/2026-03-14-pi-day.md`).
- Begin each entry with an H1 title covering the date and topic.
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
