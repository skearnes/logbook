# 2026-07-15 — RDKit maintenance planning docs

- **Date:** 2026-07-15
- **Author:** Steven Kearnes
- **Status:** final (archive)
- **Tags:** rdkit, maintenance, compile-time, boost, api-ergonomics, dead-code,
  planning
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

Archive of five planning/findings documents written while investigating RDKit
maintenance opportunities (build times, Boost footprint, API ergonomics, dead
code). They were working notes in a local `rdkit` checkout (branch
`boost-rung0-cleanup`) and were never meant to be committed there; this entry
preserves them.

## The documents

- [API_ERGONOMICS_PLAN.md](assets/API_ERGONOMICS_PLAN.md) — Python API
  ergonomics findings. Recommendation: additive and documentary fixes (a
  curated "front door" namespace and canonical-usage guide, fixing the
  `MolFromSmiles` → `None` silent-failure and property-cache footguns) rather
  than any rename-the-bindings campaign.
- [BOOST_REDUCTION_PLAN.md](assets/BOOST_REDUCTION_PLAN.md) — header-by-header
  inventory of Boost usage across `Code/` (366 of 1,271 files). Recommendation:
  pursue *less* Boost, not *zero* Boost — get the core library down to
  header-only Boost, then evaluate moving the Python bindings off Boost.Python;
  the four compiled components (`python`, `numpy`, `serialization`,
  `iostreams`) dominate the cost/benefit.
- [COMPILE_TIME_PLAN.md](assets/COMPILE_TIME_PLAN.md) — build-time reduction
  plan ranked by value ÷ effort, grounded in `-ftime-trace` profiling (~69% of
  compile time is frontend; a representative TU preprocesses to 224k lines).
  The Boost Graph Library reaching 335 TUs via `ROMol.h` is the dominant single
  contributor.
- [BGL_DECOUPLING_FINDINGS.md](assets/BGL_DECOUPLING_FINDINGS.md) — prototype
  results for removing `boost/graph` from `ROMol.h`: no measurable runtime
  regression from de-inlining the hot accessors and a real compile-time prize,
  but blocked on ROMol's public iteration accessors naming BGL iterator types.
  Outcome: handed to maintainers as the artifact for that decision; effort
  pivoted to unity builds.
- [DEAD_CODE_PLAN.md](assets/DEAD_CODE_PLAN.md) — tiered survey of deletable
  dead code (~330 lines plus 9 whole orphaned files). Recommendation: delete
  Tiers 1–2b now, hand-triage Tier 3, leave Tier 4 alone.

## Next steps

- Upstream the Tier 1–2b dead-code deletions and the doc-level API ergonomics
  fixes as RDKit PRs.
- Revisit the BGL decoupling if maintainers take up the iteration-machinery
  rewrite (coordinate with the PCH work in rdkit#9236).
