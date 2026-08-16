# Compile-Time Reduction Plan

- **Date:** 2026-07-15
- **Author:** Steven Kearnes
- **Status:** draft (working planning doc)
- **Tags:** rdkit, compile-time, build, boost-graph, profiling, planning
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

Working notes — **untracked scratch doc, do not commit.** A plan for cutting RDKit
C++ build times, ranked by value ÷ effort. Grounded in profiling of this tree
(macOS / Apple clang 21, C++20).

## Measured baseline (why this matters)

Profiled a representative core TU (`GraphMol/ROMol.cpp`) with `clang -ftime-trace`
and `-E`:

- **224,000 preprocessed lines** for a file whose own code is trivial.
- **1,590 header-parse events** per TU.
- **~69% of compile time is frontend** (parsing + template instantiation):
  Frontend 1.58s vs Backend 0.63s of a 2.3s total.
  `ParseClass` 0.47s, `InstantiateFunction` 0.36s, `InstantiateClass` 0.27s.

Root cause: pervasive header bloat, multiplied across ~1,250 TUs. The dominant
single contributor is the Boost Graph Library, pulled into **335 TUs** because
`ROMol.h:31` includes `<boost/graph/adjacency_list.hpp>` directly (and 335 `.cpp`
files include `RDKitBase.h` → `ROMol.h`).

Current infrastructure: **no unity builds, no ccache, and PCH only via the opt-in
draft PR #9236.**

## Context: PCH is already in flight (PR #9236, greglandrum)

- Adds an opt-in `Code/stdafx.h` PCH bundling STL + `dynamic_bitset` +
  `boost/algorithm/string` + `RDKitBase.h` (so it does capture BGL transitively).
- Enabled for Windows + Mac (non-SWIG) CI; **disabled for g++/Linux because it
  regresses there.**
- Does NOT cover Eigen, `boost/format`, or `boost/json`.

**Why GCC regresses:** GCC PCH is a whole-GC-heap memory-image snapshot with no
lazy deserialization and a one-PCH-per-TU / must-come-first restriction. Every TU
faults in the entire (hundreds-of-MB, esp. with `-g`) image regardless of how much
it uses; under high `-j` this causes page-cache/memory contention that serializes
parallel builds. Clang/MSVC use lazy AST deserialization and avoid this.

**Implication for this plan:** the highest-value *code* changes are the ones that
shrink real parse work for **every** compiler and for incremental builds and
external consumers — i.e. they help exactly where #9236 can't (GCC, non-PCH).
Don't duplicate #9236; stack on top of it.

---

## Ranked changes

### 1. Decouple `boost/graph` from `ROMol.h` — top structural change

**Value:** high · **Effort:** high · **Risk:** medium (public API)

`ROMol.h` includes the full BGL `adjacency_list` header, so all 335 `RDKitBase`
consumers re-parse BGL — the biggest slice of the frontend cost. PCH masks this
only for opt-in Clang/MSVC; GCC, incremental rebuilds (any `ROMol.h` edit), and
downstream library users still pay it every time.

- **Approach:** pimpl / forward-declare `MolGraph`; move BGL-using inline methods
  into `ROMol.cpp`.
- **Constraint:** `vertex_descriptor` / `edge_descriptor` / `*_ITER` are typedef'd
  in the public API (`ROMol.h:328-334`), so *full* removal is an API break. A
  partial version (move inline bodies out, keep typedefs) is feasible and still
  removes most parse cost.
- **Prereq:** maintainer discussion before touching the central header.

### 2. Unity / jumbo builds (`CMAKE_UNITY_BUILD`), per-directory

**Value:** high · **Effort:** medium · **Risk:** medium (ODR)

Parses shared headers once per batch instead of per file (typ. 2–3× on heavy
dirs). Crucially **helps GCC/Linux**, where #9236 regresses — the GCC-friendly
analog of PCH.

- **Risk:** RDKit has many anonymous-namespace / `static` helpers; concatenating
  TUs surfaces ODR collisions, symbol clashes, and macro leakage.
- **Rollout:** enable dir-by-dir, fix collisions as they surface; start with a
  heavy leaf dir (e.g. `Descriptors`, `MolDraw2D`).

### 3. Trim heavy headers the PCH doesn't cover — easy starter

**Value:** medium · **Effort:** low · **Risk:** low

Mechanical, uncontroversial, helps all builds incl. non-PCH:

- `Eigen/Dense` in a public header (1) → forward-declare / move to `.cpp`.
- `boost/format` in 6 headers → impl-only; move to `.cpp`.
- Keep `boost/json` (esp. `src.hpp`) confined to `.cpp` (mostly already is).

Validate each by diffing preprocessed line count (`clang -E | wc -l`) of a
representative consumer before/after.

### 4. `extern template` for hot templated APIs

**Value:** medium · **Effort:** medium · **Risk:** low

Currently **zero** `extern template` in the tree; `InstantiateFunction` +
`InstantiateClass` was ~0.63s/TU. Declare `extern template` in headers + one
explicit instantiation in a `.cpp` for the most-included templated code
(fingerprint generators, substruct matchers) so every TU stops re-instantiating.

### 5. Split the `RDKitBase.h` god-header

**Value:** medium · **Effort:** high · **Risk:** low-medium

Fans out to 335 TUs, dragging in Atom/Bond/RWMol/PeriodicTable/etc. Finer-grained
includes shrink incremental rebuilds. Large mechanical effort, lowest value ÷
effort of the list — do last, if at all.

### Non-code (config) — worth raising separately

- **ccache** (`CMAKE_CXX_COMPILER_LAUNCHER`) in CI/dev — not a code change, but
  eliminates recompiles; should be on.

---

## Suggested sequencing

1. **#3** first — small, low-risk PR; immediate win on all platforms; proves the
   measure-by-preprocessed-lines workflow.
2. **#2** next — biggest GCC/Linux win, complements #9236; roll out per-dir.
3. **#1** — the durable structural fix; needs maintainer buy-in given the API
   typedefs.
4. **#4 / #5** — opportunistic.

## How to measure any change

- Per-TU parse cost: `clang -E <flags> file.cpp | wc -l` before/after.
- Real hotspots: `clang -ftime-trace` → inspect `Total Frontend` / `ParseClass` /
  `InstantiateFunction`.
- Whole-build: clean build wall-clock at fixed `-j`, ideally on both Clang and GCC
  (the platform split is the whole point).
