# Decoupling boost/graph from ROMol.h — prototype findings

Working notes — **untracked scratch doc, do not commit.** Prototype + benchmark
results for plan item #1 (see COMPILE_TIME_PLAN.md). Measured on macOS / Apple
clang 21, C++20, shared-library build (RDKIT_DYN_LINK), 8 cores.

## OUTCOME (2026-06-09)

**Item #1 deferred to maintainers; pivoted to item #2 (unity builds).** The
prototype disproved the runtime objection (nil cost) and quantified the prize,
but the feasibility test revealed the real blocker: ROMol's public iteration
accessors name BGL iterator types in their signatures, so removing the header
requires rewriting the iteration machinery to be index-based (see "BLOCKER"
below). That's a maintainer-led core change. Near-term effort moved to unity
builds, which capture much of the same compile win without touching ROMol's API
or the iteration hot path. This doc is the artifact to hand maintainers for #1.

## TL;DR

The feared blocker (runtime cost of de-inlining ROMol's hot accessors) **did not
materialize** — a controlled A/B shows **no measurable regression (<1%, mixed
sign)**. The compile prize is real: **~0.51s of frontend parse per TU × 335 TUs**.
This makes the decoupling worth doing — but it remains an ABI break that touches
~42 files and the library's central type, so it needs maintainer sign-off and
should be coordinated with the PCH PR (#9236).

## The prize (compile time)

- `boost/graph/adjacency_list.hpp` = **143,286 preprocessed lines** ≈ **64% of
  ROMol.h's parse cost** (ROMol.cpp expands to 224k lines total).
- Measured wall-clock to parse it once (`clang -fsyntax-only`, best of 3):
  **0.51s** (0.52s with the header vs 0.01s empty).
- `ROMol.h` is included by **335 TUs** (everything that includes `RDKitBase.h`).
- Projected: ~335 × 0.51s ≈ **~170 CPU-seconds** of pure parse, repeated on every
  clean build and on every incremental rebuild that touches `ROMol.h`.
- Caveat: PCH (#9236) already hides this for Clang/MSVC TUs that use the PCH. The
  decoupling's *unique* value is **GCC/Linux (PCH disabled there), incremental
  rebuilds, and downstream library consumers** — exactly PCH's gaps.

## The feared cost (runtime) — measured, and it's nil

Concern: removing the BGL include forces a pimpl pointer + de-inlining ROMol's
inline graph accessors (`getNumAtoms()`, `operator[]`, `atoms()`/`bonds()`), some
of which run in tight loops.

Key de-risking fact found during scoping: **the hottest per-atom accessors,
`getAtomWithIdx()` / `getBondWithIdx()`, are ALREADY out-of-line** (defined in
ROMol.cpp). RDKit already pays call overhead for index access. In a shared-lib
build, cross-dylib calls into ROMol are already un-inlinable anyway.

Prototype: de-inlined `getNumAtoms()` and both `operator[]` overloads into
ROMol.cpp (pessimistic — `getNumAtoms()` in loop conditions is the worst case for
de-inlining). Built two variants, ran Catch2 benchmarks back-to-back on a quiet
machine (6 runs each, 40 samples):

| benchmark              | baseline | de-inlined | delta |
|------------------------|---------:|-----------:|------:|
| Descriptors spiro (us) |  616.9   |  619.1     | +0.3% |
| MolOps::addHs (ms)     |   7.214  |   7.183    | -0.4% |
| MolOps::FindSSR (ms)   |   4.023  |   4.025    | +0.0% |
| MolOps::getMolFrags(ms)|  18.060  |  17.920    | -0.8% |

All deltas < 1% with mixed signs → **no measurable runtime cost.** (Baseline
run-to-run variance was itself ~1-2% once builds stopped; early "noise" of
444→713us was concurrent-build interference, not signal.)

## What the full implementation still requires (not yet done)

The prototype only de-inlined to measure runtime; it did NOT remove the include.
To actually drop `<boost/graph/adjacency_list.hpp>` from ROMol.h:

1. **pimpl the member** — `MolGraph d_graph;` (by-value, line 988) needs the
   complete type. Must become a pointer (`std::unique_ptr<MolGraph>`); ROMol's
   destructor moves to the .cpp. Mechanically small: only **4 files** touch
   `d_graph` directly, and `getTopology()` can keep returning `const MolGraph&`
   (dereferencing the pimpl), insulating the ~335 downstream consumers.
2. **Forward-declare MolGraph** — this Boost has **no `adjacency_list_fwd.hpp`**,
   so ROMol.h needs a hand-written 7-arg `adjacency_list` template decl (specify
   all template args explicitly to avoid default-arg redefinition conflicts).
3. **Move the nested typedefs out** — `vertex_descriptor`, `edge_descriptor`,
   `EDGE_ITER`, `OEDGE_ITER`, `VERTEX_ITER`, `ADJ_ITER` (ROMol.h:328-338) require
   the complete type and so can't stay in ROMol.h. They move to a new heavy
   `MolGraph.h`. **This breaks source-compat for ~42 files** using
   `ROMol::EDGE_ITER` etc. (9 use `vertex_descriptor`/`edge_descriptor`); those
   files must `#include "MolGraph.h"`. This is the main migration cost.
4. **De-inline ~15 accessor methods** (`atoms()`, `bonds()`, `atomBonds()`,
   `operator[]`, `getNumAtoms()`, …) into ROMol.cpp. (`getNumAtoms()` could
   instead stay inline & fast by reading a cached count, sidestepping even the
   loop-condition case — optional optimization.)
5. **ABI break**: member layout + inline→out-of-line. (RDKit breaks ABI between
   releases anyway.)

## BLOCKER found during implementation (supersedes the optimism above)

A forward-declared `MolGraph` is **not sufficient**. Verified by compile test:

- `std::unique_ptr<MolGraph>` member — OK with fwd-decl.
- `const MolGraph& getTopology() const;` — OK with fwd-decl.
- **`CXXAtomIterator<MolGraph, Atom*> atoms();` — FAILS.** Naming this return type
  instantiates the default template arg `typename MolGraph::vertex_iterator`,
  which requires the *complete* `adjacency_list`. Same for `bonds()`,
  `atomBonds()`, `atomNeighbors()`, `getVertices()`, `getEdges()` and the
  `*_ITER_PAIR` typedefs.

So ROMol's **public iteration API names BGL iterator types in its member-function
signatures** — the complete header is required just to *declare the class*. The
pimpl + de-inline path cannot remove the include. The only ways to actually
decouple:

1. **Rewrite CXXAtomIterator/CXXBondIterator to be index-based** (iterate `size_t`
   0..N instead of `MolGraph::vertex_iterator`). Then the accessor return types no
   longer name MolGraph's nested types, and ROMol.h can forward-declare MolGraph.
   This is a real rewrite of the iteration machinery (bonds/edges are the tricky
   part), with its own runtime-validation burden — but it's the only clean route
   to the compile win.
2. **Type-erased ranges** for atoms()/bonds() — API break **and** reintroduces
   runtime cost in the hottest loops (iteration), i.e. it resurrects exactly the
   regression we just disproved for scalar accessors. Worse option.
3. **Don't decouple** — accept the header stays; pursue unity builds (#2) instead.

Net: the "pimpl + de-inline + migrate 42 files" scope is **insufficient**; the
real project is "make ROMol iteration index-based," which is larger and riskier.

## Recommendation

**Proceed — but socialize first.** The runtime objection is empirically dead and
the compile prize is real, so this is no longer a "probably not worth it" change.
But it's an ABI break touching the core type and ~42 files, so:

1. Open a maintainer issue with this data (prize, the nil runtime result, the
   42-file typedef migration, ABI note), and coordinate with #9236 so the two
   compile-time efforts don't overlap/conflict.
2. If greenlit, implement in the order above; keep `getNumAtoms()` inline via a
   cached count to be conservative; land the typedef migration as its own
   mechanical commit.

## Reproduce

- Prize: `printf '#include <boost/graph/adjacency_list.hpp>\nint main(){}' > x.cpp; time clang++ -I<boost> -std=c++20 -fsyntax-only x.cpp`
- Runtime: de-inline getNumAtoms()/operator[] into ROMol.cpp, `ninja Code/Bench/bench`, `./bench "[molops],[descriptors]" --benchmark-samples 40` (build both variants, run back-to-back on a quiet machine).
