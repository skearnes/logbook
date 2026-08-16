# Reducing RDKit's Dependence on Boost

- **Date:** 2026-07-15
- **Author:** Steven Kearnes
- **Status:** draft (working planning doc)
- **Tags:** rdkit, boost, dependencies, build, planning
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

Working planning doc. Grounded in a header-by-header inventory of `Code/` (1,271 C++
files; 366 include some Boost header). Counts below are file counts under `Code/`.

---

## TL;DR — honest recommendation

**Pursue *less* Boost, not *zero* Boost.** A literal no-Boost future is a person-years
trap whose payoff (drop a build dependency) doesn't reach the things RDKit users feel.
But there is a genuinely valuable, *scoped* target hiding inside the question:

> **Make the core C++ library depend only on *header-only* Boost, then evaluate moving
> the Python bindings off Boost.Python.** That, not "rip out Boost," is the prize.

Why this framing is the whole game: the only Boost components RDKit actually **compiles
and links** are `python`, `numpy`, `serialization`, and `iostreams` (CMakeLists.txt:337,
448, 458; `find_package(... COMPONENTS serialization iostreams ...)`). Everything else —
`graph`, `dynamic_bitset`, `flyweight`, `algorithm/string`, `format`, `lexical_cast`,
`optional` — is header-only `Boost::boost`. **Removing a header-only component buys you
nothing at link or packaging time** and only a little at compile time; removing a
*compiled* component shrinks the actual build/dependency footprint. So the cost/benefit is
dominated by four libraries, not by the long tail of utility headers.

The verdict, by tier:

- **Worth doing now (high value ÷ effort):** kill the two *compiled* non-Python libs —
  **`serialization`** and **`iostreams`** — from the core. This is the single highest-ROI
  move: it lets the C++ core link against header-only Boost, easing every downstream build.
- **Worth doing as hygiene (low risk, modest payoff):** the mechanical std-replaceable
  utilities (`lexical_cast`, `format`, `algorithm/string`, `random`, `optional`/`variant`,
  `noncopyable`). Co-benefit: `lexical_cast` and `format` are heavy headers, so this trims
  compile time too. But none individually moves a metric — do it as cleanup, not a project.
- **Strategic, RFC-gated, big upside:** evaluate **nanobind** to replace Boost.Python.
  This is the only path to dropping the `python`/`numpy` compiled deps, and it carries
  large *independent* wins (smaller binaries, real type stubs, faster builds). Pilot one
  module; don't frame it as "remove Boost."
- **Leave alone (stay-put verdict):** `dynamic_bitset` (84 files, no std equivalent,
  perf-critical, header-only → near-zero packaging cost) and `flyweight` (14 files,
  interning cache, no std equivalent, header-only). Keep both. **Trigger to revisit:** only
  if one blocks a specific refactor, or if a clean drop-in replacement is independently
  written and benchmarked.
- **Already covered elsewhere:** **`graph`/BGL** is the compile-time prize, but it's a
  *header-decoupling* problem already tracked in `COMPILE_TIME_PLAN.md` /
  `BGL_DECOUPLING_FINDINGS.md`, not a "remove Boost" problem. Don't duplicate it here.

---

## The mental model: Boost-in-RDKit is three unrelated problems

Conflating these is the main way this effort goes wrong. They have different owners,
risks, and payoffs:

1. **Boost.Python (the bindings).** ~125 `Wrap/*.cpp`, ~38.5K LOC, 12 `BOOST_PYTHON_MODULE`
   entry points, zero internal abstraction — raw Boost.Python throughout, plus holder
   types (`boost::shared_ptr`), custom converters, exception translators, pickle suites,
   GIL helpers (`Code/RDBoost/Wrap.h`). Replacement = **nanobind/pybind11**. Huge, partly
   mechanical, partly subtle. Compiled dep. *Strategic track of its own.*
2. **The molecular graph (BGL).** `MolGraph = boost::adjacency_list<...>` (ROMol.h:55-57);
   the blocker is public iterator typedefs (`vertex_descriptor`, `EDGE_ITER`, … ROMol.h:328-338)
   leaking BGL types into the iteration API. Header-only, but in 335 TUs transitively →
   the dominant *compile-time* cost. *This is a header-decoupling/core-API problem; tracked
   separately.*
3. **Everything else (utilities + serialization + iostreams).** Independently removable,
   mostly mechanical. *This document's main actionable scope.*

---

## Grounded inventory

`L` = header-only (link cost: none) · `C` = **compiled Boost library** (link/packaging cost).

| Boost component | Files | L/C | Replacement | Difficulty | Notes |
|---|---:|:--:|---|---|---|
| `python` + `numpy` | ~125 (Wrap) | **C** | nanobind / pybind11 | very high | The bindings. Own track. Drops 2 compiled deps. |
| `dynamic_bitset` | 84 | L | none (custom or keep) | hard | Core to fingerprints/SSS. **Keep.** No std equivalent. |
| `algorithm/string` | 72 | L | `RDKit::StringUtils` (std) | medium | split/trim/join/case. Centralize in one helper. |
| `lexical_cast` | 39 | L | `std::to_string`/`from_chars`/`sto*` | easy | Heavy header → compile co-benefit. |
| `format` | 34 | L | `std::format` (C++20) | easy–med | Heavy header. One note: a comment in FindStereo.cpp records boost::format was "WAY slower" — std::format is fine. |
| `smart_ptr`/`shared_ptr` | 22/21 | L | `std::shared_ptr` etc. | medium | **Caveat:** ~7 files use `boost::shared_ptr` as a Boost.Python *holder type* — those are coupled to track #1, do them there. |
| `serialization` + `archive` | 15 + 13 | **C** | hand-rolled save/load | med–high | ROMol = thin shim over `MolPickler` (easy). TautomerQuery/ScaffoldNetwork/FilterCatalog/Enumerate deeply embed macros (hard). Gated by `RDK_USE_BOOST_SERIALIZATION`. |
| `flyweight` | 14 | L | none (small static cache) | medium | String-keyed interning (SMARTS matchers, param tables). **Keep.** |
| `random` | 13 | L | `<random>` (std::mt19937) | easy | 1:1. A bench file already notes preferring std. |
| `graph` (direct) | 12 | L | index-based graph / pimpl | very high | Only 12 *direct* includers, but transitively in 335 TUs. Compile-time prize. See BGL docs. |
| `iostreams` | 8 | **C** | zlib directly (gzip streambuf) | medium | gzip SD/PNG streams (`RDStreams/streams.h`, `CompressedSDMolSupplier`, `PNGParser`). zlib already a dep. |
| `noncopyable` | 5 | L | `= delete` | trivial | |
| `property_map` | 4 | L | (part of BGL) | — | Falls out with BGL work. |
| `multiprecision` | 2 | L | — | — | Check if actually used; may be droppable include. |
| `optional`/`variant`/`static_assert`/`cstdint`/`foreach` | 1 each | L | std / language | trivial | Near-instant wins. |

Already gone: `boost/regex`, `boost/filesystem`, `boost/any`, `boost/tuple`, `BOOST_FOREACH`
(0 files) — RDKit has migrated these. Two dead includes spotted by the BGL pass
(`kruskal_min_spanning_tree.hpp`, `johnson_all_pairs_shortest.hpp` in MolOps.cpp) and
`boost/math` includes that are never called — drop on sight.

---

## The strategic ladder (what each rung actually buys)

Frame the effort as rungs, each a shippable milestone, each with a clear stop point:

- **Rung 0 — Hygiene.** Remove dead includes; convert the 1-file/trivial components
  (`optional`/`variant`/`static_assert`/`cstdint`/`foreach`/`noncopyable`). **Add a header
  self-containment check** (compile each public header alone in a one-line TU) — it catches
  the not-including-your-own-deps bugs that every umbrella-header removal keeps exposing
  (the iostreams removal hit two: `RDValue-taggedunion.h`, `RDLog.h`), and it's the cheap
  guard that makes every later rung safer. Low noise, no add-everything churn (unlike a
  tree-wide IWYU sweep — see Follow-ups). Days. Pure win.
- **Rung 1 — Core links only header-only Boost.** Remove the two compiled non-Python libs:
  **`iostreams`** (→ zlib streambuf) and **`serialization`** (→ hand-rolled save/load).
  *Payoff:* the core C++ library (`libRDKitGraphMol` et al., independent of the Python
  wrappers) no longer needs any *compiled* Boost — a real, measurable packaging/build
  simplification for C++ consumers and conda-forge. **This is the recommended primary goal.**
- **Rung 2 — Project links zero compiled Boost.** Add: move bindings to nanobind, dropping
  `python`/`numpy`. *Payoff:* RDKit's only remaining Boost is header-only; smaller wheels,
  real type stubs, faster wrapper builds. Big, RFC-gated, high independent value.
- **Rung 3 — Zero Boost (header-only too).** Add: replace `dynamic_bitset`, `flyweight`,
  BGL, and the string/format/lexical_cast long tail. *Payoff:* drop a build dep entirely.
  **Not recommended.** Person-years; high regression risk (fingerprint perf, graph hot path);
  the marginal benefit over Rung 2 is "no `find_package(Boost)` at all," which is cosmetic
  once nothing compiled is linked.

**Stop at Rung 1 unless the nanobind evaluation (Rung 2) pays for itself on its own merits.**
Rung 3 is a trap.

---

## Pros / cons: less-Boost vs no-Boost

**Why reduce Boost at all (the real benefits):**

- *Build & packaging:* fewer compiled Boost libs to find/build/link/version-pin. Boost
  version skew is a recurring downstream pain; every compiled component removed is one
  fewer ABI-coupled dep. (Header-only removals don't help here — hence the L/C lens.)
- *Compile time:* a few Boost headers are genuinely heavy (`graph`, `python`,
  `serialization`/`archive`, `lexical_cast`, `format`). Removing those *headers from hot
  TUs* trims parse cost — but this overlaps with the unity/PCH/BGL work and is mostly about
  BGL, not the utilities.
- *Maintainability:* `std::` idioms are known to every contributor; Boost idioms
  (converter registries, flyweight, archive macros) are tribal knowledge.
- *Binary size:* nanobind specifically produces dramatically smaller extensions than
  Boost.Python — the single biggest size lever.

**Why NOT chase zero-Boost (the costs):**

- *No std equivalent for the load-bearing pieces:* `dynamic_bitset` and `flyweight` have no
  standard replacement; `std::vector<bool>`/`std::bitset` are not substitutes. Replacing
  them means *owning* a perf-critical data structure forever.
- *Regression surface:* fingerprints, substructure search, and the graph hot path are
  exactly where RDKit's correctness and speed reputation live. Boost here is load-tested
  over ~20 years.
- *Format/ABI compatibility:* serialization and pickle touch on-disk/wire formats
  (`MolPickler` is version-locked, MolPickler.h:69-71). Changing the framing layer is fine;
  changing formats is a compatibility event.
- *Opportunity cost:* every week here is a week not on chemistry features. Boost is not on
  fire; this is optimization, not rescue.
- *`std::regex` lesson:* the standard isn't always better — a reminder to replace per-case
  on merit, not dogmatically.

---

## Recommended roadmap

1. **Rung 0 hygiene** (independent small PRs): dead includes + the trivial 1-file
   components + `noncopyable`. Uncontroversial, mergeable immediately.
2. **`random` → `<random>`, `lexical_cast` → std, `format` → `std::format`** (separate PRs,
   per-directory). Mechanical; compile co-benefit on the latter two. Low risk.
3. **`algorithm/string`:** introduce one `RDKit::StringUtils` header (std-backed
   `split`/`trim`/`join`/`to_lower`/`replace`), then migrate callers directory-by-directory.
   Centralizing first means the 72 call sites converge on one tested helper.
4. **Rung 1 — `iostreams`:** write a small zlib-backed gzip streambuf in `RDStreams/`,
   swap `boost::iostreams::filtering_stream` users over, drop the `iostreams` component.
5. **Rung 1 — `serialization`:** do the easy ROMol shim first (it already delegates to
   `MolPickler`); then decide whether the deeply-embedded optional modules
   (TautomerQuery/ScaffoldNetwork/FilterCatalog/Enumerate) are worth de-macro-ing or stay
   behind `RDK_USE_BOOST_SERIALIZATION`. Dropping the component requires all of them.
6. **Coordinate with BGL work** (separate plan) — don't re-open it here.
7. **Rung 2 — nanobind evaluation:** see below. RFC + single-module pilot before any commitment.

Each numbered item is independently shippable and independently reversible. Resist
bundling them into one mega-PR — that's how this stalls in review.

## The Boost.Python question (the strategic centerpiece)

This deserves its own RFC, not a paragraph here, but the shape:

- **Pick nanobind over pybind11.** For a from-scratch port the migration cost is similar,
  but nanobind gives much smaller binaries, faster compile, native dual `std`/`boost`
  holder support, and a cleaner type-stub story — which dovetails with the
  `Scripts/gen_rdkit_stubs/` work. pybind11's only edge is maturity/breadth, which matters
  less for a controlled, module-by-module port.
- **The hard seam is `Code/RDBoost/`** — `Wrap.h`'s converters, the exception-translator
  registry, the pickle suites, and the `boost::shared_ptr` holders. Port *that layer* first,
  behind a thin RDKit-internal shim, so individual `Wrap/*.cpp` modules can migrate against
  a stable internal API instead of raw nanobind. (The absence of such a shim today is why
  the bindings are 38.5K LOC of raw Boost.Python — building it is the enabling step.)
- **Pilot one self-contained module** (e.g. `rdGeometry` or `DataStructs`) end-to-end,
  including pickle round-trip and exception-type compatibility tests, and *measure* binary
  size + build time before proposing the full port.
- **Watch out for user-visible behavior:** exception types, pickle compatibility, and
  numpy buffer semantics are where subtle breakage hides. This is a multi-release effort
  with a deprecation/compat story, not a flag flip.

## Follow-ups (out of scope for the iostreams removal)

- **Deprecate, then delete, `CompressedSDMolSupplier` + the `RDK_BUILD_COMPRESSED_SUPPLIERS`
  option.** With gzip/bzip2 now provided as a core capability (zlib required; libbz2 via
  `RDK_USE_BZIP2`), the regular `SDMolSupplier` / `GeneralMolSupplier::getSupplier` path
  already reads `.gz` natively, so the dedicated `ForwardSDMolSupplier`-based
  `CompressedSDMolSupplier` (built only when `RDK_BUILD_COMPRESSED_SUPPLIERS=ON`, default
  OFF) is largely redundant. Plan: deprecate it (point users at the native `.gz` suppliers),
  then remove it and the option in a later release. The option is exported as
  `RDKit_HAS_COMPRESSED_SUPPLIERS` in `rdkit-config.cmake`, so removal needs a compat note.
  Kept out of the iostreams-removal change to keep that change's scope and back-compat story
  clean.

- **No tree-wide IWYU sweep.** The header self-containment check is now a Rung 0 item (it
  catches the not-including-your-own-deps bugs that umbrella-header removal keeps exposing).
  A full automated include-what-you-use pass, by contrast, is the wrong tool here: on this
  Boost/template/macro-heavy tree it is noisy, needs curated mapping files, churns hundreds
  of files, and its "add what you use" direction can *increase* includes (fighting the
  compile-time goal). The only IWYU use worth considering is *suggest-only, removal-focused*
  runs on the specific heavy headers being trimmed, to find *unused* includes worth dropping
  (compile-time positive), reviewed by hand — never a bulk apply.

## What I'm explicitly NOT recommending

- Removing `dynamic_bitset` or `flyweight` (no std equivalent; header-only; high risk, ~zero
  packaging payoff). Keep them.
- A single "de-Boost RDKit" mega-effort. Split by the three problems and by L/C.
- Chasing Rung 3 / literal zero-Boost. Stop at "zero *compiled* Boost."
- Touching pickle/serialization *formats* (vs the framing layer) — that's a compat event.

---

*Counts are checkable: `find Code -name '*.cpp' -o -name '*.h' | xargs grep -lE
'include.*<boost/<component>'`. Compiled-vs-header-only is from CMakeLists.txt:337/448/458
`find_package(Boost ... COMPONENTS ...)`. Open question for maintainers: is the goal here
build/packaging simplification (→ prioritize the compiled libs + nanobind) or compile time
(→ prioritize BGL, already in flight)? The answer reorders everything above.*
