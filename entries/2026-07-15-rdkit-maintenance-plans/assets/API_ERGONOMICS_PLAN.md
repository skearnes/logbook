# RDKit Python API Ergonomics: Findings & Improvement Plan

Status: draft for discussion. Scope: the **Python** API, especially for newer/casual users.

---

## TL;DR — honest recommendation

The pain is real and well-founded, but the highest-leverage fixes are **additive and
documentary**, not a rewrite. RDKit's API is overwhelmingly defined in C++ and exposed
verbatim through Boost.Python; the CamelCase names, getter/setter methods, and `rd*`
module split are a direct, decades-deep consequence of that. A "make it Pythonic"
campaign that renames the binding surface or collapses modules would break an enormous
installed base for cosmetic gain and would not survive review.

What *will* move the needle, in priority order:

1. **A curated "front door" + a single canonical-usage guide.** One documented import,
   one recommended way to do each common task, with the alternatives explicitly labeled
   legacy/advanced. Most of the user's confusion is "which of these three do I use?" —
   that is answerable in docs and a thin namespace, today, with zero breakage.
2. **Fix the worst *correctness* footguns**, not the cosmetics: the
   `MolFromSmiles` → `None` silent-failure model, and the property-cache/`UpdatePropertyCache`
   trap behind the hydrogen questions.
3. **Pay down the hydrogen-model confusion with documentation + unambiguous helper
   accessors**, since the underlying representation can't change.
4. **Finish what's already started**: consolidate deprecations, and improve the *new*
   type-stub pipeline's coverage (stubs already ship — see below).

A full snake_case / property-based Pythonic layer (`mol.num_atoms`) is a legitimate Tier-3
idea, but it needs an RFC and a maintainer champion before any code. Don't start there.

This document is grounded in a code read; file:line references are included so each claim
is checkable.

---

## Constraints any plan must respect

These are not opinions; they're load-bearing facts about the project.

- **Boost.Python, not pybind11.** ~125 `Wrap/*.cpp` files, 30+ `BOOST_PYTHON_MODULE`
  declarations (e.g. `Code/GraphMol/Wrap/rdchem.cpp`, `Code/GraphMol/Wrap/rdmolfiles.cpp`).
  No migration in flight. The Python names are the C++ names; changing them means touching
  bindings and breaking callers.
- **Backward compatibility is treated as sacred.** Multi-release deprecate→warn→remove
  cycle, documented in `ReleaseNotes.md` ("Deprecated code (to be removed…)" sections at
  lines 25/429/753) and `Docs/Book/BackwardsIncompatibleChanges.md`. Year.Month.Revision
  versioning (`Code/RDGeneral/versions.cpp`), no semver "major bump" escape hatch.
- **Big API changes go through GitHub Discussions first** (`Docs/Book/GettingStartedWithContributing.md`).
  Code-first proposals for cross-cutting API change will stall.
- **Additive Python-only convenience layers are an established, blessed pattern.**
  `PandasTools.py`, the IPython/Jupyter auto-integration (`rdkit/__init__.py`,
  `Chem/Draw/IPythonConsole.py`), `PropertyMol.py`, and the pure-Python convenience
  wrappers in `Chem/__init__.py` (`CanonSmiles`, `QuickSmartsMatch`, `SupplierFromFilename`)
  all show the project accepts thin Python sugar over the C++ core. **This is the lane our
  high-value work lives in.**
- **Type stubs already exist and ship.** `rdkit-stubs/CMakeLists.txt` plus the
  `Scripts/gen_rdkit_stubs/` generator (Boost.Python introspection + a `patch/` directory
  of hand-maintained `.pyi.diff` corrections) and `Scripts/patch_rdkit_docstrings/`
  (libclang-based docstring/param patching). This is recent. The opportunity is
  **coverage and correctness of the generated stubs**, not creating them from scratch.

---

## The findings (grounded)

### A. "rd\* module vs Pythonic module" duality

`Chem/__init__.py` does `from rdkit.Chem.rdchem import *`, `from rdkit.Chem.rdmolfiles
import *`, `from rdkit.Chem.rdmolops import *`, etc. So `Chem.MolFromSmiles` and
`Chem.rdmolfiles.MolFromSmiles` are the *same* object reached two ways. `AllChem.py`
re-exports a still-larger union (`Chem/AllChem.py:20-36`). Net effect: a name like
`GetMorganGenerator` is reachable as `Chem.AllChem.GetMorganGenerator`,
`Chem.rdFingerprintGenerator.GetMorganGenerator`, and via `from rdkit.Chem import *`.
Users cargo-cult `from rdkit.Chem import AllChem` because it's the only import that
"has everything," which then teaches them nothing about where things live.

The genuinely-duplicated (not just re-exported) cases:

- **Descriptors** exist as friendly wrappers (`Descriptors.MolWt`, `Crippen.MolLogP`,
  `Lipinski.NumHDonors`) over underscore-prefixed C++ originals
  (`rdMolDescriptors._CalcMolWt`, `rdMolDescriptors.CalcNumHBD`). Same value, ≥3 spellings,
  sometimes different names for the same quantity (`Descriptors.NumHDonors` vs
  `Lipinski.NumHDonors` vs `rdMolDescriptors.CalcNumHBD`). See `Chem/Descriptors.py:75`,
  `Chem/Lipinski.py:49-56`, `Chem/Crippen.py:71`, `Chem/Descriptors3D.py`.
- **Fingerprints** have a genuine old/new split: the deprecated function API
  (`GetMorganFingerprint*` in `Code/GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp:359-554`,
  marked `[[deprecated("please use MorganGenerator")]]` and emitting
  `RDLog::deprecationWarning`) vs the current generator API
  (`Code/GraphMol/Fingerprints/Wrap/FingerprintGeneratorWrapper.cpp`,
  `GetMorganGenerator(...).GetFingerprint(mol)`). Both are still reachable via `AllChem`,
  with nothing in the namespace signaling which is current.

### B. The three hydrogen "models"

These are three genuinely different things wearing overlapping words:

| What users mean | RDKit representation | Accessor | Set by |
|---|---|---|---|
| "valence" Hs (not drawn, inferred) | implicit valence count, lazily computed; `-1` = uncomputed | `Atom.GetNumImplicitHs()` | `UpdatePropertyCache()` |
| Hs pinned as a *count* on a heavy atom | `d_numExplicitHs` (uint8), sets `noImplicit=true` | `Atom.GetNumExplicitHs()` | bracket SMILES `[CH3]`, SMARTS H-queries |
| Hs as real graph **nodes** | actual H atoms in the molecule | `mol.GetNumAtoms()` after `AddHs` | `Chem.AddHs(mol)` |

Refs: `Code/GraphMol/Atom.h:185-224,400-427`, `Code/GraphMol/Atom.cpp:298-337`,
`Code/GraphMol/AddHs.cpp`, `Code/GraphMol/MolOps.h:181-350`.

The confusion is **terminological collision**, not (only) representational complexity:

- "**explicit**" means three unrelated things: explicit-H *count* (`GetNumExplicitHs`),
  explicit *valence* (`GetExplicitValence`, = bond orders + that count — *not* an H count),
  and explicit H *atoms* in the graph (post-`AddHs`).
- `GetTotalNumHs(includeNeighbors)` is the one accessor that unifies the count view, but
  it's not the obvious thing a newcomer reaches for.
- **Footgun:** `GetNumImplicitHs()` returns 0/garbage until `UpdatePropertyCache()` runs,
  because `d_implicitValence` starts at `-1`. This silently bites people who build mols
  without sanitizing (`Atom.cpp:246,303-307`).
- **Surprise:** `Chem.MolFromSmiles('[CH3]C')` yields 2 atoms, but atom 0's
  `GetNumExplicitHs()` is 3 — the explicit-H *count* survives even though no H atoms are in
  the graph and the default `removeHs=True` was applied.

Existing docs barely cover this: `Docs/Book/GettingStartedInPython.rst:552-570` only shows
`AddHs`/`RemoveHs` (the count↔graph transition), and never disentangles the vocabulary.

### C. "Two or three ways to do the same thing"

- **Parsing & error model is inconsistent.** `MolFromSmiles` swallows *all* exceptions and
  returns `None` (`rdmolfiles.cpp:75-79`); `MolFromMolBlock` catches, logs, returns `None`
  (`:148-157`); `MolFromMolFile` *raises* `IOError` on a bad file (`:164-167`). So the
  "check for None" habit users are taught is wrong for the file-based readers, and the
  silent-None readers give no reason for the failure. Plus a v1/v2 `SmilesParserParams`
  split (`Code/GraphMol/SmilesParse/SmilesParse.h`) and a `sanitize` flag bolted onto every
  one of ~15 `MolFrom*` functions.
- **Drawing** has ≥4 entry points (`Chem/Draw/__init__.py` `MolToImage`,
  `rdMolDraw2D`, `MolsToGridImage`, `IPythonConsole`) with no documented hierarchy.
- **Un-Pythonic surface** (all inherent to the C++ binding): CamelCase methods,
  getter/setter instead of properties, integer `flavor` flags
  (`MolFromSequence`, `MolFromPDBFile`), in-place `RWMol` editing.

---

## Recommendations

### Tier 1 — high leverage, low risk, do these first

**1. Write the canonical "one true way" guide.**
A single doc page (and a notebook) that, for each common task — parse, sanitize,
descriptors, fingerprints, draw, substructure, conformers — gives **one recommended
snippet** and a short "legacy / advanced alternatives" footnote. This directly answers
"which of these do I use?" without any code change. Pair each recommendation with the
import that should appear in user code (e.g. `from rdkit import Chem` +
`from rdkit.Chem import Draw`), and actively discourage reflexive `from rdkit.Chem import
AllChem` for new code.

**2. Ship a curated front-door namespace (additive, non-breaking).**
A new pure-Python module — e.g. `rdkit.Chem.recommended` (name TBD) — that imports and
re-exports *only* the current, recommended spelling of each common operation, with real
docstrings and type hints. It changes nothing about the existing namespaces; it's a
documented, narrow surface a newcomer can `from rdkit.Chem.recommended import *` and trust.
This is the same pattern as `AllChem`, but curated *down* instead of unioned *up*.

**3. Authoritative hydrogen-model documentation + a decision table.**
Adopt and document three unambiguous nicknames — **valence-Hs** (implicit),
**count-Hs** (`numExplicitHs`), **graph-Hs** (after `AddHs`) — and ship the table from
section B plus a "which do I want?" flowchart. Explicitly call out the `UpdatePropertyCache`
prerequisite and the `[CH3]` count-survival surprise. Cross-link from every H-related
docstring. **Documentation is the fix here**; the representation can't change.

**4. Make `MolFromSmiles` failure debuggable without breaking the `None` contract.**
Keep returning `None` (changing that *is* a breaking change), but route the swallowed
exception through the warning log so users get a *reason*
(`rdmolfiles.cpp:75-79` currently has a bare `catch(...)` that drops it). Optionally add an
opt-in `MolFromSmiles(..., strict=True)`-style raising path, or a documented
`SmilesParserParams` knob, so users who *want* exceptions can have them. Low risk, high
daily-pain reduction.

**5. Consolidate and complete deprecation signaling.**
The deprecation machinery is inconsistent: the old Morgan API warns
(`rdMolDescriptors.cpp`), `MCS.py` and `FastSDMolSupplier.py` warn, but many legacy spellings
are silent. Make a single tracked list (in `ReleaseNotes.md`'s deprecation section) of
"legacy spelling → recommended spelling," and ensure each legacy entry emits a
`DeprecationWarning` pointing at the replacement. This is the only thing that will *shrink*
the "three ways" problem over time instead of just documenting it.

### Tier 2 — medium effort, clear payoff, needs a little buy-in

**6. Improve the generated type stubs' coverage and accuracy.**
Build on `Scripts/gen_rdkit_stubs/` + the `patch/*.pyi.diff` mechanism. Prioritize the
high-traffic surfaces (`rdchem.Mol`/`Atom`/`Bond`, `rdmolfiles`, `rdmolops`,
`rdMolDescriptors`, `rdFingerprintGenerator`). Good stubs are what make the *existing*
CamelCase API tolerable in modern IDEs — autocomplete and signatures recover most of the
discoverability that the binding style costs. This is the pragmatic substitute for
"make it Pythonic."

**7. Unify the descriptor entry points behind one documented surface.**
Don't delete `Crippen`/`Lipinski`/`Descriptors3D` (callers depend on them), but make
`rdkit.Chem.Descriptors` the documented single source, ensure every descriptor is reachable
there under one canonical name, and add a programmatic registry/listing
(`Descriptors.descList` exists — document and lean on it) so users discover descriptors
instead of guessing module locations.

**8. Document the parser/reader error model explicitly**, and consider aligning the
file-readers' default to the same "warn + None" behavior as the string readers (behind a
deprecation cycle), so the mental model is uniform.

### Tier 3 — ambitious, RFC-gated, do not start with code

**9. A Pythonic accessor layer (`mol.num_atoms`, snake_case, properties).**
Technically feasible as an *additive* pure-Python wrapper/mixin without touching bindings,
but it doubles the surface area and risks becoming a fourth "way to do things" unless the
maintainers commit to it as *the* direction. Requires a GitHub Discussion / RFC and a
maintainer champion. Recommendation: **propose it as a design RFC, not a PR**, and gate it
on Tier 1–2 landing first so it builds on the curated namespace rather than competing with it.

**10. Long-horizon: reconsider `AllChem`'s role.** Once the curated front door exists and
docs stop recommending `AllChem` for new code, `AllChem` can be documented as "legacy
convenience union" — no removal, just demotion.

---

## Suggested sequencing

1. Open a **GitHub Discussion** summarizing this (the three pain points + the Tier-1
   proposals). Get maintainer alignment on the curated-namespace idea and the H terminology
   *before* writing the docs, so the vocabulary is blessed.
2. Land Tier 1 docs (#1, #3) — pure docs, mergeable independently, immediate user value.
3. Land #4 (warn-on-parse-failure) and #5 (deprecation consolidation) as small, focused PRs.
4. Prototype the curated namespace (#2) once the recommended-spelling list from #5 exists —
   they share the same source of truth.
5. Spin up stub-coverage work (#6) in parallel; it's independent of the rest.
6. File the Pythonic-layer RFC (#9) separately and let it bake.

## What I deliberately am **not** recommending

- Renaming the Boost.Python surface to snake_case. Breaks everyone, for cosmetics.
- Collapsing `rd*` modules or removing `AllChem`. Demote in docs; don't delete.
- Changing `MolFromSmiles` to raise by default, or changing the hydrogen representation.
  Both are large breaking changes whose cost dwarfs the benefit; address via docs + opt-ins.

---

*Findings gathered by reading the source tree; references above are checkable at the cited
file:line. Open question for maintainers: is there appetite for a curated front-door
namespace, or is the preference to keep all curation in documentation only?*
