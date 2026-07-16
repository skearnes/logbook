# Dead Code Cleanup Plan

Survey of deletable dead code across the RDKit tree. Report-first: nothing here
has been deleted yet. Candidates were found by four parallel sweeps (commented-out
code, internal unused symbols, orphaned files, unused members/locals) and the
high-value clusters were independently re-verified by grep.

**Recommendation:** delete Tiers 1, 2, and 2b now (zero risk to public API,
bindings, or builds). Triage Tier 3 by hand. Leave Tier 4 alone.
Suggested execution: one commit per tier for easy review/revert.

Estimated removal: ~330 lines of C++/Python + 9 whole files.

---

## Tier 1 — Orphaned files, never built (safest; whole-file deletes)

Each appears in **zero** `CMakeLists.txt` and is never `#include`d.

| File | What it is |
|---|---|
| `Code/MinimalLib/testMinilib.cpp` | stale test driver (`main()`) |
| `Code/SimDivPickers/pickersCLI.cpp` | unbuilt CLI |
| `Code/GraphMol/bulktest.cpp` | 2002-era test driver |
| `Code/GraphMol/sanitTest.cpp` | 2002-era test driver |
| `Code/GraphMol/Substruct/cmd_match.cpp` | 2001-era CLI |
| `Code/GraphMol/Substruct/bench.cpp` | unbuilt benchmark |
| `Code/GraphMol/MolStandardize/testPCS.cpp` | 348-line unbuilt test |
| `Code/GraphMol/GraphMol.cpp` | body is just `namespace RDKit {}` (empty) |
| `Code/GraphMol/Descriptors/datas.cpp` | 120 lines, 100% commented-out reference data |

---

## Tier 2 — Internal (file-local) unused functions (~165 lines)

Anonymous-namespace / `static` functions with no call site anywhere. Deleting them
cannot break any public API or Python/Java binding. Each verified at exactly one
repo-wide occurrence (the definition).

### `Code/GraphMol/Descriptors/GETAWAY.cpp` (~77 lines)

- `retreiveMat` (:70)
- `retreiveVect` (:75)
- `round_to_n_digits` (:92)
- `IsClose` (:97) — superseded by `IsClose2`
- `IsClose3` (:151)
- `clusterArray` (:178) — superseded by `clusterArray2`
- `AppendDouble` (:326)

### `Code/GraphMol/MolDraw2D/DrawMolMCHLasso.cpp` (~53 lines, dead since 2023 / PR #6936)

- `arcEnd` (:177)
- `adjustLineEnd` (:629)
- `calcAnglesFromXAxis` (:652)
- `getArcEnds` (:716)

### Singletons

- `SortBasedOnFirstElement` — `Code/GraphMol/SmilesParse/SmilesWrite.cpp:485`
- `canRemoveNeg` — `Code/GraphMol/MolStandardize/Charge.cpp:332`
- `getBoolDefaultValue` — `Code/GraphMol/MolInterchange/Parser.cpp:139`
- `makeReactantErrorMessage` — `Code/GraphMol/ChemReactions/SanitizeRxn.cpp:153`
- `to_std_vector` — `Code/GraphMol/ChemReactions/Wrap/Enumerate.cpp:130`

---

## Tier 2b — Dead `#if 0` / commented-out code (clearly abandoned)

- `Code/GraphMol/FileParsers/TplFileParser.cpp:298` — `#if 0` buggy duplicate overload (~18 lines)
- `Code/GraphMol/QueryBond.cpp:65` — dead constructor remnant
- `Code/GraphMol/SynthonSpaceSearch/LONG_TEST_catch_tests.cpp:406` & `:430` — two `#if 0`
  dev tests with hardcoded `/Users/david/...` paths (~65 lines)
- `rdkit/Chem/Suppliers/UnitTestSDMolSupplier.py:116` — 3 commented-out test methods
  for a removed API (~57 lines)
- `rdkit/Chem/Subshape/BuilderUtils.py:80` — Python-2 `print` loop (won't parse on Py3)

---

## Tier 3 — Triage by hand (real, but possibly wanted)

- `Code/GraphMol/Wrap/Mol.cpp:476` — `AddConformersFromTrajectory` Python binding disabled via `#if 0`
- ~10 `def _test*` methods that silently never run (several marked `FIX:`). These are
  *disabled tests* — deleting loses the record that the case was meant to work.
  Triage individually rather than bulk-delete. Locations include:
  `rdkit/Chem/UnitTestMCS.py`, `Code/GraphMol/Depictor/Wrap/testDepictor.py:115`,
  `Code/GraphMol/DistGeomHelpers/Wrap/testDistGeom.py:102`,
  `Code/GraphMol/Wrap/rough_test.py`, `Code/GraphMol/ChemReactions/Wrap/testReactionWrapper.py:310`.
- `Code/GraphMol/DetermineBonds/Wrap/testDetermineBonds.py:122` — commented `testHueckelBonds`, marked "FIX: problematic"

---

## Tier 4 — Leave alone (verified intentional)

- `#if 0` benchmark / CIP-implementation toggles (`catch_canon.cpp`, `catch_moliterators.cpp`)
- `FIX:`-marked blocks tied to open issues (`atropisomers_catch.cpp` → #8340/#8341; `hanoitest.cpp`)
- `std::cerr` debug scaffolding behind `#if 0` (MolFileParser, SmilesParseOps, vf2, MolDraw2D, etc.)
- MMFF spec-deviation `#if 0`/`#if 1` toggles (`AtomTyper.cpp`)
- `Code/GraphMol/Basement/**` and `Code/Demos/**` — RDKit's explicit attic / standalone demos
- Unused private members/locals sweep came back **clean** — well-maintained code; that
  category needs clang `-Wunused-private-field`, not grep, for authoritative coverage.

---

## Notes on false positives ruled out

- `AtomPairWrapper.cpp`, `MorganWrapper.cpp`, `RDKitFPWrapper.cpp`,
  `TopologicalTorsionWrapper.cpp` look orphaned but are `#include`d into
  `Fingerprints/Wrap/FingerprintGeneratorWrapper.cpp`.
- `AugmentedAtomData.cpp` is `#include`d by `StructChecker/StructCheckerOptions.cpp`.
- `ConnectTheDots_Small/_Medium` in `ProximityBonds.cpp` matched grep but live inside a
  `/* ... */` comment.
- `@unittest.skipIf`/`skipUnless` decorators all gate on real build/feature flags — legit.
