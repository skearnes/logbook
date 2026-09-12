# 2026-09-12 — Next-generation RDKit: scope, shape, and mechanics

- **Date:** 2026-09-12
- **Author:** Steven Kearnes
- **Status:** draft (working)
- **Tags:** rdkit, next-generation, api-design, cpp20, cpp-modules, nanobind,
  packaging, conda-forge, proposal
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

Working notes toward a proposal for a next-generation RDKit API that ships
beside the existing one in the `rdkit` package. The proposal settles scope,
shape, and mechanics; the new API itself and the chemistry under it come later.
Decisions so far:

- The new API exposes only its own types and inherits nothing implicitly.
- Its types wrap existing code, and implementations are replaced one type at a
  time.
- The new API lives in its own directory behind a CI-enforced boundary.
- It is written as C++20 headers that can gain module interfaces later.
- In Python it is `rdkit.v3`, bound with nanobind and built only with RDKit's
  nanobind wrappers. Molecules convert between the APIs by copying.
- conda-forge ships it as new outputs of the existing recipe.

Whether code below the molecule type becomes a layer both APIs share is open.

Follows [2026-07-18 — Reimplementing RDKit in Rust](../2026-07-18-rdkit-rust-reimplementation/README.md)
and the [2026-07-15 maintenance plans](../2026-07-15-rdkit-maintenance-plans/README.md).
Measured against RDKit `master` at `23378a7` (2026-09-12) and conda-forge
`rdkit-feedstock` at `6918f4e` (rdkit 2026.03.6).

## Constraints

From feedback on the July entries:

- Existing code keeps working; API changes are otherwise welcome.
- C++20 or newer, no Rust. C++ modules are welcome once support matures, but
  not required.
- One package: installing `rdkit` provides both generations, the way pydantic 2
  ships `pydantic.v1`.
- Little duplicated code.
- Scope is shape and mechanics, not the new API or its routines (the hydrogen
  model, for example).
- Splitting RDKit into several projects is an option.

## Evidence

### Where the code sits

| Layer | Contents | Lines |
| --- | --- | --- |
| Below the molecule type | Libraries that don't link `GraphMol`, even indirectly: `RDGeneral`, `DataStructs`, `Geometry`, `Numerics/Alignment`, `Numerics/EigenSolvers`, `ML/InfoTheory`, `ML/Cluster/Murtagh`, `SimDivPickers`, `Catalogs`, `ChemicalFeatures`, `RDStreams` | ~20k |
| Molecule core | Files directly in `GraphMol/`, plus `SmilesParse`, `Substruct`, `CIPLabeler` | ~62k |
| Built on the molecule type | Everything else | ~220k |
| Bindings | `Wrap/` (Boost.Python) ~41k, `nbWrap/` (nanobind) ~38k, `JavaWrappers` and `MinimalLib` ~10k | ~88k |

Counts are lines in `.cpp`, `.h`, `.hpp`, `.c`, `.cc`, `.yy`, and `.ll` files
under `Code/`. Non-binding rows skip files named `*test*`, `*catch*`, or
`*bench*` and the `Demos`, `Fuzz`, `Bench`, and `RDBoost` directories. Link
dependencies come from each `rdkit_library()` call's `LINK_LIBRARIES`.

- Most code is written against `ROMol` and `RWMol`; only ~20k lines sit below
  the molecule type.
- The boundary is blurry. The generic force-field core and `Numerics/Optimizer`
  reach `GraphMol` only through `Trajectory/Snapshot.h`, which needs only
  `Geometry`; the MMFF and UFF terms include molecule headers directly.

### RDKit's `v2` C++ API

- Since 2024.03.1, nine parser and supplier headers declare a `v2` namespace
  beside an `inline namespace v1` whose functions call v2
  (`FileParsers.h:131-139`, `SmilesParse.h:132`). The release notes keep v1
  supported "for the forseeable future".
- Nothing is duplicated because both versions share `RWMol`. v2 changes
  ownership and parameter passing but keeps v1's defaults. A namespace could
  change behavior too, just not a shared type's invariants
  ([Decision 1](#1-the-new-api-exposes-only-its-own-types)).

### RDKit's nanobind port

- [rdkit#9030](https://github.com/rdkit/rdkit/pull/9030), merged 2026-08-11,
  put an `nbWrap/` beside each of the 59 `Wrap/` directories, `rdBase`
  included (`Code/RDBoost/nbWrap/RDBase.cpp`). Both builds run the same Python
  tests.
- A build uses one: Boost.Python is on by default, nanobind is off, and enabling
  both is a CMake error.
- Modules use `nanobind_add_module(... NB_SHARED ...)` without `NB_DOMAIN`
  (`RDKitUtils.cmake:228`), so extension modules see each other's bound types.
  nanobind is pinned to 2.x ([rdkit#9535](https://github.com/rdkit/rdkit/pull/9535)).
- The discussion that started the port ([rdkit#9031](https://github.com/rdkit/rdkit/discussions/9031))
  reported two to five times lower wrapper overhead than Boost.Python in early
  `Point2D` benchmarks.
- conda-forge builds the Boost.Python wrappers (`libboost-python-devel` in
  `host`).

### How TensorFlow and pydantic ship two generations

TensorFlow is the closer analogue, a C++ core under a Python API:

- TensorFlow 2 gave `tf` to the new API and kept 1.x as `tf.compat.v1`. Release
  1.14 added `tf.compat.v2` so libraries could publish code that runs on both.
- The `tf_export` decorator records each symbol's v1 and v2 names, and both
  namespaces are generated from those records, so unchanged symbols are written
  once.
- `tf_upgrade_v2` rewrites imports to `tf.compat.v1`. Removed modules such as
  `tf.contrib` had no path.
- Mixing is limited. The generations differ in process-wide runtime behaviors,
  eager execution can be switched only once per program, and some
  `tf.compat.v1` symbols don't work under TF2 behaviors.

pydantic 2.13.4 bundles a frozen 13,166-line copy of pydantic 1.10.26 as
`pydantic.v1`, and 1.10.17 added the same namespace to the 1.x line as a bridge.
Its engine ships separately as `pydantic-core`, and its migration guide says V1
and V2 models can't be mixed.

What carries over is one package, a bridge release, and one implementation
behind both surfaces. What doesn't:

- giving existing top-level names to the new API, which the constraint on
  existing code rules out;
- process-wide switches, which Decision 1 avoids;
- a frozen copy, which for RDKit's C++ would duplicate everything and still have
  to track compilers, Boost, Python, and NumPy.

### C++ modules today

- CMake supports them from 3.28 with the Ninja (1.11+) or Visual Studio
  generators, on GCC 14+, Clang 16+, or MSVC 14.34+
  ([cmake-cxxmodules](https://cmake.org/cmake/help/latest/manual/cmake-cxxmodules.7.html)).
  AppleClang isn't supported, and Xcode's clang 21 ships no `clang-scan-deps`.
  RDKit requires CMake 3.18 and builds as C++20.
- CMake doesn't support header units. `import std;` is experimental, and Clang
  documents mixing it with `#include`d standard headers as poorly supported
  ([Clang modules documentation](https://clang.llvm.org/docs/StandardCPlusPlusModules.html)).
- Compiled module interfaces are compiler- and flag-specific, so libraries ship
  interface sources.
- Existing headers can be wrapped. Included in a module's global module fragment
  (`module;` before `export module`), their declarations stay in the global
  module, and `export using` publishes chosen names. Included after
  `export module`, they would belong to the new module and stop matching
  `librdkit`.
- Macros don't cross `import`. Of 520 headers, 53 use `PRECONDITION`, 20
  `CHECK_INVARIANT`, and 7 `BOOST_LOG`.
- Requiring CMake 3.28 enables policy CMP0155, which scans every C++20 source
  and excludes scanned sources from unity builds, the main lever in the July [compile-time plan](../2026-07-15-rdkit-maintenance-plans/assets/COMPILE_TIME_PLAN.md).

### The conda-forge recipe

- One CMake build feeds five outputs. `librdkit` (`Unspecified`, `base`,
  `data`, `runtime`), `librdkit-dev` (`dev`), and `rdkit` (`extras`, `python`,
  the `stubs` target, and a `pip install` for `.dist-info`) install by CMake
  component. `rdkit-dev` is an alias, and `rdkit-postgresql` runs the
  cartridge's install script.
- `librdkit` is pinned exactly through `run_exports`, and every dependent pins
  it exactly.
- The generators are `make` and `NMake Makefiles JOM`, and neither scans
  modules.
- The compilers are `gxx` 15, conda-forge's `clangxx` 21 (not AppleClang), and
  `vs2022`. `linux-aarch64`, `linux-ppc64le`, and `osx-arm64` are
  cross-compiled.
- `rdkit` builds the Boost.Python wrappers and generates stubs with
  `pybind11-stubgen`. `pyproject.patch` packages `rdkit`, `rdkit.*`, and
  `rdkit-stubs`.

## Decisions

### 1. The new API exposes only its own types

- It names no existing type. Anything it inherits (a type, a function, a
  default) is an explicit, recorded decision.
- A namespace can change behavior but not a shared type's invariants. Member
  functions aren't namespaced, and an `RWMol` can pass through existing
  functions and return in a state the new API forbids.
- New types hold existing objects privately, as wrappers rather than `using`
  aliases. Existing objects enter only through explicit conversion.

### 2. Implemented over existing code, replaced type by type

- A new type starts with existing code behind it. Its implementation is replaced
  when the existing representation can't hold what the API needs (an unknown
  coordinate, say) or measurement shows it is too slow.
- Duplication is limited to what each replacement forces.
- If differential tests show a replacement reproduces existing behavior exactly,
  the existing API may become an adapter over it, as v1 is over v2. This is
  optional.

### 3. A new directory with an enforced boundary

- The new API gets a top-level directory beside `Code/`, as its own CMake
  project built from the root behind an option.
- CI enforces the boundary. Public headers include only new-API headers, the
  standard library, and existing headers explicitly approved under Decision 1;
  implementation files may include any existing header.
- A directory comes first for three reasons:
  - changes to the existing code the new API wraps need both test suites in one
    CI run;
  - the new API links the existing implementation;
  - one tarball means one feedstock, with no cross-feedstock ABI pins while the
    API is changing.
- It becomes a separate project if any of these happens:
  - it needs its own release schedule;
  - a different group develops it;
  - RDKit splits first and the code the new API wraps becomes its own project;
  - it no longer links `librdkit`'s molecule code.

### 4. Headers now, modules later

CI checks from the start keep headers ready for modules:

- The public API works without macros. Checks and logging are functions built
  on `std::source_location`.
- Headers declare no `static` functions or variables and no anonymous
  namespaces, which module interfaces reject in inline code.
- Every header compiles on its own.
- Components have no dependency cycles, since modules can't import in a cycle,
  and each maps to one future module.
- The include boundary from Decision 3 holds.

Consequences:

- Without modules, a wrapped type can't be usable yet unnameable. Wrappers
  therefore use pimpl and pay a non-inlined call per access. The July prototype
  measured de-inlining `ROMol`'s hot accessors at under 1%, with mixed sign ([BGL_DECOUPLING_FINDINGS.md](../2026-07-15-rdkit-maintenance-plans/assets/BGL_DECOUPLING_FINDINGS.md)),
  and `getAtomWithIdx()` is already out of line. With modules, wrappers can
  inline.
- Ninja and CMake 3.28 aren't needed before adopting modules.
- Adding modules means interface files that wrap the same headers, shipped
  beside them. Scanning covers only the new targets, so existing unity builds
  survive.
- Adopt modules when all of these hold:
  - conda-forge's `vs2022` toolset is 14.34+ (`gxx` 15 and `clangxx` 21
    already qualify);
  - the compiler behind the macOS PyPI wheels can scan modules;
  - a wrapper module over `ROMol.h` builds on GCC, Clang, and MSVC;
  - the feedstock builds with Ninja.
- `import std;` waits until CMake's support for it isn't experimental.

### 5. Packaging

- The new API installs through its own CMake components into two new conda-forge
  outputs, a runtime library and a `-dev` package, each pinned exactly. `rdkit`
  depends on the runtime, so `conda install rdkit` provides both generations.
  `librdkit` and `rdkit-postgresql` are unchanged.
- `nanobind` joins the host requirements. `libboost-python-devel` leaves once
  the recipe builds the nanobind wrappers.
- As a separate project, the new API would get its own feedstock, as
  `pydantic-core` has.

### 6. `rdkit.v3`

- The Python subpackage is `rdkit.v3`. In C++ the namespace is `rdkit::v3`, with
  headers under `rdkit/v3/`, and any future C++ modules are `rdkit.v3.*`.
- Why `v3`:
  - it stays accurate once it is the recommended API;
  - it follows `RDKit::v1` (the original C++ API) and `RDKit::v2` (2024.03);
  - it matches `pydantic.v1` and `tf.compat.v1`;
  - it clashes with no subpackage, even ignoring case, and can't be confused
    with RDKit's year.month releases.
- The docs need a sentence on why Python has no `rdkit.v2`.
- The namespace isn't `RDKit::v3`, because names the new code doesn't define
  would fall back to the enclosing `RDKit` namespace and silently reach existing
  functions. The tree has no lowercase `rdkit` namespace.
- Python names follow PEP 8: snake_case functions, methods, and modules, and
  CapWords classes. The existing API keeps its names and modules.
- The bindings use nanobind only.
- Molecules convert by copying in both directions, because a shared molecule
  would let existing code break the new type's rules. The conversions live in
  `rdkit.v3` and are the one place it names existing types. The existing API
  never depends on `rdkit.v3`.

## Top-level imports

Importing `rdkit` or any subpackage, `rdkit.v3` included, runs
`rdkit/__init__.py`. That file:

- imports the compiled `rdBase`, whose comment says exceptions leak memory
  otherwise, and takes `__version__` from it;
- under Jupyter or Colab, imports `rdkit.Chem.Draw.IPythonConsole` (and with it
  `rdkit.Chem`, `rdchem`, `rdChemReactions`, and `rdMolDraw2D`) and sends RDKit
  logs to stderr;
- configures the `rdkit` logger;
- patches iteration on `rdBase` vector types.

### Keeping Boost.Python out of `rdkit.v3`

`rdkit.v3` is built only with the nanobind wrappers:

- In that build `rdBase` is a nanobind module, so `rdkit/__init__.py` loads no
  Boost.Python as it stands.
- Conversion goes through nanobind. An `rdkit.v3` function can take and return
  an existing `Mol`, because the extension modules share bound types.
- `rdkit.v3` ships once the nanobind wrappers are the build RDKit ships.

### Deferred

- **`from rdkit import Mol`.** Re-exporting only adds names, so this stays
  possible, with two limits:
  - a new top-level name can't reuse a subpackage's name, since importing
    `rdkit.Geometry` would overwrite a class `Geometry`;
  - a new top-level module can't be a lowercase version of an existing one,
    since `rdkit/geometry` and `rdkit/Geometry` are one directory on default
    macOS and Windows filesystems.
- **Import cost.** Until `rdkit/__init__.py` changes, importing `rdkit.v3` also
  loads `rdBase` and, in notebooks, the existing drawing stack.

Revisiting the top-level imports must preserve:

- `rdkit.__version__`;
- `rdkit.rdBase` without a separate import;
- `rdBase` loading before other extension modules;
- the `rdkit` logger setup;
- notebook rendering of existing molecules;
- the vector iteration patch.

## Revised positions

- **Rust is out.** Its July appeal was typestate. Whether C++20 can make "forgot
  to sanitize" a compile error is an API question left open.
- **No clean-room rewrite of the core.** Replacing implementations type by type,
  behind owned types, keeps that option without paying for it up front. The [hydrogen model spec](../2026-07-18-rdkit-rust-reimplementation/assets/HYDROGEN_MODEL_SPEC.md) remains input to later API work.
- **Shared types, not namespaces, limit behavior changes** (Decision 1).

## Open questions

- Whether code below the molecule type becomes a layer both APIs share, and on
  what terms. Answering it needs a survey of those ~20k lines and of the
  force-field code that reaches `GraphMol` only through `Snapshot`.
- When the nanobind build becomes the default upstream and on conda-forge, since
  `rdkit.v3` waits on it.
- Stubs and pickling for `rdkit.v3`: nanobind's stub generator or the recipe's
  `pybind11-stubgen`.
- What configures RDKit's C++ logging when only `rdkit.v3` is imported.
- Support for the existing API: fixes only or features too, and for how long.
- Whether and when Java (SWIG), JavaScript (MinimalLib), and the PostgreSQL
  cartridge get the new API.
- How PyPI wheels ship `rdkit.v3`, and what that does to wheel size.
- Which RDKit splits beyond the new directory are worth making, and when.
- Which accessors need inlining before modules arrive.
- How behavior that differs between the two APIs is documented and tested.

## Next steps

- Survey the code below the molecule type: what it contains, what the new API
  would use from it, which macros would need function forms, and what moving
  `Snapshot` would free.
- Ask RDKit's maintainers when the nanobind wrappers become the default build.
- In a nanobind build, prototype copying an existing `Mol` into a stand-in
  `rdkit.v3` type and back.
- Specify the CI checks: include boundary, self-contained headers, no internal
  linkage in headers, and no component cycles.
- Build a wrapper module over `ROMol.h`, plus a module importing it, with
  conda-forge's `gxx` 15, `clangxx` 21, and `vs2022`. This checks one of
  Decision 4's adoption conditions.
