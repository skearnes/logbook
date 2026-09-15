# 2026-09-12 — Next-generation RDKit: scope, shape, and mechanics

- **Date:** 2026-09-12
- **Author:** Steven Kearnes
- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 5)
- **Status:** draft (working)
- **Tags:** rdkit, next-generation, api-design, cpp20, cpp-modules, nanobind,
  packaging, conda-forge, proposal
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

Working notes toward a proposal for a next-generation RDKit API that ships
beside the existing one in the `rdkit` package. The proposal settles scope,
shape, and mechanics; the new API and the chemistry under it come later.
Decisions so far:

- The new API owns the types whose rules it changes, re-exports others by
  recorded decision, and inherits nothing implicitly.
- Its types wrap existing code, and implementations are replaced one type at a
  time.
- It lives in its own directory behind a CI-enforced boundary.
- It is written as C++20 headers that can gain module interfaces later.
- In Python it is `rdkit.v3`, bound with nanobind and built only with RDKit's
  nanobind wrappers. Molecules convert between the APIs by copying.
- conda-forge ships it as new outputs of the existing recipe.

Which types below the molecule it re-exports, and whether that code becomes a
layer both APIs share, is open.

Follows
[2026-07-18 — Reimplementing RDKit in Rust](../2026-07-18-rdkit-rust-reimplementation/README.md)
and the
[2026-07-15 maintenance plans](../2026-07-15-rdkit-maintenance-plans/README.md).
Measured against RDKit `master` at `23378a7` (2026-09-12) and conda-forge
`rdkit-feedstock` at `6918f4e` (rdkit 2026.03.6).

## Constraints

From feedback on the July entries:

- Existing code keeps working; API changes are otherwise welcome.
- C++20 or newer, no Rust. Rust's July appeal, typestate, becomes a C++ API
  question (Decision 1). C++ modules are welcome once support matures, but not
  required.
- One package: installing `rdkit` provides both generations, the way pydantic 2
  ships `pydantic.v1`.
- Little duplicated code, and no clean-room rewrite of the core (Decision 2).
- Scope is shape and mechanics, not the new API or its routines. The July
  [hydrogen model spec](../2026-07-18-rdkit-rust-reimplementation/assets/HYDROGEN_MODEL_SPEC.md)
  remains input to later API work.
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

- Most code is written against `ROMol` and `RWMol`.
- The boundary is blurry: the generic force-field core and `Numerics/Optimizer`
  reach `GraphMol` only through `Trajectory/Snapshot.h`, which needs only
  `Geometry`, while the MMFF and UFF terms include molecule headers directly.

### RDKit's `v2` C++ API

- Since 2024.03.1, nine parser and supplier headers declare a `v2` namespace
  beside an `inline namespace v1` whose functions call v2
  (`FileParsers.h:131-139`, `SmilesParse.h:132`). The release notes keep v1
  supported "for the forseeable future".
- Both share `RWMol`, so nothing is duplicated; v2 changes ownership and
  parameter passing but keeps v1's defaults.

### What existing code can do to a molecule

- `const` doesn't mean unchanged. Properties live in a `mutable Dict`, "a quirk
  of history" (`RDProps.h:19-20`), so `setProp`, `clearProp`, and
  `clearComputedProps` are `const`. `getRingInfo() const` returns a non-const
  `RingInfo *` (`ROMol.h:774`) that `findSSSR`, `fastFindRings`, and
  `findRingFamilies` fill from a `const ROMol &` (`MolOps.h:867-886`), and
  `MolToSmiles` stores `_smilesAtomOutputOrder` and `_smilesBondOutputOrder` on
  its `const` input (`SmilesWrite.cpp:747-750`).
- A stale cache fails when read: `getNumImplicitHs()` checks at runtime that
  `calcImplicitValence()` has run (`Atom.cpp:305-307`).
- In Python, `Chem.SanitizeMol`, `Chem.Kekulize`, `Chem.SetAromaticity`,
  `Chem.AssignStereochemistry`, and `AllChem.EmbedMolecule` change their
  argument in place.
- Behavior changes have shipped as process-wide switches. New stereo perception
  went behind a global flag in 2022.09
  ([rdkit#5309](https://github.com/rdkit/rdkit/pull/5309)), and 2023.03 removed
  the per-call `SmilesParserParams.useLegacyStereo` in its favor. The flag is
  the environment variable `RDK_USE_LEGACY_STEREO_PERCEPTION`, which the setter
  writes with `setenv` and every query reads (`Chirality.cpp:848-858`); legacy
  perception is the default (`Chirality.h:37-39`), and non-tetrahedral stereo is
  switched the same way.

### RDKit's nanobind port

- [rdkit#9030](https://github.com/rdkit/rdkit/pull/9030), merged 2026-08-11, put
  an `nbWrap/` beside each of the 59 `Wrap/` directories, `rdBase` included
  (`Code/RDBoost/nbWrap/RDBase.cpp`), and both builds run the same Python tests.
- A build uses one: Boost.Python is on by default, nanobind is off, and enabling
  both is a CMake error.
- Modules use `nanobind_add_module(... NB_SHARED ...)` without `NB_DOMAIN`
  (`RDKitUtils.cmake:228`), so they share bound types; nanobind is pinned to 2.x
  ([rdkit#9535](https://github.com/rdkit/rdkit/pull/9535)).
- A C++ type therefore has one Python class. Binding it again warns, returns the
  class registered first, sets no module attribute, and adds the second
  binding's methods to that class (nanobind 2.15.0
  `src/nb_type.cpp:1323-1339,1665`, `include/nanobind/nb_class.h:666-688`).
  `ROMol` is bound as `Mol` (`GraphMol/nbWrap/Mol.cpp:426`).
- The discussion that started the port
  ([rdkit#9031](https://github.com/rdkit/rdkit/discussions/9031)) reported two
  to five times lower wrapper overhead than Boost.Python in early `Point2D`
  benchmarks.

### How TensorFlow and pydantic ship two generations

TensorFlow, a C++ core under a Python API, is the closer analog:

- TensorFlow 2 gave `tf` to the new API and kept 1.x as `tf.compat.v1`; release
  1.14 added `tf.compat.v2` so libraries could publish code that runs on both.
- The `tf_export` decorator records each symbol's v1 and v2 names, and both
  namespaces are generated from those records, so unchanged symbols are written
  once.
- Mixing is limited: the generations differ in process-wide runtime behaviors,
  and eager execution can be switched only once per program.

pydantic 2.13.4 bundles a frozen 13,166-line copy of pydantic 1.10.26 as
`pydantic.v1`, and 1.10.17 added the same namespace to the 1.x line as a bridge.
Its engine ships separately as `pydantic-core`, and its migration guide says V1
and V2 models can't be mixed.

One package and one implementation behind both surfaces carry over. What
doesn't:

- existing top-level names for the new API, which the constraint on existing
  code rules out;
- a bridge release, unneeded because existing names never move;
- process-wide switches, which Decision 1 avoids;
- a frozen copy, which would duplicate all of RDKit's C++ and track compilers,
  Boost, Python, and NumPy.

### C++ modules today

- CMake supports modules from 3.28 with the Ninja (1.11+) or Visual Studio
  generators, on GCC 14+, Clang 16+, or MSVC 14.34+
  ([cmake-cxxmodules](https://cmake.org/cmake/help/latest/manual/cmake-cxxmodules.7.html)),
  but not header units. `import std;` is experimental, and Clang documents
  mixing it with `#include`d standard headers as poorly supported
  ([Clang modules documentation](https://clang.llvm.org/docs/StandardCPlusPlusModules.html)).
  AppleClang isn't supported: Apple's clang 21 ships no `clang-scan-deps` and
  rejects `export module` under `-std=c++20`. RDKit requires CMake 3.18 and
  builds as C++20.
- Compiled module interfaces are compiler- and flag-specific, so libraries ship
  interface sources.
- Existing headers can be wrapped. Included in a global module fragment
  (`module;` before `export module`), their declarations stay in the global
  module, and `export using` publishes chosen names; included after
  `export module`, they would belong to the new module and stop matching
  `librdkit`. The `std` modules of libc++, libstdc++, and MSVC's STL all start
  with `module;` (`std.cppm.in:14`, `std.cc.in:24`, `std.ixx:5`, fetched
  2026-09-15).
- Macros don't cross `import`. Of 520 headers, 53 use `PRECONDITION`, 20
  `CHECK_INVARIANT`, and 7 `BOOST_LOG`.
- A module interface can't expose an internal-linkage entity, such as a
  namespace-scope `static` function, from inline code or declarations; reading
  an internal constant's value is exempt, and `static` members take their
  class's linkage ([basic.link](https://eel.is/c++draft/basic.link)). GCC
  accepts such exposures from a global module fragment with a warning,
  `-Wexpose-global-module-tu-local`
  ([GCC options](https://gcc.gnu.org/onlinedocs/gcc/C_002b_002b-Dialect-Options.html)).
- Requiring CMake 3.28 at the root enables policy CMP0155, which scans every
  C++20 source and excludes scanned sources from unity builds, the main lever in
  the July
  [compile-time plan](../2026-07-15-rdkit-maintenance-plans/assets/COMPILE_TIME_PLAN.md).
  Policies are set per directory and recorded on each target, and
  `CXX_SCAN_FOR_MODULES` overrides them per target or source
  ([cmake_policy](https://cmake.org/cmake/help/latest/command/cmake_policy.html),
  [UNITY_BUILD](https://cmake.org/cmake/help/latest/prop_tgt/UNITY_BUILD.html)).
- A probe ([run.sh](assets/module_probe/run.sh)) with Homebrew Clang 23.1.1 and
  GCC 16.2.0, against conda-forge's `librdkit-dev` 2026.03.6 and Boost 1.90,
  found:
  - A module wrapping `ROMol.h` and `SmilesParse.h` in a global module fragment
    compiles on both. Its importer parses a SMILES without naming `RWMol` and,
    on Clang, links against `librdkit` and runs; naming `RDKit::ROMol` there
    fails on both.
  - GCC warns 30 times that the wrapped headers expose internal constants, 27
    from Boost and 3 from RDKit (`RDTypeTag::AnyTag`, `zero_tolerance`); Clang
    reports none.
  - A function attached to a named module gets a module-specific symbol
    (`_ZW8attached6answerv`); a wrapped one keeps its ordinary symbol.
  - Wrapped, GCC warns when an inline function uses an internal function or a
    constant's address and rejects one that uses an anonymous-namespace type,
    while Clang accepts all three. Attached, GCC rejects all three and Clang
    warns only on the function. Reading a constant's value passes everywhere.

### The conda-forge recipe

- One CMake build feeds five outputs. `librdkit` (`Unspecified`, `base`, `data`,
  `runtime`), `librdkit-dev` (`dev`), and `rdkit` (`extras`, `python`, the
  `stubs` target, and a `pip install` for `.dist-info`) install by CMake
  component; `rdkit-dev` is an alias, and `rdkit-postgresql` runs the
  cartridge's install script.
- `librdkit` is pinned exactly through `run_exports`, and every dependent pins
  it exactly.
- The generators, `make` and `NMake Makefiles JOM`, don't scan modules.
- The compilers are `gxx` 15, conda-forge's `clangxx` 21 (not AppleClang), and
  `vs2022`; `linux-aarch64`, `linux-ppc64le`, and `osx-arm64` are
  cross-compiled.
- `rdkit` builds the Boost.Python wrappers (`libboost-python-devel` in `host`)
  and generates stubs with `pybind11-stubgen`; `pyproject.patch` packages
  `rdkit`, `rdkit.*`, and `rdkit-stubs`.

## Decisions

Each decision has a strict alternative: how it would work, what it would enable,
why it isn't the choice now, and when to adopt it. Strictness that a future
capability such as C++ modules needs is worth adopting early when it is cheap.

### 1. The new API owns the types whose rules it changes

- A type whose rules, behavior, or representation the new API changes is its
  own: a wrapper holding the existing object privately, not a `using` alias. The
  molecule, atoms, and bonds are the expected cases, along with any type whose
  public members reach them, such as `Conformer` (`getOwningMol()`,
  `Conformer.h:91`) and `StereoGroup` (`Atom *` members, `StereoGroup.h:47`).
- A type used unchanged, such as a point or a bit vector, can be re-exported.
  Each re-export, like any inherited function or default, is an explicit,
  recorded decision.
- Existing objects enter owned types only through explicit conversion.

Why own a type whose rules change:

- **A shared type can't gain rules.** Existing code keeps working, so a shared
  molecule must keep allowing every state that code builds and relies on. A new
  rule, such as the hydrogen model's requirement that hydrogens be graph nodes,
  needs a type existing code can't reach.
- **Nothing guards a shared molecule.** Namespaces, like v2's, don't cover
  member functions, and `const` functions change molecules
  ([evidence](#what-existing-code-can-do-to-a-molecule)), so the new API would
  recheck its rules on every entry. An owned type checks once, at conversion,
  and only its own methods change it after.
- **Global switches don't mix.** A switch applies to every library in the
  process, the same limit on mixing TensorFlow's generations. An owned type
  carries a new behavior itself, so both behaviors meet only at conversion,
  though wrapped code reads existing switches too
  ([Open questions](#open-questions)).
- **Later changes stay possible.** Through an alias, `ROMol`'s members become
  API, down to `getTopology()` returning a Boost Graph Library graph
  (`ROMol.h:877`), and Decision 2 couldn't replace what's behind them. If C++20
  can make "forgot to sanitize" a compile error, that state lives in the type,
  and a shared `RWMol` has no room for it.

Owning costs Decision 4's non-inlined call per access and Decision 6's copy at
conversion. A re-export avoids both and lets old and new code share objects, at
a recorded price:

- its header joins the public headers (Decision 3);
- argument-dependent lookup finds its namespace's free functions for new-API
  callers;
- its implementation can't be replaced without breaking new-API code (Decision
  2);
- in Python it keeps its existing class and names unless `rdkit.v3` wraps it
  (Decision 6).

> [!NOTE]
> **Strict alternative: own every type.**
>
> - **How.** Wrap points and bit vectors too, so public headers include no
>   existing header.
> - **Enables.** No re-export costs, and no existing headers, macros, or Boost
>   in module interfaces.
> - **Not now.** Every value type would need a wrapper and a copy at each
>   boundary, against the constraint on duplicated code, and modules don't
>   require it (Decision 4).
> - **Adopt when.** A re-exported type needs a new representation, or a module
>   interface must not reach existing headers.

### 2. Implemented over existing code, replaced type by type

- A new type starts with existing code behind it and gets a new implementation
  only when the existing representation can't hold what the API needs (an
  unknown coordinate, say) or measurement shows it is too slow.
- Duplication stays limited to what each replacement forces, and a clean-room
  rewrite stays possible without being paid for up front.
- If differential tests show a replacement reproduces existing behavior exactly,
  the existing API may become an adapter over it, as v1 is over v2. This is
  optional.

> [!NOTE]
> **Strict alternative: wrapped routines read no process-wide switch.**
>
> - **How.** Each switch a wrapped routine reads, such as stereo perception's
>   environment variable (`Chirality.cpp:2658,2919`), becomes a parameter that
>   defaults to the switch for existing callers.
> - **Enables.** New-API behavior independent of environment variables and other
>   libraries' settings, closing the open question on switches, and differential
>   tests that pin both behaviors in one process.
> - **Not now.** It changes existing code routine by routine, and which routines
>   get wrapped isn't decided.
> - **Adopt when.** The first wrapped routine reads a switch.

### 3. A new directory with an enforced boundary

- The new API gets a top-level directory beside `Code/`, as its own CMake
  project built from the root behind an option that is on by default from the
  start.
- CI enforces the boundary. Public headers include only new-API headers, the
  standard library, and existing headers approved under Decision 1;
  implementation files may include any existing header.
- A directory comes first because:
  - changes to wrapped code need both test suites in one CI run;
  - the new API links the existing implementation;
  - one tarball means one feedstock, with no cross-feedstock ABI pins while the
    API is changing.
- It becomes a separate project if:
  - it needs its own release schedule;
  - a different group develops it;
  - RDKit splits first and the wrapped code becomes its own project;
  - it stops linking `librdkit`'s molecule code.

> [!NOTE]
> **Strict alternative: build against an installed RDKit.**
>
> - **How.** A CI job installs RDKit and configures the new directory on its own
>   with `find_package(rdkit)`, so it sees only installed headers and exported
>   targets.
> - **Enables.** A split whenever a trigger fires, with build-tree dependencies
>   caught in CI.
> - **Not now.** A second configure and build in CI while no split is planned.
> - **Adopt when.** A split trigger becomes likely.

### 4. Headers now, modules later

CI checks from the start keep headers ready for modules:

- The public API works without macros; checks and logging are functions built on
  `std::source_location`.
- Every header compiles on its own.
- Public headers declare no anonymous namespaces, whose types GCC rejects in
  inline code even when wrapped ([evidence](#c-modules-today)). Other internal
  linkage isn't linted.
- Components have no dependency cycles, since modules can't import in a cycle,
  and each maps to one future module.
- The include boundary from Decision 3 holds.

> [!NOTE]
> **Strict alternative: lint internal linkage.**
>
> - **How.** Public headers declare no namespace-scope `static` function or
>   variable, and their namespace-scope constants are `inline constexpr`.
> - **Enables.** Headers compile attached to a named module, as modules first
>   needs.
> - **Not now.** Wrapped headers don't need it ([evidence](#c-modules-today)),
>   and MSVC is unchecked.
> - **Adopt when.** Modules first is chosen, or a compiler rejects a wrapped
>   exposure. New headers make it nearly free to adopt early.

Consequences:

- Without modules, a wrapped type can't be usable yet unnameable, so wrappers
  use pimpl and pay a non-inlined call per access. The July prototype measured
  de-inlining `ROMol`'s hot accessors at under 1%, with mixed sign
  ([BGL_DECOUPLING_FINDINGS.md](../2026-07-15-rdkit-maintenance-plans/assets/BGL_DECOUPLING_FINDINGS.md)),
  and `getAtomWithIdx()` is already out of line. With modules, wrappers can
  inline.
- Ninja and CMake 3.28 aren't needed before adopting modules.
- Module interface files will wrap the same headers and ship beside them. The
  new directory sets CMP0155 in its own scope while the root keeps its CMake
  3.18 policies, so scanning covers only the new targets and existing unity
  builds survive.
- Adopt modules when all of these hold:
  - conda-forge's `vs2022` toolset is 14.34+ (`gxx` 15 and `clangxx` 21 already
    qualify);
  - the compiler behind the macOS PyPI wheels can scan modules;
  - a wrapper module over `ROMol.h` builds on GCC, Clang, and MSVC;
  - the feedstock builds with Ninja.
- `import std;` waits until CMake's support for it isn't experimental.

> [!NOTE]
> **Strict alternative: modules first.**
>
> - **How.** Named module interfaces attached to `rdkit.v3.*` modules from the
>   start, with existing headers only in implementation units or global module
>   fragments, CMake 3.28+ with Ninja or Visual Studio for every build, and the
>   internal-linkage lint above.
> - **Enables.** Wrapped types are reachable but unnameable, so wrappers inline
>   without pimpl, and importers see neither macros nor existing names.
> - **Not now.** conda-forge builds with `make` and `NMake Makefiles JOM`, and
>   Apple's clang 21 rejects `export module`. An entity attached to a named
>   module also can't be declared in a header
>   ([basic.link](https://eel.is/c++draft/basic.link)) and gets a different
>   symbol, so SWIG, MinimalLib, and the cartridge would need a separate
>   interface.
> - **Adopt when.** The adoption conditions above hold and every consumer of the
>   new API can import modules.

### 5. Packaging

- The new API installs through its own CMake components into two new conda-forge
  outputs, a runtime library and a `-dev` package, each pinned exactly. `rdkit`
  depends on the runtime, so `conda install rdkit` provides both generations;
  `librdkit` and `rdkit-postgresql` are unchanged.
- `nanobind` joins the host requirements, and `libboost-python-devel` leaves
  once the recipe builds the nanobind wrappers.
- As a separate project, the new API would get its own feedstock, as
  `pydantic-core` has.

> [!NOTE]
> **Strict alternative: a stable C++ ABI.**
>
> - **How.** Version the new API's shared library by ABI, check each release
>   with an ABI checker such as libabigail, and pin dependents to a compatible
>   range instead of an exact build.
> - **Enables.** C++ packages built against the new API that survive RDKit
>   releases without rebuilds, and a later release schedule of its own.
> - **Not now.** The API is changing, and exact pins keep one feedstock simple
>   (Decision 3).
> - **Adopt when.** The API is declared stable, or a split trigger makes
>   separate releases likely.

### 6. `rdkit.v3`

- The Python subpackage is `rdkit.v3`, with headers under `rdkit/v3/` and any
  future C++ modules named `rdkit.v3.*`.
- At the start it is reached only by importing `rdkit.v3` explicitly; exposing
  it from the top level comes later ([Deferred](#deferred)).
- Why `v3`:
  - stays accurate once it is the recommended API;
  - follows the parsers' `RDKit::v1` and `RDKit::v2` namespaces (2024.03);
  - matches `pydantic.v1` and `tf.compat.v1`;
  - clashes with no subpackage, even ignoring case, and can't be confused with
    RDKit's year.month releases.
- The docs need a sentence on why Python has no `rdkit.v2`.
- **Proposed C++ namespace: `rdkit::v3`.**
  - Unqualified names resolve through enclosing namespaces, even in out-of-line
    definitions, and no language feature prevents that. The enclosing `rdkit`
    declares nothing, so existing names used without `RDKit::` fail to compile.
  - Sibling generations such as `rdkit::v4` aren't searched.
  - The only lint needed bans `using namespace RDKit` in the new directory.
  - The cost is two top-level namespaces that differ only by case; the tree has
    no lowercase `rdkit` namespace.
- **Alternative: `RDKit::v3` with a linter.** It keeps one top-level namespace,
  but unqualified names fall back to `RDKit` and silently reach existing
  functions. Blocking that takes a custom Clang check (clang-query or
  clang-tidy) covering names in expressions, types, and templates, maintained
  across Clang versions.
- In either namespace, argument-dependent lookup finds `RDKit` functions for
  `RDKit` arguments, which implementation code that wraps existing objects
  relies on.
- Python names follow PEP 8: snake_case functions, methods, and modules, and
  CapWords classes. The existing API keeps its names and modules.
- The bindings use nanobind only, as their own set of wrappers in the existing
  wrappers' nanobind domain, so conversion can pass existing objects. They bind
  no C++ type the existing wrappers bind, because a second binding merges into
  the class registered first and the result depends on import order
  ([evidence](#rdkits-nanobind-port)).
- Molecules convert by copying in both directions, because a shared molecule
  would let existing code break the new type's rules. The conversions and
  re-exports are the only places `rdkit.v3` names existing types, and the
  existing API never depends on `rdkit.v3`.

> [!NOTE]
> **Strict alternative: a separate nanobind domain.**
>
> - **How.** Build `rdkit.v3` with `NB_DOMAIN`, which gives it its own nanobind
>   library and type registry (nanobind 2.15.0
>   `cmake/nanobind-config.cmake:425-437`), and convert molecules through a copy
>   that passes no bound object, such as RDKit's binary pickle.
> - **Enables.** Binding any C++ type, re-exports included, under PEP 8 names
>   with no import-order effects, and changing nanobind versions independently.
> - **Not now.** Conversion couldn't take an existing `Mol` directly, and a
>   re-exported type would get a second Python class that existing functions
>   reject.
> - **Adopt when.** `rdkit.v3` ships apart from the existing wrappers or needs
>   its own nanobind version.

## Top-level imports

Importing `rdkit` or any subpackage, `rdkit.v3` included, runs
`rdkit/__init__.py`, which:

- imports the compiled `rdBase`, without which a comment says exceptions leak
  memory, and takes `__version__` from it; loading `rdBase` sets up C++ logging
  (`Code/RDBoost/nbWrap/RDBase.cpp:285`);
- under Jupyter or Colab, imports `rdkit.Chem.Draw.IPythonConsole` (and with it
  `rdkit.Chem`, `rdchem`, `rdChemReactions`, and `rdMolDraw2D`) and sends RDKit
  logs to stderr;
- configures the `rdkit` logger;
- patches iteration on `rdBase` vector types.

A revision must keep all of this, including `rdkit.rdBase` without a separate
import, `rdBase` loading before other extension modules, and notebook rendering
of existing molecules.

`rdkit.v3` is built only with the nanobind wrappers and ships once they are the
build RDKit ships. In that build `rdBase` is a nanobind module, so
`rdkit/__init__.py` as it stands loads no Boost.Python, and an `rdkit.v3`
function can take and return an existing `Mol` because the extension modules
share bound types.

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

## Open questions

- Which types below the molecule the new API re-exports, and whether that code
  becomes a layer both APIs share. Answering it needs a survey of those ~20k
  lines and of the force-field code that reaches `GraphMol` only through
  `Snapshot`.
- When the nanobind build becomes the default upstream and on conda-forge, since
  `rdkit.v3` waits on it.
- How `rdkit.v3` keeps process-wide switches read by the code it wraps from
  changing its behavior ([evidence](#what-existing-code-can-do-to-a-molecule)).
- Stubs and pickling for `rdkit.v3`: nanobind's stub generator or the recipe's
  `pybind11-stubgen`.
- Support for the existing API: fixes only or features too, and for how long.
- Whether and when Java (SWIG), JavaScript (MinimalLib), and the PostgreSQL
  cartridge get the new API.
- How PyPI wheels ship `rdkit.v3`, and what that does to wheel size.
- Which RDKit splits beyond the new directory are worth making, and when.
- Which accessors need inlining before modules arrive.
- How behavior that differs between the two APIs is documented and tested.

## Next steps

- Survey the code below the molecule type: what it contains, which types the new
  API would re-export, which macros would need function forms, and what moving
  `Snapshot` would free.
- Ask RDKit's maintainers when the nanobind wrappers become the default build.
- In a nanobind build, prototype copying an existing `Mol` into a stand-in
  `rdkit.v3` type and back.
- Specify the CI checks: include boundary, macro-free public API, self-contained
  headers, no anonymous namespaces in public headers, no component cycles, and
  the namespace rule from Decision 6.
- Repeat the [module probe](assets/module_probe/run.sh), which passes with
  Homebrew GCC 16 and Clang 23, with conda-forge's `gxx` 15, `clangxx` 21, and
  `vs2022`, to check one of Decision 4's adoption conditions.
