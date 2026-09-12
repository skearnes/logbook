# 2026-09-12 — Next-generation RDKit: scope, shape, and mechanics

- **Date:** 2026-09-12
- **Author:** Steven Kearnes
- **Status:** draft (working; Python surface next)
- **Tags:** rdkit, next-generation, api-design, cpp20, cpp-modules, packaging,
  conda-forge, proposal
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

Working notes toward a proposal for the next generation of RDKit. Feedback on
the July work set the frame: existing code must keep working, Rust is out in
favor of C++20, both generations must be usable from the one `rdkit` package,
and the proposal decides scope, shape, and mechanics rather than the new API or
the chemistry under it.

The shape settled so far:

- The new API exposes only types it owns, and inherits nothing without an
  explicit decision.
- Those types are implemented over the existing code. Implementations are
  replaced one type at a time, when the existing representation can't carry
  what the new API needs.
- A small shared lower layer serves both generations.
- The new API lives in a new directory beside `Code/`, behind a boundary CI
  enforces, until the triggers for a separate project fire.
- It is written as C++20 headers that can gain module interfaces later.
- conda-forge ships it as new outputs of the existing recipe.

The Python surface is the next question.

Follows [2026-07-18 — Reimplementing RDKit in Rust](../2026-07-18-rdkit-rust-reimplementation/README.md),
whose Rust and clean-room framings this entry revises, and the
[2026-07-15 maintenance plans](../2026-07-15-rdkit-maintenance-plans/README.md).
Measurements are against RDKit `master` at `4fdfe06d2` (2026-07-16) and the
conda-forge `rdkit-feedstock` at `2d4244d` (rdkit 2026.03.4).

## Constraints

From feedback on the July entries:

- **Existing code keeps working.** API changes are welcome; breaking current
  callers is not.
- **C++20 or newer; no Rust.** C++ modules would be welcome once support and
  adoption improve, but they are not required.
- **One package, both generations.** `rdkit` stays the package name. Users
  install it once and use whichever generation their code needs, the way
  pydantic 2 ships `pydantic.v1`.
- **Little duplicated code.**
- **Scope is shape and mechanics.** The new API and the routines under it (the
  hydrogen model, for example) are later decisions.
- **Splitting RDKit into several projects is an option.** The codebase is
  large.

## Evidence

### Where the code sits

| Layer | Contents | Lines |
| --- | --- | --- |
| Below the molecule type | Libraries with no link dependency on `GraphMol`, direct or indirect: `RDGeneral`, `DataStructs`, `Geometry`, `Numerics/Alignment`, `Numerics/EigenSolvers`, `ML/InfoTheory`, `ML/Cluster/Murtagh`, `SimDivPickers`, `Catalogs`, `ChemicalFeatures`, `RDStreams` | ~20k |
| Molecule core | Files directly under `GraphMol/`, plus `SmilesParse`, `Substruct`, `CIPLabeler` | ~60k |
| Built on the molecule type | Everything else: file formats, drawing, fingerprints, descriptors, reactions, force fields, … | ~215k |
| Bindings | `Wrap/` directories, `JavaWrappers`, `MinimalLib` | ~50k |

Lines are counted in `.cpp`, `.h`, `.hpp`, `.c`, `.cc`, `.yy`, and `.ll` files
under `Code/`, skipping files whose names contain `test`, `catch`, or `bench`
and the `Demos`, `Fuzz`, `Bench`, and `RDBoost` directories. Link dependencies
come from the `LINK_LIBRARIES` of each `rdkit_library()` call.

Two things follow:

- Most of RDKit is written directly against `ROMol` and `RWMol`. There is little
  to share below the molecule type, and porting is expensive everywhere above
  it.
- The boundary is blurrier than the link graph suggests. The generic
  force-field core and `Numerics/Optimizer` reach `GraphMol` only through
  `Trajectory/Snapshot.h`, which itself needs only `Geometry`, while the MMFF
  and UFF terms include molecule headers directly.

### RDKit already runs two API versions side by side

RDKit 2024.03.1 introduced a `v2` C++ API beside an `inline namespace v1`, and
nine headers in the parsers and molecule suppliers now declare one. The v1
functions are adapters that build v2 parameters and call v2
(`FileParsers.h:131-139`, `SmilesParse.h:132`), so nothing is duplicated. The release notes say "The
current API is still supported and will remain so for the forseeable future".

Two features of this precedent matter here:

- It works without duplication because both versions share `RWMol`.
- Today's v2 changes ownership (`std::unique_ptr` returns) and parameter
  passing while keeping v1's defaults. That is a choice; a separate namespace
  can change behavior too. What it can't do is give a shared type new
  invariants ([Decision 1](#1-the-new-api-exposes-only-its-own-types)).

### How pydantic ships two generations

- pydantic 2.13.4 contains a frozen copy of pydantic 1.10.26 under
  `pydantic/v1/`: 26 files and 13,166 lines, beside 32,585 lines of v2 Python.
- pydantic 1.10.17 added a `pydantic.v1` namespace to the 1.x line, so code that
  imports `pydantic.v1` runs against either major version. Downstream projects
  could change their imports first and upgrade later.
- The v2 validation engine is a separately named distribution,
  `pydantic-core`.

The packaging carries over: one install, both generations, and a bridge
release. The frozen copy does not. Freezing 13k lines of Python is cheap. A
frozen copy of RDKit's C++ duplicates everything and still has to keep building
against new compilers, Boost, Python, and NumPy.

### C++ modules today

- CMake supports named modules from 3.28, with the Ninja (1.11+) and Visual
  Studio 2022/2026 generators only, on GCC 14+, Clang 16+, and MSVC toolset
  14.34+ ([cmake-cxxmodules](https://cmake.org/cmake/help/latest/manual/cmake-cxxmodules.7.html)).
  AppleClang is not on the list, and Xcode's clang 21 ships no
  `clang-scan-deps`. RDKit requires CMake 3.18 and already builds as C++20.
- Header units (`import "header.h";`) are unsupported in CMake and experimental
  in Clang.
- `import std;` sits behind an experimental CMake gate and needs Clang 18.1.2+,
  MSVC 14.36+, or GCC 15+. Clang documents mixing it with `#include`d standard
  headers as poorly supported ([Clang modules documentation](https://clang.llvm.org/docs/StandardCPlusPlusModules.html)).
- Compiled module interfaces are specific to one compiler and flag set, so
  libraries ship interface sources and each consumer compiles its own.
- Existing headers can be wrapped. Headers included in a module's global module
  fragment (`module;`, the `#include`s, then `export module`) keep their
  declarations attached to the global module, and `export using` publishes the
  chosen names. A header included after the module declaration attaches its
  declarations to the new module instead, making them different entities from
  the ones `librdkit` defines.
- Macros don't cross `import`. Of RDKit's 501 headers, 50 use `PRECONDITION`,
  19 use `CHECK_INVARIANT`, and 7 use `BOOST_LOG`.
- Under policy CMP0155, which requiring CMake 3.28 turns on, every C++20 source
  is scanned for imports, and scanned sources are left out of unity builds. Unity builds were the main lever in the July [compile-time plan](../2026-07-15-rdkit-maintenance-plans/assets/COMPILE_TIME_PLAN.md).

### The conda-forge recipe

- One CMake build feeds five outputs. Three install from CMake install
  components: `librdkit` (`Unspecified`, `base`, `data`, `runtime`),
  `librdkit-dev` (`dev`), and `rdkit` (`extras` and `python`, then the `stubs`
  target and a `pip install` for `.dist-info`). `rdkit-dev` is an alias for
  `librdkit-dev` plus `rdkit`, and `rdkit-postgresql` runs the cartridge's own
  install script. One user-facing name over several outputs is already the
  pattern.
- `librdkit` is pinned exactly through `run_exports`, and every output that
  links it pins it exactly.
- Unix builds use `make` and Windows uses `NMake Makefiles JOM`; neither
  generator can scan modules.
- Compilers are `gxx` 14 on Linux, conda-forge's `clangxx` 19 on macOS (not
  AppleClang), and `vs2022` on Windows. `linux-aarch64`, `linux-ppc64le`, and
  `osx-arm64` are cross-compiled.
- `pyproject.patch` packages `rdkit`, `rdkit.*`, and `rdkit-stubs`.

## Decisions

### 1. The new API exposes only its own types

- The new public API names no existing type. Anything it inherits (a type, a
  function, a default) is an explicit, recorded decision, never the default.
- A separate namespace can change behavior as well as signatures. The limit is
  shared types: member functions aren't namespaced, and an `RWMol` can pass
  through existing functions and come back in a state the new API forbids.
  Owning the types removes that limit.
- New types hold existing objects privately, as wrappers rather than `using`
  aliases. Existing objects cross into the new API only through explicit
  conversion.

### 2. Implemented over the existing code, replaced type by type

- A new type starts with existing code behind it. Its implementation is
  replaced when the existing representation can't hold what the new API needs
  (a coordinate that isn't known, say) or when measurement shows it is too slow.
- Duplication is limited to what a replacement forces, one type at a time.
- Where a replacement reproduces existing behavior exactly, the existing API
  may become an adapter over it, as v1 is over v2 today. That takes
  differential tests showing identical output, and it is an option rather than
  a goal.

### 3. A shared lower layer

- The ~20k lines below the molecule type serve both generations and are
  modernized in place, guarded by the existing test suite.
- Decision 1 still applies: the new API exposes a shared type only by explicit
  choice.
- Moving `Snapshot` below the molecule type would free the generic force-field
  core and `Numerics/Optimizer`. How much that adds is unmeasured.

### 4. A new directory with an enforced boundary

- The new API lives in a new top-level directory beside `Code/`, as its own
  CMake project, built from the root behind an option.
- CI enforces the boundary. Public headers include only new-API headers, the
  shared layer, and the standard library; implementation files may also include
  existing headers.
- A directory comes first for three reasons:
  - changes to the shared layer need both test suites in one CI run;
  - the new API links the existing implementation;
  - one source tarball means one feedstock, with no ABI pins across feedstocks
    while the API is still changing.
- It becomes a separate project when any of these holds:
  - it needs a release schedule independent of RDKit's;
  - a different group of contributors develops it;
  - RDKit is split first and the shared layer becomes its own project;
  - enough implementations are replaced that it no longer links `librdkit`'s
    molecule code.

### 5. Headers now, modules later

Rules checked in CI from the start keep the headers ready for module
interfaces:

- The public API works without macros. Checks and logging are functions built
  on `std::source_location`.
- Headers declare no `static` functions or variables and no anonymous
  namespaces. Module interfaces reject inline code that refers to them.
- Every header compiles on its own.
- Components have no dependency cycles, and each component maps to one future
  module. Headers can include each other in a loop; modules can't import each
  other that way.
- The include boundary from Decision 4 holds.

Consequences:

- A header can't make a wrapped type usable without also making it nameable, so
  wrappers use pimpl and pay a non-inlined call per access. The July prototype measured de-inlining `ROMol`'s hot accessors at under 1%, with mixed sign ([BGL_DECOUPLING_FINDINGS.md](../2026-07-15-rdkit-maintenance-plans/assets/BGL_DECOUPLING_FINDINGS.md)), and `getAtomWithIdx()` is already out of line. With modules, a wrapper can inline its accessors while the wrapped type stays unnameable.
- No Ninja switch or CMake 3.28 minimum is needed yet.
- Adopting modules means adding interface files that wrap the same headers,
  shipped beside them so code that doesn't use modules keeps working. Only the
  new targets are scanned, so the existing build keeps its unity builds.
- Adopt modules when all three hold:
  - conda-forge's Linux compiler pin reaches GCC 15;
  - the compiler behind the macOS PyPI wheels can scan modules;
  - a wrapper module over `ROMol.h` builds on GCC, Clang, and MSVC.

  The feedstock must also switch to Ninja.

### 6. Packaging and names

- The new generation installs through its own CMake install components.
  conda-forge gains two outputs for it, a runtime library and a `-dev` package,
  each pinned exactly. The `rdkit` output depends on the new runtime, so
  `conda install rdkit` provides both generations. `librdkit` and
  `rdkit-postgresql` are unchanged.
- If the new API becomes a separate project, it gets its own feedstock, as
  `pydantic-core` has.
- The C++ namespace can't be `v2`, because `RDKit::v2` already names part of the
  existing API.
- A Python subpackage can't differ from an existing one only by case:
  `rdkit/chem` and `rdkit/Chem` are one directory on default macOS and Windows
  filesystems.

## Revised positions

- **Rust is off the table.** The July case for Rust was typestate, not memory
  safety. Whether C++20 can make "forgot to sanitize" a compile error becomes
  part of the API design, which this proposal leaves open.
- **No clean-room rewrite of the core.** The July specs-first framing assumed a
  new implementation from the start. Replacing implementations one type at a
  time, behind types the new API owns, keeps that option without paying for it
  up front. The [hydrogen model spec](../2026-07-18-rdkit-rust-reimplementation/assets/HYDROGEN_MODEL_SPEC.md) stays input for the later API decisions.
- **Namespaces don't limit behavior; shared types do.** A new namespace over
  existing code was first dismissed as able to change signatures but not
  behavior. It can change both. The real limit is a type shared with the
  existing API, which is why Decision 1 exists.

## Open questions

- **Python surface (next).**
  - The binding library for the new API (nanobind or Boost.Python), and whether
    two binding libraries can coexist in one process.
  - The subpackage name.
  - Converting objects between generations.
  - Pickling and type stubs.
- **Support for the existing API.** Fixes only, or features too, and for how
  long.
- **Other surfaces.** Whether and when Java (SWIG), JavaScript (MinimalLib),
  and the PostgreSQL cartridge get the new API.
- **PyPI wheels.** How the wheel build ships the new generation, and what that
  does to wheel size.
- **Splitting RDKit.** Which splits beyond the new directory are worth making,
  and in what order relative to this work.
- **Pimpl overhead.** Which accessors need inlining before modules are
  available.
- **Differences between generations.** How behavior that differs between the
  two APIs is documented and tested when both are importable.

## Next steps

- Settle the Python surface.
- Specify the CI checks: the include boundary, self-contained headers, no
  internal-linkage declarations in headers, and no component cycles.
- List the shared layer's macros that need function forms.
- Measure what moving `Snapshot` below the molecule type frees.
- Build a wrapper module over `ROMol.h`, and a module importing it, on GCC 14,
  conda-forge Clang 19, and MSVC. It doesn't block anything, but it settles the
  last module-adoption condition in Decision 5.
