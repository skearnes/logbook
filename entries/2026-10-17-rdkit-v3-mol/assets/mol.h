// Copyright 2026 Steven Kearnes
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Declarations for the proposed rdkit::v3 molecule API. A design sketch: it
// compiles on its own (see check.sh) and has no implementation.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace rdkit::v3 {

// ---------------------------------------------------------------------------
// Indices and enumerations

enum class AtomIndex : std::uint32_t {};
enum class BondIndex : std::uint32_t {};

// Kekulé order. AROMATIC appears only on a molecule built with
// Relaxation::KEKULIZATION, for bonds that could not be assigned an order.
enum class BondOrder : std::uint8_t {
  ZERO,
  SINGLE,
  DOUBLE,
  TRIPLE,
  QUADRUPLE,
  DATIVE,
  AROMATIC,
};

enum class Hybridization : std::uint8_t {
  UNSPECIFIED,
  S,
  SP,
  SP2,
  SP3,
  SP3D,
  SP3D2,
  OTHER,
};

enum class AromaticityModel : std::uint8_t { RDKIT, SIMPLE, MDL };

// Checks a build may skip. A relaxed check never leaves a value undefined.
enum class Relaxation : std::uint32_t {
  NONE = 0,
  VALENCE = 1U << 0,       // Allow atoms above their allowed valence.
  KEKULIZATION = 1U << 1,  // Keep BondOrder::AROMATIC where no order fits.
  AROMATICITY = 1U << 2,   // Keep aromatic flags as given.
  CLEANUP = 1U << 3,       // Skip normalizations such as nitro groups.
  STEREO = 1U << 4,        // Keep stereo on atoms that are not stereocenters.
};

constexpr Relaxation operator|(Relaxation a, Relaxation b) {
  return static_cast<Relaxation>(static_cast<std::uint32_t>(a) |
                                 static_cast<std::uint32_t>(b));
}

constexpr bool has(Relaxation set, Relaxation flag) {
  return (static_cast<std::uint32_t>(set) & static_cast<std::uint32_t>(flag)) !=
         0;
}

// ---------------------------------------------------------------------------
// Properties: a closed set of value types that each map to one Python type.

using PropertyValue =
    std::variant<bool, std::int64_t, double, std::string,
                 std::vector<std::int64_t>, std::vector<double>,
                 std::vector<std::string>>;

using Properties = std::map<std::string, PropertyValue, std::less<>>;

// ---------------------------------------------------------------------------
// Stereo: elements name their reference atoms, so storage order of neighbors
// or bonds never changes their meaning.

enum class Chirality : std::uint8_t { CLOCKWISE, ANTICLOCKWISE };

struct TetrahedralStereo {
  AtomIndex center{};
  // Viewed from neighbors[0], the rest run in `chirality` order. An empty
  // slot is an implicit hydrogen or a lone pair.
  std::array<std::optional<AtomIndex>, 4> neighbors{};
  Chirality chirality{};
};

enum class Configuration : std::uint8_t { CIS, TRANS };

struct DoubleBondStereo {
  BondIndex bond{};
  // One neighbor of each end atom; `configuration` relates these two.
  std::array<AtomIndex, 2> reference{};
  Configuration configuration{};
};

// Atropisomers and non-tetrahedral centers join the variant as new types.
using StereoElement = std::variant<TetrahedralStereo, DoubleBondStereo>;

enum class StereoGroupKind : std::uint8_t { ABSOLUTE, AND, OR };

struct StereoGroup {
  StereoGroupKind kind{};
  std::uint32_t id{};
  std::vector<std::size_t> elements{};  // Indices into Mol::stereo().
};

// ---------------------------------------------------------------------------
// Attachments

class Conformer {
 public:
  // A 2D conformer has dimension 2 and zero z coordinates.
  Conformer(std::vector<std::array<double, 3>> positions, int dimension,
            std::uint32_t id);

  std::uint32_t id() const;
  int dimension() const;
  std::size_t num_atoms() const;
  std::span<const std::array<double, 3>> positions() const;
  std::array<double, 3> position(AtomIndex atom) const;

 private:
  std::vector<std::array<double, 3>> positions_;
  int dimension_;
  std::uint32_t id_;
};

struct SubstanceGroup {
  std::string type{};  // Molfile SGroup type, such as "SRU" or "DAT".
  std::vector<AtomIndex> atoms{};
  std::vector<BondIndex> bonds{};
  std::vector<std::size_t> parents{};  // Indices into Mol::substance_groups().
  Properties properties{};
};

// Per-atom residue data from PDB and similar formats.
struct ResidueInfo {
  std::string name{};
  std::int32_t number{};
  std::string chain{};
  std::string atom_name{};
  std::string insertion_code{};
  bool is_hetero{};
};

// ---------------------------------------------------------------------------
// Views. Atom and Bond are an index into a Mol; they hold no state of their
// own and stay valid as long as the Mol does.

class Mol;
class Atom;
class Bond;

// Random-access ranges of views. Each binds once, as a Python sequence that
// keeps its Mol alive.
class AtomRange {
 public:
  class iterator {
   public:
    using value_type = Atom;
    using difference_type = std::ptrdiff_t;
    Atom operator*() const;
    iterator &operator++();
    iterator operator++(int);
    bool operator==(const iterator &other) const;
  };
  std::size_t size() const;
  Atom operator[](std::size_t i) const;
  iterator begin() const;
  iterator end() const;
};

class BondRange {
 public:
  class iterator {
   public:
    using value_type = Bond;
    using difference_type = std::ptrdiff_t;
    Bond operator*() const;
    iterator &operator++();
    iterator operator++(int);
    bool operator==(const iterator &other) const;
  };
  std::size_t size() const;
  Bond operator[](std::size_t i) const;
  iterator begin() const;
  iterator end() const;
};

class Atom {
 public:
  AtomIndex index() const;
  const Mol &mol() const;

  std::uint8_t atomic_number() const;
  std::string_view symbol() const;
  std::int8_t formal_charge() const;
  std::optional<std::uint16_t> isotope() const;
  std::optional<std::uint32_t> map_number() const;
  std::uint8_t radical_electrons() const;

  // Hydrogens come in two kinds: graph atoms (explicit) and a count on their
  // heavy atom (implicit). total_h_count() is their sum.
  std::uint8_t implicit_h_count() const;
  std::uint8_t explicit_h_count() const;
  std::uint8_t total_h_count() const;
  bool implicit_h_count_was_given() const;

  std::uint32_t degree() const;  // Graph neighbors, explicit Hs included.
  std::uint32_t valence() const;
  bool is_aromatic() const;
  Hybridization hybridization() const;
  bool in_ring() const;

  AtomRange neighbors() const;
  BondRange bonds() const;
  std::span<const AtomIndex> neighbor_indices() const;
  const Properties &properties() const;
  const std::optional<ResidueInfo> &residue() const;

 private:
  friend class Mol;
  friend class AtomRange;
  Atom(const Mol &mol, AtomIndex index);
  const Mol *mol_;
  AtomIndex index_;
};

class Bond {
 public:
  BondIndex index() const;
  const Mol &mol() const;

  Atom begin_atom() const;
  Atom end_atom() const;
  Atom other_atom(AtomIndex atom) const;
  BondOrder order() const;
  bool is_aromatic() const;
  bool is_conjugated() const;
  bool in_ring() const;
  const Properties &properties() const;

 private:
  friend class Mol;
  friend class BondRange;
  Bond(const Mol &mol, BondIndex index);
  const Mol *mol_;
  BondIndex index_;
};

// Symmetrized SSSR, computed when the molecule is built.
class RingInfo {
 public:
  std::size_t num_rings() const;
  std::span<const std::vector<AtomIndex>> atom_rings() const;
  std::span<const std::vector<BondIndex>> bond_rings() const;
  std::uint32_t num_atom_rings(AtomIndex atom) const;
  std::uint32_t num_bond_rings(BondIndex bond) const;
  bool atom_in_ring_of_size(AtomIndex atom, std::size_t size) const;
  bool bond_in_ring_of_size(BondIndex bond, std::size_t size) const;
  std::size_t smallest_atom_ring(AtomIndex atom) const;  // 0 if acyclic.
};

// ---------------------------------------------------------------------------
// Mol: a finished, checked molecule. Immutable; copies share storage.

namespace detail {
struct MolData;
}  // namespace detail

class Mol {
 public:
  Mol();  // Empty.

  std::size_t num_atoms() const;
  std::size_t num_bonds() const;
  Atom atom(AtomIndex index) const;  // Throws std::out_of_range.
  Bond bond(BondIndex index) const;  // Throws std::out_of_range.
  std::optional<Bond> bond_between(AtomIndex a, AtomIndex b) const;
  AtomRange atoms() const;
  BondRange bonds() const;

  // Columns, one entry per atom or bond. Python sees read-only NumPy views.
  std::span<const std::uint8_t> atomic_numbers() const;
  std::span<const std::int8_t> formal_charges() const;
  std::span<const std::uint8_t> implicit_h_counts() const;
  std::span<const std::array<AtomIndex, 2>> bond_atoms() const;
  std::span<const BondOrder> bond_orders() const;

  const RingInfo &rings() const;
  std::span<const StereoElement> stereo() const;
  std::span<const StereoGroup> stereo_groups() const;
  std::span<const Conformer> conformers() const;
  std::span<const SubstanceGroup> substance_groups() const;
  const Properties &properties() const;
  Relaxation relaxations() const;

  // Edits that need no new perception return a new Mol sharing the rest.
  Mol with_properties(Properties properties) const;
  Mol with_conformers(std::vector<Conformer> conformers) const;

  std::string to_bytes() const;
  static Mol from_bytes(std::string_view bytes);

 private:
  std::shared_ptr<const detail::MolData> data_;
};

// ---------------------------------------------------------------------------
// Building

struct AtomSpec {
  std::uint8_t atomic_number = 0;
  std::int8_t formal_charge = 0;
  std::optional<std::uint16_t> isotope{};
  std::optional<std::uint32_t> map_number{};
  std::optional<std::uint8_t> implicit_h_count{};   // Empty: inferred.
  std::optional<std::uint8_t> radical_electrons{};  // Empty: inferred.
  bool aromatic = false;
  std::optional<ResidueInfo> residue{};
  Properties properties{};
};

struct BondSpec {
  BondOrder order = BondOrder::SINGLE;
  bool aromatic = false;
  Properties properties{};
};

struct BuildOptions {
  AromaticityModel aromaticity = AromaticityModel::RDKIT;
  Relaxation relax = Relaxation::NONE;
};

struct Problem {
  enum class Kind : std::uint8_t { VALENCE, KEKULIZATION, STEREO, OTHER };
  Kind kind{};
  std::vector<AtomIndex> atoms{};
  std::vector<BondIndex> bonds{};
  std::string message{};
};

class Error : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Lists every problem found, not only the first.
class BuildError : public Error {
 public:
  explicit BuildError(std::vector<Problem> problems);
  std::span<const Problem> problems() const;

 private:
  std::vector<Problem> problems_;
};

// Mutable and unchecked. Indices stay stable until build(); removed atoms and
// bonds leave gaps that build() closes.
class MolBuilder {
 public:
  MolBuilder();
  explicit MolBuilder(const Mol &mol);

  std::size_t num_atoms() const;
  std::size_t num_bonds() const;
  const AtomSpec &atom(AtomIndex index) const;
  const BondSpec &bond(BondIndex index) const;
  std::array<AtomIndex, 2> bond_atoms(BondIndex index) const;

  AtomIndex add_atom(const AtomSpec &atom);
  BondIndex add_bond(AtomIndex begin, AtomIndex end, const BondSpec &bond = {});
  void set_atom(AtomIndex index, const AtomSpec &atom);
  void set_bond(BondIndex index, const BondSpec &bond);
  void remove_atom(AtomIndex index);  // Also removes its bonds.
  void remove_bond(BondIndex index);
  AtomIndex insert(const Mol &mol);  // Returns the first new index.
  // Where an atom lands in the built Mol: build() keeps surviving atoms in
  // order, so this is the index minus the removed atoms before it.
  std::optional<AtomIndex> built_index(AtomIndex index) const;

  void add_stereo(StereoElement element);
  void add_stereo_group(StereoGroup group);
  void add_conformer(Conformer conformer);
  void add_substance_group(SubstanceGroup group);
  void set_properties(Properties properties);

  // Throws BuildError.
  Mol build(const BuildOptions &options = {}) const;
};

// ---------------------------------------------------------------------------
// Operations. Each takes a Mol and returns a new value.

struct AddHsOptions {
  bool only_on_given_atoms = false;
  std::vector<AtomIndex> atoms{};
  bool add_coordinates = true;
};

struct RemoveHsOptions {
  // By default a hydrogen stays explicit if it carries anything a count can't:
  // an isotope, a map number, properties, stereo, a substance group, a charge,
  // or a bond to anything other than one heavy atom.
  bool keep_isotopes = true;
  bool keep_mapped = true;
  bool keep_stereo_defining = true;
};

Mol add_hs(const Mol &mol, const AddHsOptions &options = {});
Mol remove_hs(const Mol &mol, const RemoveHsOptions &options = {});
Mol with_aromaticity(const Mol &mol, AromaticityModel model);
Mol renumber_atoms(const Mol &mol, std::span<const AtomIndex> new_order);
Mol combine(const Mol &a, const Mol &b);
Mol stereo_from_conformer(const Mol &mol, std::uint32_t conformer_id);
Mol without_stereo(const Mol &mol);

std::vector<Mol> fragments(const Mol &mol);
std::vector<std::vector<AtomIndex>> fragment_atoms(const Mol &mol);
std::vector<StereoElement> possible_stereo(const Mol &mol);

enum class CipLabel : std::uint8_t { R, S, r, s, M, P, m, p, E, Z, SEQ_CIS, SEQ_TRANS };
// One label per element of mol.stereo(); empty where no label applies.
std::vector<std::optional<CipLabel>> cip_labels(const Mol &mol);

std::int32_t total_formal_charge(const Mol &mol);
std::vector<std::uint8_t> adjacency_matrix(const Mol &mol);  // Row-major.
std::vector<double> distance_matrix(const Mol &mol);         // Row-major.

}  // namespace rdkit::v3
