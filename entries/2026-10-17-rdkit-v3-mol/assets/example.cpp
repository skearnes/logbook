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

// Calls against the mol.h sketch, checked for syntax only.

#include <iostream>

#include "mol.h"

namespace v3 = rdkit::v3;

v3::Mol ethanol() {
  v3::MolBuilder builder;
  const v3::AtomIndex c1 = builder.add_atom({.atomic_number = 6});
  const v3::AtomIndex c2 = builder.add_atom({.atomic_number = 6});
  const v3::AtomIndex o = builder.add_atom({.atomic_number = 8});
  builder.add_bond(c1, c2);
  builder.add_bond(c2, o);
  return builder.build();
}

void print_hydrogens(const v3::Mol &mol) {
  const v3::Mol with_hs = v3::add_hs(mol);
  for (const v3::Atom atom : with_hs.atoms()) {
    std::cout << atom.symbol() << ' ' << +atom.total_h_count() << '\n';
  }
}

// A hypervalent input fails the default build and passes a relaxed one; the
// result records the relaxation.
v3::Mol build_permissively(const v3::MolBuilder &builder) {
  try {
    return builder.build();
  } catch (const v3::BuildError &error) {
    for (const v3::Problem &problem : error.problems()) {
      std::cerr << problem.message << '\n';
    }
  }
  const v3::Mol mol = builder.build({.relax = v3::Relaxation::VALENCE});
  if (!v3::has(mol.relaxations(), v3::Relaxation::VALENCE)) {
    throw std::logic_error("relaxation not recorded");
  }
  return mol;
}

// Edits go through a builder; the original is unchanged.
v3::Mol strip_charges(const v3::Mol &mol) {
  v3::MolBuilder builder(mol);
  for (const v3::Atom atom : mol.atoms()) {
    v3::AtomSpec spec = builder.atom(atom.index());
    spec.formal_charge = 0;
    spec.implicit_h_count.reset();
    builder.set_atom(atom.index(), spec);
  }
  return builder.build();
}
