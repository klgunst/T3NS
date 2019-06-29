#!/bin/env python3
from pyscf import gto
from pyT3NS import t3ns


mol = gto.Mole()
R = 1.1208
mol.build(
    atom=f"N 0 0 -{R / 2}; N 0 0 {R / 2}",
    spin=0,
    basis='sto3g',
    symmetry=True,
    verbose=3
)

tree = t3ns.T3NS(mol)
print(tree.symmetries)
print(tree._pg_irrep)
print(mol.spin)
