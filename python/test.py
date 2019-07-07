#!/bin/env python3
from pyscf import gto
from pyT3NS import t3ns


mol = gto.Mole()
mol.build(
    atom="N 0 0 0; N 0 0 1.1208",
    spin=0,
    basis='sto3g',
    symmetry=True,
    verbose=3
)

tree = t3ns.T3NS(mol, network='DMRG')
energy = tree.kernel(D=[300, (300, 500, 1e-3)], max_sweeps=100)
entanglement = tree.disentangle()
print(f'E = {energy}, entanglement = {entanglement}')
