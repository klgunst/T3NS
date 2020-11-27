from pyscf import gto, lo
from pyT3NS import t3ns
import numpy as np
import sys
import os

"""
An examplar calculation of a T3NS-DOCI-calculation (without orbital opt).

The integrals are generated with pyscf and the eventual output file generated
on a quad-core laptop can be found in `examples/Ru_ph/doci.out`.

If you pass 'localized' as argument to the file, the orbitals will first be
localized.
"""

with open(os.path.abspath(__file__)) as f:
    print(f.read())

print()
print("######################################################################")
print()

mol = gto.M(
    atom=open('Ru_ph_ph_ph.xyz').read(),
    unit='ang',
    basis='STO3G',
    charge=2,
    verbose=3,
    symmetry='C1',
    max_memory=13000,
)

mf = mol.HF().run()

# split localizing by Pipek-Mezey
if len(sys.argv) == 2 and sys.argv[1] == 'localized':
    print("Using PM localized orbitals")
    C = np.hstack((
        lo.PM(mol, mf.mo_coeff[:, mf.mo_occ > 0]).kernel(),
        lo.PM(mol, mf.mo_coeff[:, mf.mo_occ == 0]).kernel()))
else:
    print("Using RHF orbitals")
    C = mf.mo_coeff

tree = t3ns.T3NS(mol, c=C, network='T3NS')

tree.kernel(D=50, doci=True, max_sweeps=10, verbosity=2)
tree.disentangle()
tree.kernel(D=50, doci=True, max_sweeps=10, verbosity=2)

E = tree.kernel(D=100, doci=True, max_sweeps=100,
                energy_conv=1e-5, verbosity=3)

print()
print(f"E@D=100: {E}")
print("E@D=75:", tree.kernel(D=75, doci=True, energy_conv=1e-5,
                             max_sweeps=100, verbosity=0))
print("E@D=50:", tree.kernel(D=50, doci=True, energy_conv=1e-5,
                             max_sweeps=100, verbosity=0))
