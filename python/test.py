import numpy as np
from pyscf import gto
from pyT3NS import t3ns


def renyi_entropy(svals, α=1):
    """Calculates the renyi entropy for a given set of singular values.
    """
    if svals is None:
        return 0

    omega = svals * svals
    is_zero = np.isclose(svals, np.zeros(svals.shape), atol=1e-32)
    if α == 1:
        V = -omega * np.log(omega, where=np.logical_not(is_zero))
        V[is_zero] = 0
        return np.sum(V, axis=1)
    else:
        spomega = np.sum(np.power(omega, α), axis=1)
        is_zero = np.isclose(spomega, np.zeros(spomega.shape), atol=1e-32)
        return np.log(spomega, where=np.logical_not(is_zero)) / (1 - α)


mol = gto.Mole()
mol.build(
    atom="N 0 0 0; N 0 0 1.1208",
    spin=0,
    basis='sto3g',
    symmetry='C1',
    verbose=3
)

tree = t3ns.T3NS(mol, network='DMRG')
energy = tree.kernel(D=3, sitesize=1, max_sweeps=100, verbosity=1)
print(tree._bookkeeper)
print(f'E = {energy}')
exit()

tree = t3ns.T3NS(mol)
energy = tree.kernel(D=[300, (300, 500, 1e-3)], max_sweeps=100)
svds_prev = tree.singular_values()
entanglement = tree.disentangle()

energy = tree.kernel(D=[300, (300, 500, 1e-3)], max_sweeps=100)
print(tree._bookkeeper)
svds_after = tree.singular_values()
print(f'E = {energy}, entanglement = {entanglement}')

for i in [0.25, 0.5, 1]:
    before = sum(renyi_entropy(svds_prev, i))
    after = sum(renyi_entropy(svds_after, i))
    print(f"α = {i}: {before} → {after}")
