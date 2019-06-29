import pyscf
import numpy
from pyT3NS import netw

supported_pgs = ['D2h', 'C2v', 'C2h', 'D2', 'Cs', 'C2', 'Ci', 'C1']


class T3NS:
    '''Class for the optimization of the three-legged tree tensor network
    state.

    Attributes:
        symmetries: The symmetries used for the optimization
        target: The targetstate of the optimization

        _T3NS: The wave function
        _rOps: The renormalized operators
        _netw: The network for the problem
        _c:
        _nuc:
        _h1e:
        _eri:
        _pg_irrep:
        verbose: The verbosity
    '''
    def __init__(self, mol_or_hdf5, c=None, network=None, verbose=None):
        '''Initializing the T3NS calculation.

        Args:
            mol_or_hdf5: An instance of :class:`Mole` from PySCF
            or a string to a path of the hdf5 file. In this case c and network
            should be left as defaults.

            c: The orbital coefficients.
            The default are the RHF orbitals.

            network: Network of the tensor network.
            This can be the string 'DMRG' for a linear chain or 'T3NS' for a
            T3NS that is radially extended (i.e. the maximal distance between
            two sites is minimal). An instance of :class:`Network` is also
            accepted for custom networks or a path to a network file.
            Default is 'T3NS' if Mole instance was passed, else None.
        '''
        from pyscf import symm
        if isinstance(mol_or_hdf5, pyscf.gto.Mole):
            # Initialize with RHF
            mol = mol_or_hdf5
            if network is None:
                network = 'T3NS'

            if verbose is None:
                self.verbose = mol.verbose

            self.mol = mol

            myhf = pyscf.scf.RHF(mol)
            myhf.verbose = self.verbose
            # if no orbital coefficients are given, use RHF orbitals
            if c is None:
                myhf.kernel()
                c = myhf.mo_coeff
            self._c = c
            self._eri = pyscf.ao2mo.kernel(mol, c)
            self.symmetries = ['Z2', 'U1', 'SU2']
            self.target = [int(sum(mol.nelec) % 2), sum(mol.nelec), mol.spin]
            try:
                irrep_ids = symm.label_orb_symm(mol, mol.irrep_id,
                                                mol.symm_orb, c)
                newsym = mol.groupname
                if mol.groupname == 'Dooh' or mol.groupname == 'Coov':
                    newsym = 'D2h' if mol.groupname else 'C2v'
                    if self.verbose > 1:
                        print(f'Changed point group from {mol.groupname} to '
                              f'{newsym}')

                # mod 10 needed for translation Coov and Dooh
                # See https://sunqm.github.io/pyscf/symm.html
                self._pg_irrep = [
                    symm.irrep_id2name(newsym, int(ii % 10))
                    for ii in irrep_ids
                ]
                self.symmetries.append(newsym)
                self.target.append(symm.irrep_id2name(newsym, 0))
            except ValueError:
                self._pg_irrep = None

            if isinstance(network, netw.Network):
                self._netw = network
            elif network == 'DMRG' or network == 'T3NS':
                isDMRG = network == 'DMRG'
                self._netw = netw.Network(c.shape[0], isDMRG=isDMRG)

                # exchange matrix
                Kij = numpy.zeros((c.shape[0],) * 2)
                trids = numpy.tril_indices(c.shape[0], 0)
                for el, r, c in zip(self._eri.diagonal(), trids[0], trids[1]):
                    Kij[r, c] = el
                    Kij[c, r] = el
                cost = self._netw.optimize(Kij)
                if self.verbose >= 2:
                    print(f'Optimization Exchange-cost: {cost}')

            elif isinstance(network, str):
                self._netw = netw.Network()
                self._netw.readnetworkfile(network)
            else:
                raise ValueError(f'{network} is invalid for network')

        elif isinstance(mol_or_hdf5, str):
            # h5path = mol_or_hdf5
            if verbose is None:
                self.verbose = 1

            if c is not None or network is not None:
                raise ValueError(
                    'A hdf5 file does not expect a defined c or network'
                )
            raise ValueError('Need to implement hdf5 readin still')
        else:
            raise ValueError(
                'Expects a Mole instance or a path to a hdf5 file'
            )
