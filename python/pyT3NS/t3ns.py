import pyscf
import numpy
from pyT3NS import netw, bookkeeper
from ctypes import c_int, cdll, POINTER, c_double

libt3ns = cdll.LoadLibrary("libT3NS.so")

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
        _bookkeeper: The bookkeeper of the problem
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

            # if no orbital coefficients are given, use RHF orbitals
            if c is None:
                myhf = pyscf.scf.RHF(mol)
                myhf.verbose = self.verbose
                myhf.kernel()
                c = myhf.mo_coeff

            self._c = c
            self._mol = mol
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
            self._netw.pass_network()

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

    def init_bookkeeper(self, mstates, maxD, doci=False):
        '''This initializes a bookkeeper.

        The self._network object should already been made and should have
        a cnetwork attribute (so to know that the network was already
        communicated to the C library).

        If a bookkeeper was already initialized for the T3NS isnstance, than
        it will check if symmetries or target state are changed or if a new
        mstates is given. If so, the bookkeeper will be appropriately changed
        and the T3NS itself also, if already initialized.If _rOps exists, it
        will be deleted!
        '''
        if not hasattr(self, "_netw") or \
                not hasattr(self._netw, "cnetwork"):
            raise ValueError('No network initialized yet. Make sure it has '
                             'been passed to the C library through a call to'
                             ' pass_network()')

        if hasattr(self, "_bookkeeper"):
            pbookiep = self._bookkeeper
        else:
            pbookiep = None
        if doci:
            symmetries = ['U1']
            target = [0]
            for s, t in zip(self.symmetries, self.target):
                if s == ['U1']:
                    target[0] += t
            target[0] = int(target[0] % 2)

        self._bookkeeper = bookkeeper.Bookkeeper(symmetries, target, pbookiep,
                                                 maxD, mstates)

        print(self._bookkeeper)
        if self._bookkeeper != pbookiep:
            # Bookkeeper has changed
            exit(0)

    def kernel(self, D=500, sweeps=20, e_conv=1e-6, sites=2, david_rtl=1e-5,
               david_its=100, mstates=2, doci=False):
        '''Optimization of the tensor network.

        Args:
            mstates: The minimal number of states to be kept in each symmetry
            sector at initialisation.


        For the optimization scheme:
            For any of these arguments one can also pass an iterable instead.
            The optimization scheme will than be built of different regimes.
            If multiple arguments are an iterable, they should have the same
            length. Iterables exclude tuples for D.

            D: The bond dimension of the tensor network during optimization.
            This can either be an integer or a tuple with (minD, maxD,
            trunc_err). In the second case the optimization choose the bond
            dimension that meets the truncation error but with minD and maxD as
            lower and upper boundary.

            sweeps: The maximal number of sweeps through the network

            e_conv: The minimal energy difference between consecutive sweeps
            for convergence.

            sites: Number of sites to optimize together at one steP
            Possibilities are:
                - 1 for allowing a one site optimization
                - 2 for allowing a two site optimization
                - 3 for allowing a three site optimization of one branching
                  and 2 adjacent physical tensors
                - 4 for allowing a three site optimization of one branching
                  and 3 adjacent physical tensors

            david_rtl: The tolerance for the Davidson optimization.

            david_its: The maximal number of Davidson iterations.

            noise: The amount of noise to add after each optimization.
            The level of noise is scaled as
            0.5 * noise * (discarded weight last sweep)
        '''

        # Passes the network to the static global in the C library
        self._netw.pass_network()
        self.init_hamiltonian(doci)

        if hasattr(D, '__iter__') and not isinstance(D, tuple):
            fD = D[0]
        else:
            fD = D
        if isinstance(fD, tuple):
            maxD = fD[1]
        else:
            maxD = fD
        self.init_bookkeeper(mstates, maxD, doci)

    def init_hamiltonian(self, doci):
        from pyT3NS.bookkeeper import translate_irrep
        ham = c_int.in_dll(libt3ns, 'ham')

        # only qchem or doci allowed atm. DOCI == 3, qchem == 1
        nham = 1 + doci * 2
        if nham != ham.value:
            libt3ns.destroy_hamiltonian()
        else:
            return

        ham.value = nham

        if not hasattr(self, '_nuc'):
            self._nuc = self._mol.energy_nuc()
        if not hasattr(self, '_eri'):
            self._eri = pyscf.ao2mo.kernel(self._mol, self._c)
        if not hasattr(self, '_h1e'):
            self._h1e = self._c.T @ pyscf.scf.hf.get_hcore(self._mol) @ self._c

        h1e = self._h1e.astype(numpy.float64).reshape(self._h1e.shape,
                                                      order='F')
        assert h1e.ctypes.data % 16 == 0  # Check alignment
        h1e = h1e.ctypes.data_as(POINTER(c_double))

        norb = self._c.shape[0]
        fulleri = pyscf.ao2mo.restore(1, self._eri, norb)
        fulleri = fulleri.astype(numpy.float64).reshape(fulleri.shape,
                                                        order='F')
        assert fulleri.ctypes.data % 16 == 0  # Check alignment
        fulleri = fulleri.ctypes.data_as(POINTER(c_double))

        if doci:
            libt3ns.DOCI_ham_from_integrals.argtypes = \
                [c_int, POINTER(c_double), POINTER(c_double), c_double]
            libt3ns.DOCI_ham_from_integrals(norb, h1e, fulleri, self._nuc)
        else:
            libt3ns.QC_ham_from_integrals.argtypes = \
                [c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double),
                 c_double, c_int, c_int]
            irrep = None
            if self._pg_irrep is not None:
                for s in self.symmetries:
                    if s in supported_pgs:
                        irrep = (c_int * norb)(
                            *translate_irrep(self._pg_irrep, s)
                        )

            print([i for i in irrep])
            exit(0)
            libt3ns.QC_ham_from_integrals(
                norb,
                irrep,
                h1e,
                fulleri,
                self._nuc,
                int('SU2' in self.symmetries),
                int('SENIORITY' in self.symmetries),
            )
