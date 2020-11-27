import pyscf
import numpy
from pyT3NS import netw, bookkeeper, tensors
from ctypes import c_int, cdll, POINTER, c_double, c_char, byref, c_void_p, \
    cast, Structure, c_char_p, c_bool

libt3ns = cdll.LoadLibrary("libT3NS.so")

supported_pgs = ['D2h', 'C2v', 'C2h', 'D2', 'Cs', 'C2', 'Ci', 'C1']


class SvalSelect(Structure):
    _fields_ = [
        ("minD", c_int),
        ("maxD", c_int),
        ("truncerr", c_double),
    ]

    def __init__(self, D):
        if isinstance(D, tuple):
            self.minD = D[0]
            self.maxD = D[1]
            self.truncerr = D[2]
        else:
            self.minD = D
            self.maxD = D
            self.truncerr = 0

    def __str__(self):
        return f"(min: {self.minD}, max: {self.maxD}, trunc: {self.truncerr})"


class DisentScheme(Structure):
    _fields_ = [
        ("max_sweeps", c_int),
        ("gambling", c_bool),
        ("beta", c_double),
        ("svd_sel", SvalSelect)
    ]

    def __init__(self, D, max_sweeps=30, gambling=True, beta=20):
        self.max_sweeps = max_sweeps
        self.gambling = gambling
        self.beta = beta
        self.svd_sel = SvalSelect(D)

    def __str__(self):
        if self.gambling:
            string = f"Metropolis like acceptance with beta = {self.beta}"
        else:
            string = "Accepting of best permutation"
        string += f" ({self.max_sweeps} sweeps)\n"
        return string + f"Truncation: {self.svd_sel}\n"


class Regime(Structure):
    _fields_ = [
        ("svd_sel", SvalSelect),
        ("sitesize", c_int),
        ("davidson_rtl", c_double),
        ("davidson_max_its", c_int),
        ("max_sweeps", c_int),
        ("energy_conv", c_double),
        ("noise", c_double)
    ]

    def __init__(self, D, sitesize=2, davidson_rtl=1e-5, davidson_max_its=100,
                 max_sweeps=20, energy_conv=1e-6, noise=0):
        self.svd_sel = SvalSelect(D)
        self.sitesize = sitesize
        self.davidson_rtl = davidson_rtl
        self.davidson_max_its = davidson_max_its
        self.max_sweeps = max_sweeps
        self.energy_conv = energy_conv
        self.noise = noise

    def __str__(self):
        one_two_three = {1: 'one', 2: 'two', 3: 'three', 4: 'four'}
        return f"{one_two_three[self.sitesize]} site optimization\n" + \
            f"Truncation: {self.svd_sel}\n" + \
            f"Davidson: (tol: {self.davidson_rtl}, " + \
            f"its: {self.davidson_max_its})\n" + \
            f"Maximal sweeps: {self.max_sweeps}\n" + \
            f"Energy convergence: {self.energy_conv}\n" + \
            f"Added noise: {self.noise}\n"


class OptScheme(Structure):
    _fields_ = [
        ("nrRegimes", c_int),
        ("regimes", POINTER(Regime))
    ]

    def __init__(self, D, **kwargs):
        self.nrRegimes = len(D) if isinstance(D, list) else 1

        for k, v in kwargs.items():
            if isinstance(v, list):
                assert self.nrRegimes == 1 or self.nrRegimes == len(v)
                self.nrRegimes = len(v)

        self.regimes = (Regime * self.nrRegimes)()

        for i in range(self.nrRegimes):
            reg_data = {}
            reg_data['D'] = D[i] if isinstance(D, list) else D
            for k, v in kwargs.items():
                reg_data[k] = v[i] if isinstance(v, list) else v
            self.regimes[i] = Regime(**reg_data)

    def __str__(self):
        result = ""
        for i in range(self.nrRegimes):
            result += f"Regime {i}:\n{self.regimes[i]}\n"
        return result


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
        _lastD:
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
            self.target = [int(mol.nelectron % 2), mol.nelectron, mol.spin]

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
            except (ValueError, TypeError):
                self._pg_irrep = None

            if isinstance(network, netw.Network):
                self._netw = network
            elif network == 'DMRG' or network == 'T3NS':
                isDMRG = network == 'DMRG'
                self._netw = netw.Network(c.shape[1], isDMRG=isDMRG)

                # exchange matrix
                # Kij = numpy.zeros((c.shape[1],) * 2)
                # trids = numpy.tril_indices(c.shape[1], 0)
                # for el, r, c in zip(self._eri.diagonal(), trids[0], trids[1]):
                #     Kij[r, c] = el
                #     Kij[c, r] = el
                # cost = self._netw.optimize(Kij)
                # if self.verbose >= 2:
                #     print(f'Optimization Exchange-cost: {cost}')
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

    def __del__(self):
        libt3ns.destroy_hamiltonian()
        libt3ns.clear_instructions()
        for rops in self._rOps:
            rops.delete()
        bookie = bookkeeper.Bookkeeper.in_dll(libt3ns, "bookie")
        libt3ns.destroy_bookkeeper.argtypes = [POINTER(bookkeeper.Bookkeeper)]
        libt3ns.destroy_bookkeeper(byref(bookie))

    def kernel(self, D=500, mstates=None, doci=False, **kwargs):
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
        pbookie = self._bookkeeper if hasattr(self, '_bookkeeper') else None
        if mstates is None:
            mstates = 2 if pbookie is None else 0

        # Passes the network to the static global in the C library
        self._netw.pass_network()

        # Get global bookkeeper
        self._bookkeeper = bookkeeper.Bookkeeper.in_dll(libt3ns, "bookie")
        # fillin symmetries and target state
        if doci:
            symmetries = ['U1']
            target = [0]
            for s, t in zip(self.symmetries, self.target):
                if s == 'U1':
                    target[0] += t
                if s == 'SU2' and t != 0:
                    print('Executing DOCI for non-singlet calculation')
            target[0] = target[0] // 2
        else:
            symmetries = self.symmetries
            target = self.target
        self._bookkeeper.fill_symmetry_and_target(symmetries, target)

        # Initialize the Hamiltonian
        self.init_hamiltonian(doci)

        if hasattr(D, '__iter__') and not isinstance(D, tuple):
            fD = D[0]
        else:
            fD = D
        if isinstance(fD, tuple):
            maxD = fD[1]
        else:
            maxD = fD

        self._bookkeeper.init_bookkeeper(pbookie, maxD, mstates)
        self.init_wave_function(pbookie)
        self.init_operators()
        self.energy = self.execute_optimization(D, **kwargs)
        return self.energy

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

        norb = self._c.shape[1]
        if self._eri.size == norb ** 4 or doci:
            fulleri = self._eri
        else:
            fulleri = pyscf.ao2mo.restore(1, self._eri, norb)

        if not doci:
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
                 c_double, c_int, c_int, c_int]
            irrep = None
            if self._pg_irrep is not None:
                for s in self.symmetries:
                    if s in supported_pgs:
                        irrep = (c_int * norb)(
                            *translate_irrep(self._pg_irrep, s)
                        )

            libt3ns.QC_ham_from_integrals(
                norb,
                irrep,
                h1e,
                fulleri,
                self._nuc,
                3,
                int('SU2' in self.symmetries),
                int('SENIORITY' in self.symmetries),
            )

    def init_wave_function(self, pbookie=None):
        initwav = libt3ns.init_wave_function
        initwav.argtypes = [POINTER(c_void_p), c_int,
                            POINTER(bookkeeper.Bookkeeper), c_char]

        ppbookie = None if pbookie is None else byref(pbookie)

        if not hasattr(self, '_T3NS'):
            self._T3NS = c_void_p(None)

        initwav(cast(byref(self._T3NS), POINTER(c_void_p)), 0, ppbookie,
                'r'.encode('utf8'))
        self._T3NS = \
            cast(self._T3NS, POINTER(tensors.SiteTensor))[:self._netw.nrsites]
        self._T3NS = (tensors.SiteTensor * len(self._T3NS))(*self._T3NS)

    def init_operators(self):
        initop = libt3ns.init_operators
        initop.argtypes = [POINTER(c_void_p), POINTER(tensors.SiteTensor), c_bool]

        if not hasattr(self, '_rOps') or self._rOps is None:
            self._rOps = c_void_p(None)
        else:
            self._rOps = cast(self._rOps, c_void_p)

        initop(byref(self._rOps), self._T3NS, False)
        self._rOps = \
            cast(self._rOps, POINTER(tensors.ROperators))[:self._netw.nrbonds]
        self._rOps = (tensors.ROperators * len(self._rOps))(*self._rOps)

    def execute_optimization(self, D, saveloc=None, verbosity=1, **kwargs):
        from sys import stdout
        import ctypes

        scheme = OptScheme(D, **kwargs)
        self._lastD = D[-1] if isinstance(D, list) else D

        execute = libt3ns.execute_optScheme
        execute.argtypes = [
            POINTER(tensors.SiteTensor),
            POINTER(tensors.ROperators),
            POINTER(OptScheme),
            c_char_p,
            c_int,
            POINTER(c_int),
            c_int
        ]
        execute.restype = c_double
        saveloc = saveloc if saveloc is None else saveloc.encode('utf8')
        # Flushing python
        stdout.flush()

        energy = execute(self._T3NS, self._rOps, byref(scheme), saveloc, 0,
                         POINTER(c_int)(), verbosity)
        libc = ctypes.CDLL(None)
        c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')
        # Flushing C
        libc.fflush(c_stdout)
        return energy

    def disentangle(self, D=None, verbosity=0, **kwargs):
        """Disentangles the network.
        """
        from sys import stdout
        if D is None:
            D = self._lastD
        scheme = DisentScheme(D, **kwargs)

        disent = libt3ns.disentangle_state
        disent.argtypes = [
            POINTER(tensors.SiteTensor),
            POINTER(DisentScheme),
            c_int
        ]
        disent.restype = c_double

        entanglement = disent(self._T3NS, byref(scheme), verbosity)
        stdout.flush()
        self._netw = netw.Network(get_global=True)
        self._bookkeeper = bookkeeper.Bookkeeper.in_dll(libt3ns, "bookie")
        for rops in self._rOps:
            rops.delete()
        self._rOps = None

        libt3ns.clear_instructions()
        libt3ns.reinit_hamiltonian()

        return entanglement

    def singular_values(self):
        """Returns the singular values for the different bonds in network.
        """
        from pyT3NS.c_stdout_filter import stdout_redirector
        import io
        print_sval = libt3ns.print_singular_values_wav
        print_sval.argtypes = [POINTER(tensors.SiteTensor)]

        f = io.BytesIO()
        with stdout_redirector(f):
            print_sval(self._T3NS)
        outpstring = f.getvalue().decode('utf8')
        svals = [l.split() for l in outpstring.split('\n')[:-1]]
        maxbond = max([len(s) - 1 for s in svals])
        result = numpy.zeros((self._netw.nrbonds, maxbond))

        for l in svals:
            index = int(l[0][:-1])
            sval_bonds = numpy.asarray([float(el) for el in l[1:]])
            result[index][:sval_bonds.size] = sval_bonds
        return result
