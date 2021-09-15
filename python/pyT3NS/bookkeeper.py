from ctypes import cdll, Structure, c_int, POINTER, byref, c_double, c_char_p

libt3ns = cdll.LoadLibrary("libT3NS.so")
MAX_SYMMETRIES = c_int.in_dll(libt3ns, "maxsymmetries").value


def translate_symmgroup(symmetry):
    which_symmgroup = libt3ns.which_symmgroup
    which_symmgroup.argtypes = [c_char_p, POINTER(c_int)]
    which_symmgroup.restype = c_int

    if hasattr(symmetry, '__iter__') and not isinstance(symmetry, str):
        result = []
        for s in symmetry:
            r = c_int(0)
            if which_symmgroup(s.encode('utf-8'), byref(r)) != 1:
                raise ValueError(f'{s} is not a value symmetry')
            result.append(r.value)
        return result
    else:
        r = c_int(0)
        s = symmetry
        if which_symmgroup(s.encode('utf-8'), byref(r)) != 1:
            raise ValueError(f'{s} is not a value symmetry')
        return r.value


def translate_irrep(irrep, symmetry):
    def tl_one_irrep(ir, sym):
        if not isinstance(sym, int):
            sym = translate_symmgroup(sym)
        if isinstance(ir, int):
            return ir
        else:
            val = c_int(0)
            which_irrep = libt3ns.which_irrep
            which_irrep.argtypes = [c_char_p, c_int, POINTER(c_int)]
            which_irrep.restype = c_int
            assert isinstance(ir, str)

            if which_irrep(ir.encode('utf-8'), sym, byref(val)) != 1:
                raise ValueError(f'{ir} is an invalid irrep for {sym}')
            return val.value

    if hasattr(irrep, '__iter__') and not isinstance(irrep, str):
        if hasattr(symmetry, '__iter__') and not isinstance(symmetry, str):
            return [tl_one_irrep(ir, sym) for ir, sym in zip(irrep, symmetry)]
        else:
            return [tl_one_irrep(ir, symmetry) for ir in irrep]
    else:
        return tl_one_irrep(irrep, symmetry)


class Symsecs(Structure):
    _fields_ = [
        ("bond", c_int),
        ("nrSecs", c_int),
        ("irreps", POINTER(c_int * MAX_SYMMETRIES)),
        ("fcidims", POINTER(c_double)),
        ("dims", POINTER(c_int)),
        ("totaldims", c_int)
    ]


class Bookkeeper(Structure):
    _fields_ = [
        ("nrSyms", c_int),
        ("sgs", (c_int * MAX_SYMMETRIES)),
        ("target_state", (c_int * MAX_SYMMETRIES)),
        ("nr_bonds", c_int),
        ("v_symsecs", POINTER(Symsecs)),
        ("psites", c_int),
        ("p_symsecs", POINTER(Symsecs))
    ]

    def init_bookkeeper(self, pbookie=None, maxD=500, mstates=2):
        preparebookkeeper = libt3ns.preparebookkeeper
        preparebookkeeper.argtypes = [POINTER(Bookkeeper), c_int, c_int, c_int,
                                      POINTER(c_int)]
        changedSS = c_int(0)
        preparebookkeeper(None if pbookie is None else byref(pbookie), maxD,
                          1, mstates, byref(changedSS))
        return changedSS.value

    def __str__(self):
        from io import StringIO
        from contextlib import redirect_stdout
        print_bookkeeper = libt3ns.print_bookkeeper
        print_bookkeeper.argtypes = [POINTER(Bookkeeper)]
        f = StringIO()
        with redirect_stdout(f):
            print_bookkeeper(byref(self), 0)
        return f.getvalue()

    def fill_symmetry_and_target(self, symmetries, target):
        """Fills in the passed symmetry and target
        """
        self.nrSyms = len(symmetries)
        symmetries = translate_symmgroup(symmetries)
        self.sgs = (c_int * MAX_SYMMETRIES)(*symmetries)
        self.target_state = (c_int * MAX_SYMMETRIES)(
            *translate_irrep(target, symmetries)
        )

    def shallow_copy(self):
        func = libt3ns.shallow_copy_bookkeeper
        func.argtypes = [POINTER(Bookkeeper)]
        func.restype = Bookkeeper
        return func(byref(self))
