from ctypes import cdll, Structure, c_int, POINTER, byref, c_double, c_int64


libt3ns = cdll.LoadLibrary("libT3NS.so")
STEPSPECS_MSITES = 4


class SparseBlocks(Structure):
    _fields_ = [
        ("beginblock", POINTER(c_int)),
        ("tel", POINTER(c_double))
    ]


class SiteTensor(Structure):
    _fields_ = [
        ("nrsites", c_int),
        ("sites", (c_int * STEPSPECS_MSITES)),
        ("nrblocks", c_int),
        ("qnumbers", POINTER(c_int64)),
        ("blocks", SparseBlocks)
    ]

    def __str__(self):
        from io import StringIO
        from contextlib import redirect_stdout
        print_siteTensor = libt3ns.print_siteTensor
        print_siteTensor.argtypes = [POINTER(SiteTensor)]
        f = StringIO()
        with redirect_stdout(f):
            print_siteTensor(None, byref(self))
        return f.getvalue()


class ROperators(Structure):
    _fields_ = [
        ("bond", c_int),
        ("is_left", c_int),
        ("P_operator", c_int),
        ("nrhss", c_int),
        ("begin_blocks_of_hss", POINTER(c_int)),
        ("qnumbers", POINTER(c_int64)),
        ("nrops", c_int),
        ("hss_of_ops", POINTER(c_int)),
        ("operators", POINTER(SparseBlocks))
    ]

    def __str__(self):
        from io import StringIO
        from contextlib import redirect_stdout
        print_rOperators = libt3ns.print_rOperators
        print_rOperators.argtypes = [POINTER(ROperators), c_int]
        f = StringIO()
        with redirect_stdout(f):
            print_rOperators(byref(self), 1)
        return f.getvalue()

    def delete(self):
        libt3ns.destroy_rOperators(byref(self))
