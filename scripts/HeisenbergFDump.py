import numpy


"""
2 spins M=0:
    ud  du
ud  -J  J/2
du  J/2 -J

DOCI 2 orbitalen 1 deeltje:
    10                  01
10  Jii + 2Tii + core   Kij
01  Kij                 Jii + 2Tii + core
"""


def printVal(val, i=None, j=None, k=None, ll=None):
    i = 0 if i is None else i + 1
    j = 0 if j is None else j + 1
    k = 0 if k is None else k + 1
    ll = 0 if ll is None else ll + 1
    if not numpy.isclose(val, 0):
        pre = "  " if val > 0 else " "
        print(f"{pre}{val:.16E}{i:4d}{j:4d}{k:4d}{ll:4d}")


def printDump(J, h):
    assert J.shape[0] == J.shape[1] == h.shape[0]
    assert numpy.allclose(J.diagonal(), 0)

    orbs = J.shape[0]
    print(f" &FCI NORB= {orbs},NELEC= {orbs // 2 * 2}, MS2= 0,")
    print("  ORBSYM=" + "1," * orbs)
    print("  ISYM=1,")
    print(" /")
    for i in range(orbs):
        for j in range(i):
            # print first the coulomb interactions [i i j j]
            printVal(J[i, j] * 1.25, i, i, j, j)
        for j in range(i):
            # print the exchange interactions [i j i j]
            printVal(J[i, j] * 0.5, i, j, i, j)

    # print the kinetic terms
    for i in range(orbs):
        printVal(h[i] - sum(J[i, :]), i, i)

    # print the core energy
    printVal(J.sum() / 2 - h.sum())


if __name__ == "__main__":
    # Nearest neighbour on a rectangular lattice
    PBC = True
    Jval = 1
    hval = 0
    W = 8
    H = 8
    h = numpy.array([hval] * W * H)
    J = numpy.zeros((W * H, W * H))
    for i in range(W * H):
        Lneighbour = i - 1
        Rneighbour = i + 1
        Uneighbour = i - W
        Dneighbour = i + W
        if i < W:
            # Upper boundary
            Uneighbour = W * (H - 1) + i if PBC else None
        if i >= W * (H - 1):
            # Lower boundary
            Dneighbour = int(i % (W * (H - 1))) if PBC else None
        if int(i % W) == 0:
            # Left boundary
            Lneighbour = i + W - 1 if PBC else None
        if int(i % W) == W - 1:
            # Right boundary
            Rneighbour = i - W + 1 if PBC else None

        if Lneighbour is not None:
            J[i, Lneighbour] = Jval
        if Rneighbour is not None:
            J[i, Rneighbour] = Jval
        if Uneighbour is not None:
            J[i, Uneighbour] = Jval
        if Dneighbour is not None:
            J[i, Dneighbour] = Jval
    printDump(J, h)
