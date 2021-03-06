from ctypes import cdll, c_int, POINTER, Structure, c_void_p, byref

libt3ns = cdll.LoadLibrary("libT3NS.so")


class cNetwork(Structure):
    _fields_ = [
        ("nr_bonds", c_int),
        ("psites", c_int),
        ("sites", c_int),
        ("bonds", POINTER(c_int * 2)),
        ("sitetoorb", POINTER(c_int)),
        ("nr_left_psites", POINTER(c_int)),
        ("order_psites", POINTER(c_void_p)),
        ("sweeplength", c_int),
        ("sweep", POINTER(c_int))
    ]

    def __str__(self):
        from io import StringIO
        from contextlib import redirect_stdout
        print_network = libt3ns.print_network
        print_network.argtypes = [POINTER(cNetwork)]
        f = StringIO()
        with redirect_stdout(f):
            print_network(byref(self))
        return f.getvalue()


class Site:
    def __init__(self, kind='Vacuum', nr=0):
        if kind != 'Vacuum' and kind != 'B' and kind != 'P':
            raise ValueError(f"{kind} is not a valid Site.kind")
        self.kind = kind
        self.nr = nr

    @property
    def is_physical(self):
        return self.kind == 'P'

    @property
    def is_branching(self):
        return self.kind == 'B'

    @property
    def is_vacuum(self):
        return self.kind == 'Vacuum'

    def __str__(self):
        if self.is_vacuum:
            return 'Vacuum'
        else:
            return f'{self.kind}{self.nr}'

    def __eq__(self, s):
        return self.kind == s.kind and self.nr == s.nr

    def __hash__(self):
        return hash(str(self))


class Network:
    def __init__(self, sites=None, layers=0, isMET=True, isDMRG=False,
                 get_global=False):
        if get_global:
            cnetwork = cNetwork.in_dll(libt3ns, "netw")
            self.sites = []
            for s in cnetwork.sitetoorb[:cnetwork.sites]:
                self.sites.append(Site('B' if s == -1 else 'P',
                                       self.nrB if s == -1 else self.nrP))

            self.sitemap = [
                s for s in cnetwork.sitetoorb[:cnetwork.sites] if s != -1
            ]
            self.sweep = [s for s in cnetwork.sweep[:cnetwork.sweeplength]]

            sites = [Site()] + self.sites
            self.bonds = [[sites[b + 1] for b in bond] for bond in
                          cnetwork.bonds[:cnetwork.nr_bonds]]
        else:
            self.bonds = []
            self.sites = []
            self.sitemap = []
            if sites is not None and not isDMRG:
                from math import ceil, log
                layers = int(ceil(log(ceil(sites / 3 + 1), 2)))

                part1 = Network()
                part1.make_branch(layers)
                part2 = Network()
                part2.make_branch(layers)
                part1.join(part2)
                self.make_branch(layers, True)
                self.join(part1, jointo=self.bonds[0][1])

                toremove = 3 * (2 ** layers - 1) - sites
                assert toremove >= 0

                Premove = []
                Bremove = []
                self.sitemap = list(range(self.nrP))
                for bond in self.bonds:
                    if len(Premove) == toremove:
                        break
                    if bond[0] == Site():
                        Premove.append(bond[1])
                for bond in self.bonds:
                    if len(Bremove) == toremove:
                        break
                    for psit in Premove:
                        if psit == bond[0]:
                            Bremove.append(bond[1])
                            break
                assert len(Premove) == len(Bremove) == toremove

                for i in range(toremove):
                    self.cut([Premove[i], Bremove[i]])

                self.remove_redundantbranch()
                self.sitemap = list(range(self.nrP))
            elif sites is not None and isDMRG:
                self.sites = [Site('P', i) for i in range(sites)]
                sites = [Site()] + self.sites + [Site()]
                self.bonds = [[s1, s2] for s1, s2 in zip(sites, sites[1:])]
                self.sitemap = list(range(self.nrP))
            elif layers != 0 and isMET:
                part1 = Network()
                part1.make_branch(layers)
                part2 = Network()
                part2.make_branch(layers)
                part1.join(part2)
                self.make_branch(layers, True)
                self.join(part1, jointo=self.bonds[0][1])
            elif layers != 0:
                part1 = Network()
                part1.make_branch(layers)
                part2 = Network()
                part2.make_branch(layers)
                part1.join(part2)
                part3 = Network()
                part3.make_branch(1, True)
                part3.join(part1, jointo=part3.bonds[0][1])
                part4 = Network()
                part4.make_branch(layers)
                part3.join(part4)
                self.make_branch(layers, True)
                self.join(part3, jointo=self.bonds[0][1])

    @property
    def nrP(self):
        return sum([s.is_physical for s in self.sites])

    @property
    def nrB(self):
        return sum([s.is_branching for s in self.sites])

    @property
    def nrsites(self):
        return len(self.sites)

    @property
    def nrbonds(self):
        return len(self.bonds)

    @property
    def net_bonds(self):
        toadd = {'P': 0, 'B': self.nrP, 'Vacuum': -1}
        return [[s.nr + toadd[s.kind] for s in b] for b in self.bonds]

    def __str__(self):
        assert self.nrB + self.nrP == self.nrsites
        res = ""
        res += f"NR_SITES = {self.nrsites}\n"
        res += f"NR_PHYS_SITES = {self.nrP}\n"
        res += f"NR_BONDS = {self.nrbonds}\n"
        if hasattr(self, 'sweep'):
            res += f"SWEEP_LENGTH = {len(self.sweep)}\n"
        res += "&END\n"
        for p in range(self.nrP):
            res += str(self.sitemap[p]) + " "
        res += "* " * self.nrB + "\n"
        res += "&END\n"
        if hasattr(self, 'sweep'):
            for swp in self.sweep:
                res += str(swp) + " "
            res += "\n&END\n"
        return res + "".join([f"{b[0]}\t{b[1]}\n" for b in self.net_bonds])

    def readnetworkfile(self, filename):
        with open(filename) as f:
            header_info = {}

            for line in f:
                words = line.split()
                if words[0] == "&END" or words[0] == "/END" or words[0] == "/":
                    break
                key, value = line.split('=')
                header_info[key.strip()] = int(value.strip())

            line = next(f)
            self.sitemap = []
            self.sites = []
            for sit in line.split():
                if sit == '*':
                    self.sites.append(Site('B', self.nrB))
                else:
                    self.sites.append(Site('P', self.nrP))
                    self.sitemap.append(int(sit))

            assert self.nrP == header_info['NR_PHYS_SITES']
            assert self.nrsites == header_info['NR_SITES']

            line = next(f)
            words = line.split()
            assert words[0] == "&END" or words[0] == "/END" or words[0] == "/"

            self.bonds = []
            for line in f:
                words = line.split()
                bond = [None, None]
                if words[0] == '-1':
                    bond[0] = Site()
                else:
                    bond[0] = self.sites[int(words[0])]

                if words[1] == '-1':
                    bond[1] = Site()
                else:
                    bond[1] = self.sites[int(words[1])]

                self.bonds.append(bond)

    def cut(self, bond, first=True):
        ind = self.search_bond(bond)
        site_ct = []
        if ind == -1:
            print("no such bond found.")
            return
        site_ct.append(self.bonds[ind][not first])
        self.bonds[ind][not first] = Site()
        nrP, nrB = self.nrP, self.nrB

        for s in site_ct:
            self.remove_Site(s, site_ct)
        self.adapt_numbering(nrP, nrB)

    def remove_Site(self, s, site_ct):
        self.sites.remove(s)
        for bond in self.bonds:
            if bond[0] == s:
                if bond[1] != Site():
                    site_ct.append(bond[1])
                self.bonds.remove(bond)
            if bond[1] == s:
                if bond[0] != Site():
                    site_ct.append(bond[0])
                self.bonds.remove(bond)

    def search_bond(self, bond):
        cnt = 0
        for bd in self.bonds:
            if bd[0] == bond[0] and bd[1] == bond[1]:
                return cnt
            cnt += 1
        return -1

    def join(self, tojoin, jointo=None):
        assert tojoin.bonds[-1][1] == Site()
        assert self.bonds[-1][1] == Site()

        if jointo is None:
            self.bonds.pop()
            jointo = Site('B', self.nrB)
            self.sites.append(jointo)
            self.bonds.append([self.bonds[-1][1], jointo])
            self.bonds.append([Site(), jointo])
            self.bonds.append([jointo, Site()])

        tojoin.adapt_sites(self.nrP, self.nrB)
        tojoin.bonds.pop()
        indxs = self.search_bond([Site(), jointo])
        jointo = self.bonds[indxs][1]
        tojoin.bonds.append([tojoin.bonds[-1][1], jointo])
        indxs = self.search_bond([Site(), jointo])
        self.bonds.pop(indxs)
        for bond in tojoin.bonds:
            self.bonds.insert(indxs, bond)
            indxs += 1
        self.sites.extend(tojoin.sites)
        self.adapt_numbering()

    def adapt_sites(self, addP, addB):
        toadd = {'P': addP, 'B': addB, 'Vacuum': 0}
        for sit in self.sites:
            sit.nr += toadd[sit.kind]

    def adapt_numbering(self, nrP=None, nrB=None):
        if nrP is None:
            nrP = self.nrP
            nrB = self.nrB
        newnumber = {
            'P': [-1 for i in range(nrP)],
            'B': [-1 for i in range(nrB)],
            'Vacuum': [0]
        }
        currnmbr = {
            'P': 0,
            'B': 0,
            'Vacuum': 0
        }

        for bond in self.bonds:
            if newnumber[bond[0].kind][bond[0].nr] == -1:
                newnumber[bond[0].kind][bond[0].nr] = currnmbr[bond[0].kind]
                currnmbr[bond[0].kind] += 1
            if newnumber[bond[1].kind][bond[1].nr] == -1:
                newnumber[bond[1].kind][bond[1].nr] = currnmbr[bond[1].kind]
                currnmbr[bond[1].kind] += 1
        for site in self.sites:
            site.nr = newnumber[site.kind][site.nr]

    def make_branch(self, layers, outgoing=False, firstcall=True):
        if firstcall and outgoing:
            p = Site('P', self.nrP)
            self.sites.append(p)
            self.bonds.insert(0, [Site(), p])

        if layers != 1:
            if outgoing:
                b = Site('B', self.nrB)
                self.sites.append(b)
                self.bonds += [[self.bonds[-1][1], b]]
                self.make_branch(layers - 1, firstcall=False)
                self.bonds += [[self.bonds[-1][1], b]]

                p = Site('P', self.nrP)
                self.sites.append(p)
                self.bonds += [[b, p]]
                self.make_branch(layers - 1, outgoing, firstcall=False)
            else:
                self.make_branch(layers - 1, firstcall=False)
                b = Site('B', self.nrB)
                self.sites.append(b)
                self.bonds += [[self.bonds[-1][1], b]]

                self.make_branch(layers - 1, firstcall=False)
                self.bonds += [[self.bonds[-1][1], b]]
                p = Site('P', self.nrP)
                self.sites.append(p)
                self.bonds += [[b, p]]
        else:
            if outgoing:
                self.bonds.append([self.bonds[-1][1], Site()])
            else:
                p = Site('P', self.nrP)
                self.sites.append(p)
                self.bonds.append([Site(), p])

        if firstcall and not outgoing:
            self.bonds.append([self.bonds[-1][1], Site()])

    def add_branch(self, layers, jointo=None):
        branch = Network()
        branch.make_branch(layers)
        self.join(branch, jointo)

    def remove_redundantbranch(self):
        oldnrB = self.nrB
        oldnrP = self.nrP
        for bd in self.bonds:
            if bd[0] == Site() and bd[1].is_branching:
                b = bd[1]
                self.bonds.remove(bd)
                self.sites.remove(b)
                for i in range(self.nrbonds):
                    bond = self.bonds[i]
                    if bond[1] == b:
                        assert self.bonds[i + 1][0] == b
                        self.bonds[i + 1][0] = bond[0]
                        self.bonds.remove(bond)
                        break
        self.adapt_numbering(oldnrP, oldnrB)

    def calcDistances(self):
        from numpy import ones
        distances = ones((self.nrsites,) * 2) * -1
        for b, nb in zip(self.bonds, self.net_bonds):
            if b[1].is_physical:
                # when attaching a Physical tensor
                # add distance one to all the previous bonds
                distances[nb[1], nb[1]] = 0
                add_one(distances, nb[1], nb[0])
            elif b[1].is_branching:
                # When attaching branching tensor:
                # Search two previous physical tensors.
                # distances that are not not None, in both arrays,
                # add +1 for distances[i][currsiteid] and
                # distances[currsiteid][i]
                # And once this is done, take the array distances[currsiteid]
                # Double loop over all elements i,j and do:
                #    if distances is not set jet, set distance to [i] + [j]
                distances[nb[1], nb[1]] = 0
                add_one(distances, nb[1], nb[0])

                currsitedist = distances[nb[1], :]
                for i in range(self.nrsites):
                    if currsitedist[i] == -1:
                        continue
                    for j in range(i + 1, self.nrsites):
                        if currsitedist[j] == -1 or distances[i, j] != -1:
                            continue
                        assert distances[i, j] == -1
                        assert distances[j, i] == -1
                        distances[i, j] = currsitedist[j] + currsitedist[i]
                        distances[j, i] = currsitedist[j] + currsitedist[i]
        return distances

    def optimize(self, Iij, **kwargs):
        Dij = self.calcDistances()
        psites = [site.nr for site in self.sites if site.is_physical]
        Dij = Dij[psites, :][:, psites]
        assert Iij.shape == Dij.shape
        assert Dij.size == self.nrP ** 2

        self.sitemap, cost = montecarlo(Iij, Dij, **kwargs)
        return cost

    def pass_network(self):
        '''Fills in the static global network structure in the T3NS.so library
        and puts a pointer to that global network in the network structure
        '''
        bonds = ((c_int * 2) * self.nrbonds)(
            *[(c_int * 2)(*bond) for bond in self.net_bonds]
        )
        sitetoorb = self.sitemap + [-1] * self.nrB
        sitetoorb = (c_int * len(sitetoorb))(*sitetoorb)
        if hasattr(self, 'sweep'):
            sweep = (c_int * len(self.sweep))(*self.sweep)
            sweeplength = len(self.sweep)
        else:
            sweep = None
            sweeplength = 0

        fillin_network = libt3ns.fillin_network
        fillin_network.argtypes = [c_int, c_int, c_int, POINTER((c_int * 2)),
                                   POINTER(c_int), c_int, POINTER(c_int)]

        fillin_network(self.nrbonds, self.nrP, self.nrsites, bonds, sitetoorb,
                       sweeplength, sweep)

    def plotGraph(self, grid=None):
        """Plots a graph to the current pyplot scope

        Arguments:
            grid: If a tuple (x,y), this will put the physical sites on grid
            with size (x,y). Site with order number `i` will be put on position
            (i % x, i // x).
        """
        import networkx as nx

        G = nx.Graph()
        G.add_nodes_from(self.sites)
        G.add_edges_from([b for b in self.bonds if
                          not b[0].is_vacuum and not b[1].is_vacuum])
        pos = None
        if grid is not None:
            if not isinstance(grid, tuple):
                raise ValueError('grid should be a tuple')
            if grid[0] * grid[1] != self.nrP:
                raise ValueError(
                    'Dimension of grid ({grid[0] * grid[1]}) should correspond'
                    ' with nr physical sites ({self.nrP}).'
                )
            fixed_positions = {}
            for s in self.sites:
                if s.is_physical:
                    nr = self.sitemap[s.nr]
                    fixed_positions[s] = (nr % grid[0], nr // grid[0])
            pos = nx.spring_layout(G, pos=fixed_positions,
                                   fixed=fixed_positions.keys())

        elif self.nrB == 0:
            count = 0
            pos = {}
            prev = None
            for edge in G.edges:
                assert prev is None or prev == edge[0]
                pos[edge[0]] = (count, 0)
                count += 1
                prev = edge[1]
            pos[prev] = (count, 0)

        labels = {}
        for s in self.sites:
            if s.is_physical:
                labels[s] = self.sitemap[s.nr]

        color_map = ['lightblue' if s.is_physical else 'green' for s in G]
        nx.draw(G, pos=pos, node_color=color_map,
                labels=labels, with_labels=True)


def add_one(distances, currsiteid, prevsiteid):
    from copy import deepcopy
    previousdistances = deepcopy(distances[prevsiteid, :])
    if prevsiteid == -1:
        return
    for i, val in enumerate(previousdistances):
        if val != -1:
            assert distances[i, currsiteid] == -1
            assert distances[currsiteid, i] == -1
            distances[i, currsiteid] = val + 1
            distances[currsiteid, i] = val + 1


def costfunction(state, Iij, Dij, η=2):
    """The default cost function.

    Cost is Σ Iij * |i - j|^η
    """
    from numpy import power, argsort
    nsites = len(state)
    if costfunction.Dij_η is None:
        costfunction.Dij_η = power(Dij, η)
    if costfunction.Dij_sort_id is None:
        costfunction.Dij_sort_id = argsort(Dij, axis=1)

    cost = 0
    for i in range(nsites):
        for j in range(i + 1, nsites):
            cost += 2 * Iij[state[i], state[j]] * costfunction.Dij_η[i][j]
    return cost


costfunction.Dij_η = None
costfunction.Dij_sort_id = None


def doswap(nsites, swapspaces=None):
    from random import sample
    if swapspaces is None:
        return sample(range(nsites), 2)
    else:
        swap = [sample(range(nsites), 1)[0]] * 2
        for liss in swapspaces:
            if swap[0] in liss:
                if len(liss) == 1:
                    return doswap(nsites, swapspaces)
                while swap[0] == swap[1]:
                    swap[1] = sample(liss, 1)[0]
                return swap


def montecarlo(Iij, Dij, iterations=1000, β=1, initstate=None, **kwargs):
    from random import sample, random
    from math import exp
    assert Iij.shape == Dij.shape and Iij.shape[0] == Iij.shape[1]

    nsites = Iij.shape[0]
    if initstate is None:
        initstate = sample(range(nsites), nsites)
    else:
        assert len(initstate) == nsites

    currstate = [el for el in initstate]
    beststate = [el for el in initstate]
    costfunction.Dij_η = None
    costfunction.Dij_sort_id = None
    oldcost = costfunction(currstate, Iij, Dij, **kwargs)
    bestcost = oldcost

    for it in range(iterations):
        while True:
            swapper = doswap(nsites, **kwargs)
            # swap them
            currstate[swapper[0]], currstate[swapper[1]] = \
                currstate[swapper[1]], currstate[swapper[0]]
            newcost = costfunction(currstate, Iij, Dij, **kwargs)

            if newcost < oldcost or random() < exp(-β * (newcost - oldcost)):
                # swap is accepted
                oldcost = newcost
                # new best state
                if bestcost > oldcost:
                    beststate = [i for i in currstate]
                    bestcost = oldcost
                break
            # swap is not accepted, redo swap
            currstate[swapper[0]], currstate[swapper[1]] = \
                currstate[swapper[1]], currstate[swapper[0]]
    return beststate, bestcost
