#!/bin/env python3


class Site:
    def __init__(self, kind='Vacuum', nr=0):
        if kind != 'Vacuum' and kind != 'B' and kind != 'P':
            raise ValueError(f"{kind} is not a valid Site.kind")
        self.kind = kind
        self.nr = nr

    def __str__(self):
        if self.kind == 'Vacuum':
            return self.kind
        else:
            return f'{self.kind}{self.nr}'

    def __eq__(self, s):
        return self.kind == s.kind and self.nr == s.nr


class Network:
    def __init__(self, sites=None, layers=0, isMET=True, isDMRG=False):
        self.nrP = 0
        self.nrB = 0
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
                    Premove += [bond[1]]
            for bond in self.bonds:
                if len(Bremove) == toremove:
                    break
                for psit in Premove:
                    if psit == bond[0]:
                        Bremove += [bond[1]]
                        break
            assert len(Premove) == len(Bremove) == toremove

            for i in range(toremove):
                self.cut([Premove[i], Bremove[i]])

            self.remove_redundantbranch()
            self.sitemap = list(range(self.nrP))
        elif sites is not None and isDMRG:
            self.nrP = sites
            self.sites = [Site('P', i) for i in range(sites)]
            self.bonds = [[Site(), self.sites[0]]] + \
                [[self.sites[i - 1], self.sites[i]] for i in range(1, sites)] \
                + [[self.sites[sites - 1], Site()]]
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

    def __str__(self):
        result = "nrP = {}, nrB = {}".format(self.nrP, self.nrB)
        for lel in self.bonds:
            result += "\n" + str(lel[0]) + "\t" + str(lel[1])
        return result

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
            self.nrP = 0
            self.nrB = 0
            for sit in line.split():
                if sit == '*':
                    self.sites.append(Site('B', self.nrB))
                    self.nrB += 1
                else:
                    self.sites.append(Site('P', self.nrP))
                    self.nrP += 1
                    self.sitemap.append(int(sit))
            assert self.nrP == header_info['NR_PHYS_SITES']
            assert self.nrB == header_info['NR_SITES'] - \
                header_info['NR_PHYS_SITES']

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
        nrP = self.nrP
        nrB = self.nrB
        for s in site_ct:
            self.remove_Site(s, site_ct)
        self.adapt_numbering(nrP, nrB)

    def remove_Site(self, s, site_ct):
        self.sites.remove(s)
        if s.kind == 'P':
            self.nrP -= 1
        elif s.kind == 'B':
            self.nrB -= 1
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
            self.nrB += 1

        tojoin.adapt_sites(self.nrP, self.nrB)
        self.nrP = tojoin.nrP
        self.nrB = tojoin.nrB
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
        toadd = {
            'P': addP,
            'B': addB,
            'Vacuum': 0,

        }
        self.nrP += addP
        self.nrB += addB
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

    def make_branch(self, layers, outgoing=False, firstcall=1):
        if firstcall and outgoing:
            p = Site('P', self.nrP)
            self.sites.append(p)
            self.bonds.insert(0, [Site(), p])

        if layers != 1:
            if outgoing:
                b = Site('B', self.nrB)
                self.sites.append(b)
                self.bonds += [[self.bonds[-1][1], b]]
                self.nrP += 1
                self.nrB += 1
                self.make_branch(layers - 1, firstcall=0)
                self.bonds += [[self.bonds[-1][1], b]]

                self.nrP += 1
                p = Site('P', self.nrP)
                self.sites.append(p)
                self.bonds += [[b, p]]
                self.make_branch(layers - 1, outgoing, firstcall=0)
            else:
                self.make_branch(layers - 1, firstcall=0)
                b = Site('B', self.nrB)
                self.sites.append(b)
                self.bonds += [[self.bonds[-1][1], b]]
                self.nrP += 1
                self.nrB += 1

                self.make_branch(layers - 1, firstcall=0)
                self.bonds += [[self.bonds[-1][1], b]]
                self.nrP += 1
                p = Site('P', self.nrP)
                self.sites.append(p)
                self.bonds += [[b, p]]
        else:
            if outgoing:
                self.bonds.append([self.bonds[-1][1], Site()])
                self.nrP += 1
            else:
                p = Site('P', self.nrP)
                self.sites.append(p)
                self.bonds.append([Site(), p])

        if firstcall and not outgoing:
            if self.bonds[-1][1].kind == 'P':
                self.nrP += 1
            elif self.bonds[-1][1].kind == 'B':
                self.nrB += 1
            self.bonds.append([self.bonds[-1][1], Site()])

    def add_branch(self, layers, jointo=None):
        branch = Network()
        branch.make_branch(layers)
        self.join(branch, jointo)

    def remove_redundantbranch(self):
        oldnrB = self.nrB
        oldnrP = self.nrP
        for bd in self.bonds:
            if bd[0] == Site() and bd[1].kind == 'B':
                b = bd[1]
                self.bonds.remove(bd)
                self.sites.remove(b)
                self.nrB -= 1
                for i in range(len(self.bonds)):
                    bond = self.bonds[i]
                    if bond[1] == b:
                        assert self.bonds[i + 1][0] == b
                        self.bonds[i + 1][0] = bond[0]
                        self.bonds.remove(bond)
                        break
        self.adapt_numbering(oldnrP, oldnrB)

    def printnetworkfile(self, filename=None):
        from sys import stdout
        if filename is None:
            file = stdout
        else:
            file = open(filename, 'w')
        assert self.nrB + self.nrP == len(self.sites)
        print("NR_SITES =", self.nrB + self.nrP, file=file)
        print("NR_PHYS_SITES =", self.nrP, file=file)
        print("NR_BONDS =", len(self.bonds), file=file)
        if hasattr(self, 'sweep'):
            print("SWEEP_LENGTH =", len(self.sweep), file=file)
        print("&END", file=file)
        result = ""
        for p in range(self.nrP):
            result += str(self.sitemap[p]) + " "
        print(result + "* " * self.nrB, file=file)
        print("&END", file=file)
        result = ""
        if hasattr(self, 'sweep'):
            for swp in self.sweep:
                result += str(swp) + " "
            print(result, file=file)
            print("&END", file=file)
        toadd = {
            'P': 0,
            'B': self.nrP,
            'Vacuum': -1
        }
        for bond in self.bonds:
            result = str(toadd[bond[0].kind] + bond[0].nr) + "\t"
            result += str(toadd[bond[1].kind] + bond[1].nr)
            print(result, file=file)
        if filename is not None:
            file.flush()
            file.close()

    def calcDistances(self):
        from numpy import ones
        assert self.nrB + self.nrP == len(self.sites)

        nrsites = self.nrB + self.nrP
        distances = ones((nrsites, nrsites)) * -1

        toadd = {
            'P': 0,
            'B': self.nrP,
            'Vacuum': -1
        }
        for bond in self.bonds:
            if bond[1].kind == 'P':
                # when attaching a Physical tensor
                # add distance one to all the previous bonds
                currsiteid = toadd[bond[1].kind] + bond[1].nr
                distances[currsiteid, currsiteid] = 0
                add_one(toadd, distances, currsiteid, bond[0])
            elif bond[1].kind == 'B':
                # When attaching branching tensor:
                # Search two previous physical tensors.
                # distances that are not not None, in both arrays,
                # add +1 for distances[i][currsiteid] and
                # distances[currsiteid][i]
                # And once this is done, take the array distances[currsiteid]
                # Double loop over all elements i,j and do:
                #    if distances is not set jet, set distance to [i] + [j]
                currsiteid = toadd[bond[1].kind] + bond[1].nr
                distances[currsiteid, currsiteid] = 0
                add_one(toadd, distances, currsiteid, bond[0])

                currsitedist = distances[currsiteid, :]
                for i in range(nrsites):
                    if currsitedist[i] == -1:
                        continue
                    for j in range(i + 1, nrsites):
                        if currsitedist[j] == -1 or distances[i, j] != -1:
                            continue
                        assert distances[i, j] == -1
                        assert distances[j, i] == -1
                        distances[i, j] = currsitedist[j] + currsitedist[i]
                        distances[j, i] = currsitedist[j] + currsitedist[i]
        return distances

    def optimize(self, Iij, **kwargs):
        Dij = self.calcDistances()
        psites = [site.nr for site in self.sites if site.kind == 'P']
        Dij = Dij[psites, :][:, psites]
        assert Iij.shape == Dij.shape
        self.sitemap, cost = montecarlo(Iij, Dij, **kwargs)
        return cost


def add_one(toadd, distances, currsiteid, sitetoadd):
    from copy import deepcopy
    prevsiteid = toadd[sitetoadd.kind] + sitetoadd.nr
    previousdistances = deepcopy(distances[prevsiteid, :])
    if sitetoadd == Site():
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
