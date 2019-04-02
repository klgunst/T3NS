#!/bin/env python3
from random import random, sample
import math
import numpy as np
from copy import deepcopy


# I should include an extra monte carlo step where only valid swap is between
# two sites with same irreps.
# Maybe a good strategy than is First optimize like that by putting irreps in
# trees Afterwards, I can do a (quite strict) reshuffle over whole tree.


# Print iterations progress
def printProgressBar(iteration, total, prefix='', suffix='',
                     decimals=1, length=50, fill='█'):
    """Call in a loop to create terminal progress bar
    From stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent
                                  complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration /
                                                            float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end='\r')
    # Print New Line on Complete
    if iteration == total:
        print()


def defaultcost(state, Iij, Dij, η=2, Irreps=None, γ=1):
    """The default cost function.
    Cost is Σ Iij * |i - j|^η + Σ γ * |i - closest site of same irrep|^η
    """
    nsites = len(state)
    if Iij.shape != (nsites, nsites) or Dij.shape != (nsites, nsites):
        raise Exception("Length of state is not compatible with shapes of Iij \
                        or Dij")

    if defaultcost.Dij_η is None:
        defaultcost.Dij_η = np.power(Dij, η)
    if defaultcost.Dij_sort_id is None:
        defaultcost.Dij_sort_id = np.argsort(Dij, axis=1)

    cost = 0
    for i in range(nsites):
        for j in range(i + 1, nsites):
            cost += 2 * Iij[state[i], state[j]] * defaultcost.Dij_η[i][j]
    if Irreps is None or γ == 0:
        return cost
    for i in range(nsites):
        irrepi = Irreps[state[i]]
        mindist = None
        for currsite in defaultcost.Dij_sort_id[i][1:]:
            if irrepi == Irreps[state[currsite]]:
                mindist = defaultcost.Dij_η[i][currsite]
                break
        if mindist is not None:
            cost += γ * mindist
    return cost


defaultcost.Dij_η = None
defaultcost.Dij_sort_id = None


def montecarlo(amount=1, iterations=1000, β=1, initstate=None,
               costfunction=None, nsites=None, progressShow=False,
               swapspaces=None, **kwargs):
    bestcost = 0
    for i in range(amount):
        newstate, newcost = singlemonte(iterations, β, initstate, costfunction,
                                        nsites, progressShow, swapspaces,
                                        **kwargs)
        if costfunction is None and 'Irreps' in kwargs:
            if 'η' in kwargs:
                irreplesscosts = defaultcost(newstate, kwargs['Iij'],
                                             kwargs['Dij'], kwargs['η'])
            else:
                irreplesscosts = defaultcost(newstate, kwargs['Iij'],
                                             kwargs['Dij'])
            print("iteration {} found full cost {:.4g} \
                  (without irreps {:.4g}).".format(i, newcost, irreplesscosts))
        else:
            print("iteration {} found cost {:.4g}.".format(i, newcost))
        print(newstate)

        if i == 0 or bestcost > newcost:
            beststate, bestcost = newstate, newcost

    defaultcost.Dij_η = None
    defaultcost.Dij_sort_id = None
    return beststate, bestcost


def doswap(nsites, swapspaces):
    if swapspaces is None:
        return sample(range(nsites), 2)
    else:
        swap = [] + sample(range(nsites), 1)
        swap += [swap[0]]
        for liss in swapspaces:
            if swap[0] in liss:
                while swap[0] == swap[1]:
                    swap[1] = sample(liss, 1)[0]
                return swap
        print('WHUT')


def singlemonte(iterations=1000, β=1, initstate=None, costfunction=None,
                nsites=None, progressShow=False, swapspaces=None, **kwargs):
    if initstate is None and nsites is None:
        raise Exception("Number of sites or the initial state should be \
                        initialized.")
    if initstate is not None and nsites is not None \
            and len(initstate) != nsites:
        raise Exception("The number of sites specified should correspond with \
                        the length of the initial state")

    if initstate is None:
        initstate = sample(range(nsites), nsites)
    else:
        nsites = len(initstate)

    if costfunction is None:
        costfunction = defaultcost

    currstate = deepcopy(initstate)
    beststate = deepcopy(initstate)
    oldcost = costfunction(currstate, **kwargs)
    print('First cost', oldcost)
    bestcost = oldcost
    i = 0
    modul = iterations // 500
    if modul == 0:
        modul = 1

    while i < iterations:
        # select two elements
        if progressShow and i % modul == 0:
            printProgressBar(i, iterations,
                             suffix="Cost: {0:9.4g}".format(bestcost))

        swapper = doswap(nsites, swapspaces)
        # swap them
        currstate[swapper[0]], currstate[swapper[1]] = \
            currstate[swapper[1]], currstate[swapper[0]]
        newcost = costfunction(currstate, **kwargs)

        if newcost < oldcost or random() < math.exp(-β * (newcost - oldcost)):
            i += 1
            # swap is accepted
            oldcost = newcost
            # new best state
            if bestcost > oldcost:
                beststate = deepcopy(currstate)
                bestcost = oldcost

        else:
            # swap is not accepted, redo swap
            currstate[swapper[0]], currstate[swapper[1]] = \
                currstate[swapper[1]], currstate[swapper[0]]

    # end progress bar
    printProgressBar(iterations, iterations,
                     suffix="Cost: {0:9.4g}".format(bestcost))

    return beststate, bestcost


if __name__ == "__main__":
    from sys import argv
    from os import path
    from makenetwork import network, kind
    from readFCIDUMP import getExchange, getOrbsym

    def printhelp():
        print("Usage: " + path.basename(__file__) + " NETWORK FCIDUMP")
        print("Optimizes the orbital ordering through monte-carlo.")
        print("Orbitals which are connected by large exchange matrix elements")
        print("are placed close in the network.")
        print("\tNETWORK is a path to a network file, the order of the sites")
        print("\tis overwritten in this file")
        print("\tFCIDUMP is a path to a FCIDUMP file")
        exit(0)

    if len(argv) != 3 or argv[1] == "--help" or argv[1] == "-h":
        printhelp()

    netw = network()
    netw.readnetworkfile(argv[1])

    Dij = netw.calcDistances()
    psites = [site.nr for site in netw.sites if site.kind == kind.P]
    Dij = Dij[psites, :][:, psites]

    Kij = getExchange(argv[2])
    Irreps = getOrbsym(argv[2])

    init = [netw.sitemap[i] for i in psites]

    oldbest, costold = montecarlo(1, iterations=1000, β=0.5,
                                  nsites=Kij.shape[1], progressShow=True,
                                  Iij=Kij, Dij=Dij)

    netw.sitemap = oldbest
    netw.printnetworkfile(argv[1])
