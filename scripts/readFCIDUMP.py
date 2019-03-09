from warnings import warn
import numpy as np


def getExchange(filename, norbs=None):
    with open(filename) as f:
        line = next(f)
        if not line.strip().startswith('&FCI NORB='):
            raise OSError('Error in FCIDUMP file header')

        words = line[5:].split(',')
        header_info = {}
        for word in words:
            if word.count('=') == 1:
                key, value = word.split('=')
                header_info[key.strip()] = value.strip()
        norbsread = int(header_info['NORB'])
        if norbs is not None and norbsread != norbs:
            raise OSError('Error in FCIDUMP, should be a file with {} \
            orbitals, not {}.'.format(norbsread, norbs))
        norbs = norbsread

        for line in f:
            words = line.split()
            if words[0] == "&END" or words[0] == "/END" or words[0] == "/":
                break

        # read the integrals
        Kij = np.zeros((norbs, norbs))
        for line in f:
            words = line.split()
            if len(words) != 5:
                raise OSError('Expecting 5 fields on each data line in \
                              FCIDUMP')
            value = float(words[0])

            if words[3] != '0':
                ii = int(words[1])-1
                ij = int(words[2])-1
                ik = int(words[3])-1
                il = int(words[4])-1

                # print(value, ii, ij, ik, il)
                # fcidump (i,j,k,l) with permutations
                # the exchange terms are (i,j,j,i)
                # however due to permutational invariance the terms
                # are stored in the FCIDUMP as (i,j,i,j)
                if ii == ik and ij == il and ii >= 0 \
                        and ij >= 0 and ik >= 0 and il >= 0:
                    if Kij[ii, ij] != 0 or Kij[ij, ii] != 0:
                        warn("Duplicate exchange term for orbs {},{} in {}."
                             .format(ii, ij, f.name))
                    Kij[ii, ij] = value
                    Kij[ij, ii] = value
        return Kij


def getOrbsym(filename, norbs=None):
    with open(filename, 'r') as f:
        line = next(f)
        if not line.strip().startswith('&FCI NORB='):
            raise OSError('Error in FCIDUMP file header')

        words = line[5:].split(',')
        header_info = {}
        for word in words:
            if word.count('=') == 1:
                key, value = word.split('=')
                header_info[key.strip()] = value.strip()
        norbsread = int(header_info['NORB'])
        if norbs is not None and norbsread != norbs:
            raise OSError('Error in FCIDUMP, should be a file with {} \
                          orbitals, not {}.'.format(norbsread, norbs))
        norbs = norbsread

        # skip rest of header
        line = next(f)
        if not line.strip().startswith('ORBSYM='):
            warn("No orbital symmetries found in FCIDUMP, will assume trivial \
                 irreps.", "once")
            return np.zeros(norbs, dtype=np.int)

        return np.array(line.replace('=', ',').split(',')[1:norbs + 1],
                        dtype=np.int)
