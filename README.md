T3NS: An implementation of the Three-Legged Tree Tensor Network algorithm
=========================================================================

T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
Copyright (C) 2018-2019 Klaas Gunst <Klaas.Gunst@UGent.be>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Installation
------------

T3NS can be built with CMake and depends on BLAS and LAPACK.

It is parallelized for shared memory architectures with the 
Open Multi-Processing (OpenMP) API.

In your terminal, do:

    > cd /sourcefolder
    > git clone 'https://github.com/klgunst/T3NS.git'
    > cd T3NS
    > mkdir build
    > cd build
    > cmake ..
    > make

To build with MKL and icc:

    > CC=$(which icc) cmake -DMKL=ON ..

To install:

    > sudo make install

Testing the build can be done by:

    > make test

The number of threads used by openMP can be specified by setting the 
`OMP_NUM_THREADS` variable. e.g.:

    > export OMP_NUM_THREADS=4

To see the help for T3NS:

    > T3NS --help

Running a calculation
---------------------

Examples for running calculations are given in the `examples` sub directory.
For running a T3NS calculation, one needs to specify both an input file and a
network file.

Possible options for the input file can be found through `T3NS --help`.
Mandatory options are:

* `networkfile` : The path to the defined network file for the tensor network
* `symm` : The symmetries that should be used.
* `ts` : The irreps of the targeted ground state. For each symmetry in
  `symm` you need one irrep.
* `interaction` : The type of interaction. Possibilities are:
    * If a path to a `.FCIDUMP` file, it will optimize this Hamiltonian.
    * If a path to a `.FCIDUMP` file preceded by `DOCI`, it will optimize
      within the seniority zero subspace.
    * For a nearest neighbor Hubbard calculations (nearest neighbor
      according to the network geometry) use the following format: 
      `NN_HUBBARD (t = ..., U = ...)`.

A selection of the optional options are:

* `D` : The maximal virtual bond dimension
* `SITE_SIZE` : currently supports 1 and 2 site optimization.
* `SWEEPS` : The maximal number of sweeps to be executed.
* `E_CONV` : If this energy difference between sweeps has been reached, the
  current optimization regime is stopped.

The network file is formatted as follows:
    
    NR_SITES = number of branching and physical tensors
    NR_PHYS_SITES = number of physical tensors
    NR_BONDS = number of virtual bonds in the network. At the border of the
    network there is an extra bond which connects the bordering tensor to
    nothing.
    &END
    List of sites. 
    '*' depict branching sites and numbers represent physical sites where the
    number corresponds with the orbital in the FCIDUMP (counting starts from 0).
    &END
    List of bonds.
    The bonds are specified by giving two tensors which it should connect. The
    number 'n' corresponds with the n'th tensor in the previously specified list
    of tensors. 
    To specify the bond connecting the n'th tensor with nothing (and thus
    specify that this tensor is a border), use:
    -1 n
    Every bond starts at a new line and the number of bonds should correspond
    with the previously defined amount.

Once defined the files, a T3NS calculation can be executed by:

    > T3NS inputfile

Intermediate results are stored in the current working directory in hdf5 format
and can be used to continue a calculation through:
    
    > T3NS --continue=T3NSwav.h5 inputfile

In this case, only the optimization scheme defined in the input file will be
used for the continued calculation. Other specified options will be ignored and
instead read from the hdf5 file.

Provided scripts
----------------

In the `scripts` folder, there are some python scripts provided for the
automatic generation of network files and the optimization of the ordering.

One can generate a T3NS  or DMRG network file by executing

    > ./makenetwork.py nr_sites
    > ./makenetwork.py nr_sites DMRG

This will print the network file to `stdin`. Afterwards one can optimize the
orbital ordering by providing a FCIDUMP and the previously generated network
file.
    
    > ./optimizenetwork.py networkfile fcidump

Both scripts have a rudimentary help accessed by the `-h` or `--help` arguments.
