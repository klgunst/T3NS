T3NS: An implementation of the Three-Legged Tree Tensor Network algorithm
=========================================================================

T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
Copyright (C) 2018 Klaas Gunst <Klaas.Gunst@UGent.be>

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
