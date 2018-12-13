/*
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
*/
#pragma once
#include <hdf5.h>

#include "siteTensor.h"
#include "rOperators.h"

#define H5_DEFAULT_LOCATION "./"

enum hdf5type { THDF5_INT, THDF5_DOUBLE, THDF5_EL_TYPE, THDF5_QN_TYPE };

void write_to_disk(const char * hdf5_loc, const struct siteTensor * const T3NS, 
                   const struct rOperators * const ops);

void read_from_disk(const char filename[], struct siteTensor ** const T3NS, 
                    struct rOperators ** const ops);

void write_dataset(hid_t id, const char datname[], const void * dat, 
                   hsize_t size, enum hdf5type kind);

void read_dataset(hid_t id, const char datname[], void * dat);

void write_attribute(hid_t group_id, const char atrname[], const void * atr, 
                     hsize_t size, enum hdf5type kind);

void read_attribute(hid_t id, const char atrname[], void * atr);
