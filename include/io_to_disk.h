#pragma once
#include <hdf5.h>

#include "siteTensor.h"
#include "rOperators.h"

#define H5_DEFAULT_LOCATION "./"

void write_to_disk(const char * hdf5_loc, const struct siteTensor * const T3NS, 
                   const struct rOperators * const ops);

void read_from_disk(const char filename[], struct siteTensor ** const T3NS, 
                    struct rOperators ** const ops);

void ints_write_attribute(const hid_t group_id, const char atrname[], 
                          const int * const atr, const hsize_t size);

void ints_read_attribute(const hid_t id, const char atrname[], int * const atr);

void doubles_write_attribute(const hid_t id, const char atrname[], 
                             const double * const atr, const hsize_t size);

void doubles_read_attribute(const hid_t id, const char atrname[], 
                            double * const atr);

void ints_write_dataset(const hid_t id, const char datname[],
                        const int * const dat, const hsize_t size);

void ints_read_dataset(const hid_t id, const char datname[], int * const dat);

void doubles_write_dataset(const hid_t id, const char datname[],
                           const double * const dat, const hsize_t size);

void doubles_read_dataset(const hid_t id, const char datname[],
                          double * const dat);

void EL_TYPE_write_dataset(const hid_t id, const char datname[],
                           const EL_TYPE * const dat, const hsize_t size);

void EL_TYPE_read_dataset(const hid_t id, const char datname[],
                          EL_TYPE * const dat);
