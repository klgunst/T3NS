#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hdf5.h>
#include <unistd.h>

#include "io_to_disk.h"
#include "sparseblocks.h"
#include "siteTensor.h"
#include "rOperators.h"
#include "bookkeeper.h"
#include "symsecs.h"
#include "network.h"
#include "macros.h"
#include "debug.h"
#include "hamiltonian.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void write_bookkeeper_to_disk(const hid_t file_id);

static void read_bookkeeper_from_disk(hid_t file_id);

static void write_symsec_to_disk(const hid_t id, const struct symsecs * const 
                                 ssec, const int nmbr);

static void read_symsec_from_disk(const hid_t id, struct symsecs * const ssec, 
                                  const int nmbr);

static void write_T3NS_to_disk(const hid_t file_id, 
                               const struct siteTensor * const T3NS);

static void read_T3NS_from_disk(const hid_t file_id, 
                                struct siteTensor ** const T3NS);

static void write_siteTensor_to_disk(const hid_t id, const struct siteTensor * 
                                     const tens, const int nmbr);

static void read_siteTensor_from_disk(const hid_t id, struct siteTensor * 
                                     const tens, const int nmbr);

static void write_sparseblocks_to_disk(const hid_t id, 
                                       const struct sparseblocks * block, 
                                       const int nrblocks, const int nmbr);

static void read_sparseblocks_from_disk(const hid_t id, 
                                       struct sparseblocks * block, 
                                       const int nrblocks, const int nmbr);

static void write_rOps_to_disk(const hid_t id,
                               const struct rOperators * const rOps);

static void read_rOps_from_disk(const hid_t id,
                               struct rOperators ** const rOps);

static void write_rOperator_to_disk(const hid_t id,
                                    const struct rOperators * const rOp,
                                    const int nmbr);

static void read_rOperator_from_disk(const hid_t id,
                                    struct rOperators * const rOp,
                                    const int nmbr);

static void write_network_to_disk(const hid_t id);

static void read_network_from_disk(const hid_t id);

/* ========================================================================== */

static void make_h5f_name(const char * hdf5_loc, const char hdf5file[], 
                          const int size, char hdf5_resulting[size])
{
        strncpy(hdf5_resulting, hdf5_loc, size - 1);
        hdf5_resulting[size - 1] = '\0';

        int len = strlen(hdf5_resulting);
        if (hdf5_resulting[len - 1] != '/' && hdf5file[0] != '/') {
                hdf5_resulting[len] = '/';
                hdf5_resulting[len + 1] = '\0';
                ++len;
        }
        if (hdf5_resulting[len - 1] == '/' && hdf5file[0] == '/') {
                hdf5_resulting[len - 1] = '\0';
                --len;
        }

        strncat(hdf5_resulting, hdf5file, size - len);
        hdf5_resulting[size - 1] = '\0';
}

void write_to_disk(const char * hdf5_loc, const struct siteTensor * const T3NS, 
                   const struct rOperators * const ops)
{
        if (hdf5_loc == NULL)
                return;

        hid_t file_id;
        const char hdf5nam[] = "T3NScalc.h5";
        const int size = 255;
        char hdf5file[size];
        make_h5f_name(hdf5_loc, hdf5nam, size, hdf5file);

        file_id = H5Fcreate(hdf5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        write_network_to_disk(file_id);
        write_bookkeeper_to_disk(file_id);
        write_hamiltonian_to_disk(file_id);
        write_T3NS_to_disk(file_id, T3NS);
        write_rOps_to_disk(file_id, ops);

        H5Fclose(file_id);
}

void read_from_disk(const char filename[], struct siteTensor ** const T3NS, 
                    struct rOperators ** const ops)
{
        if (access(filename, F_OK) != 0) {
                fprintf(stderr, "Error in %s: Can not read from disk.\n"
                        "%s was not found.\n", __func__, filename);
                exit(EXIT_FAILURE);
        }

        hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

        read_network_from_disk(file_id);
        read_bookkeeper_from_disk(file_id);
        read_hamiltonian_from_disk(file_id);
        read_T3NS_from_disk(file_id, T3NS);
        read_rOps_from_disk(file_id, ops);

        H5Fclose(file_id);
}

void ints_write_attribute(const hid_t group_id, const char atrname[], 
                          const int * const atr, const hsize_t size)
{
        if (atr == NULL)
                return;

        const hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        const hid_t attribute_id = H5Acreate(group_id, atrname, H5T_STD_I32LE, 
                                       dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

        H5Awrite(attribute_id, H5T_STD_I32LE, atr);

        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);
}

void ints_read_attribute(const hid_t id, const char atrname[], int * const atr)
{
        hid_t attribute_id = H5Aopen(id, atrname, H5P_DEFAULT);
        H5Aread(attribute_id, H5T_STD_I32LE, atr);
        H5Aclose(attribute_id);
}

void doubles_write_attribute(const hid_t id, const char atrname[], 
                             const double * const atr, const hsize_t size)
{
        if (atr == NULL)
                return;

        const hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        const hid_t attribute_id = H5Acreate(id, atrname, H5T_IEEE_F64LE, 
                                       dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

        H5Awrite(attribute_id, H5T_IEEE_F64LE, atr);

        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);
}

void doubles_read_attribute(const hid_t id, const char atrname[], 
                            double * const atr)
{
        hid_t attribute_id = H5Aopen(id, atrname, H5P_DEFAULT);
        H5Aread(attribute_id, H5T_IEEE_F64LE, atr);
        H5Aclose(attribute_id);
}

void ints_write_dataset(const hid_t id, const char datname[],
                        const int * const dat, const hsize_t size)
{
        if (dat == NULL)
                return;

        const hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        const hid_t dataset_id = H5Dcreate(id, datname, H5T_STD_I32LE, 
                                           dataspace_id, H5P_DEFAULT, 
                                           H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite (dataset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);

        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
}

void ints_read_dataset(const hid_t id, const char datname[], int * const dat)
{
        hid_t dataset_id = H5Dopen(id, datname, H5P_DEFAULT);
        H5Dread(dataset_id, H5T_STD_I32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

void doubles_write_dataset(const hid_t id, const char datname[],
                           const double * const dat, const hsize_t size)
{
        if (dat == NULL)
                return;

        const hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        const hid_t dataset_id = H5Dcreate(id, datname, H5T_IEEE_F64LE, 
                                           dataspace_id, H5P_DEFAULT, 
                                           H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite (dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

void doubles_read_dataset(const hid_t id, const char datname[],
                          double * const dat)
{
        hid_t dataset_id = H5Dopen(id, datname, H5P_DEFAULT);
        H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

void EL_TYPE_write_dataset(const hid_t id, const char datname[],
                           const EL_TYPE * const dat, const hsize_t size)
{
        if (dat == NULL)
                return;

        const hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        const hid_t dataset_id = H5Dcreate(id, datname, EL_TYPE_H5, 
                                           dataspace_id, H5P_DEFAULT, 
                                           H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, EL_TYPE_H5, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

void EL_TYPE_read_dataset(const hid_t id, const char datname[],
                          EL_TYPE * const dat)
{
        hid_t dataset_id = H5Dopen(id, datname, H5P_DEFAULT);
        H5Dread(dataset_id, EL_TYPE_H5, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void QN_TYPE_write_dataset(const hid_t group_id, const char datname[],
                                  const QN_TYPE * const dat, const hsize_t size)
{
        if (dat == NULL)
                return;

        const hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        const hid_t dataset_id = H5Dcreate(group_id, datname, QN_TYPE_H5, 
                                           dataspace_id, H5P_DEFAULT, 
                                           H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(dataset_id, QN_TYPE_H5, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

static void QN_TYPE_read_dataset(const hid_t id, const char datname[],
                          QN_TYPE * const dat)
{
        hid_t dataset_id = H5Dopen(id, datname, H5P_DEFAULT);
        H5Dread(dataset_id, QN_TYPE_H5, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

static void write_bookkeeper_to_disk(const hid_t file_id)
{
        const hid_t group_id = H5Gcreate(file_id, "/bookkeeper", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "nrSyms", &bookie.nrSyms, 1);
        ints_write_attribute(group_id, "sgs", (int *) bookie.sgs, bookie.nrSyms);
        ints_write_attribute(group_id, "target_state", bookie.target_state, 
                             bookie.nrSyms);
        ints_write_attribute(group_id, "nr_bonds", &bookie.nr_bonds, 1);

        for (int i = 0 ; i < bookie.nr_bonds; ++i) 
                write_symsec_to_disk(group_id, &bookie.list_of_symsecs[i], i);

        H5Gclose(group_id);
}

static void read_bookkeeper_from_disk(const hid_t file_id)
{
        const hid_t group_id = H5Gopen(file_id, "/bookkeeper", H5P_DEFAULT);

        ints_read_attribute(group_id, "nrSyms", &bookie.nrSyms);

        bookie.sgs = safe_malloc(bookie.nrSyms, enum symmetrygroup);
        ints_read_attribute(group_id, "sgs", (int *) bookie.sgs);

        bookie.target_state = safe_malloc(bookie.nrSyms, int);
        ints_read_attribute(group_id, "target_state", bookie.target_state);

        ints_read_attribute(group_id, "nr_bonds", &bookie.nr_bonds);

        bookie.list_of_symsecs = safe_malloc(bookie.nr_bonds, struct symsecs);
        for (int i = 0 ; i < bookie.nr_bonds; ++i) 
                read_symsec_from_disk(group_id, &bookie.list_of_symsecs[i], i);

        H5Gclose(group_id);
}

static void write_symsec_to_disk(const hid_t id, const struct symsecs * const 
                                 ssec, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./symsec_%d", nmbr);
        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "nrSecs", &ssec->nrSecs, 1);
        ints_write_attribute(group_id, "totaldims", &ssec->totaldims, 1);
        ints_write_dataset(group_id, "./dims", ssec->dims, ssec->nrSecs);
        ints_write_dataset(group_id, "./irreps", ssec->irreps, 
                           ssec->nrSecs * bookie.nrSyms);
        doubles_write_dataset(group_id, "./fcidims", ssec->fcidims, ssec->nrSecs);
        H5Gclose(group_id);
}

static void read_symsec_from_disk(const hid_t id, struct symsecs * const ssec, 
                                  const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./symsec_%d", nmbr);
        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        ints_read_attribute(group_id, "nrSecs", &ssec->nrSecs);
        ints_read_attribute(group_id, "totaldims", &ssec->totaldims);

        ssec->dims = safe_malloc(ssec->nrSecs, int);
        ints_read_dataset(group_id, "./dims", ssec->dims);

        ssec->irreps = safe_malloc(ssec->nrSecs, *ssec->irreps);
        ints_read_dataset(group_id, "./irreps", ssec->irreps);
        
        ssec->fcidims = safe_malloc(ssec->nrSecs, double);
        doubles_read_dataset(group_id, "./fcidims", ssec->fcidims);
        H5Gclose(group_id);
}

static void write_T3NS_to_disk(const hid_t file_id, 
                               const struct siteTensor * const T3NS)
{
        const hid_t group_id = H5Gcreate(file_id, "/T3NS", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "nrSites", &netw.sites, 1);

        for(int i = 0 ; i < netw.sites; ++i) 
                write_siteTensor_to_disk(group_id, &T3NS[i], i);

        H5Gclose(group_id);
}

static void read_T3NS_from_disk(const hid_t file_id, 
                                struct siteTensor ** const T3NS)
{
        const hid_t group_id = H5Gopen(file_id, "/T3NS", H5P_DEFAULT);
        int nrsit;

        ints_read_attribute(group_id, "nrSites", &nrsit);
        assert(nrsit == netw.sites);

        *T3NS = safe_malloc(nrsit, struct siteTensor);
        for(int i = 0 ; i < netw.sites; ++i) 
                read_siteTensor_from_disk(group_id, &(*T3NS)[i], i);

        H5Gclose(group_id);
}

static void write_siteTensor_to_disk(const hid_t id, const struct siteTensor * 
                                     const tens, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./tensor_%d", nmbr);
        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "nrsites", &tens->nrsites, 1);
        ints_write_attribute(group_id, "sites", tens->sites, tens->nrsites);
        ints_write_attribute(group_id, "nrblocks", &tens->nrblocks, 1);
        QN_TYPE_write_dataset(group_id, "./qnumbers", 
                              tens->qnumbers, tens->nrblocks * tens->nrsites);
        write_sparseblocks_to_disk(group_id, &tens->blocks, tens->nrblocks, 0);
        H5Gclose(group_id);
}

static void read_siteTensor_from_disk(const hid_t id, struct siteTensor * 
                                     const tens, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./tensor_%d", nmbr);
        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        ints_read_attribute(group_id, "nrsites", &tens->nrsites);

        tens->sites = safe_malloc(tens->nrsites, int);
        ints_read_attribute(group_id, "sites", tens->sites);
        ints_read_attribute(group_id, "nrblocks", &tens->nrblocks);

        tens->qnumbers = safe_malloc(tens->nrblocks * tens->nrsites, QN_TYPE);
        QN_TYPE_read_dataset(group_id, "./qnumbers", tens->qnumbers);
        read_sparseblocks_from_disk(group_id, &tens->blocks, tens->nrblocks, 0);
        H5Gclose(group_id);
}

static void write_sparseblocks_to_disk(const hid_t id, 
                                       const struct sparseblocks * block, 
                                       const int nrblocks, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./block_%d", nmbr);
        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "nrBlocks", &nrblocks, 1);
        ints_write_dataset(group_id, "./beginblock", 
                           block->beginblock, nrblocks + 1);
        EL_TYPE_write_dataset(group_id, "./tel", block->tel, 
                              block->beginblock[nrblocks]);

        H5Gclose(group_id);
}

static void read_sparseblocks_from_disk(const hid_t id, 
                                       struct sparseblocks * block, 
                                       const int nrblocks, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./block_%d", nmbr);
        int bloccount;

        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        ints_read_attribute(group_id, "nrBlocks", &bloccount);
        assert(bloccount == nrblocks);

        block->beginblock = safe_malloc(nrblocks + 1, int);
        ints_read_dataset(group_id, "./beginblock", block->beginblock);

        block->tel = safe_malloc(block->beginblock[nrblocks], EL_TYPE);
        EL_TYPE_read_dataset(group_id, "./tel", block->tel);

        H5Gclose(group_id);
}

static void write_rOps_to_disk(const hid_t id,
                               const struct rOperators * const rOps)
{
        const hid_t group_id = H5Gcreate(id, "/rOps", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "nrOps", &netw.nr_bonds, 1);

        for (int i = 0 ; i < netw.nr_bonds; ++i)
                write_rOperator_to_disk(group_id, &rOps[i], i);

        H5Gclose(group_id);
}

static void read_rOps_from_disk(const hid_t id,
                               struct rOperators ** const rOps)
{
        const hid_t group_id = H5Gopen(id, "/rOps", H5P_DEFAULT);
        int nrbonds;

        ints_read_attribute(group_id, "nrOps", &nrbonds);
        assert(nrbonds == netw.nr_bonds);

        *rOps = safe_malloc(nrbonds, struct rOperators);
        for (int i = 0 ; i < netw.nr_bonds; ++i)
                read_rOperator_from_disk(group_id, &(*rOps)[i], i);

        H5Gclose(group_id);
}

static void write_rOperator_to_disk(const hid_t id,
                                    const struct rOperators * const rOp,
                                    const int nmbr)
{
        char buffer[255];
        if (rOp->is_left == -1 || rOp->bond_of_operator == -1)
                return;

        sprintf(buffer, "./rOperator_%d", nmbr);

        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "bond_of_operator", &rOp->bond_of_operator, 1);
        ints_write_attribute(group_id, "is_left", &rOp->is_left, 1);
        ints_write_attribute(group_id, "P_operator", &rOp->P_operator, 1);
        ints_write_attribute(group_id, "nrhss", &rOp->nrhss, 1);

        ints_write_dataset(group_id, "./begin_blocks_of_hss",
                           rOp->begin_blocks_of_hss, rOp->nrhss + 1);
        QN_TYPE_write_dataset(group_id, "./qnumbers", rOp->qnumbers, 
                              rOp->begin_blocks_of_hss[rOp->nrhss] *
                              rOperators_give_nr_of_couplings(rOp));

        ints_write_attribute(group_id, "nrops", &rOp->nrops, 1);
        ints_write_dataset(group_id, "./hss_of_ops", rOp->hss_of_ops, rOp->nrops);

        for (int i = 0; i < rOp->nrops; ++i) {
                const int nr_blocks = rOperators_give_nr_blocks_for_operator(rOp, i);
                write_sparseblocks_to_disk(group_id, &rOp->operators[i], nr_blocks, i);
        }
        H5Gclose(group_id);
}

static void read_rOperator_from_disk(const hid_t id,
                                    struct rOperators * const rOp, 
                                    const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./rOperator_%d", nmbr);

        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        ints_read_attribute(group_id, "bond_of_operator", &rOp->bond_of_operator);
        ints_read_attribute(group_id, "is_left", &rOp->is_left);
        ints_read_attribute(group_id, "P_operator", &rOp->P_operator);
        ints_read_attribute(group_id, "nrhss", &rOp->nrhss);

        rOp->begin_blocks_of_hss = safe_malloc(rOp->nrhss + 1, int);
        ints_read_dataset(group_id, "./begin_blocks_of_hss",
                          rOp->begin_blocks_of_hss);

        rOp->qnumbers = safe_malloc(rOp->begin_blocks_of_hss[rOp->nrhss] *
                              rOperators_give_nr_of_couplings(rOp), QN_TYPE);
        QN_TYPE_read_dataset(group_id, "./qnumbers", rOp->qnumbers);

        ints_read_attribute(group_id, "nrops", &rOp->nrops);

        rOp->hss_of_ops = safe_malloc(rOp->nrops, int);
        ints_read_dataset(group_id, "./hss_of_ops", rOp->hss_of_ops);

        rOp->operators = safe_malloc(rOp->nrops, struct sparseblocks);
        for (int i = 0; i < rOp->nrops; ++i) {
                const int nr_blocks = rOperators_give_nr_blocks_for_operator(rOp, i);
                if (nr_blocks == 0) {
                        rOp->operators[i].beginblock = NULL;
                        rOp->operators[i].tel = NULL;
                }
                else
                        read_sparseblocks_from_disk(group_id, &rOp->operators[i], 
                                                    nr_blocks, i);
        }
        H5Gclose(group_id);
}

static void write_network_to_disk(const hid_t id)
{
        const hid_t group_id = H5Gcreate(id, "/network", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        ints_write_attribute(group_id, "nr_bonds", &netw.nr_bonds, 1);
        ints_write_dataset(group_id, "./bonds", (int*) netw.bonds, netw.nr_bonds * 2);

        ints_write_attribute(group_id, "psites", &netw.psites, 1);
        ints_write_attribute(group_id, "sites", &netw.sites, 1);
        ints_write_dataset(group_id, "./sitetoorb", netw.sitetoorb, netw.sites);
        ints_write_dataset(group_id, "./nr_left_psites", netw.nr_left_psites,
                           netw.nr_bonds);
        ints_write_dataset(group_id, "./order_psites", netw.order_psites,
                           netw.nr_bonds * netw.psites);

        ints_write_attribute(group_id, "sweeplength", &netw.sweeplength, 1);
        ints_write_dataset(group_id, "./sweep", netw.sweep, netw.sweeplength);

        H5Gclose(group_id);
}

static void read_network_from_disk(const hid_t id)
{
        const hid_t group_id = H5Gopen(id, "/network", H5P_DEFAULT);

        ints_read_attribute(group_id, "nr_bonds", &netw.nr_bonds);

        netw.bonds = malloc(netw.nr_bonds * sizeof netw.bonds[0]);
        ints_read_dataset(group_id, "./bonds", (int*) netw.bonds);

        ints_read_attribute(group_id, "psites", &netw.psites);
        ints_read_attribute(group_id, "sites", &netw.sites);

        netw.sitetoorb = safe_malloc(netw.sites, int);
        ints_read_dataset(group_id, "./sitetoorb", netw.sitetoorb);

        ints_read_attribute(group_id, "sweeplength", &netw.sweeplength);

        netw.sweep = safe_malloc(netw.sweeplength, int);
        ints_read_dataset(group_id, "./sweep", netw.sweep);

        H5Gclose(group_id);

        create_nr_left_psites();
        create_order_psites();
}
