#pragma once

#include "siteTensor.h"
#include "rOperators.h"
#include "symsecs.h"

struct Heffdata {
        int isdmrg;

        struct siteTensor siteObject;
        struct rOperators Operators[3];
        struct symsecs symarr[4][3]; /* nrsites * 3 */
        struct symsecs MPOsymsec;

        int rOperators_on_site[3];
        int posB; // for dmrg this is the second site
        int nr_qnB;
        QN_TYPE * qnB_arr;
        int * nr_qnBtoqnB;
        QN_TYPE ** qnBtoqnB_arr;
        int ** nrMPOcombos;
        int *** MPOs;

        int * instructions;
        int * instrbegin;
        double * prefactors;
};

void init_Heffdata(struct Heffdata * data, const struct rOperators * Operators,
                   const struct siteTensor * siteObject, int isdmrg);

void matvecT3NS(const double * vec, double * result, const void * vdata);

EL_TYPE * make_diagonal(const struct Heffdata * data);

void destroy_Heffdata(struct Heffdata * data);
