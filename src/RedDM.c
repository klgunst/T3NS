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
#include <stdlib.h>
#include <stdio.h>

#include "RedDM.h"
#include "siteTensor.h"
#include "network.h"
#include "macros.h"
#include "debug.h"

#define T3NS_REDDM_DEBUG

static EL_TYPE ** backuptel(struct siteTensor * T3NS)
{
        EL_TYPE ** backup = safe_malloc(netw.sites, *backup);
        for (int i = 0; i < netw.sites; ++i) {
                const int N = siteTensor_get_size(&T3NS[i]);
                backup[i] = safe_malloc(N, *backup[i]);
                for (int j = 0; j < N; ++j) {
                        backup[i][j] = T3NS[i].blocks.tel[j];
                }
        }
        return backup;
}

static void putback_backup(struct siteTensor * T3NS, EL_TYPE ** backup)
{
        for (int i = 0; i < netw.sites; ++i) {
                safe_free(T3NS[i].blocks.tel);
                T3NS[i].blocks.tel = backup[i];
        }
}

#ifdef T3NS_REDDM_DEBUG
static int check_orthonormality(struct siteTensor * orthocenter,
                                struct siteTensor * ortho)
{

        // Check normality of orthocenter
        double norm = 0;
        int i;
        for (i = 0; i < siteTensor_get_size(orthocenter); ++i) {
                norm += orthocenter->blocks.tel[i] * orthocenter->blocks.tel[i];
        }
        if (fabs(norm - 1) > 1e-12) {
                fprintf(stderr, "site %d : Orthocenter not normed.", 
                        orthocenter->sites[0]);
                return 1;
        }
        const int comb = get_common_bond(orthocenter->sites[0], ortho->sites[0]);
        int bonds[3];
        get_bonds_of_site(ortho->sites[0], bonds);
        for (i = 0; i < 3; ++i) { if (bonds[i] == comb) { break; } }
        assert(i != 3);

        if (!is_orthogonal(ortho, i)) {
                fprintf(stderr, "site %d : Not orthogonal.", ortho->sites[0]);
                return 1;
        }
        return 0;
}
#endif

static int initialize_rdm(struct RedDM * rdm, int mrdm)
{
        if (mrdm < 0 || mrdm > MAX_RDM) {
                fprintf(stderr, "Calculation of %d-body RDMs not supported. (Maximum: %d-body RDMs.)\n",
                        mrdm, MAX_RDM);
                return 1;
        }
        rdm->orb = netw.psites;
}

int get_RedDMs(struct siteTensor * T3NS, struct RedDM * rdm, int mrdm)
{
        printf(" >> Calculating RDMs\n");
        struct Rmatrix R;
        EL_TYPE ** telbackup = backuptel(T3NS);
        if (initialize_rdm(rdm, mrdm)) { return 1; }
        int exitcode = 0;

        int * sweep, swlength;
        if (make_simplesweep(1, &sweep, &swlength)) { return 1; }

        for (int i = 0; i < swlength; ++i) {
                // So, now going through all the sites in the sweep.
                // Assume that current site is orthogonality center, 
                // next site is orthogonal.
                struct siteTensor * orthocenter = &T3NS[sweep[i]];
                struct siteTensor * ortho = &T3NS[sweep[(i + 1) % swlength]];
                assert(orthocenter->nrsites == 1);
                assert(ortho->nrsites == 1);

#ifdef T3NS_REDDM_DEBUG
                if (check_orthonormality(orthocenter, ortho)) {
                        exitcode = 1;
                        break;
                }
#endif
                const int comb = get_common_bond(orthocenter->sites[0],
                                                 ortho->sites[0]);
                int bonds[3];
                get_bonds_of_site(orthocenter->sites[0], bonds);
                int oc_id;
                for (oc_id = 0; oc_id < 3; ++oc_id) {
                        if (bonds[oc_id] == comb) { break; }
                }
                assert(oc_id != 3);
                get_bonds_of_site(ortho->sites[0], bonds);
                int o_id;
                for (o_id = 0; o_id < 3; ++o_id) {
                        if (bonds[o_id] == comb) { break; }
                }
                assert(o_id != 3);


                struct siteTensor Q;
                if(qr(orthocenter, oc_id, &Q, &R)) {
                        exitcode = 1;
                        break;
                }
                destroy_siteTensor(orthocenter);
                *orthocenter = Q;

                struct siteTensor B;
                if(multiplyR(ortho, o_id, &R, 1, &B)) {
                        exitcode = 1;
                        break;
                }
                destroy_Rmatrix(&R);
                *ortho = B;
        }

        safe_free(sweep);
        for (int i = 0; i < netw.sites; ++i) { safe_free(telbackup[i]); }
        safe_free(telbackup);
        //putback_backup(T3NS, telbackup);
        return exitcode;
}
