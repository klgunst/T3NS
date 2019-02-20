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
#include <string.h>
#include <stdint.h>
#include "macros.h"

static void print_status(void)
{
#ifndef NDEBUG
        /* /proc/[pid]/status contains all human readible status, if i want to
         * get the memory later on and am only interested in that (e.g. to
         * write to a file every x seconds), it is better to read from
         * /proc/[pid]/statm (easier to parse) */
        FILE * status;
        char line[MY_STRING_LEN];
        if ((status = fopen("/proc/self/status", "r")) == NULL) {
                fprintf(stderr, "Error in openeing /proc/self/status\n");
                return;
        }
        while (fgets(line, sizeof line, status)) {
                fprintf(stderr, "%s", line);
        }
        fclose(status);
#endif
}

/* ========================================================================== */

void * safe_malloc_helper(long long s, size_t t, const char * typ, 
                          const char * file, int line, const char * func)
{
        void * pn = malloc(s * t);
        if (pn == NULL || s < 0) {
                const char *filname = strrchr(file, '/') + 1;
                fprintf(stderr, "%s:%d @%s :: Failed to allocate %s array of size %lld (%llu bytes)!\n"
                        "Maximal size of size_t : %lu\n", 
                        filname == NULL ? file : filname, line, func, 
                        typ, s, s*t, SIZE_MAX);
                print_status();
                exit(EXIT_FAILURE);
        }
        return pn;
}

void * safe_calloc_helper(long long s, size_t t, const char *typ, 
                          const char *file, int line, const char *func)
{
        void *pn = calloc(s, t);
        if (pn == NULL || s < 0) {
                const char *filname = strrchr(file, '/') + 1;
                fprintf(stderr, "%s:%d @%s :: Failed to reallocate %s array of size %lld (%llu bytes)!\n"
                        "Maximal size of size_t : %lu\n",
                        filname == NULL ? file : filname, line, func, 
                        typ, s, s*t, SIZE_MAX);
                print_status();
                exit(EXIT_FAILURE);
        }
        return pn;
}
