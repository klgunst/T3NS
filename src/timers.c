/*
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
*/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include <string.h>

#include "timers.h"

struct timers init_timers(const char **names, const int * keys, int n)
{
        struct timers tim = { .n = n, .timers = safe_malloc(n, *tim.timers)};
        gettimeofday(&tim.inittime, NULL);

        for (int i = 0; i < n; ++i) {
                strncpy(tim.timers[i].name, names[i], MY_STRING_LEN);
                tim.timers[i].name[MY_STRING_LEN - 1] = '\0';
                tim.timers[i].key = keys[i];

                tim.timers[i].t = 0;
                tim.timers[i].ticed = false;
                tim.timers[i].touched = false;
        }
        return tim;
}

void destroy_timers(struct timers * tim)
{
        safe_free(tim->timers);
}

static int search_key(struct timers * tim, int key)
{
        for (int i = 0; i < tim->n; ++i) {
                if (tim->timers[i].key == key) { return i; }
        }
        return -1;
}

int tic(struct timers * tim, int key)
{
        const int id = search_key(tim, key);
        if (id == -1) {
                fprintf(stderr, "Key %d not found in timer.", key);
                return 1;
        }

        if (tim->timers[id].ticed) {
                fprintf(stderr, "Timer with key %d already ticed.", key);
                return 1;
        }

        gettimeofday(&tim->timers[id].tictime, NULL);
        tim->timers[id].ticed = true;
        tim->timers[id].touched = true;

        return 0;
}

int toc(struct timers * tim, int key)
{
        const int id = search_key(tim, key);
        if (id == -1) {
                fprintf(stderr, "Key %d not found in timer.", key);
                return 1;
        }

        if (!tim->timers[id].ticed) {
                fprintf(stderr, "Timer with key %d was not ticed.", key);
                return 1;
        }

        struct timeval tv;
        gettimeofday(&tv, NULL);
        tim->timers[id].ticed = false;
        long long t_el = (tv.tv_sec - tim->timers[id].tictime.tv_sec) * 
                1000000LL + tv.tv_usec - tim->timers[id].tictime.tv_usec;
        tim->timers[id].t += t_el * 1e-6;
        return 0;
}

void print_timers(struct timers * tim, const char * prefix, bool onlytouched)
{
        for (int i = 0; i < tim->n; ++i) {
                struct timer TimTim = tim->timers[i];
                if (TimTim.ticed) {
                        fprintf(stderr, "Timer %s was ticed but not toced.\n",
                                TimTim.name);
                }
                if (!onlytouched || TimTim.touched)
                        printf("%s%-35s :: %lf sec\n",
                               prefix, TimTim.name, TimTim.t);
        }
        struct timeval tv;
        gettimeofday(&tv, NULL);
        long long t_el = (tv.tv_sec - tim->inittime.tv_sec) * 
                1000000LL + tv.tv_usec - tim->inittime.tv_usec;
        printf("%s%-35s :: %lf sec\n", prefix, "Total time", t_el * 1e-6);
}

int add_timers(struct timers * result, const struct timers * toadd)
{
        if (result->n != toadd->n) {
                fprintf(stderr, "The inputted timers do not have the same amount of timers for adding.\n");
                return 1;
        }
        for (int i = 0; i < toadd->n; ++i) {
                if (!toadd->timers[i].touched) { continue; }
                const int id = search_key(result, toadd->timers[i].key);
                if (id == -1 || strncmp(result->timers[id].name, 
                                        toadd->timers[i].name, MY_STRING_LEN) != 0) {
                        fprintf(stderr, "The inputted timers are not the same for adding them.\n");
                        return 1;
                }
                result->timers[id].t += toadd->timers[i].t;
                result->timers[id].touched = true;
        }
        return 0;
}

void reset_timers(struct timers * tim)
{
        for (int i = 0; i < tim->n; ++i) {
                tim->timers[i].ticed = false;
                tim->timers[i].t = 0;
        }
}
