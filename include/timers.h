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
#pragma once
#include <stdbool.h>
#include <time.h>
#include "macros.h"

/** 
 * @file timers.h
 *
 * The header file for timers.
 */

/// This structure defines a single timer
struct timer {
        /// The name of the timer
        char name[MY_STRING_LEN];
        /// The key of the timer
        int key;
        /// True if you performed a tic on this timer
        bool ticed;
        /// True if you at least tic-toced it once (or added)
        bool touched;
        /// The time when ticed
        struct timeval tictime;
        /// The total seconds already tictoc-ed
        double t;
};

/// A collection of timers
struct timers {
        /// The number of timers
        int n;
        /// The different timers
        struct timer * timers;
        /// The time when initialized
        struct timeval inittime;
};

/**
 * @brief Initializes timers
 *
 * @param [in] names The names of the timers, these are deep copied.
 * @param [in] keys The keys of the timers.
 * @param [in] n The number of timers, both keys and names should have n
 * elements.
 * @return The timers
 */
struct timers init_timers(const char **names, const int * keys, int n);

/// Destroys a timers structure
void destroy_timers(struct timers * tim);

/// Starts the timing of a timer with the given key.
int tic(struct timers * tim, int key);

/// Stops the timing of a timer with the given key and adds the timed time.
int toc(struct timers * tim, int key);

/** Prints the timers, with an option to add a prefix and only print the 
 * touched timers
 */
void print_timers(struct timers * tim, const char * prefix, bool onlytouched);

/// Adds two timers together if they have the same timer.name and timer.key.
int add_timers(struct timers * result, const struct timers * toadd);

/// Resets the timers
void reset_timers(struct timers * tim);
