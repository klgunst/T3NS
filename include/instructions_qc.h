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

#include "instructions.h"

void QC_fetch_pUpdate(struct instructionset * const instructions, 
                         const int bond, const int is_left);

void QC_fetch_bUpdate(struct instructionset * const instructions, 
                         const int bond, const int is_left);

void QC_fetch_merge(struct instructionset * const instructions, 
                       const int bond, int isdmrg);
