/*
 * This file is part of Load Balancing Algorithm for DSLIM
 *  Copyright (C) 2019  Muhammad Haseeb, Fahad Saeed
 *  Florida International University, Miami, FL
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef LBE_H_
#define LBE_H_

#include <chrono>
#include <ctime>
#include <dirent.h>
#include "common.h"
#include "utils.h"
#include "mods.h"
#include "slmerr.h"
#include "dslim.h"
#include "msquery.h"

#define LBE_Main               main

/*
 * FUNCTION: LBE_Initialize
 *
 * DESCRIPTION: Initialize internal peptides
 *              database from FASTA file
 *
 * INPUT:
 * @threads      : Number of parallel threads
 * @modconditions: String with mod conditions
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_Initialize(UINT threads, STRING modconditions);

/*
 * FUNCTION: LBE_Deinitialize
 *
 * DESCRIPTION: Deallocate all memory and
 *              reset variables
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_Deinitialize(VOID);

/*
 * FUNCTION: LBE_Distribute
 *
 * DESCRIPTION: Apply Distribution Policy on peptides
 *
 * INPUT:
 * @threads   : Number of parallel threads
 * @policy    : Policy ID (enum)
 * @slm_chunks: Number of distribution chunks
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_Distribute(UINT threads, DistPolicy policy, UINT& slm_chunks);

/*
 * FUNCTION: LBE_RevDist
 *
 * DESCRIPTION: Returns the actual peptide ID
 *
 * INPUT:
 * @virtID : (Distributed) Virtual Peptide ID.
 *
 * OUTPUT:
 * @realID: Actual SPI peptide ID
 */
UINT LBE_RevDist(UINT virtID);

/*
 * FUNCTION: LBE_CountPeps
 *
 * DESCRIPTION: Count peptides in FASTA and the
 *              number of mods that will be generated
 *
 * INPUT:
 * @threads      : Number of parallel threads
 * @filename     : Path to FASTA file
 * @modconditions: Mod generation conditions
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_CountPeps(UINT threads, CHAR *filename, STRING modconditions);

/*
 * FUNCTION: LBE_PrintHeader
 *
 * DESCRIPTION: Prints the LBE header
 *
 * INPUT : none
 * OUTPUT: none
 */
VOID   LBE_PrintHeader(VOID);

#endif /* LBE_H_ */
