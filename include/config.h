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

#ifndef INCLUDE_CONFIG_H_
#define INCLUDE_CONFIG_H_

/* Contains MACROS that configure SLM Index build */

/* Manually switch _OPENMP */
//#undef _OPENMP

/* Define if OS == WINDOWS*/
#define WINDOWS

/* Enable to show SLMIndex Benchmark information  */
#undef BENCHMARK

/* Number of b- and y-ions per indexed fragment   */
#define F                        40

/* Fragment mass tolerance in (Daltons * SCALE)   */
#define dF                        5 // = 0.05Da

/* Min Shared Peaks to qualify as candidate PSM   */
#define MIN_SHRD_PKS              4

/* Max fragment charge to expect                  */
#define MAXz                      3

/* Scaling factor for integer conversion          */
#define SCALE                    100

/* The resolution for ions/fragments.             */
#define RESOLUTION               (1/SCALE)

/* Min and Max length of peptide sequence         */
#define MIN_SEQ_LEN               6
#define MAX_SEQ_LEN              60

/* Min and Max precursor mass for peptide         */
#define MIN_MASS                 100
#define MAX_MASS                 5000

/* Chunk size in number of peptides (def: 1.8M)   */
#define CHUNKSIZE                1800000

/* Number of fragments (peaks) per query spectrum */
#define QALEN                    100

/* How many query spectra per chunk */
#define QCHUNK                   1000

/************** Developer Mode options: Danger Area *****************/

/* Enable DSLIM Mods construction                 */
#define VMODS

/* Enable debug checkpoints                       */
#undef DEBUG

/* Recycle and reuse the already allocated memory
 * instead of new allocations when possible       */
#undef RECYCLE

/* Future: Define if the peak intensities
 * are required for the query spectra             */
#undef REQUIRE_INTENSITY

/* TODO:  Only b- and -y-ions used in current
 *               implementation                   */

/* Define the Ion Series and add it to iSERIES    */
#define aIONS                     0
#define bIONS                     1
#define cIONS                     0
#define xIONS                     0
#define yIONS                     1
#define zIONS                     0
#define NLOSS                     0
#define IMMONIUM                  0

/* Number of ion series to generate per spectrum  */
#define iSERIES                   (aIONS + \
                                   bIONS + \
                                   cIONS + \
                                   xIONS + \
                                   yIONS + \
                                   zIONS + \
                                   NLOSS + \
                                   IMMONIUM)

/*******************************************************************/

#endif /* INCLUDE_CONFIG_H_ */
