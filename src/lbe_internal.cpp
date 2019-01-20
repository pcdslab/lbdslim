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

#include "lbe.h"
using namespace std;

/* Global Variables */
UINT       *dist;
ULONGLONG    pepCount;
ULONGLONG    modCount;
ULONGLONG  totalCount;
PepSeqs   seqPep;

vector<STRING> Seqs;

/* External Variables */
extern UINT       chunksize;
extern UINT   lastchunksize;
extern UINT         nchunks;
#ifdef VMODS
extern varEntry *modEntries;
#endif /* VMODS */

/* Static function Prototypes */
static STATUS LBE_AllocateMem(UINT N, UINT M);

/*
 * FUNCTION: LBE_AllocateMem
 *
 * DESCRIPTION: Allocate memory for Data Structures
 *
 * INPUT:
 * @N: Number of peptides
 * @M: Number of mods
 *
 * OUTPUT:
 * @status: Status of execution
 */
static STATUS LBE_AllocateMem(UINT N, UINT M)
{
    STATUS status = SLM_SUCCESS;
#ifdef VMODS
    modEntries = NULL;
#endif /* VMODS */

    /* Allocate Memory for seqPep */
    seqPep.seqs = new AA[seqPep.AAs];
    seqPep.idx = new UINT[N+1];

    /* Set the seqPep.idx[N] to AAs for consistency*/
    seqPep.idx[N] = seqPep.AAs;

    if (seqPep.seqs == NULL || seqPep.idx == NULL)
    {
        status = ERR_BAD_MEM_ALLOC;
    }

#ifdef VMODS
    /* Allocate the seqMod */
    if (status == SLM_SUCCESS)
    {
        modEntries = new varEntry[M];

        if (modEntries == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }
#endif /* VMODS */

    /* Allocate the distribution array */
    if (status == SLM_SUCCESS)
    {
        dist = new UINT[totalCount];

        if (dist == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }

    return status;
}

/*
 * FUNCTION: LBE_Initialize
 *
 * DESCRIPTION: Initialize internal peptides
 *              database from FASTA file
 *
 * INPUT:
 * @threads      : Number of parallel threads
 * @filename     : Path to FASTA file
 * @modconditions: String with mod conditions
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS LBE_Initialize(UINT threads, STRING modconditions)
{
    STATUS status = SLM_SUCCESS;
    UINT iCount = 1;
    STRING seq;

    /* Check if ">" entries are > 0 */
    if (pepCount > 0)
    {
        status = LBE_AllocateMem(pepCount, modCount);
    }
    else
    {
        status = ERR_INVLD_PARAM;
    }

    /* If Seqs was successfully filled */
    if (Seqs.size() != 0 && status == SLM_SUCCESS)
    {
        UINT idx = 0;

        /* Extract First Sequence manually */
        seq = Seqs.at(0);

        /* Increment counter */
        iCount += 2;

        memcpy((void *)&seqPep.seqs[idx], (const void *)seq.c_str(), seq.length());
        seqPep.idx[0] = idx;
        idx+=seq.length();

#ifdef DEBUG
        cout << seq << endl;
#endif /* DEBUG */

        for (UINT i = 1; i < Seqs.size(); i++)
        {
            /* Extract Sequences */
            STRING seq = Seqs.at(i);
            /* Copy into the seqPep.seqs array */
            memcpy((void *) &seqPep.seqs[idx], (const void *) seq.c_str(), seq.length());

            /* Increment the counters */
            iCount += 2;
            seqPep.idx[iCount / 2 - 1] = idx;
            idx += seq.length();

#ifdef DEBUG
            cout << seq << endl;
#endif /* DEBUG */
        }

        /* Get the peptide count */
        iCount /= 2;
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
        cout << endl << "ABORT: File could not be opened" << endl;
    }

    /* Make sure that the '>' == PEPTIDES */
    if (iCount != pepCount)
    {
        cout << endl << "pepCount != iCount - Please check the FASTA file";
        status = ERR_INVLD_SIZE;
    }

#ifdef VMODS

    /* Initialize gModInfo */
    //status = UTILS_InitializeModInfo(modconditions);

    if (modCount > 0)
    {
        /* Fill in the mods entry */
        status = MODS_GenerateMods(threads, modCount, modconditions);
    }
#endif /* VMODS */

    /* Make sure Seqs is clear anyway */
    Seqs.clear();

    if (status != SLM_SUCCESS)
    {
        (VOID) LBE_Deinitialize();
    }

    return status;
}

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
STATUS LBE_Deinitialize(VOID)
{
   (VOID) DSLIM_Deinitialize();

    /* Reset all counters */
    seqPep.AAs = 0;
    pepCount = 0;
    modCount = 0;
    totalCount = 0;

    /* Deallocate memory */
    if (seqPep.idx != NULL)
    {
        delete[] seqPep.idx;
        seqPep.idx = NULL;
    }

    if (seqPep.seqs != NULL)
    {
        delete[] seqPep.seqs;
        seqPep.seqs = NULL;
    }

    if (dist != NULL)
    {
        delete[] dist;
        dist = NULL;
    }

#ifdef VMODS
    if (modEntries != NULL)
    {
        delete[] modEntries;
        modEntries = NULL;
    }
#endif /* VMODS */

    return SLM_SUCCESS;
}

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
STATUS LBE_Distribute(UINT threads, DistPolicy policy, UINT& slm_chunks)
{
    STATUS status = 0;
    UINT N = totalCount;
    UINT p = threads;

    if ((N / p) > CHUNKSIZE)
    {
        /* Set the chunksize to max chunksize */
        chunksize = CHUNKSIZE;

        /* Calculate the number of chunks */
        nchunks = ((N % chunksize) == 0)?
                   (N / chunksize)      :
                   ((N / chunksize) + 1);
    }
    else
    {
            /* Calculate the chunksize */
        chunksize = ((N % p) == 0)?
                         (N/p)    :
                         ((N+p)/p);

        /* Set the number of chunks to p */
        nchunks = p;
    }

    /* Calculate the size of last chunk */
    UINT factor = N / chunksize;

    lastchunksize = ((N % chunksize) == 0)?
                     chunksize            :
                     N - (chunksize * factor);

    if (policy != _cyclic)
    {
        /* Initialize the array with OpenMP */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
        for (UINT k = 0; k < N; k++)
        {
            dist[k] = k;
        }
    }

    /* Apply the distribution policy */
    switch (policy)
    {
        case _random:
        {
            /* Generate a random seed */
            ULONGLONG seed = time(0);
#ifdef _OPENMP
            seed -= (seed % 997) * 72393;
            seed += (seed % 93617) * 17329;
            seed -= (seed % 2357) * 24178;
            seed += (seed % 19487) * 3175753;
#else
            seed = ((seed % 371291) + 213) * 317753;
#endif /* _OPENMP */

            std::cout << "Random Sampling Seed:\t\t" << seed << std::endl;

            for (UINT chno = 0; chno < nchunks; chno++)
            {
                UINT chsz = ((chno == nchunks - 1) && (nchunks > 1))?
                             lastchunksize                          :
                             chunksize;

                /* Shuffle the local indices */
                status = UTILS_ShuffleI(dist + (chno * chunksize), chsz, (seed << (chno + 1)));
            }

            /* Shuffle the global indices */
            status = UTILS_ShuffleI(dist, N, (seed + time(0)) << (omp_get_num_threads()));

            break;
        }

        case _chunk:
        {
            status = SLM_SUCCESS;
            break;
        }

        case _cyclic:
        {
            UINT i = 0;

            /* Initialize the array with OpenMP */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
            for (i = 0; i < lastchunksize * nchunks; i++)
            {
                /* Compute positions */
                UINT l = i % nchunks;
                UINT m = i/nchunks;
                /* Populate the cyclic sample */
                dist[(l * chunksize) + m] = i;
            }

            /* Fill the remaining peptides in (nchunks-1) chunks */
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif /* _OPENMP */
            for (i = lastchunksize * nchunks; i < totalCount; i++)
            {
                /* Compute positions */
                UINT l = i % (nchunks - 1);
                UINT m = (i- (lastchunksize * nchunks))/(nchunks -1);

                /* Populate the cyclic sample */
                dist[(l * chunksize) + lastchunksize + m] = i;
            }

            status = SLM_SUCCESS;
            break;
        }

        default:
        {
            status = ERR_INVLD_PARAM;
            break;
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Return the number of chunks created */
        slm_chunks = nchunks;
    }

    return status;
}

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
STATUS LBE_CountPeps(UINT threads, CHAR *filename, STRING modconditions)
{
    STATUS status = SLM_SUCCESS;
    pepCount = 0;
    modCount = 0;
    seqPep.AAs = 0;
    STRING line;
    FLOAT pepmass = 0.0;

#ifndef VMODS
    LBE_UNUSED_PARAM(modconditions);
#endif /* VMODS */

    /* Open file */
    ifstream file(filename);

    if (file.is_open())
    {
        while (getline(file, line))
        {
            if (line.at(0) != '>')
            {
                /* Linux has a weird \r at end of each line */
                if (line.at(line.length() - 1) == '\r')
                {
                    line = line.substr(0, line.size() - 1);
                }

                /* Transform to all upper case letters */
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);

                /* Calculate mass of peptide */
                pepmass = UTILS_CalculatePepMass((AA *)line.c_str(), line.length());

                /* Check if the peptide mass is legal */
                if (pepmass > MIN_MASS && pepmass < MAX_MASS)
                {
                    pepCount++;
                    seqPep.AAs += line.length();
                    Seqs.push_back(line);
                }
            }
        }
    }
    else
    {
        cout << endl << "FATAL: Could not read FASTA file" << endl;
        status = ERR_INVLD_PARAM;
    }

#ifdef VMODS
    /* Count the number of variable mods given
     * modification information */
    if (status == SLM_SUCCESS)
    {
        modCount = MODS_ModCounter(threads, modconditions);
    }

#endif /* VMODS */

    /* Print if everything is okay */
    if (status == SLM_SUCCESS)
    {
        /* Return the total count */
        totalCount = pepCount + modCount;

        cout << "Number of Peptides   = \t\t" << pepCount << endl;
        cout << "Number of Mods       = \t\t" << modCount << endl;
        cout << "Distributed SPI Size = \t\t" << totalCount << endl << endl;

        /* Close the file once done */
        file.close();
    }

    return status;
}

/*
 * FUNCTION: LBE_PrintHeader
 *
 * DESCRIPTION: Prints the LBE header
 *
 * INPUT : none
 * OUTPUT: none
 */
VOID LBE_PrintHeader(VOID)
{
    cout << "\n"
            "*********************************"
            "\n    Load Balancing for DSLIM\n"
            "Florida International University\n"
            "        Miami, FL, USA\n"
            "*********************************"
          << endl << endl;

    return;
}
