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

#include "dslim.h"

/* Sanity check for CHUNKSIZE */
#if (CHUNKSIZE <= 0)
#error "ABORT: The macro CHUNKSIZE must be > 0"
#endif /* (CHUNKSIZE == 0) */

/* Enternal Variables */
extern ULONGLONG    pepCount;
extern PepSeqs        seqPep;

/* Global Variables */
SLMindex          dslim; /* DSLIM Index       */
pepEntry    *pepEntries; /* SLM Peptide Index */
UINT           *SpecArr; /* Spectra Array     */

#ifdef VMODS
varEntry    *modEntries; /* SLM Mods Index    */
#endif /* VMODS */

/* Global Variables */
UINT     chunksize = 0;
UINT lastchunksize = 0;
UINT     nchunks   = 1;


/* FUNCTION: DSLIM_Construct
 *
 * DESCRIPTION: Construct DSLIM chunks
 *
 * INPUT:
 * @threads: Number of parallel threads to launch
 * @modInfo: DSLIM mods Information
 *
 * OUTPUT:
 * @status: status of execution
 */
STATUS DSLIM_Construct(UINT myid, SLM_vMods *modInfo)
{
    STATUS status = SLM_SUCCESS;

#ifdef VMODS
    /* Update gModInfo */
    status = UTILS_InitializeModInfo(modInfo);
#else
    LBE_UNUSED_PARAM(modInfo);
#endif /* VMODS */


    if (modInfo == NULL)
    {
        status = ERR_INVLD_PTR;
    }

    if (status == SLM_SUCCESS)
    {
        /* Spectra Array (SA) */
        SpecArr = new UINT[chunksize * iSERIES * F];

        /* Check if Spectra Array has been allocated */
        if (SpecArr == NULL)
        {
            status = ERR_BAD_MEM_ALLOC;
        }
    }
    if (status == SLM_SUCCESS)
    {
        /* Allocate memory for SLMChunks and SPI*/
        status = DSLIM_AllocateMemory(chunksize, 1);

        /* Construct DSLIM.iA */
        if (status == SLM_SUCCESS)
        {
            /* Distributed SLM Ions Array construction */
                /* Construct each DSLIM chunk in Parallel */
                status = DSLIM_ConstructChunk(1, myid);

                /* Apply SLM-Transform on the chunk */
                if (status == SLM_SUCCESS)
                {
                    status = DSLIM_SLMTransform(1, myid);
                }
        }
    }

    if (status == SLM_SUCCESS)
    {
        /* Construct DSLIM.bA */

            UINT *bAPtr = dslim.pepChunks.bA;

            /* Get size of the chunk */
            UINT csize = chunksize;

            UINT count = bAPtr[0];

            /* Initialize first and last bA entries */
            bAPtr[0] = 0;
            bAPtr[(MAX_MASS * SCALE)] = (csize * iSERIES * F);

            for (UINT li = 1; li <= (MAX_MASS * SCALE); li++)
            {
                UINT tmpcount = bAPtr[li];
                bAPtr[li] = bAPtr[li - 1] + count;
                count = tmpcount;
#ifdef VALIDATE_SLM
                if (bAPtr[li] < bAPtr[li - 1])
                {
                    std::cout << chunk_number << " " << li << " " << bAPtr[li - 1] << " " << count << " " << bAPtr[li] << std::endl;
                    status = ERR_INVLD_SIZE;

                    while(true);
                }
#endif /* VALIDATE_SLM */
            }

            /* Check if all correctly done */
            if (bAPtr[(MAX_MASS * SCALE)] != (csize * iSERIES * F))
            {
                status = ERR_INVLD_SIZE;
            }
    }

    /* Remove the temporary SpecArray (SA) */
    delete[] SpecArr;
    SpecArr = NULL;

    return status;
}

/*
 * FUNCTION: DSLIM_AllocateMemory
 *
 * DESCRIPTION: Allocate memory for DSLIM chunks
 *
 * INPUT:
 * @chsize: Chunk size
 * @Chunks: Number of chunks
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_AllocateMemory(UINT chsize, UINT Chunks)
{
    STATUS status = SLM_SUCCESS;

    dslim.nChunks = nchunks;

    /* Initialize direct hashing bA */
    dslim.pepChunks.bA = new UINT[(MAX_MASS * SCALE) + 1];

    if (dslim.pepChunks.bA != NULL)
    {
        /* Calculate the iA chunk size */
        INT size = chsize;

        /* Total Number of Ions = peps * #ion series * ions/ion series */
        dslim.pepChunks.iA = new UINT[(size * iSERIES * F)];

        if (dslim.pepChunks.iA == NULL)
        {
            status = ERR_INVLD_MEMORY;
        }
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }


    /* Allocate memory for SPI and Spectra Array */
    if (status == SLM_SUCCESS)
    {
        pepEntries = new pepEntry[pepCount];

        if (pepEntries == NULL)
        {
            status = ERR_INVLD_MEMORY;
        }
    }

    return status;
}

/*
 * FUNCTION: DSLIM_ConstructChunk
 *
 * INPUT:
 * @threads:      Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_ConstructChunk(UINT threads, UINT chunk_number)
{
    STATUS status = SLM_SUCCESS;

    /* Check if this chunk is the last chunk */
    BOOL lastChunk = (chunk_number == (nchunks - 1))? true: false;

    UINT SAPtr[(MAX_SEQ_LEN * iSERIES * MAXz)] = {};
    std::memset(dslim.pepChunks.bA, 0x0, (((SCALE * MAX_MASS) + 1) * sizeof(UINT)));


    if (status == SLM_SUCCESS)
    {
        UINT start_idx = chunk_number * chunksize;
        UINT interval = chunksize;

        /* Check for last chunk */
        if (lastChunk == true && nchunks > 1)
        {
            interval = lastchunksize;
        }

        for (UINT k = start_idx; k < (start_idx + interval); k++)
        {
            /* Filling point */
            UINT nfilled = (k - start_idx) * iSERIES * F;
            UINT *Spec = SAPtr;
            UINT *bAPtr = dslim.pepChunks.bA;

            /* Extract peptide Information */
            FLOAT pepMass = 0.0;
            CHAR *seq = NULL;
            INT len = 0;
            UINT pepID = LBE_RevDist(k);

#ifdef VMODS
            /* Check if pepID belongs to peps or mods */
            if (pepID >= pepCount)
            {
                /* Extract from Mods */
                varEntry *entry = modEntries + (pepID - pepCount);
                seq = &seqPep.seqs[seqPep.idx[entry->seqID]];
                len = (INT)seqPep.idx[entry->seqID + 1] - (INT)seqPep.idx[entry->seqID];

                /* Generate the Mod. Theoretical Spectrum */
                pepMass = UTILS_GenerateModSpectrum(seq, (UINT)len, Spec, entry->sites);

                /* Fill in the pepMass */
                entry->Mass = pepMass;
            }
            else
#endif /* VMODS */
            {
                /* Extract from Peps */
                pepEntry *entry = pepEntries + pepID;
                seq = &seqPep.seqs[seqPep.idx[pepID]];
                len = (INT)seqPep.idx[pepID + 1] - (INT)seqPep.idx[pepID];

                /* Generate the Theoretical Spectrum */
                pepMass = UTILS_GenerateSpectrum(seq, len, Spec);

                /* Fill in the pepMass */
                entry->Mass = pepMass;
            }

            /* If a legal peptide */
            if (pepMass >= MIN_MASS && pepMass <= MAX_MASS)
            {
                /* Calculate the length/2 of Spec */
                UINT half_len = ((iSERIES * MAX_SEQ_LEN * MAXz) / 2);

                /* Sort by ion Series and Mass */
                UTILS_Sort<UINT>(Spec, half_len, false);
                UTILS_Sort<UINT>((Spec + half_len), half_len, false);

                /* Choose SpecArr filling method */
                UINT nIons = (len - 1) * MAXz;
                UINT maxfill = F;
                UINT fill = 0;

                /* Fill the sorted ions */
                for (fill = 0; fill < nIons && fill < maxfill; fill++)
                {
                    /* Check if legal b-ion */
                    if (Spec[half_len - nIons + fill] >= (MAX_MASS * SCALE))
                    {
                        Spec[half_len - nIons + fill] = (MAX_MASS * SCALE) - 1;
                    }

                    SpecArr[nfilled + fill] = Spec[half_len - nIons + fill]; // Fill in the b-ion
                    bAPtr[SpecArr[nfilled + fill]]++; // Update the BA counter

                    /* Check for a legal y-ion */
                    if (Spec[(2 * half_len) - nIons + fill] >= (MAX_MASS * SCALE))
                    {
                        Spec[(2 * half_len) - nIons + fill] = (MAX_MASS * SCALE) - 1;
                    }

                    SpecArr[nfilled + F + fill] = Spec[(2 * half_len) - nIons + fill]; // Fill in the y-ion
                    bAPtr[SpecArr[nfilled + F + fill]]++; // Update the BA counter

                }

                /* Fill the rest with zeros */
                for (; fill < maxfill; fill++)
                {
                    SpecArr[nfilled + fill] = 0; // Fill in the b-ion
                    SpecArr[nfilled + F + fill] = 0; // Fill in the y-ion
                    bAPtr[0] += 2; // Update the BA counter
                }
            }

            /* Illegal peptide, fill in the container with zeros */
            else
            {
                /* Fill zeros for illegal peptides
                 * FIXME: Should not be filled into the chunk
                 *        and be removed from SPI as well
                 */
                std::memset(&SpecArr[nfilled], 0x0, sizeof(UINT) * iSERIES * F);
                bAPtr[0] += (iSERIES * F); // Update the BA counter
            }
        }

    }

    return status;
}

/*
 * FUNCTION: DSLIM_SLMTransform
 *
 * DESCRIPTION: Constructs SLIM Transform
 *
 * INPUT:
 * @threads     : Number of parallel threads
 * @chunk_number: Chunk Index
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_SLMTransform(UINT threads, UINT chunk_number)
{
    STATUS status = SLM_SUCCESS;

    /* Check if this chunk is the last chunk */
    UINT size = chunksize;

    UINT *iAPtr = dslim.pepChunks.iA;
    UINT iAsize = size * iSERIES * F;

    /* Construct DSLIM.iA */
    for (UINT k = 0; k < iAsize; k++)
    {
        iAPtr[k] = k;
    }

    KeyVal_Serial<INT>((INT *)SpecArr, (INT *)iAPtr, iAsize);

    /* Check integrity of SLM-Transform */
#ifdef VALIDATE_SLM
    BOOL integrity = true;

    for (UINT k = 1; k < iAsize; k++)
    {
        if (integ[iAPtr[k]] < integ[iAPtr[k - 1]])
        {
            integrity = false;
        }
    }

    if (integrity == false)
    {
        while (!integrity);
        status = ERR_INVLD_SORT;
    }
#endif /* VALIDATE_SLM */

    return status;
}

/*
 * FUNCTION: DSLIM_InitializeSC
 *
 * DESCRIPTION: Initialize Scorecard for DSLIM
 *
 * INPUT:
 * @threads: Number of parallel threads
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_InitializeSC(UINT threads)
{
    STATUS status = SLM_SUCCESS;

    threads = 1;

    /* Check if this chunk is the last chunk */
    UINT size = chunksize;

    /* Allocate memory for scorecard */
    UCHAR *SC = new UCHAR[size];
    dslim.pepChunks.sC = SC;

    if (SC != NULL)
    {
        std::memset(SC, 0x0, size);
    }
    else
    {
        status = ERR_INVLD_MEMORY;
    }

    return status;
}

/*
 * FUNCTION: DSLIM_Analyze
 *
 * DESCRIPTION: Analyze the DSLIM distribution
 *
 * INPUT:
 * @threads: Number of parallel threads
 * @avg    : Pointer to DSLIM mean load
 * @std    : Pointer to DSLIM load distribution
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_Analyze(UINT threads, DOUBLE &avg, DOUBLE &std)
{
    STATUS status = SLM_SUCCESS;
    DOUBLE sum_std = 0;
    ULONGLONG truecount = 1;

    if (nchunks <= 1)
    {
        status = ERR_INVLD_SIZE;
    }

    if (status == SLM_SUCCESS)
    {
        DOUBLE *arr = new DOUBLE[(MAX_MASS * SCALE)+1];

        /* Trivial to include zeros or last value */
        arr[(MAX_MASS * SCALE)] = 0x0;
        arr[0] = 0x0;

        /* Compute the means */

        for (UINT i = 1; i < (MAX_MASS * SCALE) -1; i++)
        {
            arr[i] = 0x0;

            for (UINT j = 0; j < nchunks; j++)
            {
                arr[i] += (dslim.pepChunks.bA[i+1] - dslim.pepChunks.bA[i]);
            }

            arr[i] /= nchunks;

            if (arr[i] > 0)
            {
                truecount++;
            }
        }


        /* Compute the Standard Deviations */
        for (UINT i = 1; i < (MAX_MASS * SCALE) - 1; i++)
        {
            DOUBLE mean = arr[i];
            arr[i] = 0x0;

            if (mean > 0)
            {
                for (UINT j = 0; j < nchunks; j++)
                {
                    arr[i] += (((DOUBLE) dslim.pepChunks.bA[i+1]
                            - (DOUBLE) dslim.pepChunks.bA[i] - mean)
                            * ((DOUBLE) dslim.pepChunks.bA[i+1]
                            - (DOUBLE) dslim.pepChunks.bA[i] - mean));
                }

                arr[i] = sqrt(arr[i] / nchunks);
            }
        }

        /* Compute mean of stdevs */
        for (UINT i = 0; i < MAX_MASS * SCALE; i++)
        {
            sum_std += arr[i];
        }

        sum_std /= truecount;
        avg = sum_std;

        sum_std = 0;

        /* Compute stdev of stdevs */
        for (UINT i = 0; i < MAX_MASS * SCALE; i++)
        {
            sum_std += ((arr[i] - avg) * (arr[i] - avg));
        }

        /* Gather the counts to truecount */
        sum_std /= truecount;
        sum_std = std::sqrt(sum_std);

        std = sum_std;

        delete[] arr;
        arr = NULL;
    }

    return status;
}

/*
 * FUNCTION: DSLIM_Deinitialize
 *
 * DESCRIPTION: Deallocate DSLIM memory
 *
 * INPUT: none
 *
 * OUTPUT:
 * @status: Status of execution
 */
STATUS DSLIM_Deinitialize()
{
    /* Deallocate all the DSLIM chunks */

    SLMchunk curr_chunk = dslim.pepChunks;

    if (curr_chunk.bA != NULL)
    {
        delete[] curr_chunk.bA;
        curr_chunk.bA = NULL;
    }

    if (curr_chunk.iA != NULL)
    {
        delete[] curr_chunk.iA;
        curr_chunk.iA = NULL;
    }

    if (curr_chunk.sC != NULL)
    {
        delete[] curr_chunk.sC;
        curr_chunk.sC = NULL;
    }


    /* Reset Global Variables */
    nchunks = 0;
    dslim.nModChunks = 0;
    dslim.nPepChunks = 0;
    dslim.nChunks = 0;
    chunksize = 0;
    lastchunksize = 0;

    return SLM_SUCCESS;
}
