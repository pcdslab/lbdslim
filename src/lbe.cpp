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

/* Define this macro to enable DSLIM's Load Balancing Quality */
#define ANALYSIS


/* FUNCTION: LBE_Main (main)
 *
 * DESCRIPTION: Driver Application
 *
 * INPUT: none
 *
 * OUTPUT
 * @status: Status of execution
 */
STATUS LBE_Main()
{
    STATUS status = SLM_SUCCESS;
    UINT slm_chunks = 0;
    UINT *QA = NULL;
    ULONGLONG *Matches = NULL;
    SLM_vMods vModInfo;
    STRING modconditions = "3 STY 2";

    /* Benchmarking */
    auto start = chrono::system_clock::now();
    auto end   = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;

    /* Database and Dataset files */
    STRING filename = "C:/work/yeast_no_red.fasta";
    STRING querypath = "C:/work";
    STRING patt = ".ms2";
    DIR*    dir;
    dirent* pdir;
    vector<STRING> queryfiles;

    /* Initialize the vModInfo */
    vModInfo.num_vars = 1;
    vModInfo.vmods_per_pep = 3;
    vModInfo.vmods[0].aa_per_peptide = 2;
    vModInfo.vmods[0].modMass = 79.97 * SCALE;
    vModInfo.vmods[0].residues[0] = 'S';
    vModInfo.vmods[0].residues[1] = 'T';
    vModInfo.vmods[0].residues[2] = 'Y';

#ifdef ANALYSIS
    /* Analysis Variables */
    ULONGLONG sum = 0;
    DOUBLE   mean = 0;
    DOUBLE   std  = 0;
#endif /* ANALYSIS */

    /* Get thread count */
    UINT threads = UTILS_GetNumProcs();

#ifndef _OPENMP
    threads = 1;
#endif /* _OPENMP */

    /* Print Header */
    LBE_PrintHeader();

    /* Print start time */
    auto start_tim = chrono::system_clock::now();
    time_t start_time = chrono::system_clock::to_time_t(start_tim);
    cout << endl << "Start Time: " << ctime(&start_time) << endl;

    /* Check for a dangling / character */
    if (querypath.at(querypath.length()- 1) == '/')
    {
        querypath = querypath.substr(0, querypath.size() - 1);
    }

    /* Open all the query files */
    dir = opendir(querypath.c_str());

    /* Check if opened */
    if (dir != NULL)
    {
        while ((pdir = readdir(dir)) != NULL)
        {
            string cfile(pdir->d_name);

            /* Only add if there is a matching file */
            if (cfile.find(patt) != std::string::npos)
            {
                queryfiles.push_back(querypath + '/' + pdir->d_name);
            }
        }
    }
    else
    {
        status = ERR_FILE_NOT_FOUND;
    }

    /* Count the number of ">" entries in FASTA */
    if (status == SLM_SUCCESS)
    {
        status = LBE_CountPeps(threads, (CHAR *) filename.c_str(), modconditions);
    }

    /* Initialize internal structures */
    if (status == SLM_SUCCESS)
    {
        start = chrono::system_clock::now();

        /* Initialize the LBE */
        status = LBE_Initialize(threads, modconditions);

        end = chrono::system_clock::now();

        /* Compute Duration */
        elapsed_seconds = end - start;
        cout << "LBE Initialized with status:\t" << status << endl;
        cout << "Elapsed Time: " << elapsed_seconds.count() << "s" <<endl << endl;
    }

    /* Distribution Algorithm */
    if (status == SLM_SUCCESS)
    {
        start = chrono::system_clock::now();

        /* Distribute peptides among cores */
        status = LBE_Distribute(threads, _chunk, slm_chunks);

        end = chrono::system_clock::now();

        /* Compute Duration */
        elapsed_seconds = end - start;
        cout << "LBE Distributed with status:\t" << status << endl;
        cout << "Elapsed Time: " << elapsed_seconds.count() << "s" <<endl<< endl;

    }

    /* DSLIM-Transform */
    if (status == SLM_SUCCESS)
    {
        start = chrono::system_clock::now();

        /* Construct DSLIM by SLM Transformation */
        status = DSLIM_Construct(threads, &vModInfo);

        end = chrono::system_clock::now();

        /* Compute Duration */
        elapsed_seconds = end - start;
        cout << "SLM-Transform with status:\t" << status << endl;
        cout << "Elapsed Time: " << elapsed_seconds.count() << "s" <<endl<< endl;
    }

#ifdef ANALYSIS
    /* Analyze DSLIM fragment distribution */
    if (status == SLM_SUCCESS && slm_chunks > 1)
    {
        status = DSLIM_Analyze(threads, mean, std);

        /* Print DSLIM stats */
        cout << endl << "Number of DSLIM Chunks =\t" << slm_chunks << endl;
        cout << "Mean (Fragment Distribution) =\t" << mean << endl;
        cout << "StDev (Fragment Distribution) =\t" << std << endl << endl;

    }
#endif /* ANALYSIS */

    /* Initialize DSLIM Scorecard Manager */
    if (status == SLM_SUCCESS)
    {
        status = DSLIM_InitializeSC(threads);
        cout << "DSLIM SC Init with status:\t" << status << endl;
    }

    /* Initialize the Matches array */
    Matches = new ULONGLONG[slm_chunks];

    /* Allocate the Query Array */
    QA = new UINT[QCHUNK * QALEN];

    if (status == SLM_SUCCESS)
    {
        /* Initialize and process Query Spectra */
        for (UINT qf = 0; qf < queryfiles.size(); qf++)
        {
            cout << endl << "Query File: " << queryfiles[qf].c_str() << endl;

            /* Initialize Query MS/MS file */
            status = MSQuery_InitializeQueryFile((CHAR *) queryfiles[qf].c_str());

            /* DSLIM Query Algorithm */
            if (status == SLM_SUCCESS)
            {
                UINT qchunk_number = 0;
                DOUBLE cum_std = 0;

                /* Extract a chunk of MS/MS spectra and
                 * query against DSLIM Index */
                for (; QA != NULL;)
                {
                    /* Extract a chunk and return the chunksize */
                    UINT ms2specs = MSQuery_ExtractQueryChunk(QA, threads);

                    /* If the chunksize is zero, all done */
                    if (ms2specs <= 0)
                    {
                        break;
                    }

                    qchunk_number++;
                    cout << endl << "Extracted Batch of size:\t" << ms2specs << endl;

                    /* Reset Matches */
                    std::memset(Matches, 0x0, slm_chunks * sizeof(ULONGLONG));

                    start = chrono::system_clock::now();

                    /* Query the chunk */
                    status = DSLIM_QuerySpectrum(QA, ms2specs, Matches, threads);
                    end = chrono::system_clock::now();

                    /* Compute Duration */
                    elapsed_seconds = end - start;
                    cout << "Queried with status:\t\t" << status
                            << endl;
                    cout << "Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;

#ifdef ANALYSIS
                    /* Analyze the hit distribution */
                    if (status == SLM_SUCCESS && slm_chunks > 1)
                    {
                        /* Reset variables */
                        sum = 0;
                        std = 0;

                        for (UINT chs = 0; chs < slm_chunks; chs++)
                        {
                            sum += Matches[chs];
                        }

                        mean = ((DOUBLE) sum) / slm_chunks;

                        for (UINT chs = 0; chs < slm_chunks; chs++)
                        {
                            std += (((DOUBLE) Matches[chs] - mean) *
                                    ((DOUBLE) Matches[chs] - mean));
                        }

                        std /= slm_chunks;
                        std = sqrt(std);

                        /* Print Results */
                        cout << "Batch Number =\t\t\t" << qchunk_number << endl;
                        cout << "Average Hits per DSLIM chunk =\t" << mean << endl;
                        cout << "Standard Deviation of Hits =\t" << std << endl;

                        cum_std += std;
                    }
#endif /* ANALYSIS */
                }
#ifdef ANALYSIS
                /* Compute average std per query chunk */
                if (qchunk_number > 1)
                {
                    cum_std /= qchunk_number;
                    cout << endl << "Mean (StDev/batch) =\t\t" << cum_std << endl;
                }
#endif /* ANALYSIS */
            }
        }
    }

    /* Deinitialize DSLIM and LBE */
    if (status == SLM_SUCCESS)
    {
        status = LBE_Deinitialize();
        cout << endl <<"LBE Deinitialized with status:\t" << status << endl;
    }

    /* Deallocate QA */
    if (QA != NULL)
    {
        delete[] QA;
    }

    /* Deallocate Matches */
    if (Matches != NULL)
    {
        delete[] Matches;
    }

    /* Print final program status */
    cout << "\n\nLBE ended with status: \t\t" << status << endl;

    /* Print end time */
    auto end_tim = chrono::system_clock::now();
    time_t end_time = chrono::system_clock::to_time_t(end_tim);
    cout << endl << "End Time: " << ctime(&end_time) << endl;
    elapsed_seconds = end_tim - start_tim;
    cout << "Total Elapsed Time: " << elapsed_seconds.count() << "s" <<endl;

    /* Make sure stdout is empty at the end */
    fflush(stdout);

    return status;
}

