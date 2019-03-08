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
#include <mpi.h>

#define MASTER               0

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
STATUS LBE_Main(int argc, char** argv)
{
    STATUS status = ERR_MPI_ERR;
    INT threads = 0;
    UINT *QA = NULL;
    ULONGLONG Matches = 0;
    INT myid = -1;
    SLM_vMods vModInfo;
    DOUBLE *tbuf = NULL;
    MPI_Comm comm = MPI_COMM_WORLD;
    DistPolicy policy = _random;

    /* Benchmarking */
    DOUBLE start = 0;
    DOUBLE end = 0;
    DOUBLE qtime = 0;
    DOUBLE mean = 0;
    DOUBLE sum = 0;
    DOUBLE stdev = 0;

    auto start_tim = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = start_tim - start_tim;

    /* Handle unused parameters to avoid compiler warnings */
    LBE_UNUSED_PARAM(argc);
    LBE_UNUSED_PARAM(argv);

    if (MPI_Init(&argc, &argv) == MPI_SUCCESS)
    {
        status = SLM_SUCCESS;
    }

    /* Database and Dataset files */
    STRING filename = "/home/mhaseeb/human_proc.fasta";
    STRING querypath = "/home/mhaseeb";
    STRING patt = ".ms2";
    DIR*    dir;
    dirent* pdir;
    vector<STRING> queryfiles;

    /* Initialize the vModInfo */
    STRING modconditions = "3 M 2 CK 1 NQ 2";

    vModInfo.num_vars = 1;
    vModInfo.vmods_per_pep = 3;
    vModInfo.vmods[0].aa_per_peptide = 2;
    vModInfo.vmods[0].modMass = 15.99 * SCALE;
    vModInfo.vmods[0].residues[0] = 'M';

    vModInfo.vmods[1].aa_per_peptide = 2;
    vModInfo.vmods[1].modMass = 114.0429 * SCALE;
    vModInfo.vmods[1].residues[0] = 'C';
    vModInfo.vmods[1].residues[1] = 'K';
	
    vModInfo.vmods[2].aa_per_peptide = 2;
    vModInfo.vmods[2].modMass = 0.98 * SCALE;
    vModInfo.vmods[2].residues[0] = 'N';
    vModInfo.vmods[2].residues[1] = 'Q';


    // Get the number of processes
    if (status == SLM_SUCCESS)
    {
        status = MPI_Comm_size(comm, &threads);

        if (status == SLM_SUCCESS)
        {
            // Get the rank of the process
            status = MPI_Comm_rank(comm, &myid);
        }
    }

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

    if (status == SLM_SUCCESS && myid == MASTER)
    {
        /* Print Header */
        LBE_PrintHeader();

        /* Print start time */
        time_t start_time = chrono::system_clock::to_time_t(start_tim);
        cout << endl << "Start Time: " << ctime(&start_time) << endl;

        cout << endl << "MPI Comm Size        =\t\t" << threads << endl;

        tbuf = new DOUBLE[threads];

        std::memset(tbuf, 0x0, sizeof(DOUBLE) * threads);

    }

    /* Count the number of ">" entries in FASTA */
    if (status == SLM_SUCCESS)
    {
        status = LBE_CountPeps(myid, (CHAR *) filename.c_str(), modconditions);
    }

    /* Initialize internal structures */
    if (status == SLM_SUCCESS)
    {
        /* Initialize the LBE */
        status = LBE_Initialize(1, modconditions);
    }

    /* Distribution Algorithm */
    if (status == SLM_SUCCESS)
    {
        DOUBLE seed = 0;

        if (myid == MASTER)
        {
            cout << "Distribution Policy  =\t\t"<< (INT) policy << endl;
        }

        if (policy == _random)
        {
            /* Generate a seed */
            if (myid == MASTER)
            {
                seed = MPI_Wtime();
                cout << "Seed Used            =\t\t" << (ULONGLONG)seed << endl;
            }

            /* Broadcast the seed */
            status = MPI_Bcast(&seed, 1, MPI_DOUBLE, MASTER, comm);

            /* Distribute peptides among cores */
            status = LBE_Distribute(myid, policy, threads, (ULONGLONG)seed);
        }
        else
        {
            /* Distribute peptides among cores */
            status = LBE_Distribute(myid, policy, threads, (ULONGLONG)seed);
        }
    }

    /* DSLIM-Transform */
    if (status == SLM_SUCCESS)
    {
        start = MPI_Wtime();

        /* Construct DSLIM by SLM Transformation */
        status = DSLIM_Construct(myid, &vModInfo);

        end = MPI_Wtime();

        /* Compute Duration */
        qtime = end - start;
    }

    /* No need of LBE strucutures beyond this point */
    if (status == SLM_SUCCESS && myid != MASTER)
    {
        status = LBE_Deinitialize();
    }

    /* Initialize DSLIM Scorecard Manager */
    if (status == SLM_SUCCESS)
    {
        status = DSLIM_InitializeSC(1);
    }

    if (myid==MASTER && status == SLM_SUCCESS)
    {
        cout << "SLM-Transform with status:\t" << status << endl;
        cout << "Elapsed Time: " << qtime << "s" <<endl<< endl;
    }
    /* Allocate the Query Array */
    QA = new UINT[QCHUNK * QALEN];

    if (status == SLM_SUCCESS)
    {
        /* Reset the qtime */
        qtime = 0;

        /* Initialize and process Query Spectra */
        for (UINT qf = 0; qf < queryfiles.size(); qf++)
        {
            /* Initialize Query MS/MS file */
            status = MSQuery_InitializeQueryFile((CHAR *) queryfiles[qf].c_str());

            /* DSLIM Query Algorithm */
            if (status == SLM_SUCCESS)
            {
                UINT qchunk_number = 0;

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

                    start = MPI_Wtime();

                    /* Query the chunk */
                    status = DSLIM_QuerySpectrum(QA, ms2specs, Matches, threads);
                    end = MPI_Wtime();

                    /* Compute Duration */
                    qtime += end - start;

                }

            }
        }
    }

    /* Deinitialize DSLIM and LBE */
    if (status == SLM_SUCCESS)
    {
        status = DSLIM_Deinitialize();

        if (myid == MASTER)
        {
            status = LBE_Deinitialize();
        }
    }

    /* Deallocate QA */
    if (QA != NULL)
    {
        delete[] QA;
    }

    if (status == SLM_SUCCESS)
    {
        status = MPI_Gather(&qtime, 1,
        MPI_DOUBLE,
        tbuf, 1,
        MPI_DOUBLE,
        MASTER, comm);
    }

    (VOID)MPI_Finalize();

    if (status == SLM_SUCCESS && myid == MASTER)
    {
        cout << endl;

        /* Analyze the load balance */
        if (status == SLM_SUCCESS && threads > 1)
        {
            DOUBLE maxqtime = 0;

            for (INT chs = 0; chs < threads; chs++)
            {
                (tbuf[chs] > maxqtime)? maxqtime=tbuf[chs]:maxqtime;
                sum += tbuf[chs];
            }

            mean = (sum) / threads;

            for (INT chs = 0; chs < threads; chs++)
            {
                stdev += fabs(tbuf[chs] - mean);
            }

            stdev /= threads;

            /* Print Results */
            cout << "Average Query Time : " << mean << "s" << endl;
            cout << "Maximum Query Time : " << maxqtime << "s" << endl;
            cout << "\nAverage Deviation  : " << stdev << "s" << endl;
            cout << "Max Load Imbalance : " << maxqtime-mean << "s" << endl;
        }

        if (threads ==1 && status == SLM_SUCCESS)
        {
            cout << "Average Query Time : " << tbuf[0] << "s" << endl;
		}

        delete[] tbuf;
        tbuf = NULL;

        /* Print end time */
        auto end_tim = chrono::system_clock::now();
        time_t end_time = chrono::system_clock::to_time_t(end_tim);
        cout << endl << "End Time: " << ctime(&end_time) << endl;
        elapsed_seconds = end_tim - start_tim;
        cout << "Total Elapsed Time: " << elapsed_seconds.count() << "s" << endl << endl;

        /* Print final program status */
        cout << "\nLBE ended with status: \t\t" << status << endl;

        /* Make sure stdout is empty at the end */
        fflush(stdout);
    }

    return status;
}
