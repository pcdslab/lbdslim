# Load Balancing Algorithm for Distributed SLIM-Transform (DSLIM)
## Authors
Muhammad Haseeb and Fahad Saeed
## Institution
Florida International University, Miami, FL, USA.

# What do you need?
1. make
2. GCC v5.4.0 or later with C++11+ support
2. Windows: MinGW (x86_64-6.3.0-win32-seh-rt_v5-rev2)

# Pre LBDSLIM Steps
1. Digest the proteome database using Protein Digestion Simulator or OpenMS.
2. Remove the redundant peptide sequences in the digested database using DBToolkit.
3. Optional: The decoy database can be generated and appended to the target database using DBToolkit.
4. Convert the MS/MS data in (mzML/mzXML/MS2) format (MS2 preferably) using msconvert.exe

# The Sample Driver Application
The sample application shows the software pipeline for successfully incorporating LBE into SLM-Transform for peptide search. The application first initializes SLM Peptide Index (SPI) which is then distributed using LBE policy among available compute units (OpenMP threads). The threads construct distributed SLM Ions Index (SII) and in turn DSLIM. The application then analyzes the constructed index to evaluate the ion distribution among DSLIM chunks. The query MS/MS data are extracted and processed by the main threads and then each threads searches the query data on only 1 OpenMP thread (to simulate distributed machine based querying). The number of matches obtained from each DSLIM chunk are returned which are further analyzed to evaluate the computational load imbalance in the system.

# Configure LBDSLIM & Sample Application
1. Configure the SLM-Transform parameters in /lbdslim/include/config.h: 
*Make sure that #define/undef WINDOWS is correct according to your Host OS else there will be Seg faults.*
2. Add the database, dataset and mods information in the ./slm.cpp (Sample Application) using the following format:

## Format

    STRING modconditions = "#1 <res> #2"; 
    /* #1 is maximum allowed modified residues per sequence, 
    <res> modified residues of 1st type, 
    #2 how many residues of this modification in peptide sequence. 
    
    For example: we want max 5 mod residues; and max STY 3 and max M 2 we will set modconditions as
    "5 STY 3 M 2"
    */
    
    <some code> 
    
    /* Database and Dataset files */
    STRING filename = "/path/to/digested/and/unduplicated/database.fasta"; // Path to processed database
    STRING querypath = "/path/to/query/dataset"; // Path to Dataset folder
    STRING patt = ".ms2/mzML/mzXML"; // Keep one, remove others
    
    <some code> 
    
    /* Initialize the vModInfo */
    vModInfo.num_vars = #3; // Types of modifications been specified (max 7) for above example, set this to: 2
    vModInfo.vmods_per_pep = #1; // #1 from the previous explanation, from above example, set this to: 5
    
    /* List of Mods Info */
    vModInfo.vmods[0].aa_per_peptide = 3// #2 from previous explanation, for above example, set this to: 3
    vModInfo.vmods[0].modMass = 79.97 * SCALE; // Mass of the specified modification, for above example: 79.97
    vModInfo.vmods[0].residues[0] = 'S'; // modified residues list (max 4 per mod type allowed), for above example S
    vModInfo.vmods[0].residues[1] = 'T'; // T
    vModInfo.vmods[0].residues[2] = 'Y'; // Y
    
    vModInfo.vmods[0].aa_per_peptide = 2; // #2 from previous explanation, for above example, set this to: 2
    vModInfo.vmods[0].modMass = 15.997 * SCALE; // Mass of the specified modification, for above example: 15.997
    vModInfo.vmods[0].residues[0] = 'M'; // modified residues list (max 4 per mod type allowed), for above example M
    
    <some code>
    /* Distribute peptides among cores */
    status = LBE_Distribute(threads, _chunk, slm_chunks); // Use either _chunk, _cyclic or _random to specify the distribution policy

# Building LBDSLIM Sample Application
1. Open the Git Bash shell (Windows) or normal Terminal (Ubuntu).
1. Navigate to LBDSLIM home directory /lbdslim/
2. Execute the following command: `make`

# Running LBDSLIM Sample Application
1. Navigate to /lbdslim
2. ./LoadBalancer.exe

# Please Note:
1. Work is being done to move LBDSLIM configuration to runtime so please bear with us.
2. Max digested peptide mass allowed: 5000Da
3. The LBDSLIM only returns the fragment-ion filtered PSMs (not the final PSMs), which can be used in any way e.g. post filtered based on precursor masses or sequence tags, formally scored for PSMs, FDRed.
4. Please contact us if you experience any trouble or bugs in the software. Thanks. :)

# Please cite our work
For queries or questions about LBDSLIM, please contact: {fsaeed, mhaseeb}@fiu.edu. Thank you.

