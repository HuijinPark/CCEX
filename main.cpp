// ReadOption
#include "../include/option.h"

// Structures
#include "../include/qubit.h"
#include "../include/bath.h"
#include "../include/defect.h"
#include "../include/cluster.h"
#include "../include/pulse.h"
#include "../include/output.h"
#include "../include/general.h"

// Simulators
#include "../include/simulator.h"

// Utilities
#include "../include/utilities.h"
#include "../include/memory.h"
#include <getopt.h>

// Read files
#include "../include/reader.h"

//global variables
int seed; // for random genearator
int rank,nprocess;
bool verbosity = false;

int main(int argc, char* argv[]){

    //Set up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get process id
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess); // Get number of processes

    // wall time check
    double time_start, time_end;
    time_start = MPI_Wtime();
    //=======================================================
    // Variables declaration
    //=======================================================

    // Input file (cce.in)
    char* fccein = NULL; 

    //Simulation objects
    Config*      cnf = Config_init();
    QubitArray*  qa  = QubitArray_init();
    Cluster*     cls = Cluster_init();
    BathArray*   ba  = BathArray_init();
    DefectArray* dfa = DefectArray_init();
    Pulse*       pls = Pulse_init();
    Output*      op  = Output_init();

    //=======================================================
    // Read options
    //=======================================================
    int c;
    int nbathfiles_fromfccein = 0;
    int nbathfiles_current = 0;
    int nstate_current = 0;

    while ((c = getopt(argc, argv, "hvf:m:q:I:s:N:S:a:B:o:")) != -1) {
        switch (c) {
            case 'f':
                
                fccein = allocChar1d(MAX_FILEPATH); 
                strcpy(fccein,optarg);

                // check if the file exists
                if (access(fccein, F_OK) == -1){
                    fprintf(stderr, "Error: fccein does not exist\n");
                    exit(EXIT_FAILURE);
                }
                
                cJSON_readOptionConfig      (cnf, fccein); // general.h
                cJSON_readOptionQubitArray  (qa,  fccein); //qubit.h
                cJSON_readOptionDefectArray (dfa, fccein); // defect.h
                cJSON_readOptionCluster     (cls, fccein); // cluster.h
                cJSON_readOptionPulse       (pls, fccein); // pulse.h
                cJSON_readOptionOutput      (op,  fccein); // output.h

                nbathfiles_fromfccein = Config_getNbathfiles(cnf);
                nbathfiles_current += nbathfiles_fromfccein;
                break;

            case 'm':
                Config_setMethod(cnf,optarg);
                Cluster_setMethod(cls,optarg);
                break;

            case 'q':
                Config_setQuantity(cnf,optarg);
                break;

            case 'I':
                nbathfiles_current++;
                if (nbathfiles_current == 1){
                    Config_freeBathfiles(cnf);
                    Config_setNbathfiles(cnf,nbathfiles_current);
                    Config_allocBathfiles(cnf);
                    Config_setBathfiles_i(cnf,optarg,nbathfiles_current-1);
                }
                else{
                    Config_setNbathfiles(cnf,nbathfiles_current);
                    Config_reallocBathfiles(cnf, nbathfiles_current-1, nbathfiles_current);
                    Config_setBathfiles_i(cnf,optarg,nbathfiles_current-1);
                }
                break;
                
            case 'N':
                Config_setNstate(cnf,atoi(optarg));
                break;

            case 'B':
                Config_setBfield_z(cnf,atof(optarg));
                break;

            case 'o':
                Output_allocOutfile(op);
                Output_setOutfile(op,optarg);
                break;

            case 'v':
                verbosity = true;
                break;

            case 'h':
                printf("Usage: %s [options]\n", argv[0]);
                printf("Options:\n");
                printf("  -h, --help\t\t\t\tPrint this message\n");
                printf("  -f, --file\t\t\t\tInput file (cce.in)\n");
                printf("  -I, --bathfiles\t\t\tInput file name (it can be more than one)\n");
                printf("  -N, --nstate\t\t\t\tNumber of bath states\n");
                printf("  -m, --method\t\t\t\tMethod (ensemble | single | cce | gcce )\n");
                printf("  -q, --quantity\t\t\tQuantity (coherence | noise | dm)\n");
                printf("  -B, --bfield\t\t\t\tMagnetic field\n");
                printf("  -v, --verbosity\t\t\tLong or short print\n");
                printf("  -o, --outputfile\t\t\tOutput file name\n");
                exit(EXIT_SUCCESS);
                break;
            default:
                fprintf(stderr, "Error: unknown option\n");
                exit(EXIT_FAILURE);
        }   
    }
    
    if (fccein == NULL){
        fprintf(stderr, "Error: input file is required\n");
        exit(EXIT_FAILURE);
    }

    if (nbathfiles_current != nbathfiles_fromfccein){
        fprintf(stderr, "Error: number of bath files is not consistent\n");
        fprintf(stderr, "nbathfiles_fromfccein = %d\n", nbathfiles_fromfccein);
        fprintf(stderr, "nbathfiles_current = %d\n", nbathfiles_current);
        exit(EXIT_FAILURE);
    }
    
    //=======================================================
    // Report options
    //=======================================================
    if (rank==0){
        Config_report(cnf);
        QubitArray_report(qa);

        if (DefectArray_getNdefect(dfa) > 0){
            for (int idf=0; idf<DefectArray_getNdefect(dfa); idf++){
                DefectArray_reportDefect_idf(dfa,idf);
            }
        }

        Cluster_report(cls);
        Pulse_report(pls);
        Output_report(op);
    }

    //=======================================================
    // Set simulators
    //=======================================================

    srand(seed);

    // Qubit
    if (rank == 0){
        printSubTitle("Qubit File");
    }
    readQubitfile(qa,cnf); // set qubit xyz position for nqubit=1

    // Bath
    if (rank == 0){
        printSubTitle("Bath related File");
    }
    readGyrofile(ba,cnf); // set gyro, spin, name (properties)
    readBathfiles(ba,qa,cnf);  // set bath xyz position and properties
    readHftensorfile(ba,qa,cnf); // set hyperfine tensor only from file
    // readQdtensorfile(ba,qa,cnf);

    // Bath State
    // readStatefile(ba,cnf); // set bath states
    
    // Extra spins
    // readAvaaxfile(ba,cnf); // set available principal axis
    // readExstatefile(ba,cnf); // set extra spin states

    // if (rank==0){
    //     printLineSection();
    //     printTitle("BathArray");
    //     BathArray_report(ba);
    //     printf("\n");
    //     printLineSection();
    // }

    // point dipole approximation if no hyperfine tensor
    int nqubit = QubitArray_getNqubit(qa);
    BathArray_setBathHypfs(ba,qa); 
    if (rank==0){
        BathArray_reportBath_hypf(ba, nqubit); //report updated hyperfine tensors
    }

    // Defect
    if (DefectArray_getNdefect(dfa) > 0){

        if (rank == 0){
            printSubTitle("Defect Paxes File");
        }

        int nspin = BathArray_getNspin(ba);
        DefectArray_allocPaxes(dfa, nspin);
        // DefectArray_setPaxesRandom(dfa,ba);
        setDefectPaxes(dfa,ba,cnf); // Read or random set
        if (rank==0){
            DefectArray_reportPaxes(dfa);
        }

        DefectArray_allocNaddspins(dfa, nspin);
        DefectArray_setNaddspins(dfa,ba);
        if (rank==0){
            DefectArray_reportNaddspins(dfa);
        }

        int nqubit = QubitArray_getNqubit(qa);
        DefectArray_allocSubbath(dfa,ba,nqubit);
        DefectArray_setSubbath(dfa,ba,qa);
        if (rank==0){
            DefectArray_reportSubbath_props(dfa);
            DefectArray_reportSubbath_hypfs(dfa, nqubit);
            DefectArray_reportSubbath_quads(dfa);
            DefectArray_reportSubbath_hypf_subs(dfa);
        }
    }
    

    // Cluster
    Cluster_clusterize(cls,ba,cnf);
    if (rank==0){
        Cluster_reportClusinfo(cls);
    }

    //=======================================================
    // Main calculation & save results
    //=======================================================
    calculate(qa,ba,dfa,cnf,pls,cls,op);
    
    //=======================================================
    // Free memory
    //=======================================================

    QubitArray_freeAll(qa);
    Cluster_freeAll(cls);
    Config_freeAll(cnf);
    BathArray_freeAll(ba);
    DefectArray_freeAll(dfa);
    Pulse_freeAll(pls);
    Output_freeAll(op);
    
    /////////////////////////////////////////////////////////
    

    //wall time check
    if (rank==0){ 
        time_end = MPI_Wtime();
        char message[500];
        sprintf(message, "JOB DONE : Wall time = %.5f s", time_end - time_start);
        printLineSection();
        printMessage(message);
        printLineSection();
        printf("\n");
    }

    MPI_Finalize();

    return 0;
}
