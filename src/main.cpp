// ReadOption
#include "../include/option.h"

// Structures
#include "../include/qubit.h"
#include "../include/bath.h"
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
    time_t start, end;
    time(&start);

    //=======================================================
    // Variables declaration
    //=======================================================

    // Input file (cce.in)
    char* fccein = allocChar1d(MAX_FILEPATH);

    //Simulation objects
    Config*     cnf = Config_init();
    QubitArray* qa  = QubitArray_init();
    Cluster*    cls = Cluster_init();
    BathArray*  ba  = BathArray_init();
    Pulse*      pls = Pulse_init();
    Output*     op  = Output_init();

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
                strcpy(fccein,optarg);
                cJSON_readOptionConfig(cnf, fccein); // general.h
                cJSON_readOptionQubitArray(qa, fccein); //qubit.h
                cJSON_readOptionCluster(cls, fccein); // cluster.h
                cJSON_readOptionPulse(pls, fccein); // pulse.h
                cJSON_readOptionOutput(op, fccein); // output.h

                nbathfiles_fromfccein = Config_getNbathfiles(cnf);
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
                Output_freeOutfilehead(op);
                Output_allocOutfilehead(op);
                Output_setOutfilehead(op,optarg);
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

    if (nbathfiles_current != nbathfiles_fromfccein){
        fprintf(stderr, "Error: number of bath files is not consistent\n");
        exit(EXIT_FAILURE);
    }

    // check if the file exists
    if (access(fccein, F_OK) == -1){
        fprintf(stderr, "Error: fccein does not exist\n");
        exit(EXIT_FAILURE);
    }
    
    //=======================================================
    // Report options
    //=======================================================
    Config_report(cnf);
    QubitArray_report(qa);
    Cluster_report(cls);
    Pulse_report(pls);
    Output_report(op);


    //=======================================================
    // Set simulators
    //=======================================================

    // Qubit
    readQubitfile(qa,cnf); // set qubit xyz position for nqubit=1

    // Bath
    readGyrofile(ba,cnf); // set gyro, spin, name (properties)
    readBathfiles(ba,qa,cnf);  // set bath xyz position and properties
    readHftensorfile(ba,qa,cnf); // set hyperfine tensor only from file
    // readQdtensorfile(ba,qa,cnf);

    BathArray_report(ba);

    // point dipole approximation if no hyperfine tensor
    int nqubit = QubitArray_getNqubit(qa);
    BathArray_setBathHypfs(ba,qa); 
    BathArray_reportBath_hypf(ba, nqubit); //report updated hyperfine tensors

    // Cluster
    Cluster_clusterize(cls,ba,cnf);
    Cluster_reportClusinfo(cls);


    //=======================================================
    // Main calculation & save results
    //=======================================================
    calculate(qa,ba,cnf,pls,cls,op);
    
    //=======================================================
    // Free memory
    //=======================================================

    QubitArray_freeAll(qa);
    Cluster_freeAll(cls);
    Config_freeAll(cnf);
    BathArray_freeAll(ba);
    Pulse_freeAll(pls);
    Output_freeAll(op);
    
    /////////////////////////////////////////////////////////
    MPI_Finalize();

    //wall time check
    time(&end);
    if (rank==0){ 
        char message[500];
        sprintf(message, "JOB DONE : Wall time = %f\n", difftime(end, start));

        printLineSection();
        printSubTitle(message);
        printLineSection();
        printf("\n");
    }

    return 0;
}