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
    double time_start_step;
    time_start = MPI_Wtime();

    // Print welcome message
    if (rank==0){
        printBanner();
    }

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

    //Cluster* cls_pcce = Cluster_init();
    //if (strcasecmp(cls->method, "pcce")==0){
    //    Cluster* cls_pcce = Cluster_init();
    //}
    //=======================================================
    // Read options
    //=======================================================
    int c;
    int nbathfiles_fromfccein = 0;
    int nbathfiles_current = 0;

    //  Declare and initialize the variables 
    // that will be used to store the arguments
    char method[MAX_CHARARRAY_LENGTH] = "";
    char quantity[MAX_CHARARRAY_LENGTH] = "";
    char** bathfiles = NULL;
    int nstate = 0;
    double bfield_z = 0.0;
    char* outfile = NULL;

    while ((c = getopt(argc, argv, "hvf:m:q:I:s:N:S:a:B:o:")) != -1) {
        switch (c) {
            case 'f':

                // set the parameter as the input file
                if (rank==0){ 
                    // print message
                    char message[500];
                    sprintf(message, "Option filename : '%s' \n", optarg);
                    printMessage(message);
                }
                fccein = optarg;

                // check if the file exists
                if (access(fccein, F_OK) == -1){
                    fprintf(stderr, "Error: Option file does NOT exist\n\n");
                    exit(EXIT_FAILURE);
                }

                // read the input file
                cJSON_readOptionConfig      (cnf, fccein); // general.h
                cJSON_readOptionQubitArray  (qa,  fccein); //qubit.h
                cJSON_readOptionDefectArray (dfa, fccein); // defect.h
                cJSON_readOptionCluster     (cls, fccein); // cluster.h
                cJSON_readOptionPulse       (pls, fccein); // pulse.h
                cJSON_readOptionOutput      (op,  fccein); // output.h

				if (rank==0){
					printf("\n	>> Read %s file successfully ..\n",optarg);
				}

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

            case 's':
                Config_freeStatefile(cnf);
                Config_allocStatefile(cnf);
                Config_setStatefile(cnf,optarg);

            case 'a':
                Config_freeAvaaxfile(cnf);
                Config_allocAvaaxfile(cnf);
                Config_setAvaaxfile(cnf,optarg);

            case 'S':
                Config_freeExstatefile(cnf);
                Config_allocExstatefile(cnf);
                Config_setExstatefile(cnf,optarg);

            case 'N':
                Config_setNstate(cnf,atoi(optarg));
                break;

            case 'B':
                Config_setBfield_z(cnf,atof(optarg));
                break;

            case 'o':
                Output_freeOutfile(op);
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
        printf("\n\n");
        printf("    ======================================================================\n");
        printf("        Report all structures \n");
        printf("    ======================================================================\n\n");

        printf("    ----------------------------------------------------------------------\n");
        Config_report(cnf); 
        printf("\n\n");

        printf("    ----------------------------------------------------------------------\n");
        QubitArray_report(qa);
        printf("\n\n");

        if (DefectArray_getNdefect(dfa) > 0){
            printf("    ----------------------------------------------------------------------\n");
            printTitle("DefectArray");
            for (int idf=0; idf<DefectArray_getNdefect(dfa); idf++){
                DefectArray_reportDefect_idf(dfa,idf);
                printf("\n");
            }
        }

        printf("    ----------------------------------------------------------------------\n");
        Cluster_report(cls);
        printf("\n\n");

        printf("    ----------------------------------------------------------------------\n");
        Pulse_report(pls);
        printf("\n\n");

        printf("    ----------------------------------------------------------------------\n");
        Output_report(op);
        printf("\n");

        printLineSection();
    }

    //=======================================================
    // Set simulators
    //=======================================================

    srand(seed);

    if (rank==0){
        printf("\n\n");
        printf("    ======================================================================\n");
        printf("        Read files \n");
        printf("    ======================================================================\n\n");
        time_start_step = MPI_Wtime();
    }

    // Qubit
    if (rank == 0){
        printSubTitle("Qubit file");
    }
    readQubitfile(qa,cnf); // set qubit xyz position for nqubit=1

    // Bath : Gyromagnetic ratio
    if (rank == 0){
        printSubTitle("Gyromagnetic ratio file");
    }
    readGyrofile(ba,cnf); // set gyro, spin, name (properties)

    // Bath : Configuration files
    if (rank == 0){
        printSubTitle("Bath configuration files");
    }
    readBathfiles(ba,qa,cnf);  // set bath xyz position and properties

    // Hyperfine tensor
    readHftensorfile(ba,qa,cnf); // set hyperfine tensor only from file

    // Quadrupole tensor
    // readQdtensorfile(ba,qa,cnf);

    // Defect
    if (DefectArray_getNdefect(dfa) > 0){

        if (rank == 0){
            printSubTitle("Defect Paxes File");
        }

        int nspin = BathArray_getNspin(ba);
        DefectArray_allocPaxes(dfa, nspin);
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

        // Update main bathspins : mainspidx, zfs(pax), detuning(pax)
        updateMainSpins_fromDefectArray(dfa,ba);
    }

    if (rank==0){
        printf("    ----------------------------------------------------------------------\n");
        printTitle("BathArray");
        BathArray_report(ba);
        printf("\n");
    }

    if (rank==0){
        // print wall clock time
        time_end = MPI_Wtime();
        printLine();
        printf("          Wall time = %.5f s\n", time_end - time_start);
        printLine();
    }
    
    //=======================================================
    // Clusterize
    //=======================================================

    if (rank==0){
        printf("\n\n");
        printf("   ======================================================================\n");
        printf("       Clusterize \n");
        printf("   ======================================================================\n\n");
        time_start_step = MPI_Wtime();
    }

    if (strcasecmp(cls->method, "pcce") ==0){
        Cluster_clusterize_pcce(cls, ba, qa, cnf, rank);
        //Cluster_clusterize_pcce(cls_pcce, cls, ba, qa, cnf, rank);
    }
    printf("In main.cpp: Finish the \" Cluster_clusterize_pcce \" function !! \n\n");
    MPI_Barrier(MPI_COMM_WORLD);
    //exit(1);
    

    if (rank==0){
        Cluster_reportClusinfo(cls);

        // print wall clock time
        time_end = MPI_Wtime();
        printLine();
        printf("          Wall time = %.5f s\n", time_end - time_start_step);
        printLine();
    }

    //=======================================================
    // Main calculation & save results
    //=======================================================

    if (rank==0){
        printf("\n\n");
        printf("   ======================================================================\n");
        printf("       Main calculation \n");
        printf("   ======================================================================\n\n");
        time_start_step = MPI_Wtime();
    }

    // MPI distribution for the clusters
    int order = Cluster_getOrder(cls);
    int*** clusters = Cluster_getClusinfo(cls);

    if (rank == 0){
        printf("rank: %d\n", rank);
        reportClusinfo(clusters,order);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //exit(1);

    int*** localclusters = MPI_getLocalClusters(order,clusters);
    // int*** localclusters = clusters; // No MPI
    // reportClusinfo(localclusters,order);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){
        printf("Done");
    }

    // Calculate the dynamics
    calculate(qa,ba,dfa,cnf,pls,cls,op,localclusters);
    MPI_Barrier(MPI_COMM_WORLD);
    exit(1);

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
