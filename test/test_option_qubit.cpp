#include <getopt.h>

// Main library
#include "../include/option.h" //json

// Common libraries
#include "../include/utilities.h"
#include "../include/memory.h"

// Structures
#include "../include/general.h"
#include "../include/qubit.h"


int verbosity = false;

int main(int argc, char* argv[]){

    /////////////////////////////////////////////////////////
    // Variables declaration
    /////////////////////////////////////////////////////////

    // Simulation objects
    QubitArray* qa = QubitArray_init();

    // Report
    printf("\n\n");
    printTitle("Before reading options");
    QubitArray_report(qa);

    /////////////////////////////////////////////////////////
    // Read options from command line
    /////////////////////////////////////////////////////////
    int c;
    int nbathfiles_current = 0;
    int nstate_current = 0;
    while ((c = getopt(argc, argv, "hvf:m:q:I:s:N:S:a:B:o:")) != -1) {
        switch (c) {
            case 'h':
                printf("Usage: %s [options]\n", argv[0]);
                printf("Options:\n");
                printf("  -h, --help\t\t\t\tPrint this message\n");

                printf("  -f, --file\t\t\t\tInput file (cce.in)\n");
                printf("  -I, --bathfiles\t\t\tInput file name (it can be more than one)\n");
                printf("  -s, --statefile\t\t\tFixed state input file\n");
                printf("  -S, --substatefile\t\t\tThe spin states of bath defect\n");
                printf("  -a, --avaaxfile\t\t\tAvailable axis of bath defect\n");
                printf("  -N, --nstate\t\t\t\tNumber of bath states\n");
                printf("  -m, --method\t\t\t\tMethod (ensemble | single | cce | gcce )\n");
                printf("  -q, --quantity\t\t\tQuantity (coherence | noise | dm)\n");
                printf("  -B, --bfield\t\t\t\tMagnetic field\n");
                printf("  -v, --verbosity\t\t\tLong or short print\n");
                printf("  -o, --outputfile\t\t\tOutput file name\n");
                exit(EXIT_SUCCESS);
                break;
            
            case 'f':
                /////////////////////////////////////////////////////////
                // Read options
                /////////////////////////////////////////////////////////
                printf("\n");
                printf("    Read options from %s file ... \n",optarg);
                cJSON_readOptionQubitArray(qa,optarg);
                break;

            case 'm':
                break;

            case 'q':
                break;

            case 'I':
                break;

            case 's':
                break;
                
            case 'N':
                break;

            case 'S':
                break;

            case 'a':
                break;

            case 'B':
                break;

            case 'o':
                break;

            case 'v':
                verbosity = true;
                break;

            default:
                fprintf(stderr, "Error: unknown option\n");
                exit(EXIT_FAILURE);
        }   
    }   

    /////////////////////////////////////////////////////////
    // Report
    /////////////////////////////////////////////////////////
    printf("\n");
    printTitle("After reading options");
    QubitArray_report(qa);

    /////////////////////////////////////////////////////////
    // Free memory
    /////////////////////////////////////////////////////////
    QubitArray_freeAll(qa);

    printf(" Success! \n\n");
    return 0;
}