#include <getopt.h>

// Main library
#include "../include/option.h" //json

// Common libraries
#include "../include/utilities.h"
#include "../include/memory.h"

// Structures
#include "../include/general.h"


bool verbosity = false;
int rank = 0;

int main(int argc, char* argv[]){

    /////////////////////////////////////////////////////////
    // Variables declaration
    /////////////////////////////////////////////////////////

    // Simulation objects
    Config* cnf = Config_init();

    // Report
    printf("\n\n");
    printTitle("Before reading options");
    Config_report(cnf);

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
                cJSON_readOptionConfig(cnf,optarg);
                Config_report(cnf);
                break;

            case 'm':
                printf("    Read \"method\" option from command line ... \n");
                Config_setMethod(cnf,optarg);
                break;

            case 'q':
                Config_setQuantity(cnf,optarg);
                break;

            case 'I':
                ;
                // nbathfiles_current = getOptionsNbathfiles(options);
                // setOptionsBathfiles(options,optarg,nbathfiles_current+1);
                break;

            case 's':
                // setOptionsStatefile(options,optarg);
                break;
                
            case 'N':
                Config_setNstate(cnf,atoi(optarg));
                // nstate_current = getOptionsNbathfiles(options);
                // setOptionsNstate(options,nstate_current+1);
                break;

            case 'S':
                // setOptionsSubstatefile(options,optarg);
                break;

            case 'a':
                // setOptionsAvaaxfile(options,optarg);
                break;

            case 'B':
                Config_setBfield_z(cnf,atof(optarg));
                break;

            case 'o':
                // setOptionsOutputfile(options,optarg);
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
    Config_report(cnf);

    /////////////////////////////////////////////////////////
    // Free memory
    /////////////////////////////////////////////////////////
    Config_freeAll(cnf);

    printf(" Success! \n\n");
    return 0;
}