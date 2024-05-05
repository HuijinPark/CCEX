#include "../include/output.h"
#include "../include/memory.h"
#include "../include/utilities.h"

bool verbosity = true;
int rank = 0;

int main(){

    Output* output = Output_init();
    Output_report(output);

    // Matrix
    int nstep = 10;
    float deltat = 0.1;
    MatrixXcd* res = new MatrixXcd[nstep];
    for (int i=0; i<nstep; i++){
        res[i] = MatrixXcd::Random(3,3);
    }

    // alloc
    Output_allocOutfilehead(output);
    Output_allocOutfile(output);

    char* outputhead = allocChar1d(MAX_FILEPATH);
    strcpy(outputhead, "./output");
    Output_setOutfilehead(output, outputhead);

    // all
    Output_setSavemode(output, "all");

    int nspin = 4;
    int* cluster = allocInt1d(nspin);
    cluster[0] = 1;
    cluster[1] = 2;
    cluster[2] = 3;
    cluster[3] = 4;

    int istate = 1;
    Output_setOutfile_all(output, cluster, nspin, istate);
    Output_report(output);

    // set
    Output_save(output, res, res, nstep, deltat);
    Output_report(output);

    // info
    Output_setSavemode(output, "info");
    double Azx = 1.0;
    double Azz = 2.0;
    char bathfile[MAX_FILEPATH] = "bathfile";
    Output_setOutfile_info(output, istate);
    Output_save_info(output, res, res, nstep, deltat, Azx, Azz, bathfile);
    Output_report(output);


    freeInt1d(cluster);
    freeChar1d(outputhead);

    Output_freeAll(output);
    delete[] res;

    return 0;
}