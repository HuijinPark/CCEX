#include "../include/output.h"
#include "../include/memory.h"

Output* Output_init(){
    Output* op = (Output*)allocArray1d(1,sizeof(Output));
    op->savemode[0] = '\0';
    op->outfile = NULL;
    return op;
}

void Output_save(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, float deltat, int istate){

    char* outfile = Output_getOutfile(op);

    // istate to string
    char istate_str[10];
    if (istate==0){
        istate_str[0] = '\0';
    }
    else{
        sprintf(istate_str, "_state%d", istate);
    }

    int dim = result_wD[0].rows();

    char* outfile_wD = allocChar1d(MAX_FILEPATH);
    char* outfile_nD = allocChar1d(MAX_FILEPATH);

    strcpy(outfile_wD, outfile);
    strcat(outfile_wD, istate_str);
    strcat(outfile_wD, "_wiDiv");

    strcpy(outfile_nD, outfile);
    strcat(outfile_nD, istate_str);
    strcat(outfile_nD, "_noDiv");

    FILE* fp_wD = fopen(outfile_wD, "w");
    FILE* fp_nD = fopen(outfile_nD, "w");
    
    if (fp_wD == NULL || fp_nD == NULL){
        printf("Error: cannot open file %s\n", outfile_wD);
        printf("Error: cannot open file %s\n", outfile_nD);
        exit(1);
    }

    // write format :
    // for all nstate , write time dm00 dm01 dm02 ... dmnn (write only)
    for (int i=0; i<nstep; i++){
        int c = 0;
        fprintf(fp_wD, "%15g\t", (float)i*deltat);
        fprintf(fp_nD, "%15g\t", (float)i*deltat);
        for (int j=0; j<dim; j++){
            for (int k=0; k<dim; k++){

                fprintf(fp_wD, "%+15.10lf", result_wD[i](j,k).real());
                fprintf(fp_nD, "%+15.10lf", result_nD[i](j,k).real());

                if (std::isnan(result_wD[i](j,k).imag())){
                    fprintf(fp_wD, "%+-4lfj%10s\t", result_wD[i](j,k).imag()," ");
                }else{
                    fprintf(fp_wD, "%+-13.10lfj\t", result_wD[i](j,k).imag());
                }

                if (std::isnan(result_nD[i](j,k).imag())){
                    fprintf(fp_nD, "%+-4lfj%10s\t", result_nD[i](j,k).imag()," ");
                }else{
                    fprintf(fp_nD, "%+-13.10lfj\t", result_nD[i](j,k).imag());
                }

                if ((c%3==2) && (j*k!=(dim-1)*(dim-1))){
                    fprintf(fp_wD, "\n");
                    fprintf(fp_wD, "%15s\t"," ");

                    fprintf(fp_nD, "\n");
                    fprintf(fp_nD, "%15s\t"," ");
                }

                c++;
            }
        }
        fprintf(fp_wD, "\n");
        fprintf(fp_nD, "\n");
    }

    fclose(fp_wD);
    fclose(fp_nD);

    freeChar1d(&outfile_wD);
    freeChar1d(&outfile_nD);
    
    return;
}

void Output_save_all(Output* op, MatrixXcd* result, int nstep, float deltat, int* cluster, int nspin, int istate){

    //////////////////////////////
    // Save the original one
    char* outfile_original = allocChar1d(MAX_FILEPATH);
    strcpy(outfile_original, Output_getOutfile(op));
    //////////////////////////////

    //////////////////////////////
    // Change outfile name

    // cluster to string
    char cluster_str[100] = "";
    for (int i=0; i<nspin; i++){
        char spin_str[100];
        sprintf(spin_str, "_%d", cluster[i]);
        strcat(cluster_str, spin_str);
    }

    // add string outfilehead, cluster, istate to outfile
    strcat(op->outfile, cluster_str);
    
    Output_save(op, result, result, nstep, deltat, istate);
    //////////////////////////////

    //////////////////////////////
    // Restore outfile name
    strcpy(op->outfile, outfile_original);
    //////////////////////////////

    freeChar1d(&outfile_original);

    return;
}

void Output_save_info(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, float deltat, int istate, double Azx, double Azz, char* bathfile){

    char* savemode = Output_getSavemode(op);

    if (strcasecmp(savemode, "info") != 0){
        printf("Error: savemode is not info\n");
        exit(1);
    }

    char* outfile = Output_getOutfile(op);

    // istate to string
    char istate_str[10];
    if (istate==0){
        istate_str[0] = '\0';
    }
    else{
        sprintf(istate_str, "_state%d", istate);
    }


    int dim = result_wD[0].rows();

    char* outfile_wD = allocChar1d(MAX_FILEPATH);
    char* outfile_nD = allocChar1d(MAX_FILEPATH);

    strcpy(outfile_wD, outfile);
    strcat(outfile_wD, istate_str);
    strcat(outfile_wD, "_wiDiv");

    strcpy(outfile_nD, outfile);
    strcat(outfile_nD, istate_str);
    strcat(outfile_nD, "_noDiv");

    FILE* fp_wD = fopen(outfile_wD, "w");
    FILE* fp_nD = fopen(outfile_nD, "w");
    
    if (fp_wD == NULL){
        printf("Error: cannot open file %s\n", outfile_wD);
        exit(1);
    }

    for (int i=0; i<nstep; i++){
        int c = 0;
        fprintf(fp_wD, "%15g\t", (float)i*deltat);
        fprintf(fp_nD, "%15g\t", (float)i*deltat);
        for (int j=0; j<dim; j++){
            for (int k=0; k<dim; k++){

                fprintf(fp_wD, "%+15.10lf", result_wD[i](j,k).real());
                fprintf(fp_nD, "%+15.10lf", result_nD[i](j,k).real());

                if (std::isnan(result_wD[i](j,k).imag())){
                    fprintf(fp_wD, "%+-4lfj%10s\t", result_wD[i](j,k).imag()," ");
                }else{
                    fprintf(fp_wD, "%+-13.10lfj\t", result_wD[i](j,k).imag());
                }

                if (std::isnan(result_nD[i](j,k).imag())){
                    fprintf(fp_nD, "%+-4lfj%10s\t", result_nD[i](j,k).imag()," ");
                }else{
                    fprintf(fp_nD, "%+-13.10lfj\t", result_nD[i](j,k).imag());
                }

                if ((c%3==2) && (j*k!=(dim-1)*(dim-1))){
                    fprintf(fp_wD, "\n");
                    fprintf(fp_wD, "%15s\t"," ");

                    fprintf(fp_nD, "\n");
                    fprintf(fp_nD, "%15s\t"," ");
                }

                c++;
            }
        }
        // write Azx, Azz, bathfile
        fprintf(fp_wD, "%+.10lf\t%+.10lf\t%s\n", Azx, Azz, bathfile);
        fprintf(fp_nD, "%+.10lf\t%+.10lf\t%s\n", Azx, Azz, bathfile);
    }

    fclose(fp_wD);
    fclose(fp_nD);
    
    freeChar1d(&outfile_wD);
    freeChar1d(&outfile_nD);
    return;
}

void Output_setSavemode(Output* op, char* savemode){

    const int opsize = 4;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"all", "normal", "info"};
    int idx = findIndexCharFix(options,0,opsize-1,savemode);
    if (idx == -1) {
        fprintf(stderr, "Error: current savemode options (%s) is not available\n",savemode);
        exit(EXIT_FAILURE);
    }
    strcpy(op->savemode, savemode);
}

char* Output_getSavemode(Output* op){
    return op->savemode;
}

void Output_allocOutfile(Output* op){
    op->outfile = allocChar1d(MAX_FILEPATH);
}

void Output_freeOutfile(Output* op){
    freeArray1d((void**)&(op->outfile));
}

void Output_setOutfile(Output* op, char* outfile){
    strcpy(op->outfile, outfile);
}

char* Output_getOutfile(Output* op){
    return op->outfile;
}

void Output_freeAll(Output* op){
    Output_freeOutfile(op);
    freeArray1d((void**)&op);
}

void Output_report(Output* op){

    printTitle("Structure Output");

    printStructElementChar("savemode", op->savemode);
    printStructElementChar("outfile", op->outfile);

}
