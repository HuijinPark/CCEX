#include "../include/general.h"
#include "../include/memory.h"
#include <string.h>

// init
Config* Config_init(){
    Config* cnf = (Config*)allocArray1d(1,sizeof(Config));
    cnf->method[0] = '\0';
    cnf->quantity[0] = '\0';
    return cnf;
}

// alloc
void Config_allocBathfiles(Config* cnf){
    cnf->bathfiles = allocChar2d(cnf->nbathfiles,MAX_FILEPATH);
}

void Config_reallocBathfiles(Config* cnf, int oldlength, int newlength){
    reallocChar2d(&(cnf->bathfiles),oldlength,newlength,MAX_FILEPATH);
}

void Config_allocBathadjust(Config* cnf){
    cnf->bathadjust = allocDouble2d(cnf->nbathfiles,3);
}

void Config_allocGyrofile(Config* cnf){
    cnf->gyrofile = allocChar1d(MAX_FILEPATH);
}

void Config_allocQubitfile(Config* cnf){
    cnf->qubitfile = allocChar1d(MAX_FILEPATH);
}

void Config_allocHf_tensorfile(Config* cnf){
    cnf->hf_tensorfile = allocChar1d(MAX_FILEPATH);
}

void Config_allocQd_tensorfile(Config* cnf){
    cnf->qd_tensorfile = allocChar1d(MAX_FILEPATH);
}

void Config_allocQd_tensorfile_woqubit(Config* cnf){
    cnf->qd_tensorfile_woqubit = allocChar1d(MAX_FILEPATH);
}

// free
void Config_freeAll(Config* cnf){
    freeArray1d((cnf));
}

void Config_freeBathfiles(Config* cnf){
    freeChar2d(cnf->bathfiles,cnf->nbathfiles);
}

void Config_freeBathadjust(Config* cnf){
    freeDouble2d(cnf->bathadjust,cnf->nbathfiles);
}

void Config_freeGyrofile(Config* cnf){
    freeChar1d(cnf->gyrofile);
}

void Config_freeQubitfile(Config* cnf){
    freeChar1d(cnf->qubitfile);
}

void Config_freeHf_tensorfile(Config* cnf){
    freeChar1d(cnf->hf_tensorfile);
}

void Config_freeQd_tensorfile(Config* cnf){
    freeChar1d(cnf->qd_tensorfile);
}

void Config_freeQd_tensorfile_woqubit(Config* cnf){
    freeChar1d(cnf->qd_tensorfile_woqubit);
}

/* Low level ---------------------------------------------------------*/

// get
char* Config_getMethod(Config* cnf){
    return cnf->method;
}

char* Config_getQuantity(Config* cnf){
    return cnf->quantity;
}

int   Config_getOrder(Config* cnf){
    return cnf->order;
}

float* Config_getBfield(Config* cnf){
    return cnf->bfield;
}

float Config_getRbath(Config* cnf){
    return cnf->rbath;
}

float Config_getRdip(Config* cnf){
    return cnf->rdip;
}

float Config_getDeltat(Config* cnf){
    return cnf->deltat;
}

int   Config_getNstep(Config* cnf){
    return cnf->nstep;
}

float Config_getRbathcut(Config* cnf){
    return cnf->rbathcut;
}

float Config_getRdipcut(Config* cnf){
    return cnf->rdipcut;
}

int   Config_getNstate(Config* cnf){
    return cnf->nstate;
}

int   Config_getSeed(Config* cnf){
    return cnf->seed;
}

char*   Config_getQubitfile(Config* cnf){
    return cnf->qubitfile;
}

char*   Config_getGyrofile(Config* cnf){
    return cnf->gyrofile;
}

int     Config_getNbathfiles(Config* cnf){
    return cnf->nbathfiles;
}

char*   Config_getBathfiles_i(Config* cnf,int i){
    return cnf->bathfiles[i];
}

double* Config_getBathadjust_i(Config* cnf,int i){
    return cnf->bathadjust[i];
}

double  Config_getDefectTotSpin(Config* cnf){
    return cnf->DefectTotSpin;
}

double  Config_getCorrTotSpin(Config* cnf){
    return cnf->CorrTotSpin;
}

char*   Config_getHf_tensorfile(Config* cnf){
    return cnf->hf_tensorfile;
}

double  Config_getHf_cutoff(Config* cnf){
    return cnf->hf_cutoff;
}

int     Config_getHf_ignore_oor(Config* cnf){
    return cnf->hf_ignore_oor;
}

int     Config_getHf_readmode(Config* cnf){
    return cnf->hf_readmode;
}

char*   Config_getQd_tensorfile(Config* cnf){
    return cnf->qd_tensorfile;
}

char*   Config_getQd_tensorfile_woqubit(Config* cnf){
    return cnf->qd_tensorfile_woqubit;
}

int     Config_getQd_readmode(Config* cnf){
    return cnf->qd_readmode;
}

// set 
void Config_setMethod(Config* cnf, char* method){

    const int opsize = 7;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"gcce","cce","dsj","itb","dsjitb","kmeans"};
    int idx = findIndexCharFix(options,0,opsize-1,method);
    if (idx == -1) {
        fprintf(stderr, "Error: current method options (%s) is not available\n",method);
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->method,method);
}

void Config_setQuantity(Config* cnf, char* quantity){

    const int opsize = 3;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"coherence","dm","noise"};
    int idx = findIndexCharFix(options,0,opsize-1,quantity);
    if (idx == -1) {
        fprintf(stderr, "Error: current quantity options (%s) is not available\n",quantity);
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->quantity,quantity);
}

void Config_setOrder(Config* cnf, int order){

    if (order < 0) {
        fprintf(stderr, "Error: current order (%d < 0) is not available\n",order);
        exit(EXIT_FAILURE);
    }

    cnf->order = order;
}

void Config_setBfield(Config* cnf, float* bfield){
    copyFloat1d(cnf->bfield,bfield,3);
}

void Config_setBfield_z(Config* cnf, float bz){
    cnf->bfield[2] = bz;
}

void Config_setRbath(Config* cnf, float rbath){

    if (rbath < 0.0) {
        fprintf(stderr, "Error: current rbath (%f < 0) is not available\n",rbath);
        exit(EXIT_FAILURE);
    }
    cnf->rbath = rbath;
}

void Config_setRdip(Config* cnf, float rdip){

    if (rdip < 0.0) {
        fprintf(stderr, "Error: current rdip (%f < 0) is not available\n",rdip);
        exit(EXIT_FAILURE);
    }
    cnf->rdip = rdip;
}

void Config_setDeltat(Config* cnf, float deltat){

    if (deltat < 0.0) {
        fprintf(stderr, "Error: current deltat (%f < 0) is not available\n",deltat);
        exit(EXIT_FAILURE);
    }
    cnf->deltat = deltat;
}

void Config_setNstep(Config* cnf, int nstep){

    if (nstep < 0) {
        fprintf(stderr, "Error: current nstep (%d < 0) is not available\n",nstep);
        exit(EXIT_FAILURE);
    }
    cnf->nstep = nstep;
}

void Config_setRbathcut(Config* cnf, float rbathcut){

    if (rbathcut < 0.0) {
        fprintf(stderr, "Error: current rbathcut (%f < 0) is not available\n",rbathcut);
        exit(EXIT_FAILURE);
    }
    cnf->rbathcut = rbathcut;
}

void Config_setRdipcut(Config* cnf, float rdipcut){

    if (rdipcut < 0.0) {
        fprintf(stderr, "Error: current rdipcut (%f < 0) is not available\n",rdipcut);
        exit(EXIT_FAILURE);
    }
    cnf->rdipcut = rdipcut;
}

void Config_setNstate(Config* cnf, int nstate){

    if (nstate < 0) {
        fprintf(stderr, "Error: current nstate (%d < 0) is not available\n",nstate);
        exit(EXIT_FAILURE);
    }
    cnf->nstate = nstate;
}

void Config_setSeed(Config* cnf, int seed){

    if (seed < 0) {
        fprintf(stderr, "Error: current seed (%d < 0) is not available\n",seed);
        exit(EXIT_FAILURE);
    }
    cnf->seed = seed;
}

void Config_setQubitfile(Config* cnf, char* qubitfile){

    if (cnf->qubitfile == NULL){
        fprintf(stderr, "Error: the memory for qubitfile should be allocated\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->qubitfile,qubitfile);
}

void Config_setGyrofile(Config* cnf, char* gyrofile){

    if (cnf->gyrofile == NULL){
        fprintf(stderr, "Error: the memory for gyrofile should be allocated\n");
        exit(EXIT_FAILURE);
    }
    strcpy(cnf->gyrofile,gyrofile);
}

void Config_setNbathfiles(Config* cnf, int nbathfiles){

    if (nbathfiles <= 0) {
        fprintf(stderr, "Error: current nbathfiles(%d) (<= 0) is not available\n",nbathfiles);
        exit(EXIT_FAILURE);
    }
    cnf->nbathfiles = nbathfiles;
}

void Config_setBathfiles_i(Config* cnf, char* bathfiles, int i){

    if (cnf->bathfiles == NULL){
        fprintf(stderr, "Error: the memory for bathfiles should be allocated\n");
        exit(EXIT_FAILURE);
    }
    strcpy(cnf->bathfiles[i],bathfiles);
}

void Config_setBathadjust_i(Config* cnf, double* bathadjust, int i){

    if (cnf->bathadjust == NULL){
        fprintf(stderr, "Error: the memory for bathadjust should be allocated\n");
        exit(EXIT_FAILURE);
    }
    copyDouble1d(cnf->bathadjust[i],bathadjust,3);
}

void Config_setDefectTotSpin(Config* cnf, double DefectTotSpin){
    cnf->DefectTotSpin = DefectTotSpin;
}

void Config_setCorrTotSpin(Config* cnf, double CorrTotSpin){
    cnf->CorrTotSpin = CorrTotSpin;
}

void Config_setHf_tensorfile(Config* cnf, char* hf_tensorfile){

    if (cnf->hf_tensorfile == NULL){
        fprintf(stderr, "Error: the memory for hf_tensorfile should be allocated\n");
        exit(EXIT_FAILURE);
    }
    strcpy(cnf->hf_tensorfile,hf_tensorfile);
}

void Config_setHf_cutoff(Config* cnf, double hf_cutoff){

    if (hf_cutoff < 0.0) {
        fprintf(stderr, "Error: current hf_cutoff (%f < 0) is not available\n",hf_cutoff);
        exit(EXIT_FAILURE);
    }
    cnf->hf_cutoff = hf_cutoff;
}

void Config_setHf_ignore_oor(Config* cnf, int hf_ignore_oor){

    if (hf_ignore_oor != 0 && hf_ignore_oor != 1) {
        fprintf(stderr, "Error: possible hf_ignore_oor (current : %d) is 0 or 1\n",hf_ignore_oor);
        exit(EXIT_FAILURE);
    }
    cnf->hf_ignore_oor = hf_ignore_oor;
}

void Config_setHf_readmode(Config* cnf, int hf_readmode){

    if (hf_readmode < 0 && hf_readmode > 3) {
        fprintf(stderr, "Error: current hf_readmode (%d < 0 or > 3) is not available\n",hf_readmode);
        exit(EXIT_FAILURE);
    }

    cnf->hf_readmode = hf_readmode;
}

void Config_setQd_tensorfile(Config* cnf, char* qd_tensorfile){

    if (cnf->qd_tensorfile == NULL){
        fprintf(stderr, "Error: the memory for qd_tensorfile should be allocated\n");
        exit(EXIT_FAILURE);
    }
    strcpy(cnf->qd_tensorfile,qd_tensorfile);
}

void Config_setQd_tensorfile_woqubit(Config* cnf, char* qd_tensorfile_woqubit){

    if (cnf->qd_tensorfile_woqubit == NULL){
        fprintf(stderr, "Error: the memory for qd_tensorfile_woqubit should be allocated\n");
        exit(EXIT_FAILURE);
    }
    strcpy(cnf->qd_tensorfile_woqubit,qd_tensorfile_woqubit);
}

void Config_setQd_readmode(Config* cnf, int qd_readmode){

    if (qd_readmode < 0 && qd_readmode > 2) {
        fprintf(stderr, "Error: possible qd_readmode is 0, 1, 2\n",qd_readmode);
        exit(EXIT_FAILURE);
    }
    cnf->qd_readmode = qd_readmode;
}



void Config_report(Config* cnf){

    printLineSection();
    printTitle("Structure Config");

    printStructElementChar("method",Config_getMethod(cnf));
    printStructElementChar("quantity",Config_getQuantity(cnf));
    printStructElementInt("order",Config_getOrder(cnf));
    printStructElementFloat1d("bfield",Config_getBfield(cnf),3);
    printStructElementFloat("rbath",Config_getRbath(cnf));
    printStructElementFloat("rdip",Config_getRdip(cnf));
    printStructElementFloat("deltat",Config_getDeltat(cnf));
    printStructElementInt("nstep",Config_getNstep(cnf));
    printStructElementFloat("rbathcut",Config_getRbathcut(cnf));
    printStructElementFloat("rdipcut",Config_getRdipcut(cnf));
    printStructElementInt("nstate",Config_getNstate(cnf));
    printStructElementInt("seed",Config_getSeed(cnf));

    printStructElementChar("qubitfile",Config_getQubitfile(cnf));
    printStructElementChar("gyrofile",Config_getGyrofile(cnf));
    printStructElementInt("nbathfiles",Config_getNbathfiles(cnf));
    for (int i = 0; i < Config_getNbathfiles(cnf); i++){
        printStructElementChar("bathfiles",Config_getBathfiles_i(cnf,i));
        printStructElementDouble1d("bathadjust",Config_getBathadjust_i(cnf,i),3);
    }

    printStructElementDouble("DefectTotSpin",Config_getDefectTotSpin(cnf));
    printStructElementDouble("CorrTotSpin",Config_getCorrTotSpin(cnf));

    printStructElementChar("hf_tensorfile",Config_getHf_tensorfile(cnf));
    printStructElementDouble("hf_cutoff",Config_getHf_cutoff(cnf));
    printStructElementInt("hf_ignore_oor",Config_getHf_ignore_oor(cnf));
    printStructElementInt("hf_readmode",Config_getHf_readmode(cnf));

    printStructElementChar("qd_tensorfile",Config_getQd_tensorfile(cnf));
    printStructElementChar("qd_tensorfile_woqubit",Config_getQd_tensorfile_woqubit(cnf));
    printStructElementInt("qd_readmode",Config_getQd_readmode(cnf));
    


    
    printLineSection();
}