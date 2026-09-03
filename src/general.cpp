#include "../include/general.h"

#include "../include/memory.h"
#include <string.h>
#include <float.h>  // DBL_EPSILON

// init
Config* Config_init(){
    Config* cnf = (Config*)allocArray1d(1,sizeof(Config));
    cnf->method[0] = '\0';
    cnf->quantity[0] = '\0';
    cnf->propagator[0] = '\0';
    cnf->evolution[0] = '\0';
    cnf->evolution_isdefault = true;
    cnf->hfmedi = false;
    cnf->knight = false;

    // Bath coordinate rotation : off. Config comes from calloc, so the two axes would
    // otherwise be [0,0,0] -- set them to +z so that "enabled" with no axis given is a
    // well-defined identity instead of a zero-norm error.
    cnf->rot_enabled = false;
    cnf->rot_bath_axis[0]  = 0.0; cnf->rot_bath_axis[1]  = 0.0; cnf->rot_bath_axis[2]  = 1.0;
    cnf->rot_qubit_axis[0] = 0.0; cnf->rot_qubit_axis[1] = 0.0; cnf->rot_qubit_axis[2] = 1.0;
    cnf->rot_reference_qubit[0] = '\0'; // empty : use Qubit[0]
    cnf->rot_qubit_position_frame[0] = '\0'; // empty : "bath"
    cnf->rot_hf_tensor_frame[0] = '\0'; // empty : not given
    cnf->rot_qd_tensor_frame[0] = '\0';

    cnf->defect_axis_reference = DEFECT_AXIS_REFERENCE_NONE;

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

void Config_reallocBathadjust(Config* cnf, int oldrow, int newrow){
    int col = 3;
    reallocDouble2d(&(cnf->bathadjust),oldrow,newrow,col);
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

void Config_allocAvaaxfile(Config* cnf){
    cnf->avaaxfile = allocChar1d(MAX_FILEPATH);
}

void Config_allocStatefile(Config* cnf){
    cnf->statefile = allocChar1d(MAX_FILEPATH);
}

void Config_allocExstatefile(Config* cnf){
    cnf->exstatefile = allocChar1d(MAX_FILEPATH);
}

void Config_alloc_flines(Config* cnf, int size){
    cnf->_flines = allocInt1d(size);

}
void Config_realloc_flines(Config* cnf, int oldsize, int newsize){
    reallocInt1d(&(cnf->_flines),newsize);
}

// free
void Config_freeAll(Config* cnf){
    Config_freeBathfiles(cnf);
    Config_freeBathadjust(cnf);
    Config_freeGyrofile(cnf);
    Config_freeQubitfile(cnf);
    Config_freeAvaaxfile(cnf);
    Config_freeStatefile(cnf);
    Config_freeExstatefile(cnf);
    Config_freeHf_tensorfile(cnf);
    Config_freeQd_tensorfile(cnf);
    Config_freeQd_tensorfile_woqubit(cnf);
    Config_free_flines(cnf);
    freeArray1d((void**)&(cnf));
}

void Config_freeBathfiles(Config* cnf){
    freeChar2d(&(cnf->bathfiles),cnf->nbathfiles);
}

void Config_freeBathadjust(Config* cnf){
    freeDouble2d(&(cnf->bathadjust),cnf->nbathfiles);
}

void Config_freeGyrofile(Config* cnf){
    freeChar1d(&(cnf->gyrofile));
}

void Config_freeQubitfile(Config* cnf){
    freeChar1d(&(cnf->qubitfile));
}

void Config_freeAvaaxfile(Config* cnf){
    freeChar1d(&(cnf->avaaxfile));
}

void Config_freeStatefile(Config* cnf){
    freeChar1d(&(cnf->statefile));
}

void Config_freeExstatefile(Config* cnf){
    freeChar1d(&(cnf->exstatefile));
}

void Config_free_flines(Config* cnf){
    freeInt1d(&(cnf->_flines));
}

void Config_freeHf_tensorfile(Config* cnf){
    freeChar1d(&(cnf->hf_tensorfile));
}

void Config_freeQd_tensorfile(Config* cnf){
    freeChar1d(&(cnf->qd_tensorfile));
}

void Config_freeQd_tensorfile_woqubit(Config* cnf){
    freeChar1d(&(cnf->qd_tensorfile_woqubit));
}

/* Low level ---------------------------------------------------------*/

// get
char* Config_getMethod(Config* cnf){
    return cnf->method;
}

bool Config_getHfmedi(Config* cnf){
    return cnf->hfmedi;
}

bool Config_getKnight(Config* cnf){
    return cnf->knight;
}


char* Config_getQuantity(Config* cnf){
    return cnf->quantity;
}
char* Config_getPropagator(Config* cnf){
    return cnf->propagator;
}
char* Config_getEvolution(Config* cnf){
    return cnf->evolution;
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

char*   Config_getAvaaxfile(Config* cnf){
    return cnf->avaaxfile;
}

char*   Config_getStatefile(Config* cnf){
    return cnf->statefile;
}

char*   Config_getExstatefile(Config* cnf){
    return cnf->exstatefile;
}

int Config_get_nflines(Config* cnf){
    return cnf->_nflines;
}

int* Config_get_flines(Config* cnf){
    return cnf->_flines;
}

int  Config_get_flines_i(Config* cnf, int i){
    return cnf->_flines[i];
}

bool Config_getRot_enabled(Config* cnf){
    return cnf->rot_enabled;
}

double* Config_getRot_bath_axis(Config* cnf){
    return cnf->rot_bath_axis;
}

double* Config_getRot_qubit_axis(Config* cnf){
    return cnf->rot_qubit_axis;
}

char* Config_getRot_reference_qubit(Config* cnf){
    return cnf->rot_reference_qubit;
}

char* Config_getRot_qubit_position_frame(Config* cnf){
    return cnf->rot_qubit_position_frame;
}

char* Config_getRot_hf_tensor_frame(Config* cnf){
    return cnf->rot_hf_tensor_frame;
}

char* Config_getRot_qd_tensor_frame(Config* cnf){
    return cnf->rot_qd_tensor_frame;
}

DefectAxisReference Config_getDefect_axis_reference(Config* cnf){
    return cnf->defect_axis_reference;
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

double*   Config_getQd_cellpara(Config* cnf){
    return cnf->qd_cellpara;
}

int     Config_getQd_readmode(Config* cnf){
    return cnf->qd_readmode;
}

// set 
void Config_setMethod(Config* cnf, char* method){

    const int opsize = 7;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"gCCE","CCE","dsj","itb","dsjitb","pCCE"};

    int idx = findIndexCharFix(options,0,opsize-1,method);
    if (idx == -1) {
        fprintf(stderr, "Error: current method options (%s) is not available\n",method);
        fprintf(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            fprintf(stderr, "%s ",options[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "Please check the input file or change the possible option set\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->method,options[idx]);
}

void Config_setQuantity(Config* cnf, char* quantity){

    const int opsize = 3;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"coherence","dm","noise"};
    int idx = findIndexCharFix(options,0,opsize-1,quantity);
    if (idx == -1) {
        fprintf(stderr, "Error: current quantity options (%s) is not available\n",quantity);
        fprintf(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            fprintf(stderr, "%s ",options[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "Please check the input file or change the possible option set\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->quantity,quantity);
}

void Config_setPropagator(Config* cnf, char* propagator){

    // gCCE free-evolution propagator.
    //   eigen : one SelfAdjointEigenSolver per cluster, exp(-i*H*tau) = M*diag(exp(-i*l*tau))*M^dag
    //   expm  : Eigen's general matrix exponential per step, the original CCEX path
    const int opsize = 2;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"eigen","expm"};
    int idx = findIndexCharFix(options,0,opsize-1,propagator);
    if (idx == -1) {
        fprintf(stderr, "Error: current propagator options (%s) is not available\n",propagator);
        fprintf(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            fprintf(stderr, "%s ",options[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "Please check the input file or change the possible option set\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->propagator,propagator);
}

void Config_setEvolution(Config* cnf, char* evolution){

    // gCCE state representation.
    //   vector : DEFAULT. rho0 is a pure state whenever nstate > 0 (QubitArray_Rho0 is
    //            always psi0*psi0^dag, and BathArray_Rho0 is too unless it is the
    //            ensemble average), so propagate |psi> instead of rho. Each step becomes
    //            matrix-VECTOR products, O(n^2) where the matrix path is O(n^3), and the
    //            partial trace collapses to reshaping |psi> into a qdim x bdim matrix Psi
    //            and forming Psi*Psi^dag.
    //   matrix : propagate the full density matrix, rho(t) = U*rho0*U^dag, then partial-
    //            trace it. What CCEX did before; works for every configuration.
    // Same physics, different arithmetic: the results agree to rounding but are not
    // guaranteed bit-identical, so a converged series should not switch mid-way.
    // vector needs nstate > 0 (an ensemble rho0 is mixed, so no state vector exists),
    // propagator=eigen (with expm the per-step matrix exponential dominates and there is
    // nothing to win) and method=gCCE. cJSON_readOptionConfig falls the DEFAULT back to
    // matrix wherever those do not hold; an EXPLICIT evolution=vector is never downgraded
    // -- calCoherenceGcce says why it cannot run and stops.
    const int opsize = 2;
    char options[opsize][MAX_CHARARRAY_LENGTH] = {"matrix","vector"};
    int idx = findIndexCharFix(options,0,opsize-1,evolution);
    if (idx == -1) {
        fprintf(stderr, "Error: current evolution options (%s) is not available\n",evolution);
        fprintf(stderr, "Available options are : ");
        for (int i = 0; i < opsize; i++){
            fprintf(stderr, "%s ",options[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "Please check the input file or change the possible option set\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->evolution,evolution);
}

/**
 * Pick the evolution path once every option source has been applied.
 *
 * This CANNOT be decided while the input file is being parsed. main.cpp reads the file
 * from inside the getopt loop (case 'f'), so -m, -N and the rest are still to come:
 * `-f ccein.json -N 1` on a file with no "nstate" key would have been resolved against
 * nstate = 0 and frozen at matrix, even though the run does have a pure rho0. Call this
 * after the loop, before Config_report, so the reported value is the one that runs.
 *
 * An explicit "evolution" is never rewritten -- only checked against the one condition
 * that is fatal outside gCCE. calCoherenceGcce enforces the other two.
 */
void Config_resolveEvolution(Config* cnf){

    bool isGCCE = (strcasecmp(cnf->method,"gcce")==0);

    if (!cnf->evolution_isdefault){
        // Only calCoherenceGcce reads evolution. Any other method would ignore it and run
        // the density-matrix path anyway, so refuse rather than let the input claim
        // something the run does not do.
        if (strcasecmp(cnf->evolution,"matrix")!=0 && !isGCCE){
            fprintf(stderr,"Error : evolution=%s is implemented for method=gCCE only "
                           "(this run has method=%s).\n",cnf->evolution,cnf->method);
            exit(EXIT_FAILURE);
        }
        return;
    }

    bool vectorFits = isGCCE
                   && (cnf->nstate > 0)
                   && (strcasecmp(cnf->propagator,"eigen")==0);
    strcpy(cnf->evolution, vectorFits ? "vector" : "matrix");
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
        fprintf(stderr, "Error: current seed(%d) ( < 0) is not available\n",seed);
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

void Config_setAvaaxfile(Config* cnf, char* avaaxfile){

    if (cnf->avaaxfile == NULL){
        fprintf(stderr, "Error: the memory for avaaxfile should be allocated\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->avaaxfile,avaaxfile);
}

void Config_setStatefile(Config* cnf, char* statefile){

    if (cnf->statefile == NULL){
        fprintf(stderr, "Error: the memory for statefile should be allocated\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->statefile,statefile);
}

void Config_setExstatefile(Config* cnf, char* exstatefile){

    if (cnf->exstatefile == NULL){
        fprintf(stderr, "Error: the memory for exstatefile should be allocated\n");
        exit(EXIT_FAILURE);
    }

    strcpy(cnf->exstatefile,exstatefile);
}

void Config_set_nflines(Config* cnf, int nflines){
    cnf->_nflines = nflines;
}

void Config_set_flines_i(Config* cnf, int fline, int i){
    cnf->_flines[i] = fline;
}

// The axes are only copied here. Whether they describe a usable frame -- finite,
// non-zero, not collinear -- is decided in one place, buildQubitAlignedRotation,
// which is also what the rotation itself is built from.
void Config_setRot_enabled(Config* cnf, bool rot_enabled){
    cnf->rot_enabled = rot_enabled;
}

void Config_setRot_bath_axis(Config* cnf, double* rot_bath_axis){
    copyDouble1d(cnf->rot_bath_axis,rot_bath_axis,3);
}

void Config_setRot_qubit_axis(Config* cnf, double* rot_qubit_axis){
    copyDouble1d(cnf->rot_qubit_axis,rot_qubit_axis,3);
}

void Config_setRot_reference_qubit(Config* cnf, char* rot_reference_qubit){

    if (strlen(rot_reference_qubit) >= MAX_CHARARRAY_LENGTH){
        fprintf(stderr,"Error: reference_qubit \"%s\" is longer than a qubit name can be (%d characters)\n",
                rot_reference_qubit,MAX_CHARARRAY_LENGTH-1);
        exit(EXIT_FAILURE);
    }
    strcpy(cnf->rot_reference_qubit,rot_reference_qubit);
}

// "bath" or "qubit" and nothing else: an unrecognized string would otherwise be read
// as "not bath", silently skipping a transformation the tensors needed.
static void Config_setTensorFrame(char* dest, char* frame, const char* key){

    if (strcasecmp(frame,"bath")!=0 && strcasecmp(frame,"qubit")!=0){
        fprintf(stderr,"Error: %s must be \"bath\" or \"qubit\" (got \"%s\")\n",key,frame);
        fprintf(stderr,"  \"bath\"  : tensor components are in the original bath Cartesian basis\n");
        fprintf(stderr,"  \"qubit\" : tensor components are already in the qubit-aligned basis\n");
        exit(EXIT_FAILURE);
    }
    strcpy(dest,frame);
}

void Config_setRot_hf_tensor_frame(Config* cnf, char* frame){
    Config_setTensorFrame(cnf->rot_hf_tensor_frame,frame,"coordinate_frame_rotation.hf_tensor_frame");
}

void Config_setRot_qd_tensor_frame(Config* cnf, char* frame){
    Config_setTensorFrame(cnf->rot_qd_tensor_frame,frame,"coordinate_frame_rotation.qd_tensor_frame");
}

// Only "bath" for now. "qubit" is refused rather than deferred: it would put the qubit
// positions in the computational frame while the bath is still in the source frame, and
// readBathfiles would then compare the two in its rbath / rbathcut selection, its mindist
// and its point-dipole tensors.
void Config_setRot_qubit_position_frame(Config* cnf, char* frame){

    if (strcasecmp(frame,"bath") != 0){
        fprintf(stderr,"Error: coordinate_frame_rotation.qubit_position_frame must be \"bath\" (got \"%s\")\n",frame);
        if (strcasecmp(frame,"qubit") == 0){
            fprintf(stderr,"\"qubit\" is not supported in this version: the bath files are read in the\n");
            fprintf(stderr,"source frame, so qubit positions already in the computational frame would be\n");
            fprintf(stderr,"compared against them before anything is rotated.\n");
        }
        fprintf(stderr,"Write every Qubit.xyz in the same Cartesian basis and origin as the bath files.\n");
        exit(EXIT_FAILURE);
    }
    strcpy(cnf->rot_qubit_position_frame,frame);
}

void Config_setDefect_axis_reference(Config* cnf, DefectAxisReference reference){
    cnf->defect_axis_reference = reference;
}

/**
 * @brief With the rotation on, require the field to lie along the computational z axis
 * @details Mirrors Config_resolveEvolution in timing: it cannot run during the file
 *          parse, because -B is applied afterwards and would change the answer.
 *
 *          Without the rotation this is a no-op, so an arbitrary field vector keeps
 *          working exactly as it always has. With the rotation on, the computational z
 *          axis IS the qubit axis, and "bfield" is read in that computational frame --
 *          so [0,0,500] means 500 G along +qubit_axis and [0,0,-500] means 500 G the
 *          other way. Both are accepted: the requirement is the LINE, not the sign.
 *          A tilted field is not supported alongside the rotation and is refused here.
 *
 *          The tolerance is tight on purpose. Nothing rotates these numbers -- they are
 *          the input as written -- so anything beyond round-off in the parse is a real
 *          transverse component that the user meant and that this version cannot honour.
*/
void Config_validateBfieldAlignment(Config* cnf){

    if (!cnf->rot_enabled){ return; }

    double bx = cnf->bfield[0], by = cnf->bfield[1], bz = cnf->bfield[2];
    double bperp = sqrt(bx*bx + by*by);
    double bnorm = sqrt(bx*bx + by*by + bz*bz);

    if (bnorm == 0.0){ return; } // no field at all : nothing to align

    double tol = 1.0e-12 * ((bnorm > 1.0) ? bnorm : 1.0);

    if (bperp > tol){
        fprintf(stderr,"Error: coordinate_frame_rotation requires bfield to be aligned\n");
        fprintf(stderr,"with the computational z axis, which represents qubit_axis.\n");
        fprintf(stderr,"Current bfield: [%g, %g, %g]\n",bx,by,bz);
        fprintf(stderr,"  transverse component |B_perp| = %g  (tolerance %g)\n",bperp,tol);
        fprintf(stderr,"Write the field as [0, 0, Bz]; Bz may be negative for the opposite\n");
        fprintf(stderr,"direction along qubit_axis. A tilted field is not supported with the\n");
        fprintf(stderr,"coordinate frame rotation.\n");
        exit(EXIT_FAILURE);
    }
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

    // Anything outside 0..3 used to fall through -- the guard read
    // "hf_readmode < 0 && hf_readmode > 3", which is never true -- and readHftensorfile
    // has no else after its `==0` / `==1||2||3` branches, so neither
    // BathArray_setBathHypfs nor the tensor file ever ran. Every bath spin kept an unset
    // hyperfine tensor and the run finished, quietly, with no hyperfine coupling at all.
    if (hf_readmode < 0 || hf_readmode > 3) {
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

void Config_setQd_cellpara(Config* cnf, double* qd_cellpara){
    copyDouble1d(cnf->qd_cellpara,qd_cellpara,3);
}

void Config_setQd_readmode(Config* cnf, int qd_readmode){

    // The supported set is not a range: 1 was removed, so 0, 2, 3 and 4 are listed one
    // by one rather than bounded. Anything else used to fall through -- the guard read
    // "qd_readmode < 0 && qd_readmode > 4", which is never true -- and then missed every
    // case of the switch in cJSON_readOptionConfig and every branch of readQdtensorfile,
    // so the run went on quietly with no quadrupole interaction at all.
    if (qd_readmode != 0 &&
        qd_readmode != 2 &&
        qd_readmode != 3 &&
        qd_readmode != 4) {
        fprintf(stderr,
                "Error: current qd_readmode (%d) is not available; "
                "supported values are 0, 2, 3 and 4\n",
                qd_readmode);
        exit(EXIT_FAILURE);
    }
    cnf->qd_readmode = qd_readmode;
}

void Config_setHfmedi(Config* cnf, bool hfmedi){

    cnf->hfmedi = hfmedi;
}

void Config_setKnight(Config* cnf, bool knight){

    cnf->knight = knight;
}


void Config_report(Config* cnf){

    printTitle("Structure Config");

    printStructElementChar("method",Config_getMethod(cnf));
    printStructElementChar("quantity",Config_getQuantity(cnf));
    printStructElementChar("propagator",Config_getPropagator(cnf));
    printStructElementChar("evolution",Config_getEvolution(cnf));
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

    printStructElementBool("frame_rotation",Config_getRot_enabled(cnf));
    if (Config_getRot_enabled(cnf)){
        printStructElementDouble1d("  bath_axis",Config_getRot_bath_axis(cnf),3);
        printStructElementDouble1d("  qubit_axis",Config_getRot_qubit_axis(cnf),3);
        printStructElementChar("  reference_qubit",
            Config_getRot_reference_qubit(cnf)[0] ? Config_getRot_reference_qubit(cnf) : (char*)"(first qubit)");
        printStructElementChar("  qubit_position_frame",
            Config_getRot_qubit_position_frame(cnf)[0] ? Config_getRot_qubit_position_frame(cnf) : (char*)"bath (default)");
    }

    printStructElementChar("defect_axis_reference",
        (Config_getDefect_axis_reference(cnf) == DEFECT_AXIS_REFERENCE_QUBIT_AXIS)
        ? (char*)"qubit_axis" : (char*)"(none)");

    printStructElementChar("avaaxfile",Config_getAvaaxfile(cnf));
    printStructElementChar("statefile",Config_getStatefile(cnf));
    printStructElementChar("exstatefile",Config_getExstatefile(cnf));

    printStructElementDouble("DefectTotSpin",Config_getDefectTotSpin(cnf));
    printStructElementDouble("CorrTotSpin",Config_getCorrTotSpin(cnf));

    printStructElementInt("hf_readmode",Config_getHf_readmode(cnf));
    if (Config_getHf_readmode(cnf)!=0){
        printStructElementChar("hf_tensorfile",Config_getHf_tensorfile(cnf));
        printStructElementDouble("hf_cutoff",Config_getHf_cutoff(cnf));
        printStructElementInt("hf_ignore_oor",Config_getHf_ignore_oor(cnf));
    }

    printStructElementInt("qd_readmode",Config_getQd_readmode(cnf));
    if (Config_getQd_readmode(cnf)==2){
        printStructElementChar("qd_tensorfile",Config_getQd_tensorfile(cnf));
    }else if (Config_getQd_readmode(cnf)==3 || Config_getQd_readmode(cnf)==4 ){
        printStructElementChar("qd_tensorfile",Config_getQd_tensorfile(cnf));
        printStructElementChar("qd_tensorfile_woqubit",Config_getQd_tensorfile_woqubit(cnf));
        printStructElementDouble1d("qd_cellpara",Config_getQd_cellpara(cnf),3);
    }

    printStructElementBool("hfmedi",Config_getHfmedi(cnf));
    printStructElementBool("knight",Config_getKnight(cnf));
}
