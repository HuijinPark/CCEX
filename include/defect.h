#ifndef __CCEX_SIMULATOR_DEFECT_H_
#define __CCEX_SIMULATOR_DEFECT_H_

#ifndef __CCEX_SIMULATOR_SIMULATOR_H_
    #include "simulator.h"
#endif

/**
 * @struct Defect
 * @brief Defect structure
*/
typedef struct{

    /* Defect information */
    char dfname[MAX_CAHRARRAY_LENGTH];
    bool apprx;
    
    /* Defect' spin information */
    int  nspin;
    char* spname[MAX_CAHRARRAY_LENGTH]; // spname[n]
    float* spin; // spin[n]
    double* gyro; // gyro[n]
    double* eq; // eq[n]

    /* Defect' spin interaction information */
    // changes for avaax, spin

    int avaax; // possible number of principal axis

    double** rxyz[3]; // rxyz[n][m]
    DoubleTensor** hypf; // hypf[n][m]
    DoubleTensor** quad; // quad[n][m]
    DoubleTensor** zfs ; // zfs[n][m]
    double** detuning; // detuning[n][m]

    // n : sub spin index
    // m : available axis

} Defect;


// /**
//  * @brief DefectArray
// */
// typedef struct {
//     int      dfname[MAX_CHARARRAY_LENGTH];
//     char     apprx[MAX_CHARARRAY_LENGTH];
//     int      naddspin   = 0;
//     char**   types      = NULL;
//     float*   spins      = NULL;
//     float*   eqs        = NULL;
//     double*  gyros      = NULL;
//     int      avaax      = 0;
//     double** rxyzs      = NULL; // (avaax*naddspin,5)
//     double** hyperfine  = NULL; // (avaax*naddspin,11)
//     double** zerofield  = NULL; // (avaax,11)
//     double** detuning   = NULL; // (avaax,3)
//     double** quadrupole = NULL; // (avaax*naddspin,11)
// }DefectProp;

// typedef struct {
//     int nspecies;
//     DefectProp* defect;
// } DefectPropArray;


typedef struct{
        
    /* Input files */
    // char* _substatefile;
    // char* _avaaxfile;

    /* Defect information */
    int ndefect;
    Defect** defect;

    /* Defect' spin information */
    int nbathspin;
    float** states; // states[b][defect.nspin]


    /* Principal axis (Essential only for &Defect) */

    /**
     * @brief Principal axis of defect (Unit : radkHz) 
     * @details Value is obtained from &Bath -> "avaaxfile" and &Defect -> "davaax" (default : 0)
     *          - avaax = 0 : The bath doesn't have defects.
     *          - avaax > 1 : The bath has defects.
    */
    int* paxes; // avaax[b] 

    // b : bath spin index 
    // n : sub spin index 

} DefectArray;

// init
DefectArray* initDefectArray();

// set
void setDefectArray(DefectArray* dfa, Options* options);

// void setDefectArray_substatefile(DefectArray* dfa, Options* options);
// void setDefectArray_avaaxfile(DefectArray* dfa, Options* options);

void setDefectArrayNdefect(DefectArray* dfa, Options* options);
void setDefectArrayDefect(DefectArray* dfa, Options* options);
void setDefectArrayDefect_idf(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_dfname(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_apprx(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_nspin(DefectArray* dfa, Options* options, int idf); // n
void setDefectArrayDefect_idf_spname(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_spin(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_gyro(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_eq(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_avaax(DefectArray* dfa, Options* options, int idf); // m
void setDefectArrayDefect_idf_rxyz(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_hypf(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_quad(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_zfs(DefectArray* dfa, Options* options, int idf);
void setDefectArrayDefect_idf_detuning(DefectArray* dfa, Options* options, int idf);

void setDefectArrayStates_random(DefectArray* dfa, BathArray* ba);
void setDefectArrayPaxes_random(DefectArray* dfa, BathArray* ba);

// find index of defectname
int findDefectIndex(DefectArray* dfa, char* dfname);

// alloc
void allocDefectArray(DefectArray* dfa);
void allocDefectArrayDefect(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_spname(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_spin(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_gyro(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_eq(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_rxyz(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_hypf(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_quad(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_zfs(DefectArray* dfa, int idf);
void allocDefectArrayDefect_idf_detuning(DefectArray* dfa, int idf);
void allocDefectArrayStates(DefectArray* dfa, BathArray* ba);
void allocDefectArrayPaxes(DefectArray* dfa, BathArray* ba);

// free
void freeDefectArray(DefectArray* dfa);
void freeDefectArrayDefect(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_spname(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_spin(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_gyro(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_eq(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_rxyz(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_hypf(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_quad(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_zfs(DefectArray* dfa, int idf);
void freeDefectArrayDefect_idf_detuning(DefectArray* dfa, int idf);
void freeDefectArrayStates(DefectArray* dfa, BathArray* ba);
void freeDefectArrayPaxes(DefectArray* dfa, BathArray* ba);

// report

// get
int getDefectArrayNdefect(DefectArray* dfa);
Defect** getDefectArrayDefect(DefectArray* dfa);
Defect* getDefectArrayDefect_idf(DefectArray* dfa, int idf);
char* getDefectArrayDefect_idf_dfname(DefectArray* dfa, int idf);
bool getDefectArrayDefect_idf_apprx(DefectArray* dfa, int idf);
int getDefectArrayDefect_idf_nspin(DefectArray* dfa, int idf);
char** getDefectArrayDefect_idf_spname(DefectArray* dfa, int idf);
float* getDefectArrayDefect_idf_spin(DefectArray* dfa, int idf);
double* getDefectArrayDefect_idf_gyro(DefectArray* dfa, int idf);
double* getDefectArrayDefect_idf_eq(DefectArray* dfa, int idf);
int getDefectArrayDefect_idf_avaax(DefectArray* dfa, int idf);
double*** getDefectArrayDefect_idf_rxyz(DefectArray* dfa, int idf);
DoubleTensor** getDefectArrayDefect_idf_hypf(DefectArray* dfa, int idf);
DoubleTensor** getDefectArrayDefect_idf_quad(DefectArray* dfa, int idf);
DoubleTensor** getDefectArrayDefect_idf_zfs(DefectArray* dfa, int idf);
double** getDefectArrayDefect_idf_detuning(DefectArray* dfa, int idf);

char* getDefectArrayDefect_idf_isp_spname(DefectArray* dfa, int idf, int isp);
float getDefectArrayDefect_idf_isp_spin(DefectArray* dfa, int idf, int isp);
double getDefectArrayDefect_idf_isp_gyro(DefectArray* dfa, int idf, int isp);
double getDefectArrayDefect_idf_isp_eq(DefectArray* dfa, int idf, int isp);

double* getDefectArrayDefect_idf_isp_iax_rxyz(DefectArray* dfa, int idf, int isp, int iax);
DoubleTensor getDefectArrayDefect_idf_isp_iax_hypf(DefectArray* dfa, int idf, int isp, int iax);
DoubleTensor getDefectArrayDefect_idf_isp_iax_quad(DefectArray* dfa, int idf, int isp, int iax);
DoubleTensor getDefectArrayDefect_idf_isp_iax_zfs(DefectArray* dfa, int idf, int isp, int iax);
double getDefectArrayDefect_idf_isp_iax_detuning(DefectArray* dfa, int idf, int isp, int iax);

// BathSpin -> Expand as DefectSpin
int makeBathSpins_fromdefect(BathSpin*** bspin_new, BathSpin* bspin, DefectArray* dfa); // return nspin
//    bathspin w paxis  ->      subspin
//         name         -> name = name_[spname]
//         spin         -> spin = [spin]
//         gyro         -> gyro = [gyro]
//         eq           -> eq   = [eq]
//         xyz          -> xyz  = xyz + [rxyz of spin, paxis]
//        state         -> state = [substate]
//       detuning       -> detuning = detuning 
//       disorder       -> disorder = disorder // only think electron spin' mf (in &Defect)
//      hypf(intqb)     -> hypf = pd apprx with qubit
//     quad(intself)    -> quad = [quad of spin, paxis] or [zfs of spin, paxis]
//      defect_hypf     -> defect_hypf(0.0) => hypf (hypf of spin, paxis)
//

void freeBathSpins_fromdefect(BathSpin*** bspin_new, int nspin);




#endif //__CCEX_SIMULATOR_DEFECT_H_