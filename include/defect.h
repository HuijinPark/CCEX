#ifndef __CCEX_DEFECT_H_
#define __CCEX_DEFECT_H_

#include "utilities.h"
#include "bath.h"

/**
 * @struct Defect
 * @brief Defect structure
*/
typedef struct{

    /* Defect information */
    char dfname[MAX_CHARARRAY_LENGTH];
    bool apprx;
    
    /* Defect' spin information */
    int  naddspin;
    char** types; // spname[n]
    float* spins; // spin[n]
    double* gyros; // gyro[n]
    double* eqs; // eq[n]  //quadrupole memomts eq in "milibarn * 10e+1" m^2 = 10e-30 m^2

    /* Defect' spin interaction information */
    // changes for avaax, spin

    int navaax; // possible number of principal axis

    double*** rxyzs; // rxyz[n+1][m][k] : n: avaax(1~navaax), m : addspin type , k : x,y,z
    MatrixXcd** hypf; // hypf[n+1][m] // (MHz)
    MatrixXcd** efg; // quad[n+1][m] // (Hartree/Bohr^2)
    MatrixXcd* zfs ; // zfs[n+1]  // (MHz)
    double* detuning; // detuning[n+1] //  (MHz)

    // n : available axis
    // m : addspin type

} Defect;

typedef struct{

    int nbathspin;
    int* naddspins; // naddspin[nbathspin]
    BathSpin*** subbath; // bspin[nbathspin][naddspin]

    int* paxes; // avaax[nbathspin] (Fix for the configuration)

    int ndefect;
    Defect** defect;

}DefectArray;



// Update BathArray (input : Defect_i, paxis(index), BathSpin)
// (1) mainspinidx = self-index
// (2) quad (zfs) = check paxis, replace the quad as zfs(axis)
// (3) detuning = check paxis, add +=detuning(axis)
void updateMainSpins_fromDefectArray(DefectArray* dfa, BathArray* ba);
void updateDisorder_main_sub(DefectArray* dfa, BathArray* ba);
void updateDisorder_sub_sub(DefectArray* dfa);
void updateOverhaus_qubit_sub(DefectArray* dfa, QubitArray* ba);

/* Low level functions ------------------------------------------------ */

// init
DefectArray* DefectArray_init();

// set defect information
void DefectArray_setNdefect(DefectArray* dfa, int ndefect);
void DefectArray_setDefect_idf_dfname(DefectArray* dfa, int idf, char* dfname);
void DefectArray_setDefect_idf_apprx(DefectArray* dfa, int idf, bool apprx);
void DefectArray_setDefect_idf_naddspin(DefectArray* dfa, int idf, int naddspin);
void DefectArray_setDefect_idf_types(DefectArray* dfa, int idf, char** types);
void DefectArray_setDefect_idf_spins(DefectArray* dfa, int idf, float* spins);
void DefectArray_setDefect_idf_gyros(DefectArray* dfa, int idf, double* gyros);
void DefectArray_setDefect_idf_eqs(DefectArray* dfa, int idf, double* eqs);
void DefectArray_setDefect_idf_navaax(DefectArray* dfa, int idf, int navaax);
void DefectArray_setDefect_idf_iax_isp_rxyz(DefectArray* dfa, int idf, int iax, int isp, double* rxyzs);
void DefectArray_setDefect_idf_iax_isp_hypf(DefectArray* dfa, int idf, int iax, int isp, MatrixXcd hypf);
void DefectArray_setDefect_idf_iax_isp_efg(DefectArray* dfa, int idf, int iax, int isp, MatrixXcd efg);
void DefectArray_setDefect_idf_iax_zfs(DefectArray* dfa, int idf, int iax, MatrixXcd zfs);
void DefectArray_setDefect_idf_iax_detuning(DefectArray* dfa, int idf, int iax, double detuning);

// find defect index from name
int DefectArray_findDefectIndex(DefectArray* dfa, char* dfname);

// get defect information
int DefectArray_getNdefect(DefectArray* dfa);
char* DefectArray_getDefect_idf_dfname(DefectArray* dfa, int idf);
bool DefectArray_getDefect_idf_apprx(DefectArray* dfa, int idf);
int DefectArray_getDefect_idf_naddspin(DefectArray* dfa, int idf);
char* DefectArray_getDefect_idf_isp_types(DefectArray* dfa, int idf, int isp);
float DefectArray_getDefect_idf_isp_spins(DefectArray* dfa, int idf, int isp);
double DefectArray_getDefect_idf_isp_gyros(DefectArray* dfa, int idf, int isp);
double DefectArray_getDefect_idf_isp_eqs(DefectArray* dfa, int idf, int isp);
int DefectArray_getDefect_idf_navaax(DefectArray* dfa, int idf);
double* DefectArray_getDefect_idf_iax_isp_rxyz(DefectArray* dfa, int idf, int iax, int isp);
MatrixXcd DefectArray_getDefect_idf_iax_isp_hypf(DefectArray* dfa, int idf, int iax, int isp);
MatrixXcd DefectArray_getDefect_idf_iax_isp_efg(DefectArray* dfa, int idf, int iax, int isp);
MatrixXcd DefectArray_getDefect_idf_iax_isp_zfs(DefectArray* dfa, int idf, int iax);
double DefectArray_getDefect_idf_iax_detuning(DefectArray* dfa, int idf, int iax);

// alloc defect information
void DefectArray_allocDefect(DefectArray* dfa); // defect : ndefect
void DefectArray_allocDefect_idf(DefectArray* dfa, int idf, int navaax, int naddspin);

// set subspin information
void DefectArray_setPaxesRandom(DefectArray* dfa, BathArray* ba); // nbathspin (// before calculate)
void DefectArray_setSubbathStatesRandom(DefectArray* dfa, BathArray* ba); // subbath->staets
void DefectArray_setNaddspins(DefectArray* dfa, BathArray* ba); // before calculate
void DefectArray_setSubbath(DefectArray* dfa, BathArray* ba, QubitArray* qa);

// alloc subspin information
void DefectArray_allocPaxes(DefectArray* dfa, int nbathspin);
void DefectArray_allocNaddspins(DefectArray* dfa, int nbathspin); // before calculate
void DefectArray_allocSubbath(DefectArray* dfa, BathArray* ba, int nqubit);

// get subspin information
int DefectArray_getNbathspin(DefectArray* dfa);
int DefectArray_getPaxes_i(DefectArray* dfa, int ibs);
int DefectArray_getNaddspins_i(DefectArray* dfa, int ibs);
BathSpin* DefectArray_getSubbath_i_isp(DefectArray* dfa, int ibs, int isp); // control this with BathSpin_get/set functions

// free 
void DefectArray_freePaxes(DefectArray* dfa);
void DefectArray_freeNaddspins(DefectArray* dfa);
void DefectArray_freeSubbath(DefectArray* dfa);
void DefectArray_freeAll(DefectArray* dfa);




// report
void DefectArray_reportAll(DefectArray* dfa);

void DefectArray_reportSubbath_props(DefectArray* dfa); // name, spin, gyro, xyz, mainspidx
void DefectArray_reportSubbath_states(DefectArray* dfa); // state 
void DefectArray_reportSubbath_hypfs(DefectArray* dfa, int nqubit); // hypf
void DefectArray_reportSubbath_quads(DefectArray* dfa); //quad
void DefectArray_reportSubbath_hypf_subs(DefectArray* dfa); // hypf_sub
void DefectArray_reportSubbath_disorders(DefectArray* dfa); // disorder
//, detuning (detuning is not implemented yet)

void DefectArray_reportDefect_idf(DefectArray* dfa, int idf);
void DefectArray_reportPaxes(DefectArray* dfa);
void DefectArray_reportNaddspins(DefectArray* dfa);

#endif //__CCEX_DEFECT_H_