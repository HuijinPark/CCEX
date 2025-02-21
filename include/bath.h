#ifndef __CCEX_BATH_H_
#define __CCEX_BATH_H_

#include "utilities.h"
#include "qubit.h"
#include "general.h"

/**
 * @struct BathSpin structure
 * @brief BathSpin structure contains each BathSpin's information 
 * @note When you make this structure, 
 *       the following action would be done automatically.
 *       (1) state value would be allocated
 *       (2) The dimension of alpha, beta would be checked if it is 2*spin+1
 *       (3) The state would be normalized
 */
typedef struct {

    /* Basic spin information (Essential) */

    /**
     * @brief the name of spin
     * @details Value can be obtained from &Bath -> bathfile (no default value)
     * @todo Add documentation for cce.in : &Bath -> bathfile 
    */ 
    char name[MAX_CHARARRAY_LENGTH]; 

    /**
     * @brief Spin quantum number
     * @details Value can be obtained from &Bath -> gyrofile (no default value)
    */
    float spin;    

    /**
     * @brief Gyromagnetic ratio (Unit : radkHz/Gauss)
     * @details Value can be obtained from &Bath -> gyrofile (no default value)
     * @todo Add documentation for cce.in : &Bath -> gyrofile 
    */
    double gyro;   

    // /**
    //  * @brief Quadrupole moment (Unit : 10^-30 [m^2])
    //  * @details Value can be obtained from &Bath -> qtensorfile (default : 0.0)
    //  *      - cf. eQ = 20.44[Q/millibarn = Q * 10^-31 m^2] , in DFT eQ = 2.044 * 10^-30 [m^2]
    //  * @note Unit in input file (cce.in) : (Unit : 10^-30 [m^2])
    //  * @todo Add documentation for cce.in : &Bath -> qtensorfile 
    //  * @todo When you read quadrupole tensor, 
    //  *       - you have to get them with unit Hatree/Bohr_radius^2 (QE) -> Think about How to do this
    //  *       - you have to get them with unit MHz (VASP) -> Think about How to do this
    //  *       - which one is better? -> Think about How to do this
    //  * @todo Now, there is no default value for eq. In future, add the information of each spin' eQ in code internally.
    // */      
    // double eq;
  
    /**
     * @brief Spin position (Unit : Angstrom)
     * @details Value can be obtained from &Bath -> bathfile (no default value)
     * @note The x,y,z position can be adjusted by &Bath -> bathadjust for each bath file
     * @todo Add documentation for cce.in : &Bath -> bathfile 
    */
    double xyz[3]; 


    /* Spin state (Optional)*/

    /**
     * @brief Spin sub-level in z-basis (if &General -> nstate != 0, i.e. single-sample approach)
     * @details Value can be obtained from the following : (no default value)
     *          - &Bath -> statefile (fixed)  // old version

     *          - &General -> nstate (randomly generated) 
     *          - &General -> seed (fixed) // new version to fix
     * @note The state would be normalized and The length of state = 2*spin+1
     * @todo Add documentation for cce.in : &Bath -> "statefile" -> "substatefile" tag
     * @todo Add documentation for cce.in : &General -> "nstate" and &General -> "seed" tag
     * @todo Replace the way to get fixed state from reading "file" to using "seed" 
    */
    float state; 


    /* Interaction Information (Optional) */

    /**
     * @brief The detuning value of a bath spin (Unit : radkHz)
     * @details The artifitial effect to shift energy level of a qubit.
     *          Value can be obtained from &Bath -> "bdetuning" (default : 0.0)
     * @note Unit in input file (cce.in) : MHz
     * @todo Add documentation for cce.in : &Bath -> "bdetuning" tag
    */
    double detuning; 

    /**
     * @brief The mean-field effect from zz-coupling with other bath spins (Unit : radkHz)
     * @details Automatically set if &General -> "nstate" tag is non-zero (default : 0.0)
     * @todo Add documentation for cce.in : &Bath -> "disorder" tag 
    */
    double disorder; 

    /**
     * @brief The hyperfine interaction between qubits and a bath spin (Unit : radkHz)
     * @details Value is obtained from &Bath -> "Atensorfile" (default : 0.0)
     *       - hypf[i] : the hyperfine interaction with i-th qubit
     * @note Unit in input file (cce.in) : MHz
     * @todo Add documentation for cce.in : &Bath -> "Atensorfile" tag
     * @todo IMPORTANT : The hyperfine interaction to be deleted if QubitArray include the one of the bath spin in A tensor file
    */
    MatrixXcd* hypf; // 2024.01.26 change naming due to &Defect.defect.hypf : hypf to intqb 보류

    /**
     * @brief The quadrupole interaction a bath spin (Unit : radkHz)
     * @details Value is obtained from &Bath -> "qtensorfile" (default : 0.0)
     * @note - Unit in input file (cce.in) &Defect -> "dquad" : MHz
     *       - Unit in input file (cce.in) &Bath -> "qtensorfile" : Hatree/Bohr_radius^2
     * @todo Above @note have to be unified with one unit (MHz seems better.. but the problem is how to get Qtensorfile..)
     * @todo Add documentation for cce.in : &Bath -> "qtensorfile" tag
     * @todo IMPORTANT : The quadrupole interaction to be deleted if QubitArray include the one of the bath spin in Q tensor file
    */
    MatrixXcd quad; // 2024.01.26 change naming due to &Defect.defect.quad : quad to intself 보류

    /* &Defect */
    MatrixXcd hypf_sub; // the hyperfine interaction with bath spin (if there is dft data)
    int mainspidx; // the connected main spin index with the coupling strength of hypf_sub

    double mindist; // minimum distance between qubit set and a bath spin
} BathSpin;

/*!
 * @struct BathSpinArray
 * @brief BathSpinArray include all Bath spins' information
 */
typedef struct {

    /**
     * @brief The number of bath spins
     * @todo Add documentation for cce.in : &Bath -> "nspin" tag
    */
    int nspin; /**< The number of BathSpin */  

    /**
     * @brief The array of BathSpin
     * @todo Add documentation for cce.in : &Bath tag
    */
    BathSpin** bath;  

    // spin properties
    int prop_nspecies; /**< The number of species */
    char** prop_names; /**< The name of species */
    double* prop_gyros; /**< The gyromagnetic ratio of species */
    float* prop_spins; /**< The spin quantum number of species */

} BathArray;


// spin pairs
void        BathArray_connectivity(int*** cmap, float*** stmap, BathArray* ba, float rdip, float rdipcut); // connectivity between bath spins
void        makeSparsemap(int*** spmap, int** cmap, int nspin); // sparsemap

// spin interactions
MatrixXcd   BathArray_int_i_j(BathArray* ba, int i, int j);
double      BathArray_dist_i_j(BathArray* ba, int i, int j);

void        BathArray_setBathHypfs(BathArray* ba, QubitArray* qa);

void        BathArray_setBathDisorders(BathArray* ba);
// ! update bath disorder if there is additional spins

double      BathArray_getBath_i_disorder_j(BathArray* ba, int isp, int jsp);

double      BathArray_getOverhaus(BathArray* ba, int iq); // overhaus field for i-th qubit
// ! update bath overhausers if there is additional spins

double      BathArray_getBath_i_overhaus_j(BathArray* ba, int isp, int iq); // overhaus field for i-th qubit

// random state generator
void        BathArray_setBathStatesRandom(BathArray* ba); // random state

// physics
int         BathArray_dim(BathArray* ba);
int         BathArray_dimBath_i(BathArray* ba, int i);
int         BathSpin_dim(BathSpin* bs);

// Hamiltonian
MatrixXcd   BathArray_ZeemanHamil(BathArray* ba, MatrixXcd** sigmas, int ib, float* bfield);
MatrixXcd   BathArray_DetuningHamil(BathArray* ba, MatrixXcd** sigmas, int ib);
MatrixXcd   BathArray_DisorderHamil(BathArray* ba, MatrixXcd** sigmas, int ib, bool rm_overlap);
MatrixXcd   BathArray_QuadHamil(BathArray* ba, MatrixXcd** sigmas, int ib);
MatrixXcd   BathArray_InteractionHamil(BathArray* ba, MatrixXcd** sigmas, int ib, int jb);
MatrixXcd** BathArray_PauliOperators(BathArray* ba);

// density matrix
MatrixXcd   BathArray_Rho0(BathArray* ba, bool isEnsemble);
MatrixXcd   BathArray_Psi0(BathArray* ba);

/* Low level functions --------------------------------------------*/

BathArray*  BathArray_init();

//report
void        BathArray_report(BathArray* ba);
void        BathArray_reportBath(BathArray* ba);
void        BathArray_reportBath_i_props(BathArray* ba, int i);
void        BathArray_reportBath_states(BathArray* ba);
void        BathArray_reportBath_detunings(BathArray* ba);
void        BathArray_reportBath_disorders(BathArray* ba);
void        BathArray_reportBath_hypf(BathArray* ba, int nqubit);
void        BathArray_reportBath_quad(BathArray* ba);
void        BathArray_reportBath_hypf_sub(BathArray* ba);
void        BathArray_reportSpinProperties(BathArray* ba);

//alloc
void        BathArray_allocBath(BathArray* ba, int nqubit);
void        BathArray_reallocBath(BathArray* ba, int nspin_old, int nspin_new, int nqubit);
void        BathArray_allocProp(BathArray* ba);
void        BathArray_reallocProp(BathArray* ba, int nspin_old, int nspin_new);
void        BathArray_allocBath_i_hypf(BathArray* ba, int i, int nqubit);

//sets
void        BathArray_setNspin(BathArray* ba, const int nspin);
void        BathArray_setBath_i(BathArray* ba, const BathSpin* bath, int i, int nqubit);
void        BathArray_setBath_i_name(BathArray* ba, const char* name, int i);
void        BathArray_setBath_i_spin(BathArray* ba, const float spin, int i);
void        BathArray_setBath_i_gyro(BathArray* ba, const double gyro, int i);
void        BathArray_setBath_i_xyz(BathArray* ba, const double* xyz, int i);
void        BathArray_setBath_i_state(BathArray* ba, const float state, int i);
void        BathArray_setBath_i_detuning(BathArray* ba, const double detuning, int i);
void        BathArray_setBath_i_disorder(BathArray* ba, const double disorder, int i);
void        BathArray_setBath_i_hypf_j(BathArray* ba, const MatrixXcd hypf, int i, int j); // j : qubit index
void        BathArray_setBath_i_quad(BathArray* ba, const MatrixXcd quad, int i);
void        BathArray_setBath_i_hypf_sub(BathArray* ba, const MatrixXcd hypf_sub, int i);
void        BathArray_setBath_i_mainspidx(BathArray* ba, const int mainspidx, int i);
void        BathArray_setProp_nspecies(BathArray* ba, const int nspecies);
void        BathArray_setProp_names_i(BathArray* ba, const char* name, const int i);
void        BathArray_setProp_gyros_i(BathArray* ba, const double gyro, const int i);
void        BathArray_setProp_spins_i(BathArray* ba, const float spin, const int i);
void        BathArray_setBath_i_mindist(BathArray* ba, const double mindist, int i);

// get
int         BathArray_getNspin(BathArray* ba);
int         BathArray_getProp_nspecies(BathArray* ba);
char**      BathArray_getProp_names(BathArray* ba);
double*     BathArray_getProp_gyros(BathArray* ba);
float*      BathArray_getProp_spins(BathArray* ba);
BathSpin*   BathArray_getBath_i(BathArray* ba, int i);
char*       BathArray_getBath_i_name(BathArray* ba, int i);
float       BathArray_getBath_i_spin(BathArray* ba, int i);
double      BathArray_getBath_i_gyro(BathArray* ba, int i);
double*     BathArray_getBath_i_xyz(BathArray* ba, int i);
float       BathArray_getBath_i_state(BathArray* ba, int i);
double      BathArray_getBath_i_detuning(BathArray* ba, int i);
double      BathArray_getBath_i_disorder(BathArray* ba, int i);
MatrixXcd   BathArray_getBath_i_hypf_j(BathArray* ba, int i, int j);
MatrixXcd   BathArray_getBath_i_quad(BathArray* ba, int i);
MatrixXcd   BathArray_getBath_i_hypf_sub(BathArray* ba, int i);
int         BathArray_getBath_i_mainspidx(BathArray* ba, int i);

// free
void        BathArray_freeAll(BathArray* ba);
void        BathArray_freeProp_names(BathArray* ba);
void        BathArray_freeProp_gyros(BathArray* ba);
void        BathArray_freeProp_spins(BathArray* ba);
void        BathArray_freeBath(BathArray* ba);
void        BathArray_freeBath_i_hypf(BathArray* ba, int i);

// BathSpin
void       BathSpin_setName(BathSpin* bs, char* name);
void       BathSpin_setName_withType(BathSpin* bs, char* name, char* type); // name : mainspin name, type : subspin
void       BathSpin_setSpin(BathSpin* bs, float spin);
void       BathSpin_setGyro(BathSpin* bs, double gyro);
void       BathSpin_setXyz(BathSpin* bs, double* xyz);
void       BathSpin_setXyz_fromRxyz(BathSpin* bs, double* xyz0, double* rxyz); // Input : xyz0 : main spin position, rxyz : sub spin relative position, Output : xyz : sub spin position
void       BathSpin_setState(BathSpin* bs, float state);
void       BathSpin_setDetuning(BathSpin* bs, double detuning);
void       BathSpin_setDisorder(BathSpin* bs, double disorder);
void       BathSpin_setHypf_i(BathSpin* bs, MatrixXcd hypf, int iq);
void       BathSpin_setQuad(BathSpin* bs, MatrixXcd quad);
void       BathSpin_setQuad_fromEFG(BathSpin* bs, MatrixXcd efg, double eq, float spin); // Input : efg : Hartree/Bohr^2, eq : 10e-30 m^2 , spin number, Output : quad : radkHz
void       BathSpin_setHypfSub(BathSpin* bs, MatrixXcd hypf_sub);
void       BathSpin_setMainspidx(BathSpin* bs, int mainspidx);

char*      BathSpin_getName(BathSpin* bs);
float      BathSpin_getSpin(BathSpin* bs);
double     BathSpin_getGyro(BathSpin* bs);
double*    BathSpin_getXyz(BathSpin* bs);
float      BathSpin_getState(BathSpin* bs);
double     BathSpin_getDetuning(BathSpin* bs);
double     BathSpin_getDisorder(BathSpin* bs);
MatrixXcd  BathSpin_getHypf_i(BathSpin* bs, int iq);
MatrixXcd  BathSpin_getQuad(BathSpin* bs);
MatrixXcd  BathSpin_getHypfSub(BathSpin* bs);
int        BathSpin_getMainspidx(BathSpin* bs);

#endif // __CCEX_BATH_H_
