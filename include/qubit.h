#ifndef __CCEX_QUBIT_H_
#define __CCEX_QUBIT_H_

#include "utilities.h"

/**
 * @struct Qubit structure
 * @brief Qubit structure contains each qubit's information 
 * @note When you make this structure, 
 *       the following action would be done automatically.
 *       (1) alpha, beta value would be allocated
 *       (2) The dimension of alpha, beta would be checked if it is 2*spin+1
 *       (3) The alpha, beta would be normalized
 */
typedef struct {
    
    /* Basic spin information (Essential) */

    /**
     * @brief The name of a qubit
     * @details Value can be obtained from &Qubit -> "qname" (default : e)
     * @todo Add documentation for cce.in : &Qubit -> "qname" tag
     */
    char name[MAX_CHARARRAY_LENGTH];

    /**
     * @brief The spin number of a qubit
     * @details Value can be obtained from &Qubit -> "qnspin" (default : 1)
     * @todo Add documentation for cce.in : &Qubit -> "qspin" tag
    */
    float spin;

    /**
     * @brief The gyromagnetic ratio of a qubit (Unit : radkHz/G)
     * @details Value is obtained from &Qubit -> "qgyro" (default : gyro_electron)
     * @note Unit in input file (cce.in) : radkHz/G
     * @todo Add documentation for cce.in : &Qubit -> "qgyro" tag
    */
    double gyro; 

    /**
     * @brief The position of a qubit (Unit : angstrom)
     * @details Value should be obtained from &Qubit -> "qxyz" or qubitfile (no default value)
     * @note Unit in input file (cce.in) : angstrom
     * @todo Add documentation for cce.in : &Qubit -> "qxyz" tag and qubitfile
    */
    double xyz[3];


    /* Interaction information (Optional) */ 

    /**
     * @brief The detuning value of a qubit (Unit : radkHz)
     * @details The artifitial effect to shift energy level of a qubit.
     *          Value can be obtained from &Qubit -> "qdetuning" (default : 0.0)
     * @note Unit in input file (cce.in) : MHz
     * @todo Add documentation for cce.in : &Qubit -> "qdetuning" tag
    */
    double detuning; 

    /**
     * @brief The energy detuning from spin bath (Unit : radkHz)
     * @details Automatically set if &Qubit -> "qoverhaus" tag is true (default : 0.0)
     * @todo Add documentation for cce.in : &Qubit -> "qoverhaus" tag 
    */
    double overhaus; 

    /* Qubit projected state (Optional) */

    /**
     * @brief The projected state of a qubit (alpha)
     * @details Value is be obtained from &Qubit -> "alphams" or -> "alphastate" tag (no default value)
     *       - Demension : 2*spin+1
     *       - alphams : The z-basis sub level of alpha (e.g. S, S-1, S-2 ... -S+1, -S)
     *       - alphastate : The spin state (matrix) of alpha (e.g. [1 0 0], [0 1 0], [0 0 1])
     * @todo Add documentation for cce.in : &Qubit -> "alphams" and -> "alphastate" tag
    */
    MatrixXcd alpha; 

    /**
     * @brief The projected state of a qubit (beta)
     * @details Value is be obtained from &Qubit -> "betams" or -> "betastate" tag (no default value)
     *      - Demension : 2*spin+1
     *      - betams : The z-basis sub level of beta (e.g. S, S-1, S-2 ... -S+1, -S)
     *      - betastate : The spin state (matrix) of beta (e.g. [1 0 0], [0 1 0], [0 0 1])
     * @todo Add documentation for cce.in : &Qubit -> "betams" and -> "betastate" tag
    */
    MatrixXcd beta;
    
} Qubit;

/**
 * @struct QubitArray
 * @brief QubitArray include information of all qubits
 * @todo Check if the psia, psib, psi0, rho0 is normalied.
 * @details Each parameter is read from the input file or options.
 *          * _parameter : Only use during reading input file and then free
 *          * parameter : Use during simulation
 * @todo alphastate, betastate : 
 *          currently the value format of this tag in input file is :
 *          e.g. alphastate = { e, 1, 0, 0, }
 *          but we need to change this format to bracket notation: 
 *          e,g. alphastate = e |-1>
 *               betastate = n1 |+1.5>  
 *               betastate = n1 i|-0.5>
 *          => generally : alphastate or betastate = qname a|ms>
 *                       : a is the constant (can be complex number)
 * @todo initstate :
 *          Like above case, we need to change the format of this tag.
 *          e.g. initstate = |+1>|1.5>|0.5> - i|+1>|0.5>|-0.5>
 *          => generally : initstate = 
 *                         a|ms_1>|ms_2>..|ms_n> + b|ms'_1>|ms'_2>..|ms'_n>
 * @todo Add doucumentation for qubitarray and cce.in "&Qubit" tag
 */
typedef struct {

    /* Use only during reading input file -----------------------------------------*/

    // QubitArray
    int* _alphaidx;      /**< Eigen state alpha index (default : NULL)             */
    int* _betaidx;       /**< Eigen state beta index (default : NULL)              */
    /*-----------------------------------------------------------------------------*/
    /* Use during simulation                                                       */
    /*-----------------------------------------------------------------------------*/
    int nqubit; /**< The number of Qubit                                           */
    Qubit** qubit; /**< Qubit information (see struct Qubit)                       */
    MatrixXcd** intmap; /**< Interaction i-j (Unit : radkHz)                    */
    MatrixXcd psia; /**< Projected state alpha for QubitArray                 */
    MatrixXcd psib; /**< Projected state beta for QubitArray                  */
    MatrixXcd psi0; /**< Initial state of qubit                               */
    bool overhaus;  /**< Overhauser field : on | off                               */
    /*-----------------------------------------------------------------------------*/
    /**<
     * About psia, psib
     * wfsize = Product_i(2*spin_i+1)
     * psia[wfsize]
     * psib[wfsize]
     * 
     * About psi0, rho0
     * psi0[wfsize]
     * rho0[wfsize][wfsize]
    */
} QubitArray;

/* High Level --------------------------------------------------------*/

// physical properties
int         QubitArray_dim(QubitArray* qa);
int         QubitArray_dimQubit_i(QubitArray* qa, int i);
double      QubitArray_mindist(double* xyz, QubitArray* qa);
void        QubitArray_setPsiaPsib_fromQubit(QubitArray* qa); // psia = kron(qubit_i_alpha,qubit_j_alpha, ...)
void        QubitArray_setPsiaPsib_fromIdx(QubitArray* qa, float* bfield); // psia, psib : eigenstate of Hq (idx)
void        QubitArray_setPsi0_fromPsiaPsib(QubitArray* qa); // psi0 = Normalize(psia + psib)

// Hamiltonian
MatrixXcd   QubitArray_TotalHamil(QubitArray* qa, MatrixXcd** sigmas, float* bfield);
MatrixXcd   QubitArray_SingleHamil(QubitArray* qa, MatrixXcd** sigmas, int iq, float* bfield);
MatrixXcd   QubitArray_ZeemanHamil(QubitArray* qa, MatrixXcd** sigmas, int iq, float* bfield);
MatrixXcd   QubitArray_DetuningHamil(QubitArray* qa, MatrixXcd** sigmas, int iq);
MatrixXcd   QubitArray_OverhausHamil(QubitArray* qa, MatrixXcd** sigmas, int iq);
MatrixXcd   QubitArray_ZFSHamil(QubitArray* qa, MatrixXcd** sigmas, int iq);
MatrixXcd   QubitArray_InteractionHamil(QubitArray* qa, MatrixXcd** sigmas, int iq, int jq);
MatrixXcd** QubitArray_PauliOperators(QubitArray* qa);
MatrixXcd*  QubitArray_PauliOperator_fromPsiaPsib(QubitArray* qa);

// density matrix
MatrixXcd   QubitArray_Rho0(QubitArray* qa);

// utils
int         QubitArray_getQubitIdx_fromName(QubitArray* qa, const char* name);

/* Low Level --------------------------------------------------------*/

QubitArray* QubitArray_init(); 

// report
void        QubitArray_report(QubitArray* qa);
void        QubitArray_reportQubit_i(QubitArray* qa, int i);
void        QubitArray_reportIntmap(QubitArray* qa);
void        QubitArray_reportPsiaPsib(QubitArray* qa);
void        QubitArray_reportPsi0(QubitArray* qa);

// alloc
void        QubitArray_allocQubit(QubitArray* qa);      // length = nqubit
void        QubitArray_allocIntmap(QubitArray* qa);     // length = (nqubit,nqubit)
void        QubitArray_alloc_alphaidx_betaidx(QubitArray* qa); // length = 1

// set
void        QubitArray_setNqubit(QubitArray* qa, const int nspin);
void        QubitArray_setOverhaus(QubitArray* qa, const bool overhaus);
void        QubitArray_setQubit(QubitArray* qa, Qubit** qubit);
void        QubitArray_set_alphaidx(QubitArray* qa, const int* alphaidx);
void        QubitArray_set_betaidx(QubitArray* qa, const int* betaidx);
void        QubitArray_setIntmap_i_j(QubitArray* qa, const MatrixXcd tensor, int i, int j);
void        QubitArray_setPsia(QubitArray* qa, const MatrixXcd psia);
void        QubitArray_setPsib(QubitArray* qa, const MatrixXcd psib);
void        QubitArray_setPsi0(QubitArray* qa, const MatrixXcd psi0);
void        QubitArray_setQubit_i_name(QubitArray* qa, const char* name, int i);
void        QubitArray_setQubit_i_spin(QubitArray* qa, const float spin, int i);
void        QubitArray_setQubit_i_gyro(QubitArray* qa, const double gyro, int i);
void        QubitArray_setQubit_i_xyz(QubitArray* qa, const double* xyz, int i);
void        QubitArray_setQubit_i_detuning(QubitArray* qa, const double detuning, int i);
void        QubitArray_setQubit_i_overhaus(QubitArray* qa, const double overhaus, int i);
void        QubitArray_setQubit_i_alpha(QubitArray* qa, const MatrixXcd alpha, int i);
void        QubitArray_setQubit_i_beta(QubitArray* qa, const MatrixXcd beta, int i);
void        QubitArray_setQubit_i_alpha_fromMs(QubitArray* qa, const float ms, int i);
void        QubitArray_setQubit_i_beta_fromMs(QubitArray* qa, const float ms, int i);

// get (This is not a copy. Watch out when you change the value with this function!)
int         QubitArray_getNqubit(const QubitArray* qa);
bool        QubitArray_getOverhaus(const QubitArray* qa);
int*        QubitArray_get_alphaidx(const QubitArray* qa);
int*        QubitArray_get_betaidx(const QubitArray* qa);
MatrixXcd** QubitArray_getIntmap(const QubitArray* qa);
MatrixXcd   QubitArray_getIntmap_i_j(const QubitArray* qa, int i, int j);
MatrixXcd   QubitArray_getPsia(const QubitArray* qa);
MatrixXcd   QubitArray_getPsib(const QubitArray* qa);
MatrixXcd   QubitArray_getPsi0(const QubitArray* qa);
char*       QubitArray_getQubit_i_name(const QubitArray* qa, int i);
float       QubitArray_getQubit_i_spin(const QubitArray* qa, int i);
double      QubitArray_getQubit_i_gyro(const QubitArray* qa, int i);
double*     QubitArray_getQubit_i_xyz(const QubitArray* qa, int i);
double      QubitArray_getQubit_i_detuning(const QubitArray* qa, int i);
double      QubitArray_getQubit_i_overhaus(const QubitArray* qa, int i);
MatrixXcd   QubitArray_getQubit_i_alpha(const QubitArray* qa, int i);
MatrixXcd   QubitArray_getQubit_i_beta(const QubitArray* qa, int i);

// free
void        QubitArray_freeAll(QubitArray* qa);
void        QubitArray_freeQubit(QubitArray* qa);
void        QubitArray_freeIntmap(QubitArray* qa);
void        QubitArray_free_alphaidx_betaidx(QubitArray* qa);

#endif // __CCEX_QUBIT_H_
