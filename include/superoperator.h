#ifndef __CCEX_SUPEROPERATOR_H_
#define __CCEX_SUPEROPERATOR_H_

#include "utilities.h"
#include "bath.h"

/**
 * @brief Jump operator definition for ME-CCE
 * @details Defines a collapse operator with its type and rate
 */
typedef struct {
    char bath_name[MAX_CHARARRAY_LENGTH]; /**< Bath spin name to apply this jump operator */
    char op_type[MAX_CHARARRAY_LENGTH];   /**< Operator type: "+-" (lowering), "-+" (raising), "z" (dephasing), "x", "y" */
    double rate;                          /**< Rate of the jump process (Unit: radkHz) */
} JumpOperator;

/**
 * @brief Collection of jump operators for ME-CCE
 */
typedef struct {
    int njump;           /**< Number of jump operators */
    JumpOperator** ops;  /**< Array of jump operator pointers */
} JumpOperatorArray;

// JumpOperatorArray management
JumpOperatorArray* JumpOperatorArray_init();
void JumpOperatorArray_alloc(JumpOperatorArray* joa, int njump);
void JumpOperatorArray_freeAll(JumpOperatorArray* joa);
void JumpOperatorArray_setOp(JumpOperatorArray* joa, int idx, const char* bath_name, const char* op_type, double rate);
int  JumpOperatorArray_getNjump(JumpOperatorArray* joa);
void JumpOperatorArray_report(JumpOperatorArray* joa);

// Superoperator utility functions

/**
 * @brief Convert operator pair to superoperator: S[rho] = A * rho * B  -->  supop = kron(A, B^T)
 * @note Uses row-major vectorization convention (matching PyCCE)
 */
MatrixXcd op_to_supop(const MatrixXcd& left, const MatrixXcd& right);

/**
 * @brief Vectorize density matrix: N x N -> N^2 x 1 (row-major order)
 */
VectorXcd mat_to_vec(const MatrixXcd& mat);

/**
 * @brief Un-vectorize: N^2 x 1 -> N x N (row-major order)
 */
MatrixXcd vec_to_mat(const VectorXcd& vec);

/**
 * @brief Coherent superoperator: L_coh = -i(H ⊗ I - I ⊗ H)
 */
MatrixXcd coherent_superoperator(const MatrixXcd& H);

/**
 * @brief Projected coherent superoperator: L_proj = -i(H_alpha ⊗ I - I ⊗ H_beta)
 * @details Used in conventional CCE where Hamiltonians are projected on qubit alpha/beta states
 */
MatrixXcd projected_coherent_superoperator(const MatrixXcd& Ha, const MatrixXcd& Hb);

/**
 * @brief Build collapse superoperator for a single jump operator C
 * @details D = C⊗C* - 1/2(I⊗(C†C)^T) - 1/2(C†C⊗I)
 */
MatrixXcd collapse_superoperator_single(const MatrixXcd& C);

/**
 * @brief Get the operator matrix from operator type string for given spin dimension
 * @param op_type Operator type string: "+-", "-+", "z", "x", "y"
 * @param dim Hilbert space dimension of the spin (2S+1)
 * @param rate Square root of the rate is multiplied to the operator
 */
MatrixXcd get_jump_operator_matrix(const char* op_type, int dim, double rate);

/**
 * @brief Build incoherent superoperator for all applicable bath spins in the cluster
 * @param ba BathArray of the cluster
 * @param bsigmas Pauli operators for each bath spin
 * @param joa JumpOperatorArray containing all jump operator definitions
 * @param bdim Total bath Hilbert space dimension
 */
MatrixXcd incoherent_superoperator(BathArray* ba, MatrixXcd** bsigmas, JumpOperatorArray* joa, int bdim);

/**
 * @brief Compute super-propagator: U(t) = exp(L * t * 2*pi)
 */
MatrixXcd simple_incoherent_propagator(double t, const MatrixXcd& L);

/**
 * @brief Propagate CPMG super-propagators
 * @param u_before_pi Superoperator before pi-pulse
 * @param u_after_pi Superoperator after pi-pulse (swapped alpha/beta)
 * @param number Number of pi-pulses
 */
MatrixXcd propagate_superpropagators(const MatrixXcd& u_before_pi, const MatrixXcd& u_after_pi, int number);

/**
 * @brief Expand a single-spin operator to the full Hilbert space
 * @param single_op Operator acting on spin at index ispin
 * @param bsigmas Pauli operators for all spins (used for dimensions)
 * @param nspin Total number of spins
 * @param ispin Index of the spin this operator acts on
 */
MatrixXcd expand_operator(const MatrixXcd& single_op, MatrixXcd** bsigmas, int nspin, int ispin);

#endif // __CCEX_SUPEROPERATOR_H_
