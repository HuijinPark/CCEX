#ifndef __CCEX_UTILITIES_UTILITIES_H_
#define __CCEX_UTILITIES_UTILITIES_H_

#include <Eigen/Dense>
#include <iostream>
#include <mpi.h>
#include <unistd.h>

#define EIGEN_USE_MKL_ALL
#define _USE_MATH_DEFINES

/* MPI Global variable ------------------------------------------*/
extern int rank; /**< Rank of the process (MPI)*/
extern int nprocess; /**< Number of processes (MPI)*/
/**--------------------------------------------------------------*/

/* Print Global variable ----------------------------------------*/
extern bool verbosity; /**< Verbosity level */
/**--------------------------------------------------------------*/

/* Data Types ---------------------------------------------------*/
/**
 * @struct DoubleComplex
 * @brief Complex number with double precision
 * @note In local PC, some complex operators from complex.h are not supported.
 * 
*/
typedef struct {
    double real; /**< Real number */
    double imag; /**< Imaginary number */
} DoubleComplex;

/**
 * @struct DoubleTensor
 * @brief Tensor with double precision 
*/
typedef struct {
    double xx, xy, xz; /**< xx, xy, xz */
    double yx, yy, yz; /**< yx, yy, yz */
    double zx, zy, zz; /**< zx, zy, zz */
} DoubleTensor;

/* Eigen Library ------------------------------------------------*/
/** 
 * @typedef doublec
 * @brief Complex number with double precision in Eigen library
*/
typedef std::complex<double>doublec;

/**
 * @typedef MatrixXcd
 * @brief Matrix with complex number with double precision in Eigen library
*/
typedef Eigen::Matrix<doublec, Eigen::Dynamic, Eigen::Dynamic> MatrixXcd;

/**
 * @namespace Eigen
 * @brief Eigen namespace
*/
using namespace Eigen;


/* Math Constants -----------------------------------------------*/
// Planck constant : 1.054 571 729 x 10^-34 [J*s]
// H_BAR is the unit conversion constant (Not the Planck constant)
#define H_BAR 1.054571729

// Unit conversion
#define MHZ_TO_RADKHZ(x) (1 * (2.0 * M_PI * 1.0e+3) * x)
#define RADKHZ_TO_MHZ(x) (1 / (2.0 * M_PI * 1.0e+3) * x)
#define KHZ_TO_RADKHZ(x) (1 * (2.0 * M_PI) * x)
#define RADKHZ_TO_KHZ(x) (1 / (2.0 * M_PI) * x)

// Gyromagnetic ratio of electron spin [radkHz/G]
#define GAMMA_ELECTRON -17608.597050;

/* File Constants -----------------------------------------------*/
#define MAX_FILEPATH 500 /**< Maximum length of file path */

/* cce.in line constants ----------------------------------------*/
#define MAX_OPTION_ARRAY_LENGTH 500 /**< Maximum length of option array */
#define MAX_OPTION_ELEMENT_LENGTH 50 /**< Maximum length of option element */
#define MAX_CHARARRAY_LENGTH 20 /**< Maximum length of structure char array */

#define MAX_MAINTAG_SUBOPTION_NUMBER 500 /**< Maximum number of suboptions for each &MAINTAG */
#define MAX_MAINTAG_SUBOPTION_LINE_LENGTH 1000 /**< Maximum length of suboption element for each &MAINTAG */


/* utils -------------------------------------------------------*/
double* MatrixXcdToDouble1d(MatrixXcd mat);
MatrixXcd Double1dToMatrixXcd(double* val, int n);

/* math functions -----------------------------------------------*/

// geometry
double dist(double spin1[],double spin2[]);
double cosTheta(double spin1[], double spin2[],double dist);
double sinTheta(double spin1[], double spin2[], double dist);
double cosPhi(double spin1[], double spin2[]);
double sinPhi(double spin1[], double spin2[]);

/* Physics functions --------------------------------------------*/

// spin state control
float* substates(float S);
MatrixXcd getSpinor(float S, float ms);
MatrixXcd kron(MatrixXcd a, MatrixXcd b);
MatrixXcd partialtrace(MatrixXcd Mij, int dimrow, int dimcol);
double calNorm(MatrixXcd m);
int normalize(MatrixXcd* m);
float findZbasisSubLevel(MatrixXcd spinor);

MatrixXcd powMatrixXcdElementWise(MatrixXcd a, int n);
MatrixXcd mulMatrixXcdElementWise(MatrixXcd a, MatrixXcd b);

/* Easy print ---------------------------------------------------*/
void printInlineMatrixXcd(char* key, MatrixXcd mat);
void printStructElementChar(char* key,char* val);
void printStructElementChar2d(char* key, char** val, int n);
void printStructElementInt(char* key, int val);
void printStructElementInt1dIdx(char* key, int* val, int n);
void printStructElementFloat(char* key, float val);
void printStructElementFloat1d(char* key, float* val, int n);
void printStructElementDouble(char* key, double val);
void printStructElementDouble1d(char* key, double* val, int n);
void printStructElementBool(char* key, bool val);
void printLine();
void printLineSection();
void printTitle(char* title);
void printSubTitle(char* title);
void printMessage(char* message);

/* Find index ---------------------------------------------------*/
int findIndexInt(int* array, int ista, int iend, int val); // find in the range of ista <= i <= iend
int findIndexCharFix(char array[][MAX_CHARARRAY_LENGTH], int ista, int iend, char* val); // find in the range of ista <= i <= iend, strcasecmp
int findIndexChar(char** array, int ista, int iend, char* val); // find in the range of ista <= i <= iend, strcasecmp

/* Type checker ---------------------------------------------------*/
int isStringDouble(char *s);

/* Quick Sort --------------------------------------------------*/
void QuickSort(int** d_Array, int left, int right);
int Partition(int** d_Array, int left, int right);
void Swap(int** d_Array, int a, int b);

/* MPI ---------------------------------------------------------*/
void para_range(int n1,int n2, int nprocs, int myrank, int*ista, int *iend);
int min(int x, int y);
int*** MPI_getLocalClusters(int order, int*** clusters);
// MatrixXcd* MPI_reduceLocalResult(int nstep, MatrixXcd* local);

#endif // __CCEX_UTILITIES_UTILITIES_H_


