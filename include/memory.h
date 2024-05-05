#ifndef __CCEX_MEMORY_H_
#define __CCEX_MEMORY_H_

#ifndef __CCEX_UTILITIES_H_
    #include "utilities.h"
#endif

// what to do : copy, alloc, realloc, free

/* Memory control for each type ----------------------------*/

// int
int* allocInt1d(size_t m);
int** allocInt2d(size_t m, size_t n);
int*** allocInt3d(size_t m, size_t n, size_t l);

void reallocInt1d(int **ptr, size_t newrow);
void reallocInt2d(int ***ptr, size_t oldrow, size_t newrow, size_t col);
void reallocInt3d(int ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2);

void freeInt1d(int *ptr);
void freeInt2d(int **ptr, size_t m);
void freeInt3d(int ***ptr, size_t m, size_t n);

// double
double* allocDouble1d(size_t m);
double** allocDouble2d(size_t m, size_t n);
double*** allocDouble3d(size_t m, size_t n, size_t l);

void reallocDouble1d(double **ptr, size_t newrow);
void reallocDouble2d(double ***ptr, size_t oldrow, size_t newrow, size_t col);
void reallocDouble3d(double ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2);

void freeDouble1d(double *ptr);
void freeDouble2d(double **ptr, size_t m);
void freeDouble3d(double ***ptr, size_t m, size_t n);

// float
float* allocFloat1d(size_t m);
float** allocFloat2d(size_t m, size_t n);
float*** allocFloat3d(size_t m, size_t n, size_t l);

void reallocFloat1d(float **ptr, size_t newrow);
void reallocFloat2d(float ***ptr, size_t oldrow, size_t newrow, size_t col);
void reallocFloat3d(float ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2);

void freeFloat1d(float *ptr);
void freeFloat2d(float **ptr, size_t m);
void freeFloat3d(float ***ptr, size_t m, size_t n);

// char
char* allocChar1d(size_t m);
char** allocChar2d(size_t m, size_t n);
char*** allocChar3d(size_t m, size_t n, size_t l);

void reallocChar1d(char **ptr, size_t newrow);
void reallocChar2d(char ***ptr, size_t oldrow, size_t newrow, size_t col);
void reallocChar3d(char ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2);

void freeChar1d(char *ptr);
void freeChar2d(char **ptr, size_t m);
void freeChar3d(char ***ptr, size_t m, size_t n);

// DoubleTensor
DoubleTensor* allocDoubleTensor1d(size_t m);
DoubleTensor** allocDoubleTensor2d(size_t m, size_t n);
DoubleTensor*** allocDoubleTensor3d(size_t m, size_t n, size_t l);

void reallocDoubleTensor1d(DoubleTensor **ptr, size_t newrow);
void reallocDoubleTensor2d(DoubleTensor ***ptr, size_t oldrow, size_t newrow, size_t col);
void reallocDoubleTensor3d(DoubleTensor ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2);

void freeDoubleTensor1d(DoubleTensor *ptr);
void freeDoubleTensor2d(DoubleTensor **ptr, size_t m);
void freeDoubleTensor3d(DoubleTensor ***ptr, size_t m, size_t n);

/* Basic Memory control ------------------------------------*/

// basic memory allocation
void* allocArray1d(size_t m, size_t size);
void** allocArray2d(size_t m, size_t n, size_t size);
void*** allocArray3d(size_t m, size_t n, size_t l, size_t size);

// basic memory reallocation
void reallocArray1d(void **ptr, size_t newrow, size_t size);
void reallocArray2d(void ***ptr, size_t oldrow, size_t newrow, size_t col, size_t size);
void reallocArray3d(void ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2, size_t size);

// basic memory free
void freeArray1d(void *ptr);
void freeArray2d(void **ptr, size_t m);
void freeArray3d(void ***ptr, size_t m, size_t n);


/* Copy dynamically allocated array ------------------------*/

// int
void copyInt1d(int* dest, const int* src, size_t m);
void copyInt1dPart(int* dest, const int* src, size_t ista, size_t iend);
void copyInt2d(int** dest, const int** src, size_t m, size_t n);
void copyInt3d(int*** dest, const int*** src, size_t m, size_t n, size_t l);

// double
void copyDouble1d(double *dest, const double *src, size_t m);
void copyDouble2d(double **dest, const double **src, size_t m, size_t n);
void copyDouble3d(double ***dest, const double ***src, size_t m, size_t n, size_t l);

// float
void copyFloat1d(float *dest, const float *src, size_t m);
void copyFloat2d(float **dest, const float **src, size_t m, size_t n);
void copyFloat3d(float ***dest, const float ***src, size_t m, size_t n, size_t l);

// DoubleComplex
void copyDoubleComplex1d(DoubleComplex *dest, const DoubleComplex *src, size_t m);
void copyDoubleComplex2d(DoubleComplex **dest, const DoubleComplex **src, size_t m, size_t n);
void copyDoubleComplex3d(DoubleComplex ***dest, const DoubleComplex ***src, size_t m, size_t n, size_t l);

// DoubleTensor
void copyDoubleTensor(DoubleTensor *dest, const DoubleTensor src);
void copyDoubleTensor1d(DoubleTensor **dest, const DoubleTensor *src, size_t m);

#endif // __CCEX_MEMORY_H_