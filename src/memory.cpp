#include <stdlib.h> 
#include "../include/memory.h"

// memory allocation for int
int* allocInt1d(size_t m) {
    return (int*)allocArray1d(m, sizeof(int));
}

int** allocInt2d(size_t m, size_t n) {
    return (int**)allocArray2d(m, n, sizeof(int));
}

int*** allocInt3d(size_t m, size_t n, size_t l) {
    return (int***)allocArray3d(m, n, l, sizeof(int));
}

void reallocInt1d(int **ptr, size_t newrow){
    reallocArray1d((void**)ptr,newrow,sizeof(int));
}

void reallocInt2d(int ***ptr, size_t oldrow, size_t newrow, size_t col){
    reallocArray2d((void***)ptr,oldrow,newrow,col,sizeof(int));
}

void reallocInt3d(int ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2){
    reallocArray3d((void****)ptr,oldrow,newrow,col1,col2,sizeof(int));
}

void freeInt1d(int **ptr) {
    freeArray1d((void**)ptr);
}

void freeInt2d(int ***ptr, size_t m) {
    freeArray2d((void***)ptr, m);
}

void freeInt3d(int ****ptr, size_t m, size_t n) {
    freeArray3d((void****)ptr, m, n);
}

// memory allocation for double
double* allocDouble1d(size_t m) {
    return (double*)allocArray1d(m, sizeof(double));
}

double** allocDouble2d(size_t m, size_t n) {
    return (double**)allocArray2d(m, n, sizeof(double));
}

double*** allocDouble3d(size_t m, size_t n, size_t l) {
    return (double***)allocArray3d(m, n, l, sizeof(double));
}

void reallocDouble1d(double **ptr, size_t newrow){
    reallocArray1d((void**)ptr,newrow,sizeof(double));
}

void reallocDouble2d(double ***ptr, size_t oldrow, size_t newrow, size_t col){
    reallocArray2d((void***)ptr,oldrow,newrow,col,sizeof(double));
}

void reallocDouble3d(double ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2){
    reallocArray3d((void****)ptr,oldrow,newrow,col1,col2,sizeof(double));
}

void freeDouble1d(double **ptr) {
    freeArray1d((void**)ptr);
}

void freeDouble2d(double ***ptr, size_t m) {
    freeArray2d((void***)ptr, m);
}

void freeDouble3d(double ****ptr, size_t m, size_t n) {
    freeArray3d((void****)ptr, m, n);
}

// memory allocation for float
float* allocFloat1d(size_t m) {
    return (float*)allocArray1d(m, sizeof(float));
}

float** allocFloat2d(size_t m, size_t n) {
    return (float**)allocArray2d(m, n, sizeof(float));
}

float*** allocFloat3d(size_t m, size_t n, size_t l) {
    return (float***)allocArray3d(m, n, l, sizeof(float));
}

void reallocFloat1d(float **ptr, size_t newrow){
    reallocArray1d((void**)ptr,newrow,sizeof(float));
}

void reallocFloat2d(float ***ptr, size_t oldrow, size_t newrow, size_t col){
    reallocArray2d((void***)ptr,oldrow,newrow,col,sizeof(float));
}

void reallocFloat3d(float ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2){
    reallocArray3d((void****)ptr,oldrow,newrow,col1,col2,sizeof(float));
}

void freeFloat1d(float **ptr) {
    freeArray1d((void**)ptr);
}

void freeFloat2d(float ***ptr, size_t m) {
    freeArray2d((void***)ptr, m);
}

void freeFloat3d(float ****ptr, size_t m, size_t n) {
    freeArray3d((void****)ptr, m, n);
}

// memory allocation for char
char* allocChar1d(size_t m) {
    return (char*)allocArray1d(m, sizeof(char));
}

char** allocChar2d(size_t m, size_t n) {
    return (char**)allocArray2d(m, n, sizeof(char));
}

char*** allocChar3d(size_t m, size_t n, size_t l) {
    return (char***)allocArray3d(m, n, l, sizeof(char));
}

void reallocChar1d(char **ptr, size_t newrow){
    reallocArray1d((void**)ptr,newrow,sizeof(char));
}

void reallocChar2d(char ***ptr, size_t oldrow, size_t newrow, size_t col){
    reallocArray2d((void***)ptr,oldrow,newrow,col,sizeof(char));
}

void reallocChar3d(char ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2){
    reallocArray3d((void****)ptr,oldrow,newrow,col1,col2,sizeof(char));
}

void freeChar1d(char **ptr) {
    freeArray1d((void**)ptr);
}

void freeChar2d(char ***ptr, size_t m) {
    freeArray2d((void***)ptr, m);
}

void freeChar3d(char ****ptr, size_t m, size_t n) {
    freeArray3d((void****)ptr, m, n);
}

// MatrixXcd
MatrixXcd* allocMatrixXcd1d(size_t m) {
    MatrixXcd* array = new MatrixXcd[m];
    return array;
}

MatrixXcd** allocMatrixXcd2d(size_t m, size_t n) {
    MatrixXcd** array = new MatrixXcd*[m];
    for (size_t i = 0; i < m; i++) {
        array[i] = new MatrixXcd[n];
    }
    return array;
}

MatrixXcd*** allocMatrixXcd3d(size_t m, size_t n, size_t l) {
    MatrixXcd*** array = new MatrixXcd**[m];
    for (size_t i = 0; i < m; i++) {
        array[i] = new MatrixXcd*[n];
        for (size_t j = 0; j < n; j++) {
            array[i][j] = new MatrixXcd[l];
        }
    }
    return array;
}

void freeMatrixXcd1d(MatrixXcd **ptr) {
    if (!ptr) return;
    delete [] *ptr;
}

void freeMatrixXcd2d(MatrixXcd ***ptr, size_t m) {
    if (!ptr) return;
    for (size_t i = 0; i < m; i++) {
        delete [] (*ptr)[i];
    }
    delete [] *ptr;
}

void freeMatrixXcd3d(MatrixXcd ****ptr, size_t m, size_t n) {
    if (!ptr) return;
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            delete [] (*ptr)[i][j];
        }
        delete [] (*ptr)[i];
    }
    delete [] *ptr;
}



// DoubleTensor
DoubleTensor* allocDoubleTensor1d(size_t m){
    return (DoubleTensor*)allocArray1d(m, sizeof(DoubleTensor));
}


DoubleTensor** allocDoubleTensor2d(size_t m, size_t n){
    return (DoubleTensor**)allocArray2d(m, n, sizeof(DoubleTensor));
}


DoubleTensor*** allocDoubleTensor3d(size_t m, size_t n, size_t l){
    return (DoubleTensor***)allocArray3d(m, n, l, sizeof(DoubleTensor));
}



void reallocDoubleTensor1d(DoubleTensor **ptr, size_t newrow){
    reallocArray1d((void**)ptr,newrow,sizeof(DoubleTensor));
}


void reallocDoubleTensor2d(DoubleTensor ***ptr, size_t oldrow, size_t newrow, size_t col){
    reallocArray2d((void***)ptr,oldrow,newrow,col,sizeof(DoubleTensor));
}


void reallocDoubleTensor3d(DoubleTensor ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2){
    reallocArray3d((void****)ptr,oldrow,newrow,col1,col2,sizeof(DoubleTensor));
}



void freeDoubleTensor1d(DoubleTensor **ptr){
    freeArray1d((void**)ptr);
}


void freeDoubleTensor2d(DoubleTensor ***ptr, size_t m){
    freeArray2d((void***)ptr, m);
}


void freeDoubleTensor3d(DoubleTensor ****ptr, size_t m, size_t n){
    freeArray3d((void****)ptr, m, n);
}


// DoubleComplex
DoubleComplex* allocDoubleComplex1d(size_t m){
    return (DoubleComplex*)allocArray1d(m, sizeof(DoubleComplex));
}


DoubleComplex** allocDoubleComplex2d(size_t m, size_t n){
    return (DoubleComplex**)allocArray2d(m, n, sizeof(DoubleComplex));
}


DoubleComplex*** allocDoubleComplex3d(size_t m, size_t n, size_t l){
    return (DoubleComplex***)allocArray3d(m, n, l, sizeof(DoubleComplex));
}



void reallocDoubleComplex1d(DoubleComplex **ptr, size_t newrow){
    reallocArray1d((void**)ptr,newrow,sizeof(DoubleComplex));
}


void reallocDoubleComplex2d(DoubleComplex ***ptr, size_t oldrow, size_t newrow, size_t col){
    reallocArray2d((void***)ptr,oldrow,newrow,col,sizeof(DoubleComplex));
}


void reallocDoubleComplex3d(DoubleComplex ****ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2){
    reallocArray3d((void****)ptr,oldrow,newrow,col1,col2,sizeof(DoubleComplex));
}



void freeDoubleComplex1d(DoubleComplex **ptr){
    freeArray1d((void**)ptr);
}


void freeDoubleComplex2d(DoubleComplex ***ptr, size_t m){
    freeArray2d((void***)ptr, m);
}


void freeDoubleComplex3d(DoubleComplex ****ptr, size_t m, size_t n){
    freeArray3d((void****)ptr, m, n);
}



// basic memory allocation

// memory allocation for 1d array
void* allocArray1d(size_t m, size_t size) {
    void *ptr = calloc(m, size);
    if (ptr == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

// memory allocation for 2d array
void** allocArray2d(size_t m, size_t n, size_t size) {
    void **ptr = (void**)calloc(m, sizeof(void*));
    if (ptr == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < m; i++) {
        ptr[i] = allocArray1d(n,size);
    }
    return ptr;
}

// memory allocation for 3d array
void*** allocArray3d(size_t m, size_t n, size_t l, size_t size) {
    void ***ptr = (void***)calloc(m, sizeof(void**));
    if (ptr == NULL) {
        perror("Memory allocation failed");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < m; i++) {
        ptr[i] = allocArray2d(n,l,size);
    }
    return ptr;
}


void reallocArray1d(void** ptr, size_t newrow, size_t size){
    *ptr = (void*)realloc(*ptr,size*newrow);
    if (*ptr == NULL){
        perror("Memory reallocation failed");
        exit(EXIT_FAILURE);
    }
}


void reallocArray2d(void*** ptr, size_t oldrow, size_t newrow, size_t col, size_t size){
    *ptr = (void**)realloc(*ptr,sizeof(void*)*newrow);
    for (int row = oldrow; row < newrow; row++){
        (*ptr)[row] = (void*)allocArray1d(col,size);
    }
    if (*ptr == NULL){
        perror("Memory reallocation failed");
        exit(EXIT_FAILURE);
    }
}

void reallocArray3d(void**** ptr, size_t oldrow, size_t newrow, size_t col1, size_t col2, size_t size){
    *ptr = (void***)realloc(*ptr,sizeof(void**)*newrow);
    for (int row = oldrow; row < newrow; row++){
        (*ptr)[row] = (void**)allocArray2d(col1,col2,size);
    }
    if (*ptr == NULL){
        perror("Memory reallocation failed");
        exit(EXIT_FAILURE);
    }
}

//free memory 1d array
void freeArray1d(void **ptr) {
    if (ptr==NULL || *ptr==NULL) {
        return;
    }
    free(*ptr);
    *ptr=NULL;
}



//free memory 2d array
void freeArray2d(void ***ptr, size_t m) {
    if (ptr == NULL || *ptr == NULL) {
        return;
    }
    for (size_t i = 0; i < m; i++) {
        freeArray1d(&((*ptr)[i]));
    }
    free(*ptr);
    *ptr = NULL;
}

//free memory 3d array
void freeArray3d(void ****ptr, size_t m, size_t n) {
    if (ptr==NULL || *ptr==NULL) {
        return;
    }
    for (size_t i = 0; i < m; i++) {
        freeArray2d(&((*ptr)[i]), n);
    }
    free(*ptr);
    *ptr=NULL;
}

/* Copy dynamically allocated array ------------------------*/

// int
void copyInt1d(int* dest, const int* src, size_t m){
    memcpy(dest,src,sizeof(int)*m);
}

void copyInt1dPart(int* dest, const int* src, size_t ista, size_t iend){
    memcpy(dest,src+ista,sizeof(int)*(iend-ista+1));
}

void copyInt2d(int** dest, const int** src, size_t m, size_t n){
    for (int i=0; i<m; i++){
        memcpy(dest[i],src[i],sizeof(int)*n);
    }
}

void copyInt3d(int*** dest, const int*** src, size_t m, size_t n, size_t l){
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            memcpy(dest[i][j],src[i][j],sizeof(int)*l);
        }
    }
}


// double
void copyDouble1d(double* dest, const double* src, size_t m){
    memcpy(dest,src,sizeof(double)*m);
}

void copyDouble2d(double** dest, const double** src, size_t m, size_t n){
    for (int i=0; i<m; i++){
        memcpy(dest[i],src[i],sizeof(double)*n);
    }
}

void copyDouble3d(double*** dest, const double*** src, size_t m, size_t n, size_t l){
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            memcpy(dest[i][j],src[i][j],sizeof(double)*l);
        }
    }
}

// float
void copyFloat1d(float* dest, const float* src, size_t m){
    memcpy(dest,src,sizeof(float)*m);
}

void copyFloat2d(float** dest, const float** src, size_t m, size_t n){
    for (int i=0; i<m; i++){
        memcpy(dest[i],src[i],sizeof(float)*n);
    }
}

void copyFloat3d(float*** dest, const float*** src, size_t m, size_t n, size_t l){
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            memcpy(dest[i][j],src[i][j],sizeof(float)*l);
        }
    }
}

// DoubleComplex
void copyDoubleComplex1d(DoubleComplex *dest, const DoubleComplex *src, size_t m){
    memcpy(dest,src,sizeof(DoubleComplex)*m);
}

void copyDoubleComplex2d(DoubleComplex **dest, const DoubleComplex **src, size_t m, size_t n){
    for (int i=0; i<m; i++){
        memcpy(dest[i],src[i],sizeof(DoubleComplex)*n);
    }
}

void copyDoubleComplex3d(DoubleComplex ***dest, const DoubleComplex ***src, size_t m, size_t n, size_t l){
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            memcpy(dest[i][j],src[i][j],sizeof(DoubleComplex)*l);
        }
    }
}

// DoubleTensor
void copyDoubleTensor(DoubleTensor *dest, const DoubleTensor src){
    dest->xx = src.xx;    dest->xy = src.xy;    dest->xz = src.xz;
    dest->yx = src.yx;    dest->yy = src.yy;    dest->yz = src.yz;
    dest->zx = src.zx;    dest->zy = src.zy;    dest->zz = src.zz;
}

void copyDoubleTensor1d(DoubleTensor **dest, const DoubleTensor *src, size_t m){
    memcpy(dest,src,sizeof(DoubleTensor)*m);
}