#ifndef __CCEX_READER_H_
#define __CCEX_READER_H_

#include "qubit.h"
#include "bath.h"
#include "defect.h"
#include "general.h"


void readQubitfile(QubitArray* qa, Config* cnf);
void readBathfiles(BathArray* ba, QubitArray* qa, Config* cnf);

// Pure validation of the coordinate_frame_rotation input against the QubitArray.
// Opens no file and changes no state; call it once every option source is in and before
// any reader runs, so that an unsupported combination fails before a file is touched.
void validateCoordinateFrameRotationInputs(QubitArray* qa, Config* cnf);

/**
 * @enum HypfProvenanceState
 * @brief Where one stored hyperfine tensor came from, and so which frame it is in
 * @details HYPF_UNSET is not a fallback: readHftensorfile fills the whole array with it
 *          and then has to overwrite every entry, so a storage path that forgets to say
 *          where its tensor came from is a fatal error rather than a silent guess.
 *          It cannot be 0 for that reason -- allocInt2d zeroes, and 0 is a real state.
*/
typedef enum {
    HYPF_UNSET           = -1, /**< not yet decided; must not survive readHftensorfile */
    HYPF_FILE_TENSOR     =  0, /**< read from the tensor file, so in the hf_tensor_frame basis */
    HYPF_SOURCE_GEOMETRY =  1  /**< computed here from the source geometry, so bath frame */
} HypfProvenanceState;

/**
 * @brief Which frame each stored hyperfine tensor is in, per [bath spin][qubit]
 * @details Allocated and filled by readHftensorfile. applyCoordinateFrameRotation needs
 *          the distinction because both states can occur in the SAME run: a spin the
 *          tensor file does not cover falls back to a point-dipole tensor built from the
 *          source geometry even when every other spin was matched.
 * @ref HypfProvenanceState, readHftensorfile, applyCoordinateFrameRotation
*/
typedef int** HypfProvenance;

void applyCoordinateFrameRotation(BathArray* ba, QubitArray* qa, DefectArray* dfa, Config* cnf,
                                  HypfProvenance hypf_bathframe);
void readGyrofile(BathArray* ba, Config* cnf);
void readHftensorfile(BathArray* ba, QubitArray* qa, Config* config, HypfProvenance* hypf_bathframe);
// Kept for callers that do not need the provenance; allocates and releases it itself.
void readHftensorfile(BathArray* ba, QubitArray* qa, Config* config);
void readQdtensorfile(BathArray* ba, QubitArray* qa, Config* config);

void setBathStates(BathArray* ba, Config* cnf, int i); // read or random set
void setDefectPaxes(DefectArray* dfa, BathArray* ba, Config* cnf); // read or random set
void setSubbathStates(DefectArray* dfa, BathArray* ba, Config* cnf, int i); // read or random set

int compare_dist(const void *a, const void *b);

// Read tensor files
bool READ_BD_vertex(const char* inputfile, double*** vertex, double*** center, double*** normal, char** vertex_condi);
void READ_Tensor(const char* inputfile, double*** Tensor,int numCol);
void READ_Tensor_etc(const char* inputfile, char** names, int nspecies, double*** Tensor,char* condition, int numCol);
void READ_Tensor_const(const char* inputfile, char** names, int nspecies, double** Array,char* condition);
bool READ_Tensor_array(const char* inputfile, double(*Array)[3],char* condition,int numCol);
int READ_Tensor_ver(const char* inputfile,double* SpinFactor, double DefectTotSpin, double* CorrTotSpin);

void printHfInfo_version(int version, double DefectTotSpin, double CorrTotSpin, double SpinFactor);
void printHfInfo_BD(double** vertex, double** center, double** normal, double MinDif[3], double MaxDif[3], bool Usingvertex);
void printHfInfo_etc(double** A_Etc, double* A_Gfactor, char** names, int nspecies, int mode);
void printHfInfo_tensor(double** AtensorArray, int mode);


bool CheckPosition(const double Posi[3],double* refPosi);
bool FIND_AtomPosi(const double refPosi[3], double** TensorValue, int* num);
bool CheckBD_Range(const double Posi[3],const double minRange[3],const double maxRange[3],double err);
bool CheckBD_vertex(const double Posi[3],double** vertex, double** center, double** normal, double err);
bool CheckingUniqeBD(const double Posi[3], double** center, double** normal);
void CreatePlaneInfo(double** vertex, double*** center, double*** normal);
double ReDefinediff(double difXYZ,const  double minRange,const double maxRange,const double Copy_Length);
void CalCenter(double** center, double p1[], double p2[], double p3[], double p4[]);
void CalNormal(double** normal, double u[], double v[]);
void vector_diff(double (*result)[3], double p1[], double p2[]);
double vector_dot(double u[], double v[]);
void vector_cross(double** cross_P, double u[], double v[], bool norm);


#endif // __CCEX_READER_H_
