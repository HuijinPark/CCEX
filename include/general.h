#ifndef __CCEX_SIMULATOR_GENERAL_H_
#define __CCEX_SIMULATOR_GENERAL_H_

#include "utilities.h"

/**
 * @struct CalConfig
 * @brief This structure contains the parameters for the simulation.
 * @details Each parameter is read from the input file or options.
 *          and those parameters are used in the simulation directly.
 */
typedef struct {

    // Method
    char method[MAX_CHARARRAY_LENGTH];    /**< CCE method : ensemble  | single | gcce */
    char quantity[MAX_CHARARRAY_LENGTH];  /**< Measurment : coherence | noise  | dm   */

    // General options
    int   order;        /**< Maximum nuclear spins included in a cluster ( >= 0 )*/
    float bfield[3];    /**< Magnetic field (Bx,Bx,Bz) (Unit G)                  */
    float rbath;        /**< Bath radius (Unit angstrom)                         */
    float rdip;         /**< Dipole radius (Unit angstrom)                       */
    float deltat;       /**< Time interval (Unit ms)                             */
    int   nstep;        /**< Time step                                           */
    float rbathcut;     /**< Cutoff to remove spins within rbathcut              */
    float rdipcut;      /**< Cutoff to remove spin pairs within rdipcut          */
    int   nstate;       /**< Number of bath states (single and gcce method only) */
    int   seed;         /**< Fixed random value (if none, time(null) )           */

    // Qubit and Bath file options
    char* qubitfile; /**< Qubit file name */
    char* gyrofile; /**< Gyro file name */
    int nbathfiles; /**< Number of bath files */
    char** bathfiles; /**< Bath file names */
    double** bathadjust; /**< Bath adjust values */
    char* avaaxfile; /**< Available principal axis */
    char* statefile; /**< State file */
    char* exstatefile ; /**< extra spin' state file */
    int* _flines; /**< file line number of each bath spins */

    // tensorfile-related
    double DefectTotSpin; /**< Total spin of defect : default : 1 */
                        // Central spin number in reading Atensor file (default S = 1)
                        // It is nesessary for calulation of Hyperfin tensor with Atensor file
                        // <ex.1>
                        // If you use the Q.E., you will set this value which you want (S=5/2 --> DefectTotSpin = 5/2)
                        // <ex.2>
                        // If using VASP version, you don't consider the DefectTotSpin value (DefectTotSpin = 1)
                        // <ex.3>
                        // But if you wat to use other central spin which is different from A-file in VASP version,
                        // you have to set the central spin which you want (S=3/2 --> DefectTotSpin =3/2)
                        // , and then you will change CorrTotSpin !!!

    double CorrTotSpin; /**< Total spin of correlation : default : 0 */
                        // Correlation value of total spin in A-tensor (default : NULL)
                        // <ex.1>
                        // If A-tensor file of Q.E., automatically it adapt the central spin S = 1/2 into A-file
                        // --> CorrTotSpin = 1/2
                        // <ex.2>
                        // If using VASP version, you don't consider the CorrTotSpin value (CorrTotSpin = NULL)
                        // <ex.3>
                        // But if you want to use other central spin which is different from A-file in VASP version,
                        // you have to set the central spin of A-file (VASP; Half of Total magnetic moment (S/2))

    // DFT Hyperfine tensor
    int hf_readmode; /**< Read mode for hyperfine tensor file : 0 : off , 1~3 : fermi & DFT */ 
    char* hf_tensorfile; /**< Hyperfine tensor file name */
    double hf_cutoff ; /**< Cutoff for hyperfine tensor file */
    int hf_ignore_oor; /**< Ignore nuclear spins out of range the hyperfine tensor file */
    

    // DFT Quadrupole tensor
    int qd_readmode; /**< Read mode for quadrupole tensor file : 0 : off , 1 : exp , 2: dft */
    char* qd_tensorfile; /**< Quadrupole tensor file name */
    char* qd_tensorfile_woqubit ; /**< Quadrupole tensor file name without qubit */
    

    // semi-classical options
    // int Interval_filter; // Default : 500      //gsl interval in integration
    // int Step_filter; // Default : 100          //filter function time table step
    // double epsabs_filter; // Default : 1e-10
    // double tolerance_filter; // Default : 1e-10
    // double RoundOff_err_filter; // Default : 1e+6
    // int SaveCorr; // Default : 1    //1 : save the correlation function

    //filter function time table
    // std::vector<double> FilterTable;    //[0 ~ 1], step : Step_filter
    
} Config;

Config* Config_init();
void Config_freeAll(Config* cnf);
void Config_report(Config* cnf);

/* Low level ----------------------------------------------------*/

// get
char*   Config_getMethod(Config* cnf);
char*   Config_getQuantity(Config* cnf);
int     Config_getOrder(Config* cnf);
float*  Config_getBfield(Config* cnf);
float   Config_getRbath(Config* cnf);
float   Config_getRdip(Config* cnf);
float   Config_getDeltat(Config* cnf);
int     Config_getNstep(Config* cnf);
float   Config_getRbathcut(Config* cnf);
float   Config_getRdipcut(Config* cnf);
int     Config_getNstate(Config* cnf);
int     Config_getSeed(Config* cnf);
char*   Config_getQubitfile(Config* cnf);
char*   Config_getGyrofile(Config* cnf);
int     Config_getNbathfiles(Config* cnf);
char*   Config_getBathfiles_i(Config* cnf,int i);
double* Config_getBathadjust_i(Config* cnf,int i);
char*   Config_getAvaaxfile(Config* cnf);
char*   Config_getStatefile(Config* cnf);
char*   Config_getExstatefile(Config* cnf);
double  Config_getDefectTotSpin(Config* cnf);
double  Config_getCorrTotSpin(Config* cnf);
char*   Config_getHf_tensorfile(Config* cnf);
double  Config_getHf_cutoff(Config* cnf);
int     Config_getHf_ignore_oor(Config* cnf);
int     Config_getHf_readmode(Config* cnf);
char*   Config_getQd_tensorfile(Config* cnf);
char*   Config_getQd_tensorfile_woqubit(Config* cnf);
int     Config_getQd_readmode(Config* cnf);

// set 
void Config_setMethod(Config* cnf, char* method);
void Config_setQuantity(Config* cnf, char* quantity);
void Config_setOrder(Config* cnf, int order);
void Config_setBfield(Config* cnf, float* bfield);
void Config_setBfield_z(Config* cnf, float bz);
void Config_setRbath(Config* cnf, float rbath);
void Config_setRdip(Config* cnf, float rdip);
void Config_setDeltat(Config* cnf, float deltat);
void Config_setNstep(Config* cnf, int nstep);
void Config_setRbathcut(Config* cnf, float rbathcut);
void Config_setRdipcut(Config* cnf, float rdipcut);
void Config_setNstate(Config* cnf, int nstate);
void Config_setSeed(Config* cnf, int seed);
void Config_setQubitfile(Config* cnf, char* qubitfile);
void Config_setGyrofile(Config* cnf, char* gyrofile);
void Config_setNbathfiles(Config* cnf, int nbathfiles);
void Config_setBathfiles_i(Config* cnf, char* bathfiles, int i);
void Config_setBathadjust_i(Config* cnf, double* bathadjust, int i);
void Config_setAvaaxfile(Config* cnf, char* avaaxfile);
void Config_setStatefile(Config* cnf, char* statefile);
void Config_setExstatefile(Config* cnf, char* exstatefile);

void Config_setDefectTotSpin(Config* cnf, double DefectTotSpin);
void Config_setCorrTotSpin(Config* cnf, double CorrTotSpin);
void Config_setHf_tensorfile(Config* cnf, char* hf_tensorfile);
void Config_setHf_cutoff(Config* cnf, double hf_cutoff);
void Config_setHf_ignore_oor(Config* cnf, int hf_ignore_oor);
void Config_setHf_readmode(Config* cnf, int hf_readmode);
void Config_setQd_tensorfile(Config* cnf, char* qd_tensorfile);
void Config_setQd_tensorfile_woqubit(Config* cnf, char* qd_tensorfile_woqubit);
void Config_setQd_readmode(Config* cnf, int qd_readmode);

// alloc
void Config_allocBathfiles(Config* cnf);
void Config_reallocBathfiles(Config* cnf, int oldsize, int newsize);
void Config_allocBathadjust(Config* cnf);
void Config_allocGyrofile(Config* cnf);
void Config_allocQubitfile(Config* cnf);
void Config_allocAvaaxfile(Config* cnf);
void Config_allocStatefile(Config* cnf);
void Config_allocExstatefile(Config* cnf);
void Config_allocHf_tensorfile(Config* cnf);
void Config_allocQd_tensorfile(Config* cnf);
void Config_allocQd_tensorfile_woqubit(Config* cnf);


// free
void Config_freeBathfiles(Config* cnf);
void Config_freeBathadjust(Config* cnf);
void Config_freeGyrofile(Config* cnf);
void Config_freeQubitfile(Config* cnf);
void Config_freeAvaaxfile(Config* cnf);
void Config_freeStatefile(Config* cnf);
void Config_freeExstatefile(Config* cnf);
void Config_freeHf_tensorfile(Config* cnf);
void Config_freeQd_tensorfile(Config* cnf);
void Config_freeQd_tensorfile_woqubit(Config* cnf);





#endif // __CCEX_SIMULATOR_GENERAL_H_