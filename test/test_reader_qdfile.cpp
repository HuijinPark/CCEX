#include "../include/bath.h"
#include "../include/memory.h"
#include "../include/utilities.h"
#include "../include/reader.h"
#include "../include/qubit.h"
#include "../include/general.h"

bool verbosity = true;
int rank = 0;
int nprocess = 1;

int main(int argc, char* argv[]){

    //////////////////////////////////////////////////////////////////
    // Required variables
    Config* cnf = Config_init();
    QubitArray* qa = QubitArray_init();

    // Set Config : bathfile, qubitfile, nbathfiles
    Config_setNbathfiles(cnf, 1);
    Config_allocBathfiles(cnf);
    
    Config_setBathfiles_i(cnf, "./../example/code_verification/CCE_Reprod/Bath_Data/4.hexagonal_data/bath_1", 0);

    Config_allocQubitfile(cnf);
    Config_setQubitfile(cnf, "./../example/code_verification/CCE_Reprod/Bath_Data/4.hexagonal_data/defect");

    Config_setQd_readmode(cnf, 1);    
    Config_allocQd_tensorfile(cnf);
    Config_setQd_tensorfile(cnf, "./../example/code_verification/CCE_Reprod/Bath_Data/4.hexagonal_data/Qfile_vertex");
    
    // Config_allocQd_tensorfile_woqubit(cnf, 1);
    // Config_setQd_tensorfile_woqubit(cnf, "./../example/code_verification/CCE_Reprod/Bath_Data/4.hexagonal_data/Qfile_vertex");

    double rbath = 5.0;
    Config_setRbath(cnf, rbath);
    
    // Report
    Config_report(cnf);
    //////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////
    // Read qubit file
    QubitArray_setNqubit(qa, 1);
    QubitArray_allocQubit(qa);
    
    // Qubit properties
    const double gamma_ele = GAMMA_ELECTRON;
    QubitArray_setQubit_i_gyro(qa,gamma_ele, 0);
    QubitArray_setQubit_i_name(qa, "qubit", 0);
    QubitArray_setQubit_i_spin(qa, 1.0, 0);

    // File read
    readQubitfile(qa, cnf);

    // Report
    QubitArray_report(qa);
    //////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////
    // BathArray : initialze
    BathArray* ba = BathArray_init();

    // Bath spin properties
    BathArray_setProp_nspecies(ba, 6); // B,N,C
    BathArray_allocProp(ba);

    // Boron
    BathArray_setProp_i(ba, 0, "10B", 2.8746786, 3  );
    BathArray_setProp_i(ba, 1, "11B", 8.5847044, 1.5);
    
    // Nitrogen
    BathArray_setProp_i(ba, 2, "14N", 1.933778, 1  );
    BathArray_setProp_i(ba, 3, "15N", -2.712622, 0.5);

    // Carbon
    BathArray_setProp_i(ba, 4, "12C", 0.0, 0.0);
    BathArray_setProp_i(ba, 5, "13C", 6.728284, 0.5);

    BathArray_report(ba);

    // Read bath file
    readBathfiles(ba, qa, cnf);
    
    // Relative distance difference
    double eta_x = 0.0;
    double eta_y = 0.0;
    double eta_z = 0.0;

    // Read hyperfine tensorfile
    Config_setHf_readmode(cnf, 0);
    readHftensorfile(ba, qa, cnf);

    // Read quadrupole tensorfile
    Config_setQd_readmode(cnf, 1);

    // Free
    Config_freeAll(cnf);
    QubitArray_freeAll(qa);
    BathArray_freeAll(ba);

    return 0;
}