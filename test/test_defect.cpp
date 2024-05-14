#include "../include/qubit.h"
#include "../include/defect.h"
#include "../include/bath.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include <iostream>

bool verbosity = true;
int rank = 0;
int nprocess = 1;
int seed = time(NULL);

int main(int argc, char* argv[]){

    srand(seed);
    printStructElementInt("seed",seed);

    // Diamond NV center in P1 bath

    QubitArray* qa = QubitArray_init();
    BathArray* ba = BathArray_init();
    DefectArray* dfa = DefectArray_init();

    ////////////////////////////////////////////
    // Qubit
    ////////////////////////////////////////////
    // Allocate Qubit memory
    int nqubit = 1;
    QubitArray_setNqubit(qa,nqubit);
    QubitArray_allocQubit(qa);
    QubitArray_allocIntmap(qa);

    // Qubit spin properties
    int iq = 0;
    char name[MAX_CHARARRAY_LENGTH] = "NV";
    float spin = 1.0;
    double xyz[3] = {0.0,0.0,0.0};
    double gyro = GAMMA_ELECTRON;

    QubitArray_setQubit_i_name(qa,name,iq); // name
    QubitArray_setQubit_i_spin(qa,spin,iq); // spin
    QubitArray_setQubit_i_xyz(qa,xyz,iq); // xyz
    QubitArray_setQubit_i_gyro(qa,gyro,iq); // gyro

    // Overhaus on/off
    bool io_overhaus = true;
    QubitArray_setOverhaus(qa,io_overhaus);

    // Qubit detuning and overhaus
    double detuning = 0.0;
    double overhaus = 0.0;
    QubitArray_setQubit_i_detuning(qa,detuning,iq); // detuning
    QubitArray_setQubit_i_overhaus(qa,overhaus,iq); // overhaus

    // Qubit alpha and beta
    float alpha_ms = 1.0;
    float beta_ms = 0.0;
    QubitArray_setQubit_i_alpha_fromMs(qa,alpha_ms,iq); //alpha
    QubitArray_setQubit_i_beta_fromMs(qa,beta_ms,iq); //beta

    // Report
    QubitArray_reportQubit_i(qa,iq);

    // Qubit psi0, psia, psib
    QubitArray_setPsiaPsib_fromQubit(qa); // psia = kron(qubit_i_alpha,qubit_j_alpha, ...)
    QubitArray_setPsi0_fromPsiaPsib(qa); // psi0 = Normalize(psia + psib)

    QubitArray_report(qa);

    ////////////////////////////////////////////
    // BathArray
    ////////////////////////////////////////////

    // Allocate Bath memory
    int nbathspin = 5;
    BathArray_setNspin(ba,nbathspin);
    BathArray_allocBath(ba,nqubit);
    
    // Bath spin properties
    for (int ibs=0;ibs<nbathspin;ibs++){
        double gyro = GAMMA_ELECTRON;
        BathArray_setBath_i_name(ba,"P1",ibs); // name
        BathArray_setBath_i_spin(ba,0.5,ibs); // spin
        BathArray_setBath_i_gyro(ba,gyro,ibs); // gyro
        double xyz[3] = {10.0,10.0,10.0};
        xyz[0] = 10.0 + 10.0*ibs;
        xyz[1] = 10.0 + 10.0*ibs;
        xyz[2] = 10.0 + 10.0*ibs;
        BathArray_setBath_i_xyz(ba,xyz,ibs); // xyz
    }
    
    BathArray_reportBath(ba);

    // Bath detuning and disorder
    for (int ibs=0;ibs<nbathspin;ibs++){
        BathArray_setBath_i_detuning(ba,0.0,ibs); // detuning
        BathArray_setBath_i_disorder(ba,0.0,ibs); // disorder
    }

    // Bath hyperfine tensors
    BathArray_setBathHypfs(ba,qa); // set hyperfine tensors
    BathArray_reportBath_hypf(ba,nqubit);

    ////////////////////////////////////////////
    // DefectArray
    ////////////////////////////////////////////

    // Allocate Defect memory
    int ndefect = 1;
    DefectArray_setNdefect(dfa,ndefect);
    DefectArray_allocDefect(dfa);

    // Defect properties
    char dfname[MAX_CHARARRAY_LENGTH] = "N";
    bool apprx = true;
    int naddspin = 1;
    char** types = allocChar2d(1,MAX_CHARARRAY_LENGTH);    strcpy(types[0],"N");
    float* spins = allocFloat1d(1);                        spins[0] = 1.0;
    double* gyros = allocDouble1d(1);                      gyros[0] = GAMMA_ELECTRON;
    double* eqs = allocDouble1d(1);                        eqs[0] = 2.044;
    int navaax = 4 + 1;
    
    DefectArray_allocDefect_idf(dfa,0,navaax,naddspin);
    DefectArray_setDefect_idf_dfname(dfa,0,dfname); // name
    DefectArray_setDefect_idf_apprx(dfa,0,apprx); // approximation
    DefectArray_setDefect_idf_naddspin(dfa,0,naddspin); // naddspin
    DefectArray_setDefect_idf_types(dfa,0,types); // types
    DefectArray_setDefect_idf_spins(dfa,0,spins); // spins
    DefectArray_setDefect_idf_gyros(dfa,0,gyros); // gyros
    DefectArray_setDefect_idf_eqs(dfa,0,eqs); // eqs
    DefectArray_setDefect_idf_navaax(dfa,0,navaax); // navaax

    // DefectArray_reportDefect_idf(dfa,0);

    freeChar2d(&types,1);
    freeFloat1d(&spins);
    freeDouble1d(&gyros);
    freeDouble1d(&eqs);

    // Defect principal axis
    double rxyz[3] = {0.0,0.0,0.0};
    rxyz[0] = 0.1; rxyz[1] = 0.1; rxyz[2] = 0.1;
    DefectArray_setDefect_idf_iax_isp_rxyz(dfa,0,1,0,rxyz); // rxyz

    rxyz[0] = 0.2; rxyz[1] = 0.2; rxyz[2] = 0.2;
    DefectArray_setDefect_idf_iax_isp_rxyz(dfa,0,2,0,rxyz); // rxyz

    rxyz[0] = 0.3; rxyz[1] = 0.3; rxyz[2] = 0.3;
    DefectArray_setDefect_idf_iax_isp_rxyz(dfa,0,3,0,rxyz); // rxyz

    rxyz[0] = 0.4; rxyz[1] = 0.4; rxyz[2] = 0.4;
    DefectArray_setDefect_idf_iax_isp_rxyz(dfa,0,4,0,rxyz); // rxyz

    // DefectArray_reportDefect_idf(dfa,0);

    // Defect hyperfine tensors
    for (int iax=1; iax<navaax;iax++){
        DefectArray_setDefect_idf_iax_isp_hypf(dfa,0,iax,0,MatrixXcd::Random(3,3)); // hypf
        DefectArray_setDefect_idf_iax_isp_efg(dfa,0,iax,0,MatrixXcd::Random(3,3)); // quad
        DefectArray_setDefect_idf_iax_zfs(dfa,0,iax,MatrixXcd::Random(3,3)*1e6); // zfs
        
        // DefectArray_setDefect_idf_iax_detuning(dfa, int idf, int iax, double detuning);
    }

    // Report
    DefectArray_reportDefect_idf(dfa,0);

    ////////////////////////////////////////////
    // Generate Principal axis
    ////////////////////////////////////////////
    DefectArray_allocPaxes(dfa, nbathspin);
    DefectArray_setPaxesRandom(dfa,ba);
    DefectArray_reportPaxes(dfa);

    ////////////////////////////////////////////
    // Naddspins
    ////////////////////////////////////////////
    DefectArray_allocNaddspins(dfa, nbathspin);
    DefectArray_setNaddspins(dfa,ba);
    DefectArray_reportNaddspins(dfa);

    ////////////////////////////////////////////
    // DefectArray -> subbath
    ////////////////////////////////////////////
    /** What we do here :
     * 1. Set subbath name, xyz, spin, gyro, (eq,) hypf, mainspidx, 
     * 2. Set quad(zfs),hypf_sub (depending on principal axis) 
     * 2. Set detuning = 0.0 (No detuning)
     * 3. Set state = 0.0 and disorder = 0.0 (will be updated later)
     */
    DefectArray_allocSubbath(dfa,ba, nqubit);
    DefectArray_setSubbath(dfa,ba,qa); 
    DefectArray_reportSubbath_props(dfa);
    DefectArray_reportSubbath_hypfs(dfa, nqubit);
    DefectArray_reportSubbath_quads(dfa);
    DefectArray_reportSubbath_hypf_subs(dfa);
    
    ////////////////////////////////////////////
    // BathState, Disorder, Overhauser
    ////////////////////////////////////////////

    // Randomize bath states
    printMessage("\nRandomize bath states\n");
    BathArray_setBathStatesRandom(ba);
    BathArray_reportBath_states(ba);
    printf("\n");

    // Calculate bath disorders
    printMessage("\nCalculate bath disorders\n");
    BathArray_setBathDisorders(ba);
    BathArray_reportBath_disorders(ba);
    printf("\n");

    // Calculate overhauser fields
    bool isOvh = QubitArray_getOverhaus(qa);
    if (isOvh){
        for (int _iq=0; _iq<nqubit; _iq++){
            double overhaus = BathArray_getOverhaus(ba,_iq);
            QubitArray_setQubit_i_overhaus(qa,overhaus,_iq);
        }
    }
    QubitArray_reportQubit_overhaus(qa);


    // Randomize subbath states
    printMessage("\nRandomize subbath states\n");
    DefectArray_setSubbathStatesRandom(dfa,ba);

    // Update BathArray from DefectArray
    /**
     * What we do here:
     * 1. mainspinidx = self-index
     * 2. quad (zfs) = check paxis, replace the quad as zfs(axis)
     * 3. detuning = check paxis, add +=detuning(axis)
     */
    updateMainSpins_fromDefectArray(dfa,ba);

    // Update Disorder of mainspin, subspin from DefectArray
    updateDisorder_main_sub(dfa,ba);

    // Update Disorder of subspin from DefectArray
    updateDisorder_sub_sub(dfa);

    // Update Overhaus qubit from DefectArray
    updateOverhaus_qubit_sub(dfa,qa);

    // Report
    BathArray_reportBath_disorders(ba);
    DefectArray_reportSubbath_disorders(dfa);
    QubitArray_reportQubit_overhaus(qa);

    ////////////////////////////////////////////
    // Free
    ////////////////////////////////////////////
    QubitArray_freeAll(qa);
    BathArray_freeAll(ba);
    DefectArray_freeAll(dfa);
    return 0;
}