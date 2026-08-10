#include "../include/pulse.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include <iostream>
#include <fstream>
#include <string>

#define MAX_SEQUENCE 100

/* Low-level functions ------------------------------------------------------*/

Pulse* Pulse_init(){
    Pulse* pulse = (Pulse*)allocArray1d(1, sizeof(Pulse));
    pulse->npulse       = 0;
    pulse->pulseiter    = false;

    // Qubit spin pulse operator 
    pulse->pulse_fracs      = NULL;
    pulse->pulse_axes       = NULL;
    pulse->pulse_angles     = NULL;
    pulse->sequence         = NULL;
    pulse->sequence_indices = NULL; // Used in conventional CCE method!!

    // For "time-dependent Qubit spins' pulse operator" 
    pulse->pulsename[0] = '\0';
    pulse->pulse_time   = 0.0;
    pulse->detuning_factor = 1.0;

    // Bath spin pulse operator
    pulse->bnpulse             =    0;
    pulse->bpulse_defect[0]    = '\0';
    pulse->bspin               =  0.0;
    pulse->balphams            =  0.0;
    pulse->bbetams             =  0.0;
    pulse->bpulse_energy_shift =  0.0;
    pulse->bpulse_energy_shifts = NULL;
    pulse->n_bpulse_shift      =    0;
    // MatrixXcd bpsia; do not need to initialize (like qubit)
    // MatrixXcd bpsib; do not need to initialize (like qubit)

    pulse->bpulse_fracs        = NULL;
    pulse->bpulse_axes         = NULL;
    pulse->bpulse_angles       = NULL;
    pulse->bsequence           = NULL;
    //pulse->bsequence_indices   = NULL;

    // Total sequence (qubit's pulse + bath's pulse)
    pulse->total_npulse           = 0;
    pulse->total_pulse_fracs      = NULL;
    pulse->total_sequence         = NULL;
    pulse->total_sequence_indices = NULL;
    pulse->total_pulse_axes       = NULL;
    pulse->total_pulse_angles     = NULL;
    pulse->total_bpulse_axes      = NULL;
    pulse->total_bpulse_angles    = NULL;
    //pulse->qpulseidx = NULL;
    //pulse->bpulseidx = NULL;

    return pulse;
}

//////////////////////////////
//     Setting Function     // 
//////////////////////////////
// 1) For Qubit spin pulse operator !!
void Pulse_setNpulse(Pulse* pulse, int npulse){
    pulse->npulse = npulse;
}

void Pulse_setPulseiter(Pulse* pulse, bool pulseiter){
    pulse->pulseiter = pulseiter;
}

void Pulse_setPulsename(Pulse* pulse, char* pulsename){
    strcpy(pulse->pulsename, pulsename);
}

void Pulse_setPulseTime(Pulse* pulse, double pulse_time){
    pulse->pulse_time = pulse_time;
}

void Pulse_setPulseDetuningFactor(Pulse* pulse, double detuning_factor){
    pulse->detuning_factor = detuning_factor;
}

// 2) Bath spin pulse operator 
void Pulse_setBNpulse(Pulse* pulse, int bnpulse){
    pulse->bnpulse = bnpulse;
}

void Pulse_setBPulse_Defect(Pulse* pulse, char* bpulse_defect){
    strcpy(pulse->bpulse_defect, bpulse_defect);
}

void Pulse_setBPulse_bspin(Pulse* pulse, double bspin){
    pulse->bspin = bspin;
}

void Pulse_setBPulse_alphams(Pulse* pulse, double balphams){
    pulse->balphams = balphams;
}

void Pulse_setBPulse_betams(Pulse* pulse, double bbetams){
    pulse->bbetams = bbetams;
}

void Pulse_setBPulse_bpsia_fromMs(Pulse* pulse){
    double S = Pulse_getBPulse_bspin(pulse);
    const double balphams = pulse->balphams;
    MatrixXcd bpsia = getSpinor(S, balphams);
    pulse->bpsia = bpsia;
}

void Pulse_setBPulse_bpsib_fromMs(Pulse* pulse){
    double S = Pulse_getBPulse_bspin(pulse);
    const double bbetams = pulse->bbetams;
    MatrixXcd bpsib = getSpinor(S, bbetams);
    pulse->bpsib = bpsib;
}

void Pulse_setBPulse_EnergyShift(Pulse* pulse, double bpulse_energy_shift){
    pulse->bpulse_energy_shift = MHZ_TO_RADKHZ(bpulse_energy_shift);
    Pulse_setBPulse_EnergyShifts(pulse, &bpulse_energy_shift, 1);
}

void Pulse_setBPulse_EnergyShifts(Pulse* pulse, double* shifts, int n){
    if (pulse->bpulse_energy_shifts != NULL){
        freeDouble1d(&(pulse->bpulse_energy_shifts));
    }
    pulse->n_bpulse_shift = n;
    pulse->bpulse_energy_shifts = allocDouble1d(n > 0 ? n : 1);
    for (int i=0; i<n; i++){
        pulse->bpulse_energy_shifts[i] = MHZ_TO_RADKHZ(shifts[i]);
    }
    if (n > 0){ pulse->bpulse_energy_shift = pulse->bpulse_energy_shifts[0]; }
}

// 3) Total pulse sequence (Qubit + Bath)
void Pulse_setTotNPulse(Pulse* pulse, int total_npulse){
    pulse->total_npulse = total_npulse;
}

void Pulse_setTotSequence(Pulse* pulse, double** total_sequence){
    pulse->total_sequence = total_sequence;
}

void Pulse_setTotPulseFracs(Pulse* pulse, double* total_pulse_fracs){
    pulse->total_pulse_fracs = total_pulse_fracs;
}

void Pulse_setTotPulseAngles(Pulse* pulse, double* total_pulse_angles){
    pulse->total_pulse_angles = total_pulse_angles;
}

void Pulse_setTotPulseAxes(Pulse* pulse, char* total_pulse_axes){
    pulse->total_pulse_axes = total_pulse_axes;
}

void Pulse_setTotBPulseAngles(Pulse* pulse, double* total_bpulse_angles){
    pulse->total_bpulse_angles = total_bpulse_angles;
}

void Pulse_setTotBPulseAxes(Pulse* pulse, char* total_bpulse_axes){
    pulse->total_bpulse_axes = total_bpulse_axes;
}

void assign_sequence_indices(Pulse* pulse) {
    for (int i = 0; i < (pulse->npulse)+1; ++i) {
        int found = 0;
        for (int j = 0; j < i; ++j) {
            if ((pulse->sequence)[i][2] == (pulse->sequence)[j][2]) {
                pulse->sequence_indices[i] = pulse->sequence_indices[j];
                found = 1;
                break;
            }
        }
        if (!found) {
            (pulse->sequence_indices)[i] = i;
        }
    }
}

void assign_total_sequence_indices(Pulse* pulse) {
    for (int i = 0; i < (pulse->total_npulse)+1; ++i) {
        int found = 0;
        for (int j = 0; j < i; ++j) {
            if ((pulse->total_sequence)[i][2] == (pulse->total_sequence)[j][2]) {
                pulse->total_sequence_indices[i] = pulse->total_sequence_indices[j];
                found = 1;
                break;
            }
        }
        if (!found) {
            (pulse->total_sequence_indices)[i] = i;
        }
    }
}

//////////////////////////////
//     Getting Function     // 
//////////////////////////////
// 1) Qubit pulse
int Pulse_getNpulse(Pulse* pulse){
    return pulse->npulse;
}

bool Pulse_getPulseiter(Pulse* pulse){
    return pulse->pulseiter;
}

char* Pulse_getPulsename(Pulse* pulse){
    return pulse->pulsename;
}

double Pulse_getPulseTime(Pulse* pulse){
    return pulse->pulse_time;
}

double Pulse_getPulseDetuningFactor(Pulse* pulse){
    return pulse->detuning_factor;
}

double** Pulse_getSequence(Pulse* pulse){
    return pulse->sequence;
}

double* Pulse_getPulseFracs(Pulse* pulse){
    return pulse->pulse_angles;
}

char* Pulse_getPulseAxes(Pulse* pulse){
    return pulse->pulse_axes;
}

double* Pulse_getPulseAngles(Pulse* pulse){
    return pulse->pulse_angles;
}

int* Pulse_getSequenceIndices(Pulse* pulse){
    return pulse->sequence_indices;
}

//////////////////////////
// 2) Bath pulse
int Pulse_getBNpulse(Pulse* pulse){
    return pulse->bnpulse;
}

char* Pulse_getBPulse_Defect(Pulse* pulse){
    return pulse->bpulse_defect;
}

double Pulse_getBPulse_bspin(Pulse* pulse){
    return pulse->bspin;
}

double Pulse_getBPulse_balphams(Pulse* pulse){
    return pulse->balphams;
}

double Pulse_getBPulse_bbetams(Pulse* pulse){
    return pulse->bbetams;
}

MatrixXcd Pulse_getBPulse_bpsia(Pulse* pulse){
    return pulse->bpsia;
}

MatrixXcd Pulse_getBPulse_bpsib(Pulse* pulse){
    return pulse->bpsib;
}

double Pulse_getBPulse_EnergyShift(Pulse* pulse){
    return pulse->bpulse_energy_shift;
}


double** Pulse_getBSequence(Pulse* pulse){
    return pulse->bsequence;
}

double* Pulse_getBPulseFracs(Pulse* pulse){
    return pulse->bpulse_angles;
}

char* Pulse_getBPulseAxes(Pulse* pulse){
    return pulse->bpulse_axes;
}

double* Pulse_getBPulseAngles(Pulse* pulse){
    return pulse->bpulse_angles;
}

///////////////////////////////////
// 3) Total pulse 
double** Pulse_getTotSequence(Pulse* pulse){
    return pulse->total_sequence;
}

double* Pulse_getTotPulseFracs(Pulse* pulse){
    return pulse->total_pulse_fracs;
}

char* Pulse_getTotPulseAxes(Pulse* pulse){
    return pulse->total_pulse_axes;
}

double* Pulse_getTotPulseAngles(Pulse* pulse){
    return pulse->total_pulse_angles;
}

char* Pulse_getTotBPulseAxes(Pulse* pulse){
    return pulse->total_bpulse_axes;
}

double* Pulse_getTotBPulseAngles(Pulse* pulse){
    return pulse->total_bpulse_angles;
}

int Pulse_getTotNpulse(Pulse* pulse){
    return pulse->total_npulse;
}

int* Pulse_getTotSequenceIndices(Pulse* pulse){
    return pulse->total_sequence_indices;
}

//////////////////////////
//      Allocation      //
//////////////////////////
void Pulse_allocSequence(Pulse* pulse){
    pulse->sequence = allocDouble2d((pulse->npulse)+1, 3);
}

void Pulse_allocAxes(Pulse* pulse){
    pulse->pulse_axes = allocChar1d((pulse->npulse)+1);
}

void Pulse_allocFracs(Pulse* pulse){
    pulse->pulse_fracs = allocDouble1d((pulse->npulse)+1);
}

void Pulse_allocAngles(Pulse* pulse){
    pulse->pulse_angles = allocDouble1d((pulse->npulse)+1);
}

void Pulse_allocSequenceIndices(Pulse* pulse){
    pulse->sequence_indices = allocInt1d((pulse->npulse)+1);
}

void Pulse_allocBSequence(Pulse* pulse){
    pulse->bsequence = allocDouble2d((pulse->bnpulse)+1, 3);
}

void Pulse_allocBFracs(Pulse* pulse){
    pulse->bpulse_fracs = allocDouble1d((pulse->bnpulse)+1);
}

void Pulse_allocBAxes(Pulse* pulse){
    pulse->bpulse_axes = allocChar1d((pulse->bnpulse)+1);
}

void Pulse_allocBAngles(Pulse* pulse){
    pulse->bpulse_angles = allocDouble1d((pulse->bnpulse)+1);
}

//void Pulse_allocBSequenceIndices(Pulse* pulse){
//    pulse->bsequence_indices = allocInt1d((pulse->bnpulse)+1);
//}

void Pulse_allocTotalSequenceIndices(Pulse* pulse){
    pulse->total_sequence_indices = allocInt1d((pulse->total_npulse)+1);
}
// 3) total pulse-related varaibles -> "build_total_sequence function" 

//////////////////////////
//         Free         //
//////////////////////////
void Pulse_freeSequence(Pulse* pulse){
    int npulse = Pulse_getNpulse(pulse);
    freeDouble2d(&(pulse->sequence), npulse+1);
}

void Pulse_freeFracs(Pulse* pulse){
    freeDouble1d(&(pulse->pulse_fracs));
}

void Pulse_freeAngles(Pulse* pulse){
    freeDouble1d(&(pulse->pulse_angles));
}

void Pulse_freeAxes(Pulse* pulse){
    freeChar1d(&(pulse->pulse_axes));
}

void Pulse_freeSequenceIndices(Pulse* pulse){
    int npulse = Pulse_getNpulse(pulse);
    freeInt1d(&(pulse->sequence_indices));
}

void Pulse_freeBSequence(Pulse* pulse){
    int bnpulse = Pulse_getBNpulse(pulse);
    freeDouble2d(&(pulse->bsequence), bnpulse+1);
}

void Pulse_freeBFracs(Pulse* pulse){
    freeDouble1d(&(pulse->bpulse_fracs));
}

void Pulse_freeBAngles(Pulse* pulse){
    freeDouble1d(&(pulse->bpulse_angles));
}

void Pulse_freeBAxes(Pulse* pulse){
    freeChar1d(&(pulse->bpulse_axes));
}

////////////
/// 3) Total pulse sequence
void Pulse_freeTotSequence(Pulse* pulse){
    int total_npulse = Pulse_getTotNpulse(pulse);
    freeDouble2d(&(pulse->total_sequence), total_npulse+1);
}

void Pulse_freeTotFracs(Pulse* pulse){
    freeDouble1d(&(pulse->total_pulse_fracs));
}

void Pulse_freeTotAngles(Pulse* pulse){
    freeDouble1d(&(pulse->total_pulse_angles));
}

void Pulse_freeTotAxes(Pulse* pulse){
    freeChar1d(&(pulse->total_pulse_axes));
}

void Pulse_freeTotBAngles(Pulse* pulse){
    freeDouble1d(&(pulse->total_bpulse_angles));
}

void Pulse_freeTotBAxes(Pulse* pulse){
    freeChar1d(&(pulse->total_bpulse_axes));
}

void Pulse_freeTotSequenceIndices(Pulse* pulse){
    int total_npulse = Pulse_getTotNpulse(pulse);
    freeInt1d(&(pulse->total_sequence_indices));
}

void Pulse_freeAll(Pulse* pulse){
    //Pulse_freeTauFractions(pulse);
    Pulse_freeSequence(pulse);
    Pulse_freeFracs(pulse);
    Pulse_freeAngles(pulse);
    Pulse_freeAxes(pulse);
    Pulse_freeSequenceIndices(pulse);

    Pulse_freeBSequence(pulse);
    Pulse_freeBFracs(pulse);
    Pulse_freeBAngles(pulse);
    Pulse_freeBAxes(pulse);

    Pulse_freeTotSequence(pulse);
    Pulse_freeTotFracs(pulse);
    Pulse_freeTotAngles(pulse);
    Pulse_freeTotAxes(pulse);
    Pulse_freeTotBAngles(pulse);
    Pulse_freeTotBAxes(pulse);
    Pulse_freeTotSequenceIndices(pulse);

    freeArray1d((void**)&pulse);
}






////////////////////////////////////////////
// Pulse sequence 
////////////////////////////////////////////
void allocateDefaultSequence(Pulse* pulse) {
    if (strcasecmp(pulse->pulsename, "WAHUHA") == 0) {
        Pulse_setPulsename(pulse,"WAHUHA");
        generateSequenceWAHUHA(pulse);
    } else{
        if (pulse->npulse == 0 || strcasecmp(pulse->pulsename, "Ramsey") == 0) {
            Pulse_setPulsename(pulse,"Ramsey");
            generateSequenceRamsey(pulse);
        }
        else if (pulse->npulse == 1 || strcasecmp(pulse->pulsename, "HahnEcho") == 0) {
            Pulse_setPulsename(pulse,"HahnEcho");
            generateSequenceHE(pulse);
        }
        else if (pulse->npulse > 1 || strcasecmp(pulse->pulsename, "CPMG") == 0) {
            Pulse_setPulsename(pulse,"CPMG");
            generateSequenceCPMG(pulse);
        }
    }
}

void generateSequenceRamsey(Pulse* pulse){
    
    (pulse->sequence)[0][0]      = 0.0; 
    (pulse->sequence)[0][1]      = 1.0; 
    (pulse->sequence)[0][2]      = 1.0; 
    (pulse->pulse_fracs)[0]      = 0.0; 
    (pulse->pulse_axes)[0]       = 'I';
    (pulse->pulse_angles)[0]     = 0.0; 
    (pulse->sequence_indices)[0] =   0;

}

void generateSequenceHE(Pulse* pulse){
    (pulse->sequence)[0][0]      =   0.0; 
    (pulse->sequence)[0][1]      =   0.5; 
    (pulse->sequence)[0][2]      =   0.5; 
    (pulse->pulse_fracs)[0]      =   0.5; 
    (pulse->pulse_axes)[0]       =   'X';
    (pulse->pulse_angles)[0]     = 180.0; 
    (pulse->sequence_indices)[0] =     0;

    (pulse->sequence)[1][0]      = 0.5; 
    (pulse->sequence)[1][1]      = 1.0; 
    (pulse->sequence)[1][2]      = 0.5; 
    (pulse->pulse_fracs)[1]      = 0.5; 
    (pulse->pulse_axes)[1]       = 'I';
    (pulse->pulse_angles)[1]     =   0; 
    (pulse->sequence_indices)[1] =   0;
}

void generateSequenceWAHUHA(Pulse* pulse){

    (pulse->sequence)[0][0]      =       0.0; 
    (pulse->sequence)[0][1]      =   1.0/6.0; 
    (pulse->sequence)[0][2]      =   1.0/6.0; 
    (pulse->pulse_fracs)[0]      =   1.0/6.0; 
    (pulse->pulse_axes)[0]       =       'X';
    (pulse->pulse_angles)[0]     =      90.0; 
    (pulse->sequence_indices)[0] =         0;

    (pulse->sequence)[1][0]      =   1.0/6.0; 
    (pulse->sequence)[1][1]      =   2.0/6.0; 
    (pulse->sequence)[1][2]      =   1.0/6.0; 
    (pulse->pulse_fracs)[1]      =   1.0/6.0; 
    (pulse->pulse_axes)[1]       =       'Y';
    (pulse->pulse_angles)[1]     =     270.0; 
    (pulse->sequence_indices)[1] =         0;

    (pulse->sequence)[2][0]      =   2.0/6.0; 
    (pulse->sequence)[2][1]      =   4.0/6.0; 
    (pulse->sequence)[2][2]      =   2.0/6.0; 
    (pulse->pulse_fracs)[2]      =   2.0/6.0; 
    (pulse->pulse_axes)[2]       =       'Y';
    (pulse->pulse_angles)[2]     =      90.0; 
    (pulse->sequence_indices)[2] =         2;

    (pulse->sequence)[3][0]      =   4.0/6.0; 
    (pulse->sequence)[3][1]      =   5.0/6.0; 
    (pulse->sequence)[3][2]      =   1.0/6.0; 
    (pulse->pulse_fracs)[3]      =   1.0/6.0; 
    (pulse->pulse_axes)[3]       =       'X';
    (pulse->pulse_angles)[3]     =     270.0; 
    (pulse->sequence_indices)[3] =         0;

    (pulse->sequence)[4][0]      =   5.0/6.0; 
    (pulse->sequence)[4][1]      =       1.0; 
    (pulse->sequence)[4][2]      =   1.0/6.0; 
    (pulse->pulse_fracs)[4]      =   1.0/6.0; 
    (pulse->pulse_axes)[4]       =       'I';
    (pulse->pulse_angles)[4]     =       0.0; 
    (pulse->sequence_indices)[4] =         0;
}

void generateSequenceCPMG(Pulse* pulse){

    // n : pulse
    // c : sequence i
    // sequence (2c-1)/2n * T

    int npulse = pulse -> npulse;
    double start = 0.0;
    double end = 1.0;

    for (int i=0; i<(npulse+1); i++){

        if (i==0 && i!=npulse){
            start = 0.0;
            end = double(2*(i+1)-1)/double(2*npulse); 
        }
        else if (i!=0 && i==npulse){
            start = double(2*i-1)/double(2*npulse);
            end = 1.0;
        }
        else if (i==0 && i==npulse){
            start = 0.0;
            end = 1.0;
        }
        else{
            start = double(2*i-1)/double(2*npulse);
            end = double(2*(i+1)-1)/double(2*npulse); 
        }

        (pulse->sequence)[i][0]   =       start; 
        (pulse->sequence)[i][1]   =         end;
        (pulse->sequence)[i][2]   = end - start; 
        (pulse->pulse_fracs)[i]   = end - start; 
        (pulse->pulse_axes)[i]    =         'X';
        (pulse->pulse_angles)[i]  =       180.0; 

        // Find 3-rd value
        if (i==0 || i==npulse){
            pulse->sequence_indices[i] = 0; 
        }
        else{
            pulse->sequence_indices[i] = 1; 
        }
    }
}


//void generateSequenceEqual(double*** sequence, int npulse){
void generateSequenceEqual(Pulse* pulse){

    // n : pulse
    // c : sequence i
    // sequence c/n * T
    int npulse = pulse->npulse;
    double start = 0.0;
    double end = 1.0;
    double diff = 1.0/double(npulse+1);

    for (int i=0; i<(npulse+1); i++){

        (pulse->sequence)[i][0] = diff * i; 
        (pulse->sequence)[i][1] = diff * (i+1);
        (pulse->sequence)[i][2] = diff;

        (pulse->sequence_indices)[i]= 0; 

    }
}

/////////////////////////////////////////////////////// 
void Pulse_report(Pulse* pulse){

    int      npulse    = Pulse_getNpulse(pulse);
    char*    pulsename = Pulse_getPulsename(pulse);
    bool     pulseiter = Pulse_getPulseiter(pulse);
    double** sequence  = Pulse_getSequence(pulse);

    printTitle("Structure Pulse");

    printStructElementChar("pulsename", Pulse_getPulsename(pulse));
    printStructElementInt("npulse", Pulse_getNpulse(pulse));
    printStructElementDouble("detuning_factor", Pulse_getPulseDetuningFactor(pulse));
    printStructElementDouble("pulse_time", Pulse_getPulseTime(pulse));
    printStructElementChar("bpulse_defect", pulse->bpulse_defect);
    printStructElementInt("bnpulse", Pulse_getBNpulse(pulse));
    printStructElementInt("Total npulse", Pulse_getTotNpulse(pulse));
    printStructElementDouble("detuning_factor", Pulse_getPulseDetuningFactor(pulse));
    printStructElementDouble("pulse_time", Pulse_getPulseTime(pulse));

    printStructElementBool("pulseiter", Pulse_getPulseiter(pulse));
    printf("%27s * true  - pulse is applied to each qubit\n"," ");
    printf("%27s * false - pulse is applied to qubit-array\n"," ");
    
    if (sequence != NULL){
        printSubTitle("Sequence");
        printLine();
        //printf("%6s               %8s %8s %8s %8s \n"," ","start","end","diff\n");
        //for (int i=0; i<npulse+1; i++){
        //    printf("%6s sequence[%d] : %8.3lf %8.3lf %8.3lf\n"
        //    ," ",i,sequence[i][0],sequence[i][1],sequence[i][2]);
        //}
        //printf("%9s   %9s   %9s   %9s   %9s   %9s\n", "Start", "End", "Diff", "Axis", "Angle", "Index");
        printf("     [Qubit Spin Pulse Sequence]\n");
        printf("%9s   %9s   %9s   %9s   %9s   \n", "Start", "End", "Diff", "Q-Axis", "Q-Angle");
        for (int i=0; i<(pulse->npulse)+1; i++){
            printf("%9.3f   ", pulse->sequence[i][0]);
            printf("%9.3f   ", pulse->sequence[i][1]);
            printf("%9.3f   ", pulse->sequence[i][2]);
            printf("%9c   "  , pulse->pulse_axes[i]);
            printf("%9.2f   \n", pulse->pulse_angles[i]);
        }
        printLine();
        printf("     [Bath Spin Pulse Sequence]\n");
        printf("%9s   %9s   %9s   %9s   %9s   \n", "Start", "End", "Diff", "B-Axis", "B-Angle");
        for (int i=0; i<(pulse->bnpulse)+1; i++){
            printf("%9.3f   ", pulse->bsequence[i][0]);
            printf("%9.3f   ", pulse->bsequence[i][1]);
            printf("%9.3f   ", pulse->bsequence[i][2]);
            printf("%9c   "  , pulse->bpulse_axes[i]);
            printf("%9.2f  \n", pulse->bpulse_angles[i]);
            //printf("%9d\n", pulse->bsequence_indices[i]);
        }
        printLine();
        printf("     [Qubit & Bath Spin Pulse Sequence]\n");
        printf("%9s   %9s   %9s   %9s   %9s   %9s   %9s   %9s\n", 
                "Start", "End", "Diff", "Q-Axis", "Q-Angle", "B-Axis", "B-Angle", "Indices");
        for (int i=0; i<(pulse->total_npulse)+1; i++){
            printf("%9.3f   ", pulse->total_sequence[i][0]);
            printf("%9.3f   ", pulse->total_sequence[i][1]);
            printf("%9.3f   ", pulse->total_sequence[i][2]);
            printf("%9c   "  , pulse->total_pulse_axes[i]);
            printf("%9.2f   ", pulse->total_pulse_angles[i]);
            printf("%9c   "  , pulse->total_bpulse_axes[i]);
            printf("%9.2f   ", pulse->total_bpulse_angles[i]);
            printf("%9d\n", pulse->total_sequence_indices[i]);
        }
        printLine();
    }
}
