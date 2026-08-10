#ifndef __CCEX_PULSE_H_
#define __CCEX_PULSE_H_

#include "./../include/utilities.h"

/**
 * @struct Pulse
 * @brief  This structure contains the pulse-related parameters.
 */
typedef struct {

    int npulse; /**< The number of pulse (default : -)*/
    bool pulseiter; /**< The pulse iteration (default : false)*/
    /**<
     * About npulse (Details)
     * npulse = 0 : Ramsey (default)
     * npulse = 1 : HahnEcho (default)
     * npulse > 1 : CPMG (default)
     * npulse = 4 & pulsename = "WAHUHA" : WAHUHA (optional)
    */
    char pulsename[100]; /**< The pulsename (default : "None")*/
    /**<
     * About pulsename (Details)
     * pulsename = HahnEcho | Ramsey | CPMG | WAHUHA | Manual
     * above pulsename option without Maunal  doesn't need a sequence 
     * "Manual" option must need a sequence
     * If you want "Manual" option, then you have to choose the "gCCE" method!!!
    */

    /**
     *  @brief : Related with pulse sequence.
     */
    double** sequence; /**< pulse sequence (default : NULL)*/
    double*  pulse_fracs; 
    double*  pulse_angles; 
    /**<
     * puulse angles: Array of rotation angle of pulse.
     * Input file unit: Degree  ==> In code, radian.
     *   len(sequence) = len(pulse_angles)
     */
    char*    pulse_axes;
    /**<
     * pulse axes: Array of rotation axis of pulse.
     * you can choose the axis of rotation angle. 
     * Input file parameters: "X", "Y", "Z", "I"(Identity) 
     *   len(sequence) = len(pulse_axes)
     */
    int*     sequence_indices;
    /**<
     *  sequence_indices: The Index that have the same difference value
     *                        if there is no the same difference value 
     *                        then give the current index
     *                        ( This would reduce the calculational time cost ) 
     *   len(sequence) = len(sequence_indices)
     */
    double pulse_time;
    double detuning_factor;

    int      bnpulse; /**< The number of pulse (default : -)*/
    char     bpulse_defect[50];
    double   bspin;
    double   balphams;
    double   bbetams;
    MatrixXcd bpsia;
    MatrixXcd bpsib;
    double   bpulse_energy_shift;
    double*  bpulse_energy_shifts; /**< pump frequencies, rad*kHz (default : one entry) */
    int      n_bpulse_shift;       /**< how many pump frequencies (default : 1) */
    double** bsequence; /**< pulse sequence for specific bath spin (default : NULL)*/
    double*  bpulse_fracs; 
    double*  bpulse_angles; 
    char*    bpulse_axes;
    //int*     bsequence_indices;


    // total sequence 
    int      total_npulse;
    double** total_sequence; 
    double*  total_pulse_fracs;
    double*  total_pulse_angles; 
    char*    total_pulse_axes;
    double*  total_bpulse_angles; 
    char*    total_bpulse_axes;
    int*     total_sequence_indices;
    //int*     qpulseidx;
    //int*     bpulseidx;


    //
    /**<
     * About sequence (Details)
     * 
     * If you want to adjust the pulse sequence, that you have to set the pulse name : "Manual", 
     * and you have to give the "sequence" variable.
     *
     * -------------------------------------------------------------------
     * ex)          !                                            !(t_free)
     * Pulse Index  0    1    2    3    4    5    6    7    8    |
     * Pulse (8#)   |____|____|____|____|____|____|____|____|____|
     *              |tau1|tau2|           .......           |tau9|
     *              |seq1|seq2|           .......           |seq9|
     * it's okay for tau1 != tau2.
     * -------------------------------------------------------------------
     *
     * Sequence format:
     *  "pulsename" : "Manual",
     *  "sequence"  : [["frac_1",     "axis_1", float(angle_1) (deg)],
     *                 ["frac_2",     "axis_2", float(angle_2) (deg)],
     *                 ["frac_3",     "axis_3", float(angle_3) (deg)],
     *                                     ...
     *                 ["frac_n",          "I",                  0.0]]
     *
     *     where Summation_n(frac_n) must be float(1.0).
     *  And the rotation axis and angle of last sequence(sequence[-1]) is 0
     * 
     * ex1) Ramsey sequence 
     *    inputfile tag:
     *       - "npulse": 0
     *           or
     *       - "pulsename": "Manual",
     *       - "sequence" : [["1.0", "I", 0.0]]
     *
     *
     * ex2) Hahn-Echo sequence 
     *    inputfile tag:
     *       - "npulse" : 1
     *           or
     *       - "pulsename" : "Manual",
     *       - "sequence"  : [["0.5", "X", 180.0],
     *                        ["0.5", "I",   0.0]]
     *
     * ex3) CPMG sequence 
     *    inputfile tag:
     *       - "npulse" > 1
     *           or
     *       - "pulsename": "Manual",
     *       - "sequence" : [["1/{1tau}", "X", 180.0],
     *                       ["1/{2tau}", "X", 180.0],
     *                       ["1/{2tau}", "X", 180.0],
     *                                  ...
     *                       ["1/{2tau}", "X", 180.0],
     *                       ["1/{1tau}", "I",   0.0]]
     *
     * ex3) WAHUHA sequence 
     *    inputfile tag:
     *       - "npulse"   : 4
     *       - "pulsename": "WAHUHA",
     *           or
     *       - "npulse"   : 4
     *       - "pulsename": "Manual",
     *       - "sequence" : [["1/6", "X",  90.0],
     *                       ["1/6", "Y", 270.0],
     *                       ["2/6", "Y",  90.0],
     *                       ["1/6", "X", 270.0],
     *                       ["1/6", "I",   0.0]]
     *
     * | ------------------------------------------------------------- |
     * | ----------------- !        WARNING       ! ------------------ | 
     * | ------------------------------------------------------------- |
     * | - In CCEX CODE same variable name(sequence), but differnt!  - |
     *
     * In code, Array format:
     *  sequence is started with 0.0 and end with 1.0
     *  sequence[npulse+1][3] (defulat : NULL)
     *  sequence[ipulse][0] = Fraction of Previous pulse
     *  sequence[ipulse][1] = Fraction of current pulse
     *  sequence[ipulse][2] = Difference of Fractions "[i][1] - [i][0]"
     *
     * e.g. sequence parameter : Ramsey
     *      sequence parameter :
     *          sequence[0]      = {0.0 , 1.0 , 1.0}
     *          pulse_angles     = {0.0}
     *          pulse_axes       = {"I"}
     *          sequence_indices = {0}
     *
     * e.g. sequence parameter : Hahn-echo
     *          sequence[0]      = {"0.0", "0.5" , 0.5 }
     *          sequence[1]      = {"0.5", "1.0" , 0.5 }
     *          pulse_angles     = {180.0}
     *          pulse_axes       = {"X"}
     *          sequence_indices = {0, 0}
     *
     * e.g. sequence parameter : CPMG
     *         sequence[0]     = 0.0 , 0.25 , 0.25 
     *         sequence[1]     = 0.25 , 0.75 , 0.5 
     *         sequence[2]     = 0.75 , 1.0 , 0.25 
     *         pulse_angles    = {180.0,  180.0}
     *         pulse_axes      = {"X", "X"}
     *         sequence_indices = {0, 1, 1, 0}
     *
     * e.g. sequence parameter : WAHUHA
     *         sequence[0]     = 0.0 , 0.25 , 0.25 
     *         sequence[1]     = 0.25 , 0.75 , 0.5 
     *         sequence[2]     = 0.75 , 1.0 , 0.25 
     *         sequence[3]     = 0.75 , 1.0 , 0.25 
     *         sequence[4]     = 0.75 , 1.0 , 0.25 
     *         pulse_angles    = {180.0,  180.0}
     *         pulse_axes      = {"X", "X"}
     *         sequence_indices = {0, 1, 1, 0}
    */
} Pulse;

Pulse* Pulse_init();

void Pulse_report(Pulse* pulse);

//////////////////////
// Setting function //
//////////////////////
void Pulse_setNpulse(Pulse* pulse, int npulse);
void Pulse_setPulseiter(Pulse* pulse, bool pulseiter);
void Pulse_setPulsename(Pulse* pulse, char* pulsename);
void Pulse_setPulseTime(Pulse* pulse, double pulse_time);
void Pulse_setPulseDetuningFactor(Pulse* pulse, double detuning_factor);

void Pulse_setBNpulse(Pulse* pulse, int bnpulse);
void Pulse_setBPulse_Defect(Pulse* pulse, char* bpulse_defect);
void Pulse_setBPulse_bspin(Pulse* pulse, double bspin);
void Pulse_setBPulse_alphams(Pulse* pulse, double balphams);
void Pulse_setBPulse_betams(Pulse* pulse, double bbetams);
void Pulse_setBPulse_bpsia_fromMs(Pulse* pulse);
void Pulse_setBPulse_bpsib_fromMs(Pulse* pulse);
void Pulse_setBPulse_EnergyShift(Pulse* pulse, double bpulse_energy_shift);
void Pulse_setBPulse_EnergyShifts(Pulse* pulse, double* shifts, int n);

void Pulse_setTotNPulse(Pulse* pulse, int total_npulse);
void Pulse_setTotSequence(Pulse* pulse, double** total_sequence);
void Pulse_setTotPulseFracs(Pulse* pulse, double* total_pulse_fracs);
void Pulse_setTotPulseAngles(Pulse* pulse, double* total_pulse_angles);
void Pulse_setTotPulseAxes(Pulse* pulse, char* total_pulse_axes);
void Pulse_setTotBPulseAngles(Pulse* pulse, double* total_bpulse_angles);
void Pulse_setTotBPulseAxes(Pulse* pulse, char* total_bpulse_axes);

void assign_sequence_indices(Pulse* pulse);
void assign_total_sequence_indices(Pulse* pulse);

// Get function
int Pulse_getNpulse(Pulse* pulse);
bool Pulse_getPulseiter(Pulse* pulse);
char* Pulse_getPulsename(Pulse* pulse);
double** Pulse_getSequence(Pulse* pulse);
double* Pulse_getPulseFracs(Pulse* pulse);
double* Pulse_getPulseAngles(Pulse* pulse);
char* Pulse_getPulseAxes(Pulse* pulse);
int* Pulse_getSequenceIndices(Pulse* pulse);

double Pulse_getPulseTime(Pulse* pulse);
double Pulse_getPulseDetuningFactor(Pulse* pulse);

int Pulse_getBNpulse(Pulse* pulse);
char* Pulse_getBPulse_Defect(Pulse* pulse);
double Pulse_getBPulse_bspin(Pulse* pulse);
double Pulse_getBPulse_balphams(Pulse* pulse);
double Pulse_getBPulse_bbetams(Pulse* pulse);
MatrixXcd Pulse_getBPulse_bpsia(Pulse* pulse);
MatrixXcd Pulse_getBPulse_bpsib(Pulse* pulse);
double Pulse_getBPulse_EnergyShift(Pulse* pulse);
double** Pulse_getBSequence(Pulse* pulse);
double* Pulse_getBPulseFracs(Pulse* pulse);
double* Pulse_getBPulseAngles(Pulse* pulse);
char* Pulse_getBPulseAxes(Pulse* pulse);

int Pulse_getTotNpulse(Pulse* pulse);
double** Pulse_getTotSequence(Pulse* pulse);
double* Pulse_getTotPulseFracs(Pulse* pulse);
double* Pulse_getTotPulseAngles(Pulse* pulse);
double* Pulse_getTotBPulseAngles(Pulse* pulse);
char* Pulse_getTotPulseAxes(Pulse* pulse);
char* Pulse_getTotBPulseAxes(Pulse* pulse);
int* Pulse_getTotSequenceIndices(Pulse* pulse);

// allocation
void Pulse_allocSequence(Pulse* pulse);
void Pulse_allocAxes(Pulse* pulse);
void Pulse_allocAngles(Pulse* pulse);
void Pulse_allocSequenceIndices(Pulse* pulse);
void Pulse_allocBSequence(Pulse* pulse);
void Pulse_allocBAxes(Pulse* pulse);
void Pulse_allocBAngles(Pulse* pulse);
void Pulse_allocFracs(Pulse* pulse);
void Pulse_allocBFracs(Pulse* pulse);
void Pulse_allocTotalSequenceIndices(Pulse* pulse);

// Free function
void Pulse_freeSequence(Pulse* pulse);
void Pulse_freeFracs(Pulse* pulse);
void Pulse_freeAngles(Pulse* pulse);
void Pulse_freeAxes(Pulse* pulse);
void Pulse_freeSequenceIndices(Pulse* pulse);

void Pulse_freeBSequence(Pulse* pulse);
void Pulse_freeBFracs(Pulse* pulse);
void Pulse_freeBAngles(Pulse* pulse);
void Pulse_freeBAxes(Pulse* pulse);

void Pulse_freeTotSequence(Pulse* pulse);
void Pulse_freeTotFracs(Pulse* pulse);
void Pulse_freeTotAngles(Pulse* pulse);
void Pulse_freeTotAxes(Pulse* pulse);
void Pulse_freeTotBAngles(Pulse* pulse);
void Pulse_freeTotBAxes(Pulse* pulse);
void Pulse_freeTotSequenceIndices(Pulse* pulse);

void Pulse_freeAll(Pulse* pulse);
// generate pulse sequence
void generateSequenceRamsey(Pulse* pulse);
void generateSequenceHE(Pulse* sequence);
void generateSequenceWAHUHA(Pulse* sequence);
void generateSequenceCPMG(Pulse* pulse);
void allocateDefaultSequence(Pulse* pulse);



#endif // __CCEX_PULSE_H_
