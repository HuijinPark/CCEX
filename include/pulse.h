#ifndef __CCEX_PULSE_H_
#define __CCEX_PULSE_H_

#include "./../include/utilities.h"

/**
 * @struct Pulse
 * @brief  This structure contains the pulse-related parameters.
 */
typedef struct {

    int npulse; /**< The number of pulse (default : 0)*/
    bool pulseiter; /**< The pulse iteration (default : false)*/
    /**<
     * About npulse (Details)
     * npulse = 0 : Ramsey (default)
     * npulse = 1 : HahnEcho (default)
     * npulse > 1 : CPMG (default)
    */
    char pulsename[100]; /**< The pulsename (default : "None")*/
    /**<
     * About pulsename (Details)
     * pulsename = HahnEcho | Ramsey | CPMG
     * above pulsename option doesn't need a sequence 
    */
    double** sequence; /**< pulse sequence (default : NULL)*/
    /**<
     * About sequence (Details)
     * By using this parameter, we can adjust/delay pulse timing. (No unit)
     * 
     * Input file format:
     *  sequence = { t1, t2, t3, ... , tn } (tn < 1.0, n = npulse) (only pulse timing)
     *  ti is the time fraction of i-th pulse
     * 
     * Array format:
     *  sequence is started with 0.0 and end with 1.0
     *  sequence[npulse+1][4] (defulat : NULL)
     *  sequence[ipulse][0] = Fraction of Previous pulse
     *  sequence[ipulse][1] = Fraction of current pulse
     *  sequence[ipulse][2] = Difference of Fractions "[i][1] - [i][0]"
     *  sequence[ipulse][3] = The Index that have the same difference value
     *                        if there is no the same difference value 
     *                        then give the current index
     *                        ( This would reduce the calculational time cost )
     * e.g. Input file format : 
     *          npulse = 0 (Ramasey)
     *      sequence parameter :
     *          sequence[0] = 0.0 , 1.0 , 1.0 , 0
     * e.g. Input file format : 
     *          sequence = { 0.5 } (HahnEcho)
     *      sequence parameter :
     *          sequence[0] = 0.0 , 0.5 , 0.5 , 0
     *          sequence[1] = 0.5 , 1.0 , 0.5 , 0
     * e.g. Input file format :
     *         sequence = { 0.25, 0.75 } (CPMG)
     *      sequence parameter :
     *         sequence[0] = 0.0 , 0.25 , 0.25 , 0
     *         sequence[1] = 0.25 , 0.75 , 0.5 , 1
     *         sequence[2] = 0.75 , 1.0 , 0.25 , 0
    */
} Pulse;

Pulse* Pulse_init();

void Pulse_report(Pulse* pulse);

void Pulse_setNpulse(Pulse* pulse, int npulse);
void Pulse_setPulsename(Pulse* pulse, char* pulsename);
void Pulse_setSequence_fromName(Pulse* pulse);
void Pulse_setSequence_fromInput(Pulse* pulse, double* sequenceinput);
void Pulse_setPulseiter(Pulse* pulse, bool pulseiter);

int Pulse_getNpulse(Pulse* pulse);
char* Pulse_getPulsename(Pulse* pulse);
double** Pulse_getSequence(Pulse* pulse);
bool Pulse_getPulseiter(Pulse* pulse);

// allocation
void Pulse_allocSequence(Pulse* pulse);
void Pulse_freeSequence(Pulse* pulse);
void Pulse_freeAll(Pulse* pulse);

// generate pulse sequence
void generateSequenceRamsey(double*** sequence);
void generateSequenceHE(double*** sequence);
void generateSequenceCPMG(double*** sequence, int npulse);
void generateSequenceEqual(double*** sequence, int npulse);



#endif // __CCEX_PULSE_H_
