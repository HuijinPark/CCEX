#include "../include/pulse.h"
#include "../include/memory.h"

bool verbosity = true;
int rank = 0;

int main(){

    Pulse* pulse = Pulse_init();
    Pulse_report(pulse);

    // Ramsey
    Pulse_setNpulse(pulse, 0);
    Pulse_setPulsename(pulse, "Ramsey");

    Pulse_allocSequence(pulse);
    Pulse_setSequence_fromName(pulse);
    Pulse_report(pulse);

    // HahnEcho
    Pulse_setNpulse(pulse, 1);
    Pulse_setPulsename(pulse, "HahnEcho");

    Pulse_allocSequence(pulse);
    Pulse_setSequence_fromName(pulse);
    Pulse_report(pulse);

    // CPMG
    Pulse_setNpulse(pulse, 10);
    Pulse_setPulsename(pulse, "CPMG");

    Pulse_allocSequence(pulse);
    Pulse_setSequence_fromName(pulse);
    Pulse_report(pulse);

    // Equal
    Pulse_setNpulse(pulse, 20);
    Pulse_setPulsename(pulse, "Equal");

    Pulse_allocSequence(pulse);
    Pulse_setSequence_fromName(pulse);
    Pulse_report(pulse);

    // Arbitrary
    double* seqinput = allocDouble1d(5);
    seqinput[0] = 0.1323;
    seqinput[1] = 0.3;
    seqinput[2] = 0.35;
    seqinput[3] = 0.7;
    seqinput[4] = 0.99;

    Pulse_setNpulse(pulse, 5);
    Pulse_setPulsename(pulse, "Arbitrary");
    Pulse_setSequence_fromInput(pulse, seqinput);
    freeDouble1d(seqinput);
    Pulse_report(pulse);
    
    Pulse_freeAll(pulse);
    return 0;
}