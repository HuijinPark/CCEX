#include "../../include/calculation/simulator.h"
#include "../../include/calculation/gpu/eigen_refinement.h"

// A18: gCCE evolution=matrix with evolution_update = recursive. calPropagatorGcce (propagator=eigen) with the phases of
// each distinct interval carried over from the previous step: computed (refinedPhases) at steps 0, 16, 32, ... and in
// between the previous step's times exp(-i lambda f deltat), as propagateStateGcce does for evolution=vector (A16 b).
// The calls come for istep = 0, 1, 2, ... with tfree = istep deltat (accumulated). Everything else - the free
// propagators M diag(phase) M^dag and their product with the pulses - is the same product in the same order.
// A file of its own so that the exact calPropagatorGcce in simulator_gcce.cpp compiles as before.
MatrixXcd calPropagatorGcceRecursive(Pulse* pulse, double tfree, GcceWork* w, int istep, double deltat){

    double** sequence        = Pulse_getSequence(pulse);
    int*    sequence_indices = Pulse_getSequenceIndices(pulse);
    int     npulse           = Pulse_getNpulse(pulse);

    const int gccePhaseAnchor = 16;
    MatrixXcd Utotal;
    if ((int)w->Ufrees.size() < npulse+1){ w->Ufrees.resize(npulse+1); }
    if ((int)w->phases.size() < npulse+1){ w->phases.resize(npulse+1); w->phaseSteps.resize(npulse+1); }

    for (int ipulse=0; ipulse<npulse+1; ipulse++){
        int sameTauIndex = sequence_indices[ipulse];
        if (sameTauIndex == ipulse){
            if (istep % gccePhaseAnchor == 0){
                refinedPhases(w->evals,w->evalCorrections,tfree*sequence[ipulse][2],w->phases[ipulse]);
                if (istep == 0)
                    refinedPhases(w->evals,w->evalCorrections,deltat*sequence[ipulse][2],w->phaseSteps[ipulse]);
            } else {
                w->phases[ipulse].array() *= w->phaseSteps[ipulse].array();
            }
            w->scratch.noalias()        = w->M * w->phases[ipulse].asDiagonal();
            w->Ufrees[ipulse].noalias() = w->scratch * w->Madj;
        }
        else{
            w->Ufrees[ipulse] = w->Ufrees[sameTauIndex];
        }
        const MatrixXcd& Ufree = w->Ufrees[ipulse];

        if (ipulse==0 && ipulse != npulse){
            Utotal = w->UpulsesExpanded[ipulse] * Ufree;
        }
        else if (ipulse==0 && ipulse == npulse){
            Utotal = Ufree;
        }
        else if (ipulse==npulse){
            Utotal = Ufree * Utotal;
        }
        else {
            Utotal = w->UpulsesExpanded[ipulse] * Ufree * Utotal;
        }
    }
    return Utotal;
}
