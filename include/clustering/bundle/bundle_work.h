#ifndef CCEX_ADAPTIVE_WORK_H
#define CCEX_ADAPTIVE_WORK_H
#include "bundle_partition.h"
#include "../../calculation/simulator.h"
#include "../../calculation/scaled_product.h"
#include <fstream>
#include <map>
#include <memory>
#include <string>

// Internal, explicit file protocol used by the refinement driver (src/clustering/algo_adaptive). Runs keep
// the entire bath, its state sampling and external mean fields intact.
struct AdaptiveWork {
    bool active=false;
    bool physics=false;
    std::string mode,response;
    std::vector<std::vector<int>> bundles,probes,extra_clusters;
    // mode "factors" ((internal notes)): the one state and time step whose bundle and
    // connected-pair values are wanted. -1 in every other mode.
    int state=-1,step=-1;
    // mode "except" ((internal notes)): the expansion terms kept out of the product P_K,
    // and their position in that list (the probe number of their curves).
    std::vector<std::vector<int>> except;
    std::map<std::vector<int>,int> exceptIndex;
};
const AdaptiveWork& adaptiveWork(const Config*);
void resetAdaptiveWork();
void adaptiveWritePlan(const Config*,const AdaptivePlan&,const std::vector<int>&,
                       const std::vector<AdaptiveEdge>&,int qubitDimension,int bundleCap);
void adaptiveProbe(QubitArray*,BathArray*,DefectArray*,Config*,Pulse*,Output*,int state);
// after the state loop of a probe step, on every rank: frees the node memory budget the probes shared (A1)
void adaptiveProbeFinish();
// The sites of one probe as a bath of their own, and the Hilbert dimension its calculation works in (the qubit's
// times the sites' for gCCE, the sites' alone for cCCE); throws when a site is outside the loaded bath or the
// dimension leaves the long long range. Shared by adaptiveProbe and the connection check (2026-10-08).
struct ProbeSubset {
    std::unique_ptr<BathArray,decltype(&BathArray_freeAll)> bath{nullptr,BathArray_freeAll};
    long long dim=1;
};
ProbeSubset adaptiveProbeSubset(const std::vector<int>& sites,BathArray*,DefectArray*,QubitArray*,bool gcce);
void adaptiveWriteCurve(const Config*,const MatrixXcd*,int state);
// mode "factors", in the state the request names, on every rank: the value at the requested step of every connected
// bundle pair and of the bundle of each such pair, one row per cluster in response.rank-R.factors.tsv
// (adaptive-regroup-stage3.md 1절). Calculated as a probe is (adaptiveProbeSubset + calCoherenceCce), under the
// same node memory budget; the caller frees the budget with adaptiveProbeFinish.
void adaptiveFactors(QubitArray*,BathArray*,DefectArray*,Config*,Pulse*,Output*,int state);
// mode "except" (adaptive-regroup-stage3b.md 1절). The position of a term in the request's except list, or -1.
int adaptiveExceptIndex(const int* sites,int n);
// This rank's file of the excepted terms' raw curves for a state, in the probe format (rank 0 adds L0 = 1).
std::unique_ptr<std::ofstream> adaptiveExceptCurves(const Config*,int state);
void adaptiveExceptValues(std::ostream&,int index,const MatrixXcd* values,int nstep);
// Every rank: the excepted terms this state met, summed over the ranks, have to be the whole list - once each.
void adaptiveExceptCheck(long long matched);
// Rank 0: P_K of a state, response.except.state-S.tsv (step, mantissa real and imaginary as %La, exponent).
void adaptiveWriteExcept(const Config*,const std::vector<ScaledComplex>& kept,int state);
void adaptiveWritePhysics(const Config*,BathArray*,QubitArray*,DefectArray*);
#endif
