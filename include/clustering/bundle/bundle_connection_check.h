#ifndef CCEX_ADAPTIVE_CONNECTION_CHECK_H
#define CCEX_ADAPTIVE_CONNECTION_CHECK_H
#include "bundle_partition.h"
#include "../../calculation/simulator.h"

// Root prepares a conservative candidate. The initial installed plan is the
// complete graph; only a successful, bounded local audit can replace it.
void prepareAdaptiveConnectionCheck(const Config*,const AdaptivePlan&,
    const std::vector<int>&,const std::vector<AdaptiveEdge>&,int qubitDimension,int bundleCap);
void checkAdaptiveConnections(QubitArray*,BathArray*,DefectArray*,Config*,Pulse*,Cluster*,Output*);
bool restoreAdaptiveCheckState(BathArray*,DefectArray*,int state);
bool adaptiveCheckState(QubitArray*,BathArray*,DefectArray*,Config*,Pulse*,Output*,int state);
void clearAdaptiveConnectionCheck();
#endif
