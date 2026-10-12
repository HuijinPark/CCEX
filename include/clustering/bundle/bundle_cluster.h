#ifndef CCEX_CLUSTER_ADAPTIVE_H
#define CCEX_CLUSTER_ADAPTIVE_H

#include "../cluster.h"
#include "../../physical/defect.h"
#include "bundle_partition.h"

void validateAdaptiveConfig(const Config*,const Cluster*);
void clusterizeAdaptive(Cluster*,BathArray*,QubitArray*,DefectArray*,const Config*);
void installAdaptivePlan(Cluster*,BathArray*,const Config*,const AdaptivePlan& rootPlan);
const std::vector<std::vector<int>>& adaptiveInstalledBundles();
// The connected bundle pairs of the installed plan (indices into adaptiveInstalledBundles), which the factors work mode
// calculates ((internal notes) 5절 2).
const std::vector<std::pair<int,int>>& adaptiveInstalledEdges();
void clearAdaptiveInstalledBundles();

#endif
