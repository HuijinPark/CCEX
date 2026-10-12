#ifndef __CCEX_CLUSTER_PAIRS_H_
#define __CCEX_CLUSTER_PAIRS_H_

#include "clustering/cluster.h"

// CCE/gCCE at order 2 (no addsubclus, no verbosity): the cluster list that
// clusterizeHash makes, made without the hash and the dense nspin x nspin maps,
// with the pair search shared out between the ranks. Collective: every rank calls it.
void clusterizePairs(Cluster* cls, BathArray* ba, float rdip, float rdipcut);

#endif
