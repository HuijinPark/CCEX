/**
 * @file cluster_pairs.h
 * @brief The order-2 cluster list of CCE/gCCE made without the hash, shared out between the ranks.
 */
#ifndef __CCEX_CLUSTER_PAIRS_H_
#define __CCEX_CLUSTER_PAIRS_H_

#include "clustering/cluster.h"

// CCE/gCCE at order 2 (no addsubclus, no verbosity): the cluster list that
// clusterizeHash makes, made without the hash and the dense nspin x nspin maps,
// with the pair search shared out between the ranks. Collective: every rank calls it.
// A38: with a cost (Cluster_costFunction), only this rank's pair rows are allocated (Cluster::pairsLocal).
void clusterizePairs(Cluster* cls, BathArray* ba, float rdip, float rdipcut,
                     const std::function<double(int,const int*)>& cost = nullptr);

#endif
