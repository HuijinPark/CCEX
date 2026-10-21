/**
 * @file pcce_partition_cache.h
 * @brief The optional on-disk cache of a pCCE partition.
 */
#ifndef CCEX_PCCE_PARTITION_CACHE_H
#define CCEX_PCCE_PARTITION_CACHE_H

#include "clustering/cluster_pcce.h"
#include <string>
#include <vector>

// Capture the selected bath before shrink_restspins removes the remainder.
std::string pccePartitionCacheIdentity(BathArray*, const Partition_info*, const Cluster*, const Config*);
void pccePartitionCached(BathSpin**, Partition_info*, Cluster*, Point*, int*, const Config*, const std::string&);

// The partition pccePartitionCached gave last (searched or read), for the rdip ladder (rdip_ladder.h): the centres, and the
// centre of each spin kept.
struct PccePartition { std::vector<Point> centers; std::vector<int> assigned; };
const PccePartition& pcceLastPartition();

#endif
