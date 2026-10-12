#ifndef CCEX_CLUSTER_CONTRIBUTIONS_H
#define CCEX_CLUSTER_CONTRIBUTIONS_H
#include "../../calculation/simulator.h"
#include <map>
#include <vector>

// Adds missing proper BUNDLE subsets at power zero for diagnostics only.
// Production powers and Hamiltonians are preserved. Disabled by default.
void prepareContributionClusters(const Config*,Cluster*,BathArray*);

class ClusterContributions {
    const Config* config;
    BathArray* bath;
    std::vector<std::vector<int>> sites,bundleIds;
    std::map<std::vector<int>,size_t> index;
    std::vector<int> powers,dimensions;
    std::vector<doublec> raw;
    double captureSeconds=0;
public:
    ClusterContributions(const Config*,Cluster*,BathArray*);
    void capture(const int*,int,int,const MatrixXcd*);
    void finish(const MatrixXcd*,int,const Output*);
};
#endif
