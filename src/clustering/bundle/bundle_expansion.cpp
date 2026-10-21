/**
 * @file bundle_expansion.cpp
 * @brief Bundle clustering: checking a plan and turning its bundle expansion into site clusters.
 */
#include "../../../include/clustering/bundle/bundle_partition.h"
#include "../../../include/clustering/bundle/bundle_text.h"
#include "../../../include/clustering/bundle/bundle_powers.h"
#include <algorithm>
#include <climits>
#include <map>
#include <stdexcept>

// The expansion itself is bundle_powers.cpp, shared with the algorithm half (2026-10-08); this file checks the plan,
// turns a failure into this half's message and the bundle sets into site terms.
namespace {
using adaptive_bundle::require;

const char* message(bundle_expansion::Failure kind) {
    using bundle_expansion::Failure;
    switch(kind) {
        case Failure::BadOrder: return "adaptive requires positive bundle order";
        case Failure::MaxClusters: return "expansion_max_clusters exceeded before coherence calculation";
        case Failure::MaxWork: return "expansion_max_work exceeded before coherence calculation";
        case Failure::BadEdge: return "invalid adaptive bundle edge";
        case Failure::DuplicateEdge: return "duplicate adaptive bundle edge";
        case Failure::BadExtra: return "invalid selective bundle cluster";
        case Failure::ExtraOutside: return "selective bundle outside partition";
        case Failure::ExtraDisconnected: return "selective cluster must be connected in the retained graph";
        case Failure::PowerTooLarge: return "adaptive inclusion-exclusion power exceeds 32-bit storage";
        case Failure::PowerOverflow: return "adaptive inclusion-exclusion intermediate overflow";
    }
    return "adaptive expansion failed";
}

bundle_expansion::Limits limitsOf(AdaptiveLimits limits) {
    require(limits.max_clusters>0 && limits.max_work>0,"adaptive expansion limits must be positive");
    auto clamp=[](size_t v) {return v>size_t(LLONG_MAX) ? LLONG_MAX : (long long)v;};
    return {clamp(limits.max_clusters),clamp(limits.max_work)};
}
}

std::vector<std::vector<int>> adaptiveBundleClusters(const AdaptivePlan& p,int order,AdaptiveLimits limits) {
    const bundle_expansion::Limits bounds=limitsOf(limits);
    require(p.bundles.size()<size_t(INT_MAX),"too many adaptive bundles");
    try {return bundle_expansion::family(int(p.bundles.size()),p.edges,order,bounds);}
    catch(const bundle_expansion::Error& e) {throw std::runtime_error(message(e.kind));}
}


std::vector<AdaptiveTerm> expandAdaptivePlan(const AdaptivePlan& p,int nsites,int order,AdaptiveLimits limits,
                                          const std::vector<std::vector<int>>& extra_bundles) {
    require(nsites>=0 && p.bundles.size()==p.dimensions.size(),"invalid adaptive site/bundle dimensions");
    std::vector<bool> seen(nsites,false);
    for(size_t g=0;g<p.bundles.size();++g) {
        require(!p.bundles[g].empty() && p.dimensions[g]>0,"empty adaptive bundle");
        for(int s:p.bundles[g]) {
            require(s>=0 && s<nsites && !seen[s],"duplicate or invalid adaptive site");seen[s]=true;
        }
    }
    require(std::find(seen.begin(),seen.end(),false)==seen.end(),"adaptive partition omitted sites");
    const bundle_expansion::Limits bounds=limitsOf(limits);
    require(p.bundles.size()<size_t(INT_MAX),"too many adaptive bundles");
    std::map<std::vector<int>,int> powers;
    try {powers=bundle_expansion::powers(int(p.bundles.size()),p.edges,order,bounds,extra_bundles);}
    catch(const bundle_expansion::Error& e) {throw std::runtime_error(message(e.kind));}
    std::vector<AdaptiveTerm> result;
    for(const auto& entry:powers) {
        AdaptiveTerm term{{},entry.second,int(entry.first.size())};
        for(int g:entry.first)term.sites.insert(term.sites.end(),p.bundles[g].begin(),p.bundles[g].end());
        std::sort(term.sites.begin(),term.sites.end());result.push_back(std::move(term));
    }
    std::sort(result.begin(),result.end(),[](const AdaptiveTerm& a,const AdaptiveTerm& b) {
        return a.sites.size()!=b.sites.size() ? a.sites.size()<b.sites.size() : a.sites<b.sites;
    });
    return result;
}
