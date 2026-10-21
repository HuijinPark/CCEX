#include "clustering/bundle/bundle_partition.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

namespace {
int checks=0;
template<class F> void rejects(F f) {
    bool failed=false;try {f();}catch(const std::exception&){failed=true;}
    assert(failed);++checks;
}
void audit(const AdaptivePlan& p,const std::vector<int>& dims,const std::vector<AdaptiveEdge>& input,
           int order,int bundleCap,double duration,double budget,bool exhaustive) {
    std::vector<int> labels(dims.size(),-1);
    for(size_t g=0;g<p.bundles.size();++g) {
        long long d=1;
        for(int i:p.bundles[g]) {assert(labels[i]==-1);labels[i]=g;d*=dims[i];}
        assert(d==p.dimensions[g] && d<=bundleCap);
    }
    assert(std::find(labels.begin(),labels.end(),-1)==labels.end());
    struct Weight {long double strength=0;bool allowed=false;};
    std::map<std::pair<int,int>,Weight> pairs;
    // Reconstruct every cross-group contribution from the original site pairs.
    for(const auto& e:input) {
        int a=labels[e.a],b=labels[e.b];if(a==b)continue;if(a>b)std::swap(a,b);
        auto& w=pairs[{a,b}];w.strength+=e.strength;w.allowed|=e.allowed;
    }
    const std::set<std::pair<int,int>> kept(p.edges.begin(),p.edges.end());
    assert(kept.size()==p.edges.size());
    long double omitted=0,total=0;size_t nonzero=0,dropped=0;
    std::vector<std::pair<std::pair<int,int>,Weight>> positive;
    for(const auto& entry:pairs) {
        if(entry.second.strength==0)continue;
        positive.push_back(entry);total+=entry.second.strength;++nonzero;
        if(!kept.count(entry.first)) {omitted+=entry.second.strength;++dropped;}
        else {
            assert(order==2 && entry.second.allowed);
        }
    }
    assert(std::abs(p.total_strength-total)<1e-12L);
    assert(std::abs(p.omitted_strength-omitted)<1e-12L);
    assert(p.positive_pairs==nonzero && p.omitted_pairs==dropped);
    assert(omitted*duration<=budget+1e-12L);
    expandAdaptivePlan(p,dims.size(),order);
    if(exhaustive) {
        assert(positive.size()<=10);size_t minimum=positive.size()+1;
        for(unsigned mask=0;mask<(1U<<positive.size());++mask) {
            long double missing=0;size_t count=0;bool valid=true;
            for(size_t i=0;i<positive.size();++i) {
                const auto& x=positive[i];
                if(mask&(1U<<i)) {
                    ++count;valid &= order==2 && x.second.allowed;
                } else missing+=x.second.strength;
            }
            if(valid && missing*duration<=budget)minimum=std::min(minimum,count);
        }
        assert(kept.size()==minimum);
    }
    ++checks;
}
}

int main() {
    // Ten individually weak links exceed the TOTAL budget: only three may drop.
    std::vector<AdaptiveEdge> weak;
    for(int i=0;i<5;++i)for(int j=i+1;j<5;++j)weak.push_back({i,j,.125});
    auto p=makeAdaptiveAutoPlan({2,2,2,2,2},weak,2,2,100,1,.375);
    assert(p.edges.size()==7 && p.omitted_strength==.375);
    audit(p,{2,2,2,2,2},weak,2,2,1,.375,true);
    // Longer total evolution time requires more connections at a fixed budget.
    auto longer=makeAdaptiveAutoPlan({2,2,2,2,2},weak,2,2,100,4,.375);
    assert(longer.edges.size()==10);
    // A strong merge used to be able to force a loss: the merged bundle's link to
    // its neighbour could exceed the cluster cap and have to be dropped, so the
    // budget decided. With one cap the link is simply kept - the merge costs
    // nothing at either budget.
    const std::vector<AdaptiveEdge> tri={{0,1,10},{1,2,.5},{0,2,.25}};
    // What the budget still decides is the remaining link, of strength .75: kept
    // when dropping it would overspend, dropped when it would not.
    auto kept=makeAdaptiveAutoPlan({2,2,2},tri,2,4,1,1,.5);
    assert(kept.bundles.size()==2 && kept.edges.size()==1 && kept.omitted_strength==0);
    auto dropped=makeAdaptiveAutoPlan({2,2,2},tri,2,4,1,1,.75);
    assert(dropped.bundles.size()==2 && dropped.edges.empty() && dropped.omitted_strength==.75L);
    audit(kept,{2,2,2},tri,2,4,1,.5,true);
    audit(dropped,{2,2,2},tri,2,4,1,.75,true);
    // A distance-rejected site pair can be captured by a retained group union.
    const std::vector<AdaptiveEdge> bridge={{0,1,100,true},{0,2,10,false},{1,2,.5,true}};
    auto bridged=makeAdaptiveAutoPlan({2,2,2},bridge,2,4,50,1,0);
    assert(bridged.bundles.size()==2 && bridged.edges.size()==1 && bridged.omitted_strength==0);
    audit(bridged,{2,2,2},bridge,2,4,1,0,true);
    // Distance-rejected pairs are audited, never silently lost. (An over-cap pair
    // used to be its own category; with one cap a link is never refused for its
    // dimension, so there is nothing of that kind left to audit.)
    rejects([]{makeAdaptiveAutoPlan({2,2},{{0,1,1,false}},2,2,1,1,.5);});
    auto cut=makeAdaptiveAutoPlan({2,2},{{0,1,1,false}},2,2,1,1,1);
    assert(cut.radius_excluded_strength==1);
    auto zero=makeAdaptiveAutoPlan({2,2},{{0,1,1}},2,2,1,0,0);
    assert(zero.edges.empty());
    auto one=makeAdaptiveAutoPlan({2,2},{{0,1,1}},1,4,0,1,0);
    assert(one.bundles.size()==1);
    rejects([]{makeAdaptiveAutoPlan({2,2},{{0,1,1}},1,2,0,1,0);});
    rejects([]{makeAdaptiveAutoPlan({2,2},{{0,1,1}},2,2,0,1,-1);});
    rejects([]{makeAdaptiveAutoPlan({2,2},{{0,1,1}},2,2,0,NAN,1);});
    rejects([]{makeAdaptiveAutoPlan({2,2},{{0,1,INFINITY}},2,2,0,1,1);});
    rejects([]{makeAdaptiveAutoPlan({2,2},{{0,1,1},{1,0,2}},2,2,0,1,1);});

    std::mt19937 rng(75191);
    for(int trial=0;trial<240;++trial) {
        const int n=2+rng()%4;std::vector<int> dims(n,2);std::vector<AdaptiveEdge> edges;
        for(int a=0;a<n;++a)for(int b=a+1;b<n;++b)edges.push_back({a,b,double(rng()%17)/8,bool(rng()%4)});
        const double duration=trial%4==0 ? 0 : trial%4==1 ? .25 : 1;
        // Budget includes compulsory distance omissions, so feasibility is known.
        double budget=double(rng()%9)/8;
        for(const auto& e:edges)if(!e.allowed)budget+=e.strength*duration;
        auto plan=makeAdaptiveAutoPlan(dims,edges,2,2,100,duration,budget);
        audit(plan,dims,edges,2,2,duration,budget,true);
        auto shuffled=edges;std::shuffle(shuffled.begin(),shuffled.end(),rng);
        for(auto& e:shuffled)if(rng()%2)std::swap(e.a,e.b);
        auto again=makeAdaptiveAutoPlan(dims,shuffled,2,2,100,duration,budget);
        assert(again.bundles==plan.bundles && again.edges==plan.edges && again.omitted_strength==plan.omitted_strength);
    }
    for(int trial=0;trial<300;++trial) {
        const int n=2+rng()%15;std::vector<int> dims(n);for(auto& d:dims)d=2+rng()%2;
        std::vector<AdaptiveEdge> edges;
        for(int a=0;a<n;++a)for(int b=a+1;b<n;++b)edges.push_back({a,b,double(rng()%17)/8});
        const double duration=.25,budget=double(trial%11)/4;
        auto plan=makeAdaptiveAutoPlan(dims,edges,2,12,1,duration,budget);
        audit(plan,dims,edges,2,12,duration,budget,false);
        auto shuffled=edges;std::shuffle(shuffled.begin(),shuffled.end(),rng);
        auto again=makeAdaptiveAutoPlan(dims,shuffled,2,12,1,duration,budget);
        assert(again.bundles==plan.bundles && again.edges==plan.edges && again.omitted_strength==plan.omitted_strength);
        if(budget==0) {
            // Nonzero complete graph: no pruning must reproduce the old policy.
            for(auto& e:edges)if(e.strength==0)e.strength=.0625;
            auto old=makeAdaptivePlan(dims,edges,2,12,1);
            auto fresh=makeAdaptiveAutoPlan(dims,edges,2,12,1,duration,0);
            assert(old.bundles==fresh.bundles && old.edges==fresh.edges);
        }
    }
    std::cout<<"PASS "<<checks<<" auto checks: exhaustive link subsets, aggregate budget, time horizon, mixed dimensions, distance/cap exclusions, merges and determinism\n";
}
