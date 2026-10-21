#include "clustering/bundle/bundle_partition.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <complex>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>

namespace {
using Z=std::complex<double>;
int mask(const std::vector<int>& sites) { int m=0;for(int i:sites)m|=1<<i;return m; }
Z signal(int m) { return std::polar(0.85+0.001*(m%101),0.007*(m%37)); }
int checks=0;
void checkExpansion(const std::vector<std::vector<int>>& groups,int n) {
    for(int graph=0;graph<4;++graph)for(int order:{1,2}) {
        AdaptivePlan p;p.bundles=groups;p.dimensions.assign(groups.size(),1);
        for(int a=0;a<int(groups.size());++a)for(int b=a+1;b<int(groups.size());++b)
            if(graph==0 || (graph==1 && b==a+1) || (graph==2 && a==0))p.edges.emplace_back(a,b);
        const Z zero=signal(0);Z expected=zero;
        for(const auto& g:groups)expected*=signal(mask(g))/zero;
        if(order==2)for(auto e:p.edges) {
            const int a=mask(groups[e.first]),b=mask(groups[e.second]);
            expected*=signal(a|b)*zero/(signal(a)*signal(b));
        }
        Z actual=zero;
        for(const auto& t:expandAdaptivePlan(p,n,order))actual*=std::pow(signal(mask(t.sites))/zero,t.power);
        assert(std::abs(expected-actual)<1e-12);++checks;
    }
}
void partitions(int at,int n,std::vector<std::vector<int>>& groups) {
    if(at==n) { checkExpansion(groups,n);return; }
    for(size_t i=0;i<groups.size();++i) {
        groups[i].push_back(at);partitions(at+1,n,groups);groups[i].pop_back();
    }
    groups.push_back({at});partitions(at+1,n,groups);groups.pop_back();
}
template<class F> void rejects(F f) {
    bool thrown=false;try{f();}catch(const std::exception&){thrown=true;}assert(thrown);++checks;
}
// One cap: a bundle may not exceed it. A cluster is `order` bundles and its
// dimension follows from that (terminology 4절, 2026-09-18) - the pair and term
// caps this used to assert went away with the cluster cap.
void verify(const AdaptivePlan& p,const std::vector<int>& dims,const std::vector<AdaptiveEdge>& edges,int bundleCap,int order) {
    expandAdaptivePlan(p,int(dims.size()),order);
    std::vector<int> label(dims.size(),-1);
    for(int g=0;g<int(p.bundles.size());++g) {
        long long d=1;for(int i:p.bundles[g]){label[i]=g;d*=dims[i];}
        assert(d==p.dimensions[g] && d<=bundleCap);
    }
    std::set<std::pair<int,int>> expected;
    if(order==2)for(auto e:edges) {
        int a=label[e.a],b=label[e.b];if(a>b)std::swap(a,b);
        if(a!=b)expected.emplace(a,b);
    }
    assert((std::set<std::pair<int,int>>(p.edges.begin(),p.edges.end())==expected));
    ++checks;
}
}

int main() {
    std::vector<std::vector<int>> groups;
    for(int n=0;n<=7;++n)partitions(0,n,groups);
    std::mt19937 rng(1938);
    int blocked=0,unequal=0;
    for(int repeat=0;repeat<500;++repeat) {
        const int n=2+rng()%18;
        std::vector<int> dims(n);
        for(auto& d:dims)d=std::vector<int>{2,3,6}[rng()%3];
        std::vector<AdaptiveEdge> edges;
        for(int a=0;a<n;++a)for(int b=a+1;b<n;++b)
            if(rng()%3)edges.push_back({a,b,double(rng()%11)});
        for(int order:{1,2})for(int bundleCap:{6,12,36,128}) {
            const double threshold=repeat%4;
            auto p=makeAdaptivePlan(dims,edges,order,bundleCap,threshold);
            verify(p,dims,edges,bundleCap,order);
            blocked+=p.blocked_bundle>0;
            if(p.bundles.size()>1)unequal+=p.bundles.front().size()!=p.bundles.back().size();
            auto shuffled=edges;std::shuffle(shuffled.begin(),shuffled.end(),rng);
            for(auto& e:shuffled)if(rng()%2)std::swap(e.a,e.b);
            auto other=makeAdaptivePlan(dims,shuffled,order,bundleCap,threshold);
            assert(p.bundles==other.bundles && p.dimensions==other.dimensions && p.edges==other.edges);
        }
    }
    assert(blocked && unequal);
    // The strongest pair merges to the cap; the third site would take it past.
    auto p=makeAdaptivePlan({2,2,2},{{0,1,10},{1,2,1}},2,4,0);
    assert(p.bundles.size()==2 && p.blocked_bundle==1 && p.edges.size()==1);
    auto noMerge=makeAdaptivePlan({2,3,6},{{0,1,0},{1,2,2}},2,36,1e20);
    assert(noMerge.bundles.size()==3 && noMerge.edges.size()==2);
    auto large=makeAdaptivePlan({46340,46340},{{0,1,1}},2,INT_MAX,0);
    assert(large.bundles.size()==1 && large.dimensions[0]==2147395600);
    // 46341^2 does not fit an int. fits() never forms the product, so the merge is
    // blocked rather than overflowed - it was the cluster cap that used to throw.
    auto tooLarge=makeAdaptivePlan({46341,46341},{{0,1,1}},2,INT_MAX,0);
    assert(tooLarge.bundles.size()==2 && tooLarge.blocked_bundle==1);++checks;
    rejects([]{makeAdaptivePlan({2,3},{{0,1,1}},2,2,0);});
    rejects([]{makeAdaptivePlan({2},{{0,0,1}},1,2,0);});
    rejects([]{makeAdaptivePlan({2,2},{{0,1,NAN}},1,2,0);});
    rejects([]{makeAdaptivePlan({2,2},{{0,1,1},{1,0,2}},1,2,0);});
    rejects([]{makeAdaptivePlan({2}, {},0,2,0);});
    rejects([]{AdaptivePlan p;p.bundles={{0},{0}};p.dimensions={2,2};expandAdaptivePlan(p,2,2);});
    rejects([]{AdaptivePlan p;p.bundles={{0}};p.dimensions={2};expandAdaptivePlan(p,2,2);});
    std::cout<<"PASS "<<checks<<" checks: independent complex inclusion-exclusion, coverage, "
             <<"mixed dimensions, deterministic ties, the bundle cap and overflow\n";
}
