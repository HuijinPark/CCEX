#include "clustering/bundle/bundle_partition.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

int bits(int x){return __builtin_popcount(unsigned(x));}
int mask(const std::vector<int>& ids){int x=0;for(int i:ids)x|=1<<i;return x;}
bool connected(int m,const AdaptivePlan& p) {
    int seen=m&-m,previous=0;
    while(seen!=previous) {
        previous=seen;
        for(auto e:p.edges)if(m&(1<<e.first) && m&(1<<e.second)) {
            if(seen&(1<<e.first))seen|=1<<e.second;
            if(seen&(1<<e.second))seen|=1<<e.first;
        }
    }
    return seen==m;
}
int checks=0;
void check(const AdaptivePlan& p,int order,const std::vector<std::vector<int>>& extras={}) {
    const int n=p.bundles.size();std::set<int> family;
    // Independent exhaustive oracle: inspect all masks and their subsets.
    for(int m=1;m<(1<<n);++m)if(bits(m)<=order && connected(m,p))
        for(int s=m;s;s=(s-1)&m)family.insert(s);
    for(const auto& ids:extras)for(int s=mask(ids);s;s=(s-1)&mask(ids))family.insert(s);
    std::map<int,int> expected,actual;
    for(int s:family) {
        int exponent=0;
        for(int t:family)if((t&s)==s)exponent+=(bits(t^s)%2 ? -1:1);
        int physical=0;for(int g=0;g<n;++g)if(s&(1<<g))physical|=mask(p.bundles[g]);
        if(exponent)expected[physical]=exponent;
    }
    int ns=0;for(auto& g:p.bundles)ns+=g.size();
    for(const auto& t:expandAdaptivePlan(p,ns,order,{},extras))actual[mask(t.sites)]=t.power;
    assert(actual==expected);
    if(order>=n && connected((1<<n)-1,p))assert(actual.size()==1 && actual.begin()->second==1);
    ++checks;
}
template<class F> void rejects(F f){bool threw=false;try{f();}catch(const std::runtime_error&){threw=true;}assert(threw);++checks;}
int main() {
    for(int n=1;n<=5;++n)for(int graph=0;graph<(1<<(n*(n-1)/2));++graph) {
        AdaptivePlan p;for(int i=0;i<n;++i){p.bundles.push_back({i});p.dimensions.push_back(2);}
        int bit=0;for(int a=0;a<n;++a)for(int b=a+1;b<n;++b,++bit)if(graph&(1<<bit))p.edges.emplace_back(a,b);
        for(int order=1;order<=n+1;++order)check(p,order);
        std::vector<std::vector<int>> extras;
        for(int m=1;m<(1<<n);++m)if(bits(m)==3 && connected(m,p)) {
            std::vector<int> ids;for(int i=0;i<n;++i)if(m&(1<<i))ids.push_back(i);
            extras.push_back(ids);
        }
        check(p,2,extras);
        if(extras.size()>1) {extras.resize(2);check(p,2,extras);}
    }
    // A cycle needs disconnected diagonal subsets in its order-3 correction.
    AdaptivePlan loop;loop.bundles={{0,3},{1},{2},{4}};loop.dimensions={6,2,3,2};
    loop.edges={{0,1},{0,3},{1,2},{2,3}};
    auto family=adaptiveBundleClusters(loop,3);
    assert(std::find(family.begin(),family.end(),std::vector<int>{0,2})!=family.end());
    check(loop,3);check(loop,4);
    rejects([&]{expandAdaptivePlan(loop,5,3,{5,10000});});
    rejects([&]{expandAdaptivePlan(loop,5,3,{100,10});});
    // One cap (terminology 4절): the two strongest pairs merge up to it and the
    // link between them is what gets blocked. The cluster cap that used to refuse
    // the first merge outright - and adaptiveDimensionsFit, which asked about it -
    // no longer exist.
    auto blocked=makeAdaptivePlan({2,2,2,2},{{0,1,10},{1,2,2},{2,3,1}},3,4,0);
    assert(blocked.bundles.size()==2 && blocked.blocked_bundle>0);
    std::mt19937 rng(912);
    for(int rep=0;rep<100;++rep) {
        std::vector<int> dims={2,3,2,2,3,2};std::vector<AdaptiveEdge> e;
        for(int a=0;a<6;++a)for(int b=a+1;b<6;++b)if(rng()%2)e.push_back({a,b,double(1+rng()%20)});
        for(int order:{3,4,5,6}) {
            auto p=makeAdaptivePlan(dims,e,order,12,0);check(p,order);
            auto other=makeAdaptiveAutoPlan(dims,e,order,12,0,.01,10);
            check(other,order);
            std::reverse(e.begin(),e.end());
            auto same=makeAdaptivePlan(dims,e,order,12,0);
            assert(same.bundles==p.bundles && same.edges==p.edges);
        }
    }
    std::cout<<"PASS "<<checks<<" high-order oracle/limit checks\n";
}
