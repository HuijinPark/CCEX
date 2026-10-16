/**
 * @file rdip_ladder.h
 * @brief The rdip ladder: the coherence of every smaller rdip assembled from one order-2 run (input key rdip_ladder).
 */
#ifndef CCEX_AUTO_CONVERGENCE_RDIP_LADDER_H
#define CCEX_AUTO_CONVERGENCE_RDIP_LADDER_H
// rdip ladder ((internal notes)): one order-2 run (cCCE, gCCE or pCCE) also assembles the product of
// every rung r_k < rdip. The expansion is over units - the spins (cCCE, gCCE) or the partitions (pCCE, whose partition does not
// depend on rdip). Two units are a pair cluster when their connection distance D lies in (rdipcut, rdip): the spin distance
// (BathArray_connectivity), the centre distance (pcce_connections centers) or the smallest spin distance above rdipcut between
// the two partitions (pcce_connections spins). So rung k keeps the pairs with D < r_k (power 1) and every single unit with power
// 1 - deg_k(u); every factor is L_C / L_0 (L_0 = 1 for cCCE and pCCE, the qubit alone for gCCE), as in the run's own product.
// solveTilde keeps only clusters of nonzero power, so the single units the run drops (deg(u) = 1) are re-added at power 0: the
// main product skips them and the rungs use them.
#include <vector>
#include "calculation/scaled_product.h"
#include "clustering/cluster.h"
#include "clustering/cluster_pcce.h"
#include "common/memory.h"
#include "common/output.h"
#include "common/utilities.h"
#include "input/general.h"
#include "physical/bath.h"

namespace rdipladder {

[[noreturn]] void refuse(const char* why);
// After clusterize, on every rank: refuse what the ladder does not cover and re-add the dropped single units at power 0.
void prepare(Config* c,Cluster* cls,BathArray* ba,Output* op);
// Rank 0: write rung k as <outfile>_rdip<r_k>_state<s>_{wiDiv,noDiv}.
void save(Output* op,float rung,MatrixXcd* wD,MatrixXcd* nD,int nstep,double deltat,int istate);

// The units of the expansion and the distance that connects two of them.
struct Units {
    int size=1;                                    // spins per unit: 1, or sK (pCCE)
    bool centres=false;                            // pCCE pcce_connections centers
    double cut=0;
    BathArray* ba=nullptr;
    std::vector<int> of;                           // unit of each spin (pCCE; empty: the spin itself)
    std::vector<std::vector<int>> members;         // spins of each unit, ascending (pCCE)
    std::vector<Point> centre;                     // centre of each unit (pCCE)
    Units(Config* c,Cluster* cls,BathArray* ba);
    int count() const { return of.empty() ? BathArray_getNspin(ba) : int(members.size()); }
    int unit(int spin) const { return of.empty() ? spin : of[spin]; }
    // connected at rdip r exactly when distance(a, b) < r (the comparisons of the run's own connectivity); +inf: never
    double distance(int a,int b) const;
};

// Per-rank products of the rungs for one state.
struct Ladder {
    bool on=false;
    int K=0,nstep=0;
    std::vector<float> r;                          // rungs, then rdip
    std::vector<std::vector<int>> deg;             // deg[k][unit] for k = 0..K (K = rdip)
    std::vector<std::vector<ScaledComplex>> with,without;
    Units units;
    Ladder(Config* c,Cluster* cls,BathArray* ba,int steps):on(c->nrdip_ladder>0),K(c->nrdip_ladder),nstep(steps),units(c,cls,ba){
        if(!on)return;
        r.assign(c->rdip_ladder,c->rdip_ladder+K);r.push_back(Config_getRdip(c));
        const int n=units.count();
        deg.assign(K+1,std::vector<int>(n,0));
        for(int a=0;a<n;++a)for(int b=a+1;b<n;++b){
            const double d=units.distance(a,b);
            for(int k=0;k<=K;++k)if(d<r[k]){++deg[k][a];++deg[k][b];}
        }
    }
    void reset(){
        if(!on)return;
        with.assign(K,std::vector<ScaledComplex>(nstep));without.assign(K,std::vector<ScaledComplex>(nstep));
    }
    void add(const int* cl,int n,int power,const MatrixXcd* res,const std::vector<ScaledComplex>& inverseZero,bool unitZero){
        if(!on)return;
        if(n==units.size){
            const int u=units.unit(cl[0]);
            if(power!=1-deg[K][u])refuse("found a single-unit power that is not 1 - (pairs of the unit), the cluster set is not the plain rdip graph");
            for(int k=0;k<K;++k){
                const int p=1-deg[k][u];
                if(p==0)continue;
                for(int t=0;t<nstep;++t)assembleFactor(with[k][t],without[k][t],res[t](0,0),inverseZero[t],p,unitZero,nullptr);
            }
            return;
        }
        if(n!=2*units.size)refuse("found a cluster that is neither one unit nor two");
        if(power!=1)refuse("found a pair of power other than 1");
        const double d=units.distance(units.unit(cl[0]),units.unit(cl[units.size]));
        for(int k=0;k<K;++k){
            if(!(d<r[k]))continue;
            for(int t=0;t<nstep;++t)assembleFactor(with[k][t],without[k][t],res[t](0,0),inverseZero[t],1,unitZero,nullptr);
        }
    }
    void normalize(){
        for(auto& v:with)for(auto& x:v)x.normalize();
        for(auto& v:without)for(auto& x:v)x.normalize();
    }
};

}  // namespace rdipladder
#endif
