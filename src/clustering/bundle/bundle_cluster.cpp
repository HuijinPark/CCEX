/**
 * @file bundle_cluster.cpp
 * @brief Bundle clustering (cluster.mode bundle): checking the settings and making the clusters of bundles.
 */
#include "../../../include/clustering/bundle/bundle_cluster.h"
#include "../../../include/clustering/bundle/bundle_text.h"
#include "../../../include/clustering/bundle/bundle_partition.h"
#include "../../../include/clustering/bundle/bundle_work.h"
#include "../../../include/clustering/bundle/bundle_connection_check.h"
#include "../../../include/calculation/hamiltonian.h"
#include "../../../include/common/memory.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace {
[[noreturn]] void fail(const std::string& message) {
    // Any rank may reach MPI_Abort first. Print on that rank so input errors
    // remain visible even when another rank has not reached its diagnostic yet.
    fprintf(stderr,"Error: adaptive clustering (rank %d): %s\n",rank,message.c_str());
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,1);std::abort();
}
using adaptive_bundle::require;
int spinDimension(double spin) {
    const double dim=2*spin+1;
    require(std::isfinite(dim) && dim>=1 && dim<INT_MAX && std::floor(dim)==dim,
            "spin quantum numbers must be finite nonnegative half-integers");
    return int(dim);
}
int product(int a,int b,int cap,const char* option="the int range") {
    if(a<=0 || b<=0 || a>cap/b)
        throw std::runtime_error(std::string("dimension product ")+std::to_string(a)+" * "+std::to_string(b)+
            " exceeds "+option+" (available dimension "+std::to_string(cap)+")");
    return a*b;
}
AdaptivePlan buildPlan(BathArray* ba,QubitArray* qa,DefectArray* dfa,const Config* c) {
    require(ba && qa && dfa,"adaptive clustering requires bath, qubit and defect metadata");
    const int n=ba->nspin;
    const bool automatic=!strcasecmp(c->bundle_connections,"auto");
    const bool checked=!strcasecmp(c->bundle_connections,"checked");
    require(n>=0 && n<INT_MAX,"invalid bath size");
    int qubitFactor=1;
    if(!strcasecmp(c->method,"gcce")) for(int i=0;i<qa->nqubit;++i)
        qubitFactor=product(qubitFactor,spinDimension(qa->qubit[i]->spin),INT_MAX);
    // The only dimension cap. A cluster is `order` bundles, so its dimension
    // follows from this and the order; it is reported, not capped (terminology
    // 4절, 2026-09-18).
    const int bundleCap=c->bundle_max_dim;
    std::vector<int> dimensions;
    std::vector<std::vector<BathSpin*>> sites(n);
    for(int i=0;i<n;++i) {
        sites[i].push_back(ba->bath[i]);
        const int id=DefectArray_findDefectIndex(dfa,ba->bath[i]->name);
        if(id>=0) {
            const int count=DefectArray_getDefect_idf_naddspin(dfa,id);
            require(count>=0 && dfa->naddspins && dfa->naddspins[i]==count,
                    "defect subspins were not initialized before clustering");
            for(int j=0;j<count;++j)sites[i].push_back(DefectArray_getSubbath_i_isp(dfa,i,j));
        }
        int d=1;
        for(auto* s:sites[i]) {
            require(s!=nullptr,"missing adaptive site spin");
            require(std::isfinite(s->gyro),"non-finite gyromagnetic ratio");
            for(double x:s->xyz)require(std::isfinite(x),"non-finite site coordinate");
            d=product(d,spinDimension(s->spin),bundleCap,"bundle_max_dim");
        }
        dimensions.push_back(d);
    }
    std::vector<AdaptiveEdge> edges;
    for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) {
        bool connected=false;
        double strength=0;
        for(size_t ai=0;ai<sites[i].size();++ai)for(size_t bi=0;bi<sites[j].size();++bi) {
            auto* a=sites[i][ai];auto* b=sites[j][bi];
            const double r=std::hypot(a->xyz[0]-b->xyz[0],a->xyz[1]-b->xyz[1],a->xyz[2]-b->xyz[2]);
            if(!std::isfinite(r) || r<=0)
                throw std::runtime_error("distinct bath sites have coincident/invalid spin positions: site="+
                    std::to_string(i)+" component="+std::to_string(ai)+" and site="+
                    std::to_string(j)+" component="+std::to_string(bi)+" (component 0 is the host)");
            if(checked || (r>c->rdipcut && ((automatic && !c->adaptive_rdip_is_set) || r<c->rdip)))connected=true;
            // Point-dipole eigenvalues are DS*(1,1,-2). H_BAR is CCEX's
            // unit conversion constant, not the SI numerical value of hbar.
            const double norm=std::sqrt(6.0)*std::abs(H_BAR*a->gyro*b->gyro)/std::pow(r,3);
            strength+=3.0*double(a->spin)*double(b->spin)*norm;
        }
        if(!automatic && !connected)continue;
        // Different indivisible defect sites interact through point dipoles,
        // as in BathArray_int_i_j. Sum ALL cross-component bounds, including
        // subspins. Auto mode also audits pairs excluded by the distance cuts.
        // ||sum T_ab I_a J_b|| <= 3 ||T||_F I J. This is a conservative coupling
        // norm bound in rad/ms, not a predicted coherence error or resonance test.
        require(std::isfinite(strength),"non-finite adaptive interaction strength");
        edges.push_back({i,j,strength,connected});
    }
    const double duration=double(static_cast<long double>(c->deltat)*(c->nstep-1));
    const auto& work=adaptiveWork(c);
    const auto physical=(work.active || checked) ? edges : std::vector<AdaptiveEdge>();
    auto planDims=dimensions;
    if(!work.bundles.empty()) {
        std::vector<int> labels(n,-1);planDims.clear();
        for(size_t g=0;g<work.bundles.size();++g) {
            int d=1;
            for(int s:work.bundles[g]) {
                require(s<n && labels[s]==-1,"requested bundles repeat or exceed loaded sites");
                labels[s]=int(g);d=product(d,dimensions[s],bundleCap,"bundle_max_dim");
            }
            planDims.push_back(d);
        }
        require(std::find(labels.begin(),labels.end(),-1)==labels.end(),"requested bundles must cover every loaded site");
        std::map<std::pair<int,int>,AdaptiveEdge> contracted;
        for(const auto& e:edges) {
            int a=labels[e.a],b=labels[e.b];if(a==b)continue;if(a>b)std::swap(a,b);
            auto key=std::make_pair(a,b);
            auto found=contracted.find(key);
            if(found==contracted.end())contracted[key]={a,b,e.strength,e.allowed};
            else {found->second.strength+=e.strength;found->second.allowed|=e.allowed;}
        }
        edges.clear();for(const auto& e:contracted)edges.push_back(e.second);
    }
    const bool merge=work.bundles.empty();
    const AdaptiveLimits limits{size_t(c->expansion_max_clusters),size_t(c->expansion_max_work)};
    AdaptivePlan plan=automatic ? makeAdaptiveAutoPlan(planDims,std::move(edges),c->order,
        bundleCap,c->bundle_merge_cutoff,duration,c->bundle_connection_budget,merge,limits) :
        makeAdaptivePlan(planDims,std::move(edges),c->order,bundleCap,c->bundle_merge_cutoff,merge,limits);
    if(!merge) {
        for(auto& bundle:plan.bundles) {
            std::vector<int> expanded;
            for(int g:bundle)expanded.insert(expanded.end(),work.bundles[g].begin(),work.bundles[g].end());
            std::sort(expanded.begin(),expanded.end());bundle=std::move(expanded);
        }
    }
    for(const auto& sites:work.extra_clusters) {
        int dimension=qubitFactor;
        for(int s:sites) {
            require(s<int(dimensions.size()),"selective site outside loaded bath");
            dimension=product(dimension,dimensions[s],INT_MAX);
        }
    }
    adaptiveWritePlan(c,plan,dimensions,physical,qubitFactor,bundleCap);
    adaptiveWritePhysics(c,ba,qa,dfa);
    if(checked)prepareAdaptiveConnectionCheck(c,plan,dimensions,physical,qubitFactor,bundleCap);
    printf("[adaptive] algorithm=%s method=%s bundle_order=%d sites=%d bundles=%zu pairs=%zu "
           "eligible_edges=%zu merges=%zu blocked_bundle=%zu qubit_dim=%d "
           "bundle_max_dim=%d merge_cutoff_rad_ms=%.17g\n",
           checked ? "coherence-checked-links-v1" : (automatic ? "coupling-budget-v1" : "coupling-bound-v1"),c->method,c->order,n,plan.bundles.size(),plan.edges.size(),plan.eligible,plan.merged,
           plan.blocked_bundle,qubitFactor,bundleCap,c->bundle_merge_cutoff);
    for(size_t g=0;g<plan.bundles.size();++g) {
        printf("[adaptive-bundle] id=%zu dim=%d sites=",g,plan.dimensions[g]);
        for(size_t j=0;j<plan.bundles[g].size();++j)
            printf("%s%d",j ? "," : "",plan.bundles[g][j]);
        printf("\n");
    }
    if(automatic) {
        printf("[adaptive-connections] mode=auto duration_ms=%.17g budget=%.17g positive_pairs=%zu kept_pairs=%zu omitted_pairs=%zu total_strength_rad_ms=%.17Lg omitted_strength_rad_ms=%.17Lg omitted_score=%.17Lg radius_excluded_pairs=%zu radius_excluded_score=%.17Lg rdip_limited=%d\n",
               duration,c->bundle_connection_budget,plan.positive_pairs,plan.edges.size(),plan.omitted_pairs,
               plan.total_strength,plan.omitted_strength,plan.omitted_strength*duration,
               plan.radius_excluded_pairs,plan.radius_excluded_strength*duration,
               int(c->adaptive_rdip_is_set));
        if(verbosity)for(const auto& e:plan.selected_connections)
            printf("[adaptive-link] bundles=%d,%d strength_rad_ms=%.17g score=%.17Lg\n",
                   e.a,e.b,e.strength,static_cast<long double>(e.strength)*duration);
    }
    return plan;
}

std::vector<int> pack(const AdaptivePlan& p) {
    size_t count=2+3*p.bundles.size()+2*p.edges.size();
    for(const auto& g:p.bundles)count+=g.size();
    require(count<size_t(INT_MAX),"adaptive plan exceeds MPI message size");
    std::vector<int> data;data.reserve(count);
    data.push_back(int(p.bundles.size()));
    for(size_t g=0;g<p.bundles.size();++g) {
        data.push_back(int(p.bundles[g].size()));data.push_back(p.dimensions[g]);
        data.insert(data.end(),p.bundles[g].begin(),p.bundles[g].end());
    }
    require(p.edges.size()<size_t(INT_MAX),"too many adaptive bundle pairs");
    data.push_back(int(p.edges.size()));
    for(auto e:p.edges) { data.push_back(e.first);data.push_back(e.second); }
    return data;
}
AdaptivePlan unpack(const std::vector<int>& data) {
    AdaptivePlan p;size_t at=0;
    auto next=[&]() { require(at<data.size(),"truncated adaptive plan");return data[at++]; };
    const int bundles=next();require(bundles>=0,"invalid bundle count");
    for(int g=0;g<bundles;++g) {
        const int n=next(),d=next();require(n>0 && size_t(n)<=data.size()-at,"invalid bundle length");
        p.dimensions.push_back(d);p.bundles.emplace_back(data.begin()+at,data.begin()+at+n);at+=n;
    }
    const int pairs=next();require(pairs>=0 && size_t(pairs)<=(data.size()-at)/2,"invalid pair count");
    for(int i=0;i<pairs;++i) { int a=next(),b=next();p.edges.emplace_back(a,b); }
    require(at==data.size(),"extra data in adaptive plan");return p;
}
}

void validateAdaptiveConfig(const Config* c,const Cluster* cls) {
    if(strcasecmp(c->cluster_mode,"adaptive")) {
        if(c->adaptive_options_present)fail("the cluster.bundle block and cluster.adaptive.work_file require cluster.mode bundle");
        return;
    }
    if(strcasecmp(c->method,"cce") && strcasecmp(c->method,"gcce"))
        fail("cluster.mode bundle requires method cCCE or gCCE");
    if(strcasecmp(c->quantity,"coherence"))fail("cluster.mode bundle currently supports quantity=coherence");
    if(c->order<1)fail("cluster.mode bundle requires positive bundle order");
    if(c->expansion_max_clusters<1 || c->expansion_max_work<1)fail("cluster.bundle.max_clusters and cluster.bundle.max_work must be positive");
    if(c->bundle_max_dim<1)
        fail("cluster.mode bundle requires cluster.bundle.max_dim");
    if(!std::isfinite(c->bundle_merge_cutoff) || c->bundle_merge_cutoff<0)
        fail("cluster.mode bundle requires a nonnegative cluster.bundle.merge_cutoff");
    const bool automatic=!strcasecmp(c->bundle_connections,"auto");
    const bool checked=!strcasecmp(c->bundle_connections,"checked");
    if(!automatic && !checked && strcasecmp(c->bundle_connections,"distance"))fail("cluster.bundle.connections must be distance or auto");
    if(automatic) {
        if(!std::isfinite(c->bundle_connection_budget) || c->bundle_connection_budget<0)
            fail("cluster.bundle.connections auto requires a finite nonnegative cluster.bundle.connection_budget");
        if(!std::isfinite(c->deltat) || c->deltat<0 || c->nstep<1)
            fail("cluster.bundle.connections auto requires finite nonnegative deltat and nstep >= 1");
    } else if(c->bundle_connection_budget>=0)fail("cluster.bundle.connection_budget requires cluster.bundle.connections auto");
    if(checked) {
        if(!std::isfinite(c->deltat) || c->deltat<0 || c->nstep<1 || c->nstate<0)
            fail("rdip omission requires finite nonnegative deltat, nstep >= 1 and nstate >= 0");
        if(c->rdipcut!=0)fail("an inner rdip radius requires the outer one (rdip [inner, outer]) in the simplified bundle input");
        // Explicit plan/probe requests are an advanced protocol, with their own
        // validation policy. Avoid silently mixing two independent optimizers.
        const auto& work=adaptiveWork(c);
        if(work.active && (work.mode!="run" || !work.bundles.empty()))
            fail("automatic checked connections require a normal run, use explicit cluster.bundle.connections for plan/probe or fixed bundles");
    }
    if(!std::isfinite(c->rdip) || !std::isfinite(c->rdipcut) || c->rdipcut<0 ||
       (((!automatic && !checked) || c->adaptive_rdip_is_set) && (c->rdip<=0 || c->rdipcut>=c->rdip)))
        fail("cluster.mode bundle requires 0 <= inner < outer rdip");
    if(strcasecmp(c->pcce_cache_mode,"off"))fail("cluster.mode bundle does not use cluster.pcce.cache.mode");
    if(cls && cls->nk)for(int i=0;i<=cls->order;++i)
        if(cls->nk[i])fail("cluster.mode bundle currently requires cluster.selection.max_per_order all 0 (all retained bundle clusters)");
}

void clusterizeAdaptive(Cluster* cls,BathArray* ba,QubitArray* qa,DefectArray* dfa,const Config* c) {
    validateAdaptiveConfig(c,cls);
    try {
        AdaptivePlan plan;
        if(rank==0)plan=buildPlan(ba,qa,dfa,c);
        installAdaptivePlan(cls,ba,c,plan);
    } catch(const std::exception& e) {fail(e.what());}
}

static std::vector<std::vector<int>> installedBundles;
static std::vector<std::pair<int,int>> installedEdges;
const std::vector<std::vector<int>>& adaptiveInstalledBundles(){return installedBundles;}
const std::vector<std::pair<int,int>>& adaptiveInstalledEdges(){return installedEdges;}
void clearAdaptiveInstalledBundles(){installedBundles.clear();installedEdges.clear();}

void installAdaptivePlan(Cluster* cls,BathArray* ba,const Config* c,const AdaptivePlan& rootPlan) {
    try {
        std::vector<int> data;
        if(rank==0)data=pack(rootPlan);
        int length=int(data.size());MPI_Bcast(&length,1,MPI_INT,0,MPI_COMM_WORLD);
        require(length>=2,"empty adaptive MPI plan");
        if(rank!=0)data.resize(length);
        MPI_Bcast(data.data(),length,MPI_INT,0,MPI_COMM_WORLD);
        const AdaptivePlan plan=unpack(data);
        installedBundles=plan.bundles;
        installedEdges=plan.edges;
        std::vector<std::vector<int>> extraBundles;
        std::vector<int> owner(ba->nspin,-1);
        for(size_t g=0;g<plan.bundles.size();++g)for(int s:plan.bundles[g])owner[s]=int(g);
        for(const auto& sites:adaptiveWork(c).extra_clusters) {
            std::set<int> ids;
            for(int s:sites) {require(s<ba->nspin,"selective site outside loaded bath");ids.insert(owner[s]);}
            std::vector<int> expanded;
            long long dimension=1;
            for(int g:ids) {
                require(plan.dimensions[g]>0 && dimension<=INT_MAX/plan.dimensions[g],
                        "selective cluster dimension leaves the int range");
                dimension*=plan.dimensions[g];
                expanded.insert(expanded.end(),plan.bundles[g].begin(),plan.bundles[g].end());
            }
            std::sort(expanded.begin(),expanded.end());
            require(expanded==sites,"selective cluster splits an explicit bundle");
            extraBundles.emplace_back(ids.begin(),ids.end());
        }
        const auto terms=expandAdaptivePlan(plan,ba->nspin,c->order,
            {size_t(c->expansion_max_clusters),size_t(c->expansion_max_work)},extraBundles);
        if(adaptiveWork(c).mode=="except") {
            // Every excepted set has to be a term the run calculates (a nonzero power): the driver works X out of the
            // same exponent table, so one that is not there means the two tables differ (adaptive-regroup-stage3b.md
            // 1절). Checked here, on every rank, before any state is calculated.
            std::set<std::vector<int>> calculated;
            for(const auto& t:terms)if(t.power!=0)calculated.insert(t.sites);
            for(const auto& sites:adaptiveWork(c).except)
                require(calculated.count(sites)!=0,"an except entry is not a term of this partition's expansion");
        }
        int maximum=0;
        for(const auto& t:terms)maximum=std::max(maximum,int(t.sites.size()));
        std::vector<std::vector<const AdaptiveTerm*>> bySize(size_t(maximum)+1);
        for(const auto& t:terms)bySize[t.sites.size()].push_back(&t);
        Cluster_freeClusinfo(cls);Cluster_freeNk(cls);Cluster_setOrder(cls,maximum);
        Cluster_allocNk(cls);Cluster_allocClusinfo(cls,maximum);
        // calculate() uses L0 * product((L_C/L0)^power), so its 0th flag is 1.
        // It must NOT be replaced by 1-sum(power): L0 normalization is already
        // applied to every factor in both the matrix and scaled gCCE paths.
        cls->clusinfo[0][0][0]=1;
        for(int n=1;n<=maximum;++n) {
            const auto& bucket=bySize[n];
            require(bucket.size()<size_t(INT_MAX-1),"too many adaptive propagation terms");
            const int count=int(bucket.size());
            free(cls->clusinfo[n][0]);free(cls->clusinfo[n]);
            cls->clusinfo[n]=allocInt2d(count+1,n+1);
            cls->clusinfo[n][0][0]=count+1;
            for(int i=0;i<count;++i) {
                cls->clusinfo[n][i+1][0]=bucket[i]->power;
                std::copy(bucket[i]->sites.begin(),bucket[i]->sites.end(),cls->clusinfo[n][i+1]+1);
            }
        }
        if(rank==0) {
            std::map<int,size_t> counts{{1,0},{2,0}};
            for(const auto& t:terms)++counts[t.bundle_order];
            printf("[adaptive-expansion] bundle_order=%d max_sites=%d terms=%zu",c->order,maximum,terms.size());
            for(const auto& item:counts)printf(" bundle%d_terms=%zu",item.first,item.second);
            printf("\n");
        }
    } catch(const std::exception& e) { fail(e.what()); }
}
