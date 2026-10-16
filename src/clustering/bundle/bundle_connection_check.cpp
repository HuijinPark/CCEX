/**
 * @file bundle_connection_check.cpp
 * @brief Bundle clustering: the bounded check of which bundle connections may be omitted.
 */
#include "../../../include/clustering/bundle/bundle_connection_check.h"
#include "../../../include/clustering/bundle/bundle_text.h"
#include <climits>
#include "../../../include/clustering/bundle/bundle_work.h"
#include "../../../include/clustering/bundle/bundle_cluster.h"
#include "../../../include/calculation/scaled_product.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>

namespace {
using Sites=std::vector<int>;
using Pair=std::pair<int,int>;
struct Domain {
    Sites bundles,sites;
    std::vector<AdaptiveTerm> before,after;
    double base_error=0,candidate_error=0,change=0;
    int worst_state=0,worst_step=0;
};
struct Check {
    bool auditing=false;
    AdaptivePlan full,candidate;
    Sites dimensions;
    std::vector<AdaptiveEdge> physical;
    std::vector<Domain> domains;
    std::vector<Sites> probes;
    std::map<Sites,size_t> probe_index;
    std::vector<Pair> removed;
    std::vector<std::vector<float>> states;
    int qubit_dim=1,bundle_cap=0,states_checked=0;
    size_t eligible=0,dimension_limited=0;
    std::string status="no_weak_candidate";
    double seconds=0;
} check;

using adaptive_bundle::require;
[[noreturn]] void die(const char* why) {
    fprintf(stderr,"Error: adaptive connection check (rank %d): %s\n",rank,why);
    MPI_Abort(MPI_COMM_WORLD,1);std::abort();
}
bool fits(const Sites& bundles,const Config*) {
    // Every bundle set is checkable now: the cluster dimension follows from the
    // bundle cap and the order instead of being capped on its own.
    int d=check.qubit_dim;
    for(int g:bundles) {
        if(d>INT_MAX/check.full.dimensions[g])return false;   // int range only
        d*=check.full.dimensions[g];
    }
    return true;
}
Sites sitesOf(const Sites& bundles) {
    Sites sites;
    for(int g:bundles)sites.insert(sites.end(),check.full.bundles[g].begin(),check.full.bundles[g].end());
    std::sort(sites.begin(),sites.end());return sites;
}
std::vector<AdaptiveTerm> localTerms(const AdaptivePlan& plan,const Sites& bundles,const Config* c) {
    // The expander requires a complete cover of LOCAL site IDs. Translate both
    // directions; global IDs are not generally contiguous in a selected domain.
    const auto sites=sitesOf(bundles);
    AdaptivePlan local;
    std::map<int,int> bundle_id;
    for(int g:bundles) {
        bundle_id[g]=int(local.bundles.size());Sites row;
        for(int s:plan.bundles[g])row.push_back(int(std::lower_bound(sites.begin(),sites.end(),s)-sites.begin()));
        local.bundles.push_back(row);local.dimensions.push_back(plan.dimensions[g]);
    }
    for(auto e:plan.edges)if(bundle_id.count(e.first) && bundle_id.count(e.second))
        local.edges.emplace_back(bundle_id.at(e.first),bundle_id.at(e.second));
    auto terms=expandAdaptivePlan(local,int(sites.size()),c->order,
        {size_t(c->expansion_max_clusters),size_t(c->expansion_max_work)});
    for(auto& t:terms)for(int& s:t.sites)s=sites[s];
    return terms;
}
std::vector<BathSpin*> stateSpins(BathArray* ba,DefectArray* dfa) {
    std::vector<BathSpin*> out;
    for(int s=0;s<ba->nspin;++s) {
        out.push_back(ba->bath[s]);
        if(dfa->naddspins)for(int j=0;j<dfa->naddspins[s];++j)out.push_back(dfa->subbath[s][j]);
    }
    return out;
}
void broadcastProbes() {
    int count=int(check.probes.size());MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank!=0)check.probes.resize(count);
    for(auto& ids:check.probes) {
        int n=int(ids.size());MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        if(rank!=0)ids.resize(n);
        MPI_Bcast(ids.data(),n,MPI_INT,0,MPI_COMM_WORLD);
    }
}
ScaledComplex reconstruct(const std::vector<AdaptiveTerm>& terms,int t,int steps,
                          const std::vector<std::complex<double>>& raw) {
    const ScaledComplex zero(raw[t]);
    ScaledComplex out=zero;
    for(const auto& term:terms) {
        const auto index=check.probe_index.at(term.sites);
        out=out*(ScaledComplex(raw[(index+1)*steps+t])*zero.inverse()).power(term.power);
    }
    return out;
}
double error(std::complex<double> a,std::complex<double> b) {
    if(!std::isfinite(a.real()) || !std::isfinite(a.imag()) ||
       !std::isfinite(b.real()) || !std::isfinite(b.imag()))return INFINITY;
    return std::abs(a-b)/(ccex::adaptiveCheckAtol+ccex::adaptiveCheckRtol*std::abs(b));
}
void compareDomains(const Config* c,const std::vector<std::complex<double>>& raw,int state) {
    const bool gcce=!strcasecmp(c->method,"gcce");
    for(auto& d:check.domains) {
        const auto offset=(check.probe_index.at(d.sites)+1)*c->nstep;
        const auto before0=gcce ? reconstruct(d.before,0,c->nstep,raw).inverse() : ScaledComplex();
        const auto after0=gcce ? reconstruct(d.after,0,c->nstep,raw).inverse() : ScaledComplex();
        const auto exact0=gcce ? ScaledComplex(raw[offset]).inverse() : ScaledComplex();
        for(int t=0;t<c->nstep;++t) {
            const auto exact=(ScaledComplex(raw[offset+t])*exact0).value();
            const auto before=(reconstruct(d.before,t,c->nstep,raw)*before0).value();
            const auto after=(reconstruct(d.after,t,c->nstep,raw)*after0).value();
            d.base_error=std::max(d.base_error,error(before,exact));
            const double e=error(after,exact);
            if(e>d.candidate_error) {d.candidate_error=e;d.worst_state=state;d.worst_step=t;}
            d.change=std::max(d.change,error(after,before));
        }
    }
    ++check.states_checked;
}
using adaptive_bundle::jsonArray;
void number(std::ostream& f,double x) {if(std::isfinite(x))f<<x;else f<<"null";}
void report(const Config* c,const Output* op,bool accepted) {
    auto& selected=accepted ? check.candidate : check.full;
    adaptiveWritePlan(c,selected,check.dimensions,check.physical,check.qubit_dim,check.bundle_cap);
    const std::string path=std::string(op->outfile)+".adaptive.json";
    std::ofstream f(path);require(bool(f),"cannot write adaptive connection report");
    f<<std::setprecision(17)<<"{\"mode\":\"rdip-omitted-checked\",\"status\":\""<<check.status
     <<"\",\"accuracy_certified\":false,\"scope\":\"selected domains; full-bath order convergence unresolved\","
     <<"\"order\":"<<c->order
     <<",\"qubit_dimension\":"<<check.qubit_dim<<",\"atol\":"<<ccex::adaptiveCheckAtol
     <<",\"rtol\":"<<ccex::adaptiveCheckRtol<<",\"domain_limit\":"<<ccex::adaptiveCheckDomains
     <<",\"probe_limit\":"<<ccex::adaptiveCheckProbes<<",\"memory_limit_bytes\":"<<ccex::adaptiveCheckBytes
     <<",\"weak_candidates\":"<<check.eligible<<",\"proposed_omissions\":"<<check.removed.size()
     <<",\"accepted_omissions\":"<<(accepted ? check.removed.size() : 0)
     <<",\"dimension_limited_domains\":"<<check.dimension_limited
     <<",\"unique_probes_per_state\":"<<check.probes.size()
     <<",\"states_checked\":"<<check.states_checked<<",\"steps_per_state\":"<<c->nstep
     <<",\"check_seconds\":"<<check.seconds<<",\"bundles\":[";
    for(size_t i=0;i<selected.bundles.size();++i) {if(i)f<<',';jsonArray(f,selected.bundles[i]);}
    f<<"],\"edges\":[";
    for(size_t i=0;i<selected.edges.size();++i) {auto e=selected.edges[i];f<<(i ? "," : "")<<'['<<e.first<<','<<e.second<<']';}
    f<<"],\"proposed_removed_edges\":[";
    for(size_t i=0;i<check.removed.size();++i) {auto e=check.removed[i];f<<(i ? "," : "")<<'['<<e.first<<','<<e.second<<']';}
    f<<"],\"domains\":[";
    for(size_t i=0;i<check.domains.size();++i) {
        const auto& d=check.domains[i];f<<(i ? "," : "")<<"{\"sites\":";jsonArray(f,d.sites);
        f<<",\"bundles\":";jsonArray(f,d.bundles);
        f<<",\"baseline_max_scaled_error\":";number(f,d.base_error);
        f<<",\"candidate_max_scaled_error\":";number(f,d.candidate_error);
        f<<",\"pruning_max_scaled_change\":";number(f,d.change);
        f<<",\"passed\":"<<(d.base_error<=1 && d.candidate_error<=1 && d.change<=1 && check.states_checked>0 ? "true" : "false")
         <<",\"worst_state\":"<<d.worst_state<<",\"worst_step\":"<<d.worst_step<<'}';
    }
    f<<"]}\n";f.flush();require(bool(f),"writing adaptive connection report failed");
    printf("[adaptive-check] status=%s proposed=%zu omitted=%zu domains=%zu probes_per_state=%zu states=%d seconds=%.6f report=%s full_bath_accuracy=unresolved\n",
        check.status.c_str(),check.removed.size(),accepted ? check.removed.size() : 0,check.domains.size(),
        check.probes.size(),check.states_checked,check.seconds,path.c_str());
}
}

void prepareAdaptiveConnectionCheck(const Config* c,const AdaptivePlan& plan,
    const std::vector<int>& dimensions,const std::vector<AdaptiveEdge>& physical,int qdim,int bundleCap) {
    check=Check();check.full=check.candidate=plan;check.dimensions=dimensions;
    check.physical=physical;check.qubit_dim=qdim;check.bundle_cap=bundleCap;
    if(plan.bundles.size()<=1) {check.status="single_bundle_no_links";return;}
    if(c->order==1) {check.status="order_one_no_link_expansion_unchecked";return;}
    const double duration=double(c->deltat)*(c->nstep-1);
    if(duration==0) {check.status="zero_duration_no_pruning";return;}
    Sites labels(dimensions.size());
    for(size_t g=0;g<plan.bundles.size();++g)for(int s:plan.bundles[g])labels[s]=int(g);
    std::map<Pair,long double> weights;
    for(const auto& e:physical) {
        int a=labels[e.a],b=labels[e.b];if(a>b)std::swap(a,b);
        if(a!=b)weights[{a,b}]+=e.strength;
    }
    std::vector<Pair> weak=plan.edges;
    std::sort(weak.begin(),weak.end(),[&](Pair a,Pair b) {return weights[a]!=weights[b] ? weights[a]<weights[b] : a<b;});
    // A proposal scale, NOT an error bound. Only coherence comparisons admit
    // omissions. Untested links stay present, including zero-weight links that
    // may still transmit central-spin-mediated correlations in gCCE.
    const long double proposal=std::sqrt(ccex::adaptiveCheckAtol+ccex::adaptiveCheckRtol);
    long double score=0;
    std::set<Sites> domains;
    for(auto e:weak) {
        if(weights[e]<=0)continue;
        if((score+weights[e]*duration)>proposal)break;
        ++check.eligible;
        if(check.removed.size()>=ccex::adaptiveCheckDomains)continue;
        Sites bundle{e.first,e.second};
        if(!fits(bundle,c)) {++check.dimension_limited;continue;}
        // Include the strongest neighboring bundle that fits, so a check sees
        // more than an isolated weak pair when the dimension cap permits it.
        int neighbor=-1;long double strongest=-1;
        for(size_t g=0;g<plan.bundles.size();++g)if(int(g)!=e.first && int(g)!=e.second) {
            const auto key=[&](int a) {return Pair(std::min(a,int(g)),std::max(a,int(g)));};
            const long double w=weights[key(e.first)]+weights[key(e.second)];
            auto expanded=bundle;expanded.push_back(int(g));std::sort(expanded.begin(),expanded.end());
            if(w>strongest && fits(expanded,c)) {neighbor=int(g);strongest=w;}
        }
        if(neighbor>=0) {bundle.push_back(neighbor);std::sort(bundle.begin(),bundle.end());}
        else if(plan.bundles.size()>2)++check.dimension_limited;
        domains.insert(bundle);check.removed.push_back(e);score+=weights[e]*duration;
    }
    if(check.removed.empty())return;
    std::set<Pair> removed(check.removed.begin(),check.removed.end());
    auto& edges=check.candidate.edges;
    edges.erase(std::remove_if(edges.begin(),edges.end(),[&](Pair e){return removed.count(e);}),edges.end());
    check.candidate.omitted_strength=score/duration;check.candidate.omitted_pairs=removed.size();
    Sites all;for(size_t g=0;g<plan.bundles.size();++g)all.push_back(int(g));
    if(fits(all,c)) {domains.clear();domains.insert(all);}
    std::set<Sites> probes;
    for(const auto& bundles:domains) {
        Domain d;d.bundles=bundles;d.sites=sitesOf(bundles);
        d.before=localTerms(plan,bundles,c);d.after=localTerms(check.candidate,bundles,c);
        probes.insert(d.sites);
        for(const auto& t:d.before)probes.insert(t.sites);
        for(const auto& t:d.after)probes.insert(t.sites);
        check.domains.push_back(std::move(d));
    }
    if(probes.size()>ccex::adaptiveCheckProbes) {check.status="probe_budget_unverified_links_kept";return;}
    check.probes.assign(probes.begin(),probes.end());
    for(size_t i=0;i<check.probes.size();++i)check.probe_index[check.probes[i]]=i;
    check.status="pending";
}

bool restoreAdaptiveCheckState(BathArray* ba,DefectArray* dfa,int state) {
    if(check.auditing || check.states.empty())return false;
    const auto spins=stateSpins(ba,dfa);
    require(size_t(state)<check.states.size() && check.states[state].size()==spins.size(),"missing checked spin state");
    for(size_t i=0;i<spins.size();++i)spins[i]->state=check.states[state][i];
    return true;
}

bool adaptiveCheckState(QubitArray* qa,BathArray* ba,DefectArray* dfa,Config* c,Pulse* pulse,Output* op,int state) {
    if(!check.auditing)return false;
    try {
        if(state>0) {
            auto& saved=check.states[state];
            for(auto* s:stateSpins(ba,dfa))saved.push_back(s->state);
        }
        const bool gcce=!strcasecmp(c->method,"gcce");
        const size_t count=(check.probes.size()+1)*c->nstep;
        std::vector<std::complex<double>> raw(count),combined(rank==0 ? count : 0);
        Output quiet=*op;strcpy(quiet.savemode,"coherence");
        if(rank==0) {
            std::unique_ptr<MatrixXcd[]> zero(gcce ? calCoherenceGcce(qa,nullptr,c,pulse,&quiet) : nullptr);
            for(int t=0;t<c->nstep;++t)raw[t]=gcce ? zero[t](0,0) : std::complex<double>(1,0);
        }
        for(size_t i=rank;i<check.probes.size();i+=nprocess) {
            const ProbeSubset probe=adaptiveProbeSubset(check.probes[i],ba,dfa,qa,gcce);
            require(probe.dim<=INT_MAX,"checked domain dimension leaves the int range");
            std::unique_ptr<MatrixXcd[]> values(gcce ? calCoherenceGcce(qa,probe.bath.get(),c,pulse,&quiet) :
                                                       calCoherenceCce(qa,probe.bath.get(),c,pulse));
            for(int t=0;t<c->nstep;++t)raw[(i+1)*c->nstep+t]=values[t](0,0);
        }
        MPI_Reduce(raw.data(),combined.data(),int(count),MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);
        if(rank==0)compareDomains(c,combined,state);
    } catch(const std::exception& e) {die(e.what());}
    return true;
}

void checkAdaptiveConnections(QubitArray* qa,BathArray* ba,DefectArray* dfa,Config* c,Pulse* pulse,Cluster* cls,Output* op) {
    if(strcasecmp(c->cluster_mode,"adaptive") || strcasecmp(c->bundle_connections,"checked"))return;
    try {
        const double start=MPI_Wtime();
        int run=0;
        if(rank==0 && check.status=="pending") {
            // Curves are discarded after each state; only spin states are kept
            // for replay. No srand/reseeding, and no unbounded curve cache.
            const long double bytes=static_cast<long double>(stateSpins(ba,dfa).size())*c->nstate*sizeof(float)+
                (static_cast<long double>(c->nstate)+1)*sizeof(std::vector<float>)+
                static_cast<long double>(check.probes.size()+1)*c->nstep*sizeof(std::complex<double>)*2;
            if(bytes>ccex::adaptiveCheckBytes)check.status="memory_budget_unverified_links_kept";
            else run=1;
        }
        MPI_Bcast(&run,1,MPI_INT,0,MPI_COMM_WORLD);
        if(run) {
            broadcastProbes();check.states.resize(size_t(c->nstate)+1);check.auditing=true;
            // calculate performs the same state/defect/Overhauser preparation
            // as production, then its audit hook bypasses the ordinary product.
            calculate(qa,ba,dfa,c,pulse,cls,op,nullptr);
            check.auditing=false;
        }
        int accepted=0;
        if(rank==0) {
            if(run) {
                bool base=true,candidate=true;
                for(const auto& d:check.domains) {
                    base &= d.base_error<=1;
                    candidate &= d.candidate_error<=1 && d.change<=1;
                }
                accepted=base && candidate;
                check.status=accepted ? "local_checks_passed" : (base ? "pruning_failed_links_restored" : "base_unconverged_links_restored");
            }
            check.seconds=MPI_Wtime()-start;report(c,op,accepted);
        }
        MPI_Bcast(&accepted,1,MPI_INT,0,MPI_COMM_WORLD);
        if(accepted)installAdaptivePlan(cls,ba,c,check.candidate);
        // Retain only the sampled states until production completes.
        auto saved=std::move(check.states);check=Check();check.states=std::move(saved);
    } catch(const std::exception& e) {die(e.what());}
}

void clearAdaptiveConnectionCheck() {check=Check();}
