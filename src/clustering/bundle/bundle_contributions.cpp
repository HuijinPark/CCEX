#include "../../../include/clustering/bundle/bundle_contributions.h"
#include "../../../include/clustering/bundle/cluster_log.h"
#include "../../../include/clustering/bundle/bundle_cluster.h"
#include "../../../include/clustering/bundle/bundle_work.h"
#include "../../../include/common/memory.h"
#include <algorithm>
#include <climits>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <set>
#include <stdexcept>

namespace {
using Sites=std::vector<int>;
[[noreturn]] void fail(const std::string& message) {
    fprintf(stderr,"[cluster-contributions] rank=%d %s\n",rank,message.c_str());
    MPI_Abort(MPI_COMM_WORLD,1);std::abort();
}
std::vector<Sites> bundles(const Config* c,int n) {
    if(!strcasecmp(c->cluster_mode,"adaptive"))return adaptiveInstalledBundles();
    std::vector<Sites> out;for(int i=0;i<n;++i)out.push_back({i});return out;
}
Sites owners(const Sites& sites,const std::vector<Sites>& bundles,const Sites& owner) {
    std::set<int> ids;for(int s:sites)ids.insert(owner.at(s));
    Sites expanded;for(int g:ids)expanded.insert(expanded.end(),bundles.at(g).begin(),bundles.at(g).end());
    std::sort(expanded.begin(),expanded.end());
    if(expanded!=sites)throw std::runtime_error("cluster splits an installed bundle");
    return Sites(ids.begin(),ids.end());
}
Sites ownerMap(const std::vector<Sites>& bundles,int n) {
    Sites out(n,-1);for(size_t g=0;g<bundles.size();++g)for(int s:bundles[g])out.at(s)=int(g);
    return out;
}
template<class F> void subsets(const Sites& ids,F visit) {
    Sites sub;
    std::function<void(size_t)> walk=[&](size_t first){
        for(size_t i=first;i<ids.size();++i) {
            sub.push_back(ids[i]);visit(sub);walk(i+1);sub.pop_back();
        }
    };walk(0);
}
size_t rowBytes(const Config* c,size_t siteCount) {
    // Raw curves, map/vector bookkeeping, site/bundle indices, output scratch.
    const long double bytes=sizeof(doublec)*static_cast<long double>(c->nstep)+256+16*siteCount;
    if(bytes>std::numeric_limits<size_t>::max())throw std::runtime_error("diagnostic table size overflow");
    return size_t(bytes);
}
void list(std::ostream& out,const Sites& values) {
    for(size_t j=0;j<values.size();++j){if(j)out<<',';out<<values[j];}
}
}

void prepareContributionClusters(const Config* c,Cluster* cls,BathArray* ba) {
    if(!c->cluster_contributions)return;
    try {
        if(strcasecmp(c->quantity,"coherence") || (strcasecmp(c->method,"cce") && strcasecmp(c->method,"gcce")))
            throw std::runtime_error("recording supports cCCE/gCCE coherence (standard or adaptive)");
        if(adaptiveWork(c).mode=="probe" || adaptiveWork(c).mode=="plan")return;
        const auto gs=bundles(c,ba->nspin);const auto own=ownerMap(gs,ba->nspin);
        std::map<Sites,int> family;
        size_t bytes=0,work=0;const long double cap=1048576.L*c->cluster_contribution_max_mb;
        const int order=Cluster_getOrder(cls);
        auto add=[&](const Sites& s,int power){
            if(!family.count(s)) {
                bytes+=rowBytes(c,s.size());
                if(bytes>cap)throw std::runtime_error("raw table exceeds cluster_contributions.max_mb before propagation");
                family.emplace(s,power);
            }
        };
        for(int n=1;n<=order;++n)for(int i=1;i<cls->clusinfo[n][0][0];++i) {
            Sites s(cls->clusinfo[n][i]+1,cls->clusinfo[n][i]+1+n);std::sort(s.begin(),s.end());
            add(s,cls->clusinfo[n][i][0]);
        }
        const auto production=family;
        for(const auto& entry:production) {
            const auto ids=owners(entry.first,gs,own);
            subsets(ids,[&](const Sites& sub){
                if(++work>size_t(c->expansion_max_work))throw std::runtime_error("diagnostic subset traversal exceeds bundle.max_work");
                Sites s;for(int g:sub)s.insert(s.end(),gs[g].begin(),gs[g].end());
                std::sort(s.begin(),s.end());add(s,0);
            });
        }
        // Append only missing rows. Retain existing production order/powers.
        for(int n=1;n<=order;++n) {
            std::vector<Sites> missing;
            for(const auto& entry:family)if(int(entry.first.size())==n && !production.count(entry.first))missing.push_back(entry.first);
            if(missing.empty())continue;
            const int old=cls->clusinfo[n][0][0];
            if(missing.size()>size_t(INT_MAX-old))throw std::runtime_error("diagnostic cluster count overflow");
            int** rows=allocInt2d(old+int(missing.size()),n+1);
            // Legacy row zero contains ONE int, not an (n+1)-wide cluster.
            for(int i=1;i<old;++i)std::copy(cls->clusinfo[n][i],cls->clusinfo[n][i]+n+1,rows[i]);
            rows[0][0]=old+int(missing.size());
            for(size_t i=0;i<missing.size();++i) {
                rows[old+i][0]=0;std::copy(missing[i].begin(),missing[i].end(),rows[old+i]+1);
            }
            freeInt2d(&cls->clusinfo[n],old);cls->clusinfo[n]=rows;
        }
        if(rank==0)printf("[cluster-contributions] existing_rows=%zu diagnostic_added=%zu table_budget_estimate_bytes=%zu per_rank=1\n",
                           production.size(),family.size()-production.size(),bytes);
    }catch(const std::exception& e){fail(e.what());}
}

ClusterContributions::ClusterContributions(const Config* c,Cluster* cls,BathArray* ba):config(c),bath(ba) {
    if(!c->cluster_contributions)return;
    const auto gs=bundles(c,ba->nspin);const auto own=ownerMap(gs,ba->nspin);
    for(int n=1;n<=Cluster_getOrder(cls);++n)for(int i=1;i<cls->clusinfo[n][0][0];++i) {
        Sites s(cls->clusinfo[n][i]+1,cls->clusinfo[n][i]+n+1);std::sort(s.begin(),s.end());
        index.emplace(s,sites.size());sites.push_back(s);bundleIds.push_back(owners(s,gs,own));
        powers.push_back(cls->clusinfo[n][i][0]);
    }
    raw.resize(sites.size()*size_t(c->nstep));dimensions.resize(sites.size());
}

void ClusterContributions::capture(const int* cluster,int n,int dimension,const MatrixXcd* values) {
    if(!config->cluster_contributions)return;
    const double start=MPI_Wtime();
    Sites s(cluster,cluster+n);std::sort(s.begin(),s.end());const size_t i=index.at(s);
    for(int t=0;t<config->nstep;++t)raw[i*config->nstep+t]=values[t](0,0);
    dimensions[i]=dimension;captureSeconds+=MPI_Wtime()-start;
}

void ClusterContributions::finish(const MatrixXcd* zero,int state,const Output* output) {
    if(!config->cluster_contributions)return;
    const double start=MPI_Wtime();const size_t nt=config->nstep;
    for(size_t offset=0;offset<raw.size();) {
        const int count=int(std::min(raw.size()-offset,size_t(INT_MAX)));
        MPI_Reduce(rank==0 ? MPI_IN_PLACE : raw.data()+offset,raw.data()+offset,count,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);
        offset+=count;
    }
    MPI_Reduce(rank==0 ? MPI_IN_PLACE : dimensions.data(),dimensions.data(),int(dimensions.size()),MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0)try {
        std::ofstream out(std::string(output->outfile)+".contributions.state-"+std::to_string(state)+".tsv");
        out<<std::setprecision(17);
        out<<"# schema=ccex-contributions-v1 full_time_window=1 nstep="<<nt<<" deltat_ms="<<config->deltat
           <<" method="<<config->method<<" seed="<<config->seed<<" normalization="
           <<(!strcasecmp(config->method,"gcce") ? "factor_t0" : "none")<<"\n";
        out<<"sites\tbundles\tproduction_power\tdimension\tmax_abs_log_tilde\tundefined_points\tworst_step\n";
        std::map<Sites,size_t> byBundle;
        for(size_t i=0;i<sites.size();++i)byBundle.emplace(bundleIds[i],i);
        // Every proper bundle subset has fewer sites and occurs earlier.
        for(size_t i=0;i<sites.size();++i) {
            const auto atZero=clusterLog(raw[i*nt])-clusterLog(zero[0](0,0));
            for(size_t t=0;t<nt;++t) {
                auto value=clusterLog(raw[i*nt+t])-clusterLog(zero[t](0,0));
                if(!strcasecmp(config->method,"gcce"))value-=atZero;
                raw[i*nt+t]=value;
            }
            subsets(bundleIds[i],[&](const Sites& sub){
                if(sub.size()==bundleIds[i].size())return;
                const size_t j=byBundle.at(sub);
                for(size_t t=0;t<nt;++t)raw[i*nt+t]-=raw[j*nt+t];
            });
            double maximum=0;size_t worst=0,undefined=0;
            for(size_t t=0;t<nt;++t) {
                auto& z=raw[i*nt+t];z.imag(std::remainder(z.imag(),2*std::acos(-1.)));
                if(!std::isfinite(z.real()) || !std::isfinite(z.imag())){++undefined;continue;}
                const double magnitude=std::abs(z);
                if(magnitude>maximum){maximum=magnitude;worst=t;}
            }
            list(out,sites[i]);out<<'\t';list(out,bundleIds[i]);out<<'\t'<<powers[i]<<'\t'<<dimensions[i]<<'\t';
            if(undefined)out<<"nan";else out<<maximum;
            out<<'\t'<<undefined<<'\t'<<worst<<'\n';
        }
        out.close();if(!out)throw std::runtime_error("cannot write complete contribution table");
        if(state==(config->nstate ? 1 : 0)) {
            std::ofstream metadata(std::string(output->outfile)+".contributions.sites.tsv");
            metadata<<std::setprecision(17)<<"site\tname\tspin\tgyro\tx_A\ty_A\tz_A\n";
            for(int s=0;s<bath->nspin;++s) {
                const auto xyz=BathArray_getBath_i_xyz(bath,s);
                metadata<<s<<'\t'<<BathArray_getBath_i_name(bath,s)<<'\t'<<BathArray_getBath_i_spin(bath,s)<<'\t'
                        <<BathArray_getBath_i_gyro(bath,s)<<'\t'<<xyz[0]<<'\t'<<xyz[1]<<'\t'<<xyz[2]<<'\n';
            }
            metadata.close();if(!metadata)throw std::runtime_error("cannot write contribution site identities");
        }
    }catch(const std::exception& e){fail(e.what());}
    double local=MPI_Wtime()-start+captureSeconds,maximum=0;
    MPI_Reduce(&local,&maximum,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if(rank==0)printf("[cluster-contributions] state=%d rows=%zu capture_reduce_write_max_s=%.9f raw_bytes_per_rank=%zu extra_solver_time_excluded=1\n",
                       state,sites.size(),maximum,raw.size()*sizeof(doublec));
}
