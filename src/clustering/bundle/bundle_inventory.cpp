#include "../../../include/clustering/bundle/bundle_cluster.h"
#include "../../../include/clustering/bundle/bundle_text.h"
#include "../../../include/input/json.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

using adaptive_bundle::jsonString;

void writeClusterInventory(Cluster* cls,BathArray* ba,QubitArray* qa,DefectArray* dfa,
                           const Config* c,const char* outfile,int loadedSites) {
    try {
        const std::string path=std::string(outfile)+".clusters.json";
        std::ofstream f(path);
        if(!f)throw std::runtime_error("cannot open inventory output");
        f<<std::setprecision(17)<<"{\"schema\":1,\"dynamics_performed\":false,\"method\":";
        jsonString(f,c->method);f<<",\"cluster_mode\":";jsonString(f,c->cluster_mode);
        f<<",\"input_order\":"<<c->order<<",\"max_site_order\":"<<cls->order
         <<",\"loaded_sites\":"<<loadedSites<<",\"retained_sites\":"<<ba->nspin
         <<",\"qubit_dimension\":"<<(!strcasecmp(c->method,"gcce") ? QubitArray_dim(qa) : 1)
         <<",\"zero_power\":"<<cls->clusinfo[0][0][0]<<",\"sites\":[";
        for(int s=0;s<ba->nspin;++s) {
            auto* spin=ba->bath[s];
            f<<(s ? "," : "")<<"{\"id\":"<<s<<",\"name\":";jsonString(f,spin->name);
            f<<",\"xyz\":["<<spin->xyz[0]<<','<<spin->xyz[1]<<','<<spin->xyz[2]
             <<"],\"spin\":"<<spin->spin<<",\"component_dimensions\":["<<BathSpin_dim(spin);
            if(dfa->naddspins)for(int j=0;j<dfa->naddspins[s];++j)
                f<<','<<BathSpin_dim(dfa->subbath[s][j]);
            f<<']';
            if(dfa->paxes)f<<",\"defect_axis\":"<<dfa->paxes[s];
            f<<'}';
        }
        f<<"],\"bundles\":";
        // pCCE bundles are available in its partition cache; zero-power rows alone
        // cannot reconstruct all bundles. Never invent a partition from these rows.
        if(!strcasecmp(c->method,"pcce"))f<<"null";
        else {
            auto bundles=adaptiveInstalledBundles();
            if(strcasecmp(c->cluster_mode,"adaptive")) {
                bundles.clear();for(int s=0;s<ba->nspin;++s)bundles.push_back({s});
            }
            f<<'[';
            for(size_t g=0;g<bundles.size();++g) {f<<(g ? "," : "");adaptive_bundle::jsonArray(f,bundles[g]);}
            f<<']';
        }
        f<<",\"clusters\":[";size_t count=0;
        for(int n=1;n<=cls->order;++n)for(int i=1;i<cls->clusinfo[n][0][0];++i) {
            const int* row=cls->clusinfo[n][i];
            f<<(count++ ? "," : "")<<"{\"power\":"<<row[0]<<",\"sites\":[";
            for(int j=1;j<=n;++j)f<<(j>1 ? "," : "")<<row[j];
            f<<"]}";
        }
        f<<"]}\n";f.close();
        if(!f)throw std::runtime_error("cannot finish inventory output");
        printf("[cluster-plan] loaded_sites=%d retained_sites=%d clusters=%zu dynamics=0 path=%s\n",
               loadedSites,ba->nspin,count,path.c_str());
    } catch(const std::exception& e) {
        fprintf(stderr,"Error: cluster inventory: %s\n",e.what());MPI_Abort(MPI_COMM_WORLD,1);
    }
}
