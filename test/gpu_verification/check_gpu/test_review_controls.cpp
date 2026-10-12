#include "calculation/gpu/eigen_refinement.h"
#include "calculation/gpu/gpu_batch_memory.h"
#include <iostream>
#include <stdexcept>
int main(){
    try {
        // A 256-dimensional state fits eight works into 128 MiB. Treating it
        // as a density matrix would reject even two works at this time batch.
        for(bool gcce:{false,true}){
            auto vector=gpuClusterMemoryEstimate(256,2,64,4,gcce,true);
            auto matrix=gpuClusterMemoryEstimate(256,2,64,4,gcce,false);
            if(vector*8>128.L*1024*1024 || matrix*2<128.L*1024*1024)
                throw std::runtime_error("state-aware batch admission regressed");
        }
        bool invalid=false;
        try {gpuClusterMemoryEstimate(0,2,64,4,true,true);}catch(const std::invalid_argument&){invalid=true;}
        if(!invalid)throw std::runtime_error("invalid estimate accepted");
        Eigen::MatrixXcd h(3,3);h<<1.,.2,.3,.2,2.,.4,.3,.4,3.;
        Eigen::MatrixXcd vectors=Eigen::MatrixXcd::Identity(3,3);
        Eigen::VectorXd values=h.diagonal().real(),corrections;
        bool failed=false;
        try {refineEigensystem(h,vectors,values,corrections,0);}catch(const std::runtime_error&){failed=true;}
        if(!failed)throw std::runtime_error("exhausted refinement silently accepted");
        vectors.setIdentity();values=h.diagonal().real();
        refineEigensystem(h,vectors,values,corrections,64);
        Eigen::MatrixXcd diagonal=(values+corrections).asDiagonal();
        if((h*vectors-vectors*diagonal).norm()>1e-12 ||
           (vectors.adjoint()*vectors-Eigen::MatrixXcd::Identity(3,3)).norm()>1e-12)
            throw std::runtime_error("refinement failed eigenpair/orthogonality oracle");
        std::cout<<"PASS vector memory admission, invalid inputs, strict refinement failure/convergence\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
