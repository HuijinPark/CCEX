#include "calculation/gpu/gpu_backend.h"
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
int main(){
    char description[512];if(!gpuInitialize(0,1,description,sizeof(description))){std::cerr<<description<<'\n';return 77;}
    try {
        std::complex<double> h[]={1.,.2,.2,2.},psi[]={1.,0.};
        bool rejected=false;
        try {
            std::unique_ptr<GpuGcce,void(*)(GpuGcce*)> w(
                gpuGcceCreate(2,2,true,1,h,psi,0,nullptr,true,0),gpuGcceDestroy);
        }catch(const std::runtime_error& e){
            rejected=std::string(e.what()).find("refinement did not converge")!=std::string::npos;
        }
        if(!rejected)throw std::runtime_error("exhausted GPU refinement did not reject");
        {
            std::unique_ptr<GpuGcce,void(*)(GpuGcce*)> w(
                gpuGcceCreate(2,2,true,1,h,psi,0,nullptr,false,0),gpuGcceDestroy);
            std::complex<double> reduced[4];double time=.3,delay=1;
            gpuGcceEvolve(w.get(),1,&time,&delay,reduced,nullptr);
            if(std::abs(reduced[0]+reduced[3]-std::complex<double>(1,0))>1e-12)
                throw std::runtime_error("explicit opt-out after failure returned wrong trace");
        }
        std::cout<<"PASS strict GPU refinement failure and explicit opt-out recovery\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';gpuFinalize();return 1;}
    gpuFinalize();
}
