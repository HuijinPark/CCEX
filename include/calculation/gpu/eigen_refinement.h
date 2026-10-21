/**
 * @file eigen_refinement.h
 * @brief The compensated refinement of eigensystems and the phases computed from them.
 */
#ifndef CCEX_EIGEN_REFINEMENT_H
#define CCEX_EIGEN_REFINEMENT_H

#include <Eigen/Dense>
#include <cmath>
#include "compensated.h"
#include "accelerator_defaults.h"
#include <limits>
#include <stdexcept>

inline void refinedPhases(const Eigen::VectorXd& values,
        const Eigen::VectorXd& corrections, double tau, Eigen::VectorXcd& phase){
    phase.resize(values.size());
    for (int k=0;k<values.size();++k){
        const double angle=values[k]*tau;
        const double tail=std::fma(values[k],tau,-angle)+corrections[k]*tau;
        // |tail| < 2^-27 (always, unless the corrections are large): cos(tail) rounds to 1 and sin(tail) to tail, which
        // is what the library returns there, so the two calls are skipped and the product below is unchanged (A16)
        double c=1.0,s=tail;
        if(!(std::fabs(tail)<0x1p-27)){c=std::cos(tail);s=std::sin(tail);}
        phase[k]=std::complex<double>(std::cos(angle),-std::sin(angle))*
                 std::complex<double>(c,-s);
    }
}

// Refine the complete eigensystem in its almost-diagonal basis. Correcting only
// eigenvalues leaves small off-diagonal residuals in nearly degenerate bath
// subspaces. The residual is formed with compensated products and diagonal
// subtraction; Jacobi rotations retain two-double diagonal energies.
inline void refineEigensystem(const Eigen::MatrixXcd& h,Eigen::MatrixXcd& vectors,
        Eigen::VectorXd& values,Eigen::VectorXd& corrections,int maxSweeps=ccex::refinementSweeps){
    using Z=std::complex<double>;
    const int d=values.size();
    if(d<1 || maxSweeps<0 || h.rows()!=d || h.cols()!=d || vectors.rows()!=d || vectors.cols()!=d)
        throw std::invalid_argument("invalid gCCE refinement dimensions or sweep limit");
    Eigen::MatrixXcd residual(d,d);
    for(int col=0;col<d;++col)for(int row=0;row<d;++row){
        DoublePair real{0,0},imag{0,0};
        for(int k=0;k<d;++k){
            const Z v=vectors(k,col);
            if(k==row){
                const auto diagonal=pairDifference({h(row,row).real(),0},{values[col],0});
                real=pairAdd(real,pairScale(diagonal,v.real()));
                imag=pairAdd(imag,pairScale(diagonal,v.imag()));
            }else{
                const Z a=row>k ? h(row,k) : std::conj(h(k,row));
                if(a==Z(0,0))continue;
                real=pairAdd(real,pairAdd(pairProduct(a.real(),v.real()),pairProduct(-a.imag(),v.imag())));
                imag=pairAdd(imag,pairAdd(pairProduct(a.real(),v.imag()),pairProduct(a.imag(),v.real())));
            }
        }
        residual(row,col)={real.high+real.low,imag.high+imag.low};
    }
    Eigen::MatrixXcd projected=vectors.adjoint()*residual;
    corrections.resize(d);
    for(int i=0;i<d;++i){
        const auto value=pairAdd({values[i],0},{projected(i,i).real(),0});
        values[i]=value.high;corrections[i]=value.low;projected(i,i)=0;
        for(int j=0;j<i;++j){
            const Z entry=.5*(projected(i,j)+std::conj(projected(j,i)));
            projected(i,j)=entry;projected(j,i)=std::conj(entry);
        }
    }
    const double epsilon=std::numeric_limits<double>::epsilon();
    const double tolerance=64*epsilon*epsilon*values.cwiseAbs().maxCoeff();
    for(int sweep=0;sweep<maxSweeps;++sweep){
        for(int p=0;p<d;++p)for(int q=p+1;q<d;++q){
            const Z off=projected(p,q);const double magnitude=std::abs(off);
            if(magnitude<=tolerance)continue;
            const auto delta=pairDifference({values[q],corrections[q]},{values[p],corrections[p]});
            const double gap=delta.high+delta.low;
            const double t=gap==0 ? 1 : std::copysign(2*magnitude/(std::abs(gap)+std::hypot(gap,2*magnitude)),gap);
            const double c=1/std::sqrt(1+t*t);const Z s=off*(t*c/magnitude);
            const auto shift=pairProduct(t,magnitude);
            const auto vp=pairDifference({values[p],corrections[p]},shift);
            const auto vq=pairAdd({values[q],corrections[q]},shift);
            values[p]=vp.high;corrections[p]=vp.low;values[q]=vq.high;corrections[q]=vq.low;
            for(int k=0;k<d;++k){
                if(k!=p && k!=q){
                    const Z a=projected(k,p),b=projected(k,q);
                    const Z ap=c*a-std::conj(s)*b,aq=s*a+c*b;
                    projected(k,p)=ap;projected(p,k)=std::conj(ap);
                    projected(k,q)=aq;projected(q,k)=std::conj(aq);
                }
                const Z a=vectors(k,p),b=vectors(k,q);
                vectors(k,p)=c*a-std::conj(s)*b;vectors(k,q)=s*a+c*b;
            }
            projected(p,q)=projected(q,p)=0;
        }
        if(projected.cwiseAbs().maxCoeff()<=tolerance)return;
    }
    throw std::runtime_error("gCCE eigensystem refinement did not converge");
}
#endif
