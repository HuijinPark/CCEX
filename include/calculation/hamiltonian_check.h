/**
 * @file hamiltonian_check.h
 * @brief The finiteness and hermiticity check of a Hamiltonian before it goes to the GPU.
 */
#ifndef CCEX_HAMILTONIAN_CHECK_H
#define CCEX_HAMILTONIAN_CHECK_H
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <complex>

// What the CUDA cCCE path asks of a Hamiltonian before it goes to the card: every
// entry finite, and max|h - h^dagger| <= 1e-12 * max|h|. It is the verdict of
//     h.allFinite() && !((h.adjoint()-h).cwiseAbs().maxCoeff() > 1e-12*h.cwiseAbs().maxCoeff())
// entry for entry - the same std::abs on the same numbers - without the three
// full-size temporaries and the four walks over the matrix that expression costs
// (25 s of a 143 s setup at dimension 16384). |h(i,j) - conj(h(j,i))| is the same
// number for (i,j) and (j,i), so the upper triangle is enough; it is walked in
// tiles because h(j,i) runs across the columns; and an entry pair that is zero, or
// a difference that is zero, leaves both maxima where they are and is not measured.
inline bool finiteAndHermitian(const Eigen::MatrixXcd& h) {
    const Eigen::Index n=h.rows(),tile=64;
    double skew=0,scale=0;
    for(Eigen::Index jt=0;jt<n;jt+=tile)for(Eigen::Index it=0;it<=jt;it+=tile) {
        const Eigen::Index jend=std::min(n,jt+tile),iend=std::min(n,it+tile);
        for(Eigen::Index j=jt;j<jend;++j)for(Eigen::Index i=it;i<iend && i<=j;++i) {
            const std::complex<double> upper=h(i,j),lower=h(j,i);
            if(upper==0.0 && lower==0.0)continue;
            if(!std::isfinite(upper.real()) || !std::isfinite(upper.imag()) ||
               !std::isfinite(lower.real()) || !std::isfinite(lower.imag()))return false;
            scale=std::max(scale,std::max(std::abs(upper),std::abs(lower)));
            const std::complex<double> difference=std::conj(lower)-upper;
            if(difference!=0.0)skew=std::max(skew,std::abs(difference));
        }
    }
    return !(skew>1e-12*scale);
}
#endif
