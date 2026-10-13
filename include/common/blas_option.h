#pragma once
// A46: Intel MKL for the large clusters (do_compile.sh and CMake turn it on when they find MKL; CCEX_BLAS=off turns it off).
// Without CCEX_BLAS_DIRECT nothing here changes the build: CcexEigensolver is Eigen's solver, as before.
// With it, the large clusters go to the library (src/common/blas_option.cpp): the cCCE ensemble Hahn step's
// X = (M diag(d)) M^dagger by zgemm from dimension ccexBlas::gemmMin on, and the cCCE/gCCE cluster eigensolve by zheevd
// from ccexBlas::eigMin on. The thresholds are where MKL was faster and the h-BN runs gained (A46 step 0, idm job 17337:
// zgemm pays from 12, zheevd from 27; 21 and 63 ran as fast and changed fewer outputs). Below them Eigen runs, as before.
#include <Eigen/Dense>

#ifndef CCEX_BLAS_DIRECT
using CcexEigensolver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>;
#else
#include <complex>
#include <vector>

namespace ccexBlas {
constexpr int gemmMin = 21;
constexpr int eigMin = 63;
// X = (M diag(d)) M^dagger by zgemm; scratch holds M diag(d)
void hahnProduct(const Eigen::MatrixXcd& M, const Eigen::VectorXcd& d, Eigen::MatrixXcd& scratch, Eigen::MatrixXcd& X);
}

// Eigen's SelfAdjointEigenSolver below ccexBlas::eigMin, LAPACK zheevd (lower triangle, as Eigen reads) from it on
class CcexEigensolver {
public:
    void compute(const Eigen::MatrixXcd& H);
    const Eigen::MatrixXcd& eigenvectors() const { return lapack_ ? V_ : es_.eigenvectors(); }
    const Eigen::VectorXd& eigenvalues() const { return lapack_ ? w_ : es_.eigenvalues(); }
    Eigen::ComputationInfo info() const { return lapack_ ? info_ : es_.info(); }
private:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es_;
    bool lapack_ = false;
    Eigen::MatrixXcd V_;
    Eigen::VectorXd w_;
    Eigen::ComputationInfo info_ = Eigen::Success;
    int wsDim_ = -1;                                   // the dimension the workspace below was sized for
    std::vector<std::complex<double>> work_;
    std::vector<double> rwork_;
    std::vector<int> iwork_;
};
#endif
