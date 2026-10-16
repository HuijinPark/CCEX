/**
 * @file cluster_log.h
 * @brief The complex logarithm of a cluster coherence, NaN for zero or non-finite values.
 */
#ifndef CCEX_CLUSTER_LOG_H
#define CCEX_CLUSTER_LOG_H
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
// Zero/nonfinite raw coherence has no finite log. Never floor or clamp it.
inline std::complex<double> clusterLog(std::complex<double> z) {
    const double scale=std::max(std::abs(z.real()),std::abs(z.imag()));
    if(scale==0 || !std::isfinite(z.real()) || !std::isfinite(z.imag()))
        return {std::numeric_limits<double>::quiet_NaN(),0};
    return {std::log(scale)+std::log(std::hypot(z.real()/scale,z.imag()/scale)),std::arg(z)};
}
#endif
