#include "clustering/bundle/cluster_log.h"
#include <cassert>
#include <iostream>
int main() {
    const auto hi=clusterLog({1e308,1e308}),lo=clusterLog({1e-308,-1e-308});
    assert(std::isfinite(hi.real()) && std::isfinite(lo.real()));
    assert(std::abs(hi.real()-(std::log(1e308)+std::log(2.)/2))<1e-12);
    assert(std::abs(lo.real()-(std::log(1e-308)+std::log(2.)/2))<1e-12);
    assert(std::abs(hi.imag()-std::acos(-1.)/4)<1e-15);
    assert(std::abs(lo.imag()+std::acos(-1.)/4)<1e-15);
    assert(!std::isfinite(clusterLog({0,0}).real()));
    assert(!std::isfinite(clusterLog({std::numeric_limits<double>::infinity(),1}).real()));
    assert(!std::isfinite(clusterLog({1,std::numeric_limits<double>::quiet_NaN()}).real()));
    assert(std::abs(clusterLog(std::polar(.25,2.5))-std::complex<double>(std::log(.25),2.5))<1e-15);
    std::cout<<"PASS cluster log extremes, zero, nonfinite and phase\n";
}
