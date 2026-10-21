// Differential checks against the previous production distance/assignment path.
// The old routines remain compiled in cluster_pcce.cpp; no MPI runtime is needed.
#include "clustering/cluster_pcce.h"
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>

int rank=0, nprocess=1, seed=19;
bool verbosity=false;
static int comparisons=0;
using Clock=std::chrono::steady_clock;
static double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now()-start).count();
}

static void check(std::vector<Point> spins, std::vector<Point> centers, bool iterate) {
    const int n=spins.size(), k=centers.size();
    std::vector<double> storage(size_t(n)*k);
    std::vector<double*> rows(n);
    for(int i=0;i<n;++i)rows[i]=storage.data()+size_t(i)*k;
    std::vector<int> old(n), now(n), grid(n), counts(k), active(k);
    PcceTrialWorkspace workspace(n,k);
    auto other=centers;
    for(int it=0;it<(iterate?80:1);++it) {
        calc_dist_to_centers(spins.data(),centers.data(),rows.data(),n,k);
        assignCentroids(old.data(),rows.data(),n,k);
        assignCentroidsFromPoints(now.data(),spins.data(),other.data(),n,k,counts.data(),active.data());
        assignCentroidsGrid(grid.data(),spins.data(),other.data(),n,k,workspace);
        assert(old==now && grid==now); ++comparisons;
        for(int c=0;c<k;++c)assert(counts[c]==n/k);
        if(!iterate)break;
        const bool a=updateCentroids(spins.data(),old.data(),centers.data(),n,k);
        const bool b=updateCentroids(spins.data(),now.data(),other.data(),n,k);
        assert(a==b && !memcmp(centers.data(),other.data(),k*sizeof(Point)));
        double x=calculateInertia(spins.data(),centers.data(),old.data(),n,k);
        double y=calculateInertia(spins.data(),other.data(),now.data(),n,k);
        assert(!memcmp(&x,&y,sizeof(double)));
        if(a)break;
    }
}

static void benchmark(const char* path, int capacity, int trials) {
    std::ifstream input(path); assert(input);
    std::string line; std::getline(input,line);
    std::vector<Point> original;
    while(std::getline(input,line)) {
        std::istringstream row(line); Point p;
        if(row>>p.x>>p.y>>p.z)original.push_back(p);
    }
    const int k=original.size()/capacity, n=k*capacity;
    original.resize(n);assert(k>0);
    std::vector<double> storage(size_t(n)*k);
    std::vector<double*> rows(n);
    for(int i=0;i<n;++i)rows[i]=storage.data()+size_t(i)*k;
    std::vector<int> assignment(n), shuffle(n), counts(k), active(k);
    std::vector<Point> centers(k);
    std::vector<int> reference;
    PcceTrialWorkspace workspace(n,k);
    for(int mode=0;mode<3;++mode) {
        double init=0,assign=0,update=0,score=0;int iterations=0;
        std::vector<int> result;
        auto total=Clock::now();
        for(int t=0;t<trials;++t) {
            PcceRandom rng(73,t);auto spins=original;
            for(int i=0;i<n;++i)shuffle[i]=i;
            auto start=Clock::now();
            shuffle_spins(spins.data(),shuffle.data(),n,&rng);
            choose_initial_method(true,spins.data(),centers.data(),n,k,&rng);
            init+=seconds(start);
            for(int it=0;it<1000;++it) {
                start=Clock::now();
                if(mode==0) {
                    calc_dist_to_centers(spins.data(),centers.data(),rows.data(),n,k);
                    assignCentroids(assignment.data(),rows.data(),n,k);
                } else if(mode==1) assignCentroidsFromPoints(assignment.data(),spins.data(),centers.data(),n,k,counts.data(),active.data());
                else assignCentroidsGrid(assignment.data(),spins.data(),centers.data(),n,k,workspace);
                assign+=seconds(start); start=Clock::now();
                bool done=updateCentroids(spins.data(),assignment.data(),centers.data(),n,k);
                update+=seconds(start);++iterations;
                if(done)break;
            }
            start=Clock::now();
            // Force scores to be evaluated without folding them into timing output.
            volatile double inertia=calculateInertia(spins.data(),centers.data(),assignment.data(),n,k);
            volatile double sil=silhouetteCoefficient(spins.data(),assignment.data(),centers.data(),n,k);
            (void)inertia;(void)sil; score+=seconds(start);
            result.insert(result.end(),assignment.begin(),assignment.end());
        }
        if(mode==0)reference=result;else assert(reference==result);
        printf("mode=%s sites=%d sK=%d trials=%d iterations=%d total_s=%.6f init_s=%.6f assign_s=%.6f update_s=%.6f score_s=%.6f\n",
               mode==2?"grid":mode?"active":"original",n,capacity,trials,iterations,seconds(total),init,assign,update,score);
        fflush(stdout);
    }
}

int main(int argc,char** argv) {
    if(argc==4) { benchmark(argv[1],atoi(argv[2]),atoi(argv[3]));return 0; }
    for(int k:{1,2,3,7,32,64,131,700})for(int size:{1,2,3,7})for(int seed:{1,19,73,123}) {
        const int n=k*size;
        for(int shape=0;shape<5;++shape) {
            PcceRandom rng(seed,shape);
            std::vector<Point> spins(n),centers(k);
            for(int i=0;i<n;++i) {
                if(shape==0)spins[i]={double(i%5),double(i%3),0}; // ties, duplicate sites
                else if(shape==1)spins[i]={rng.unit()*80-40,rng.unit()*80-40,rng.unit()*80-40};
                else if(shape==2)spins[i]={double(i),std::nextafter(double(i),INFINITY),double(-i)};
                else if(shape==3)spins[i]={double(i%2)*1e150,double(i%3)*1e-150,0};
                else spins[i]={double(rng.index(9))*0.25,double(rng.index(9))*0.25,double(rng.index(9))*0.25}; // lattice: ties
            }
            for(bool plus:{false,true}) {
                choose_initial_method(plus,spins.data(),centers.data(),n,k,&rng);
                check(spins,centers,true);
            }
        }
    }
    for(double extreme:{0.0,-0.0,1e308,std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()}) {
        std::vector<Point> spins(24,Point{extreme,0,0}),centers(8,Point{0,0,0});
        check(spins,centers,false);
    }
    // P2-01 of review 20261012-233058-452278: extent > 0 but extent/side rounds to a cell width of 0 (the grid path
    // converted Inf to int); 64 centers, 128 spins, x alternating 0 and the smallest subnormal
    {
        const double tiny=std::numeric_limits<double>::denorm_min();
        std::vector<Point> spins(128),centers(64);
        for(int i=0;i<128;++i)spins[i]={i%2?tiny:0.0,0,0};
        for(int c=0;c<64;++c)centers[c]={c%2?tiny:0.0,0,0};
        check(spins,centers,false);
    }
    printf("PASS %d assignment/iteration comparisons; exact centers, inertia and capacity\n",comparisons);
}
