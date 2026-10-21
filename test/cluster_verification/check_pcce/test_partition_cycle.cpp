// Compare the optimized trial against the previous, unskipped production loop.
#include "clustering/cluster_pcce.h"
#include <array>
#include <cassert>
#include <climits>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

int rank=0,nprocess=1,seed=19;
bool verbosity=false;
static int comparisons=0,skipped_cases=0,tracker_cases=0;

static void trackerChecks() {
    // An independent finite-state oracle exercises transients, long/odd periods,
    // boundary remainders and INT_MAX without storing the entire trajectory.
    for(int transient:{0,1,2,7,31,257})for(int period:{1,2,3,6,8,37,1031})
    for(int limit:{1,2,3,7,8,31,1000,1001,8192,INT_MAX})for(int encoding:{0,1,2}) {
        std::array<double,2> centers{};std::array<int,2> assigned{};
        PcceCycleTracker tracker(sizeof(centers),sizeof(assigned));
        for(int repeat=0;repeat<2;++repeat) {
            tracker.reset();int iteration=0,detected=0,skipped=0,executed=0;
            auto state=[&](int i){return i<transient?i:transient+(i-transient)%period;};
            int current=0;
            while(iteration<limit) {
                current=state(iteration);++iteration;++executed;
                // Vary only one component at a time. Adjacent binary64 values
                // and identical centers with different assignments must differ.
                centers[0]=encoding==1?0.0:1.0+double(current)*0x1p-52;
                centers[1]=encoding==2 && current%2 ? -0.0 : 0.0;
                assigned[0]=encoding==1?current:0;
                if(!detected) {
                    detected=int(tracker.observe(centers.data(),assigned.data()));
                    if(detected) {
                        assert(detected==period);
                        skipped=((limit-iteration)/detected)*detected;
                        iteration+=skipped;
                    }
                }
            }
            assert(current==state(limit-1));assert(executed+skipped==limit);
            assert(skipped==0 || skipped%detected==0);
            ++tracker_cases;
        }
    }
    // Equal floating-point values with distinct signs still have different bits.
    double center=0.0;int assigned=0;PcceCycleTracker sign(sizeof(center),sizeof(assigned));
    assert(!sign.observe(&center,&assigned));center=-0.0;
    assert(!sign.observe(&center,&assigned));center=0.0;
    assert(!sign.observe(&center,&assigned));center=-0.0;
    assert(sign.observe(&center,&assigned)==2);
}

static PcceTrialResult reference(Point* spins,Point* centers,int* assigned,
                                 int n,int k,int limit) {
    std::vector<int> counts(k),active(k);PcceTrialResult result;
    while(!result.converged && result.iterations<limit) {
        assignCentroidsFromPoints(assigned,spins,centers,n,k,counts.data(),active.data());
        result.converged=updateCentroids(spins,assigned,centers,n,k);++result.iterations;
    }
    return result;
}

static void compare(std::vector<Point> spins,std::vector<Point> initial,int limit,
                    PcceTrialWorkspace& workspace,bool large=false) {
    const int n=spins.size(),k=initial.size();auto old=initial,now=initial;
    std::vector<int> before(n),after(n);
    const auto expected=reference(spins.data(),old.data(),before.data(),n,k,limit);
    const auto result=pcceOptimizeTrial(spins.data(),now.data(),after.data(),n,k,limit,workspace);
    assert(expected.iterations==result.iterations && expected.converged==result.converged);
    assert(result.executed+result.skipped==result.iterations);
    assert(result.skipped==0 || (!result.converged && result.cycle_period>0 && result.skipped%result.cycle_period==0));
    assert(before==after && !memcmp(old.data(),now.data(),k*sizeof(Point)));
    double a=calculateInertia(spins.data(),old.data(),before.data(),n,k);
    double b=calculateInertia(spins.data(),now.data(),after.data(),n,k);
    assert(!memcmp(&a,&b,sizeof(double)));
    a=silhouetteCoefficient(spins.data(),before.data(),old.data(),n,k);
    b=silhouetteCoefficient(spins.data(),after.data(),now.data(),n,k);
    assert(!memcmp(&a,&b,sizeof(double)));
    ++comparisons;if(result.skipped)++skipped_cases;
    if(large && result.cycle_period && !result.converged) {
        // Once the proven cycle has been traversed by the unskipped reference,
        // use the matching phase to check the production arithmetic at INT_MAX.
        int short_limit=limit+(INT_MAX-limit)%result.cycle_period;
        old=initial;now=initial;
        auto end=reference(spins.data(),old.data(),before.data(),n,k,short_limit);
        auto big=pcceOptimizeTrial(spins.data(),now.data(),after.data(),n,k,INT_MAX,workspace);
        assert(!end.converged && !big.converged && big.iterations==INT_MAX);
        assert(big.executed+big.skipped==INT_MAX);
        assert(before==after && !memcmp(old.data(),now.data(),k*sizeof(Point)));
        ++comparisons;
    }
}

static void fixture(const char* path,int size,int initial_seed,int trials,int limit) {
    std::ifstream file(path);assert(file);std::string line;std::getline(file,line);
    std::vector<Point> original;
    while(std::getline(file,line)){std::istringstream row(line);Point p;if(row>>p.x>>p.y>>p.z)original.push_back(p);}
    const int k=original.size()/size,n=k*size;original.resize(n);assert(k>0);
    PcceTrialWorkspace workspace(n,k);
    for(int trial=0;trial<trials;++trial) {
        PcceRandom rng(initial_seed,trial);auto spins=original;
        std::vector<Point> centers(k);std::vector<int> shuffle(n);
        for(int i=0;i<n;++i)shuffle[i]=i;
        shuffle_spins(spins.data(),shuffle.data(),n,&rng);
        choose_initial_method(true,spins.data(),centers.data(),n,k,&rng);
        for(int end:{limit,limit+1})compare(spins,centers,end,workspace,trial<3);
    }
}

int main(int argc,char** argv) {
    if(argc==6)fixture(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));
    else {
        trackerChecks();
        for(int k:{1,2,3,7,31})for(int size:{1,2,3,7}) {
            const int n=k*size;PcceTrialWorkspace workspace(n,k);
            for(int seed:{1,73,123})for(int shape:{0,1,2})for(bool plus:{false,true}) {
                PcceRandom rng(seed,shape);std::vector<Point> spins(n),centers(k);
                for(int i=0;i<n;++i) {
                    if(shape==0)spins[i]={double(i%3),double(i%5),0};
                    else if(shape==1)spins[i]={rng.unit()*80-40,rng.unit()*80-40,rng.unit()*80-40};
                    else spins[i]={double(i%2)*1e150,double(i%3)*1e-150,0};
                }
                choose_initial_method(plus,spins.data(),centers.data(),n,k,&rng);
                for(int limit:{1,2,3,7,31,32,99,100})compare(spins,centers,limit,workspace);
            }
        }
        assert(skipped_cases>0);
    }
    printf("PASS trial comparisons=%d skipped_cases=%d finite_state_checks=%d; exact centers, assignments, scores and final phase\n",comparisons,skipped_cases,tracker_cases);
}
