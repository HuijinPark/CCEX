// A62: the column-block one-state cCCE evolution (cceEvolveColumns, input block column_evolution) against the path it
// replaces (evolveSingleState through cceEvolve), on random Hermitian Hamiltonians of bath dimension 300 and 600 above the
// threshold: exact and recursive phases; Ramsey, Hahn echo, CPMG-8, a Manual sequence with four different intervals and
// one with repeated ones; a block size that does not divide nstep, one larger than nstep and 1. The two sum the same
// products in another order, so every output is held to 1e-12 relative to the larger magnitude (1e-15 absolute where both
// are below 1e-3, which these cases do not reach; the count is printed). Also: with the block on but the dimension below
// min_dim, and with it off, the outputs are the bits of the path before A62. Links the objects of main.out but main.o.
#include "calculation/simulator.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>

int rank=0,nprocess=1,seed=19;
bool verbosity=false;

namespace {
struct Sequence { const char* name; std::vector<double> fractions; };
// fractions as Pulse builds them; sequence_indices as assign_sequence_indices sets them
struct Built {
    std::vector<std::vector<double>> rows; std::vector<double*> seq; std::vector<int> idx; Pulse pulse{};
    explicit Built(const Sequence& s){
        const int n = int(s.fractions.size());
        rows.assign(size_t(n),std::vector<double>(3,0.0)); seq.resize(size_t(n)); idx.resize(size_t(n));
        for (int i=0;i<n;i++){ rows[i][2] = s.fractions[i]; seq[i] = rows[i].data(); }
        for (int i=0;i<n;i++){ idx[i] = i; for (int j=0;j<i;j++) if (rows[i][2]==rows[j][2]){ idx[i] = idx[j]; break; } }
        pulse.npulse = n-1; pulse.sequence = seq.data(); pulse.sequence_indices = idx.data();
        strcpy(pulse.pulsename, n==1 ? "Ramsey" : !strcmp(s.name,"Manual") || !strcmp(s.name,"Manual repeated") ? "Manual" : n==2 ? "HahnEcho" : "CPMG");
    }
};
std::vector<doublec> evolve(Config* cnf, Pulse* pulse, CceWorkspace& ws){
    std::vector<MatrixXcd> res(size_t(cnf->nstep));
    cceEvolve(cnf,pulse,res.data(),ws);
    std::vector<doublec> out(res.size());
    for (size_t i=0;i<res.size();i++) out[i] = res[i](0,0);
    return out;
}
}

int main(){
    std::vector<Sequence> sequences = {
        {"Ramsey",{1.0}}, {"Hahn",{0.5,0.5}},
        {"CPMG-8",{1.0/16,1.0/8,1.0/8,1.0/8,1.0/8,1.0/8,1.0/8,1.0/8,1.0/16}},
        {"Manual",{0.1,0.3,0.25,0.35}}, {"Manual repeated",{0.2,0.3,0.2,0.3}}};
    long outputs=0, failed=0, tiny=0, otherBits=0;
    double worst=0;
    for (int D : {300,600}){
        // Hermitian Hamiltonians (rad/us scale of a bath), Hbeta a little away from Halpha, a superposition state
        unsigned long long x = 0x9e3779b97f4a7c15ULL*(unsigned long long)D;
        auto uniform = [&x](double a){ x ^= x << 13; x ^= x >> 7; x ^= x << 17; return a*(2.0*double(x >> 11)/9007199254740992.0-1.0); };
        CceWorkspace ws;
        ws.Halpha.resize(D,D); ws.Hbeta.resize(D,D);
        for (int i=0;i<D;i++){
            ws.Halpha(i,i) = doublec(uniform(30.0),0.0); ws.Hbeta(i,i) = ws.Halpha(i,i)+doublec(uniform(1.0),0.0);
            for (int k=i+1;k<D;k++){
                ws.Halpha(i,k) = doublec(uniform(3.0),uniform(3.0)); ws.Halpha(k,i) = std::conj(ws.Halpha(i,k));
                ws.Hbeta(i,k) = ws.Halpha(i,k)+doublec(uniform(0.1),uniform(0.1)); ws.Hbeta(k,i) = std::conj(ws.Hbeta(i,k));
            }
        }
        MatrixXcd psi(D,1);
        for (int i=0;i<D;i++) psi(i,0) = doublec(uniform(1.0),uniform(1.0));
        psi /= psi.norm();
        cceSetup(ws,&psi);
        Config* cnf = Config_init();
        strcpy(cnf->method,"CCE"); strcpy(cnf->quantity,"coherence");
        cnf->nstate = 1; cnf->deltat = 0.003; cnf->nstep = D==300 ? 77 : 41;
        for (const Sequence& s : sequences){
            if (D==600 && strcmp(s.name,"Hahn") && strcmp(s.name,"CPMG-8")) continue;
            Built b(s);
            std::vector<doublec> exactBefore;
            for (bool recursive : {false,true}){
                strcpy(cnf->evolution_update, recursive ? "recursive" : "exact");
                cnf->column_evolution = false;
                const std::vector<doublec> before = evolve(cnf,&b.pulse,ws);
                if (!recursive) exactBefore = before;
                else {
                    // how far the two phase rules lie apart here: a column path with the other rule would be this far off
                    double d = 0;
                    for (size_t i=0;i<before.size();i++) d = std::max(d,std::abs(before[i]-exactBefore[i])/std::abs(before[i]));
                    printf("  D %d %-15s recursive against exact (the path before): largest relative difference %.2e\n",D,s.name,d);
                }
                // on, but the dimension below min_dim: the bits of the path before
                cnf->column_evolution = true; cnf->column_evolution_min_dim = D+1; cnf->column_evolution_columns = 32;
                const std::vector<doublec> below = evolve(cnf,&b.pulse,ws);
                if (std::memcmp(before.data(),below.data(),before.size()*sizeof(doublec))){
                    printf("FAIL D %d %s %s: below min_dim the outputs are not the bits of the path before\n",
                           D,s.name,recursive ? "recursive" : "exact");
                    ++failed;
                }
                cnf->column_evolution_min_dim = D;
                for (int columns : {32,1,cnf->nstep+5}){
                    cnf->column_evolution_columns = columns;
                    const std::vector<doublec> after = evolve(cnf,&b.pulse,ws);
                    double caseWorst = 0; long bad = 0;
                    for (size_t i=0;i<after.size();i++){
                        ++outputs;
                        if (std::memcmp(&before[i],&after[i],sizeof(doublec))) ++otherBits;
                        const double d = std::abs(before[i]-after[i]), m = std::max(std::abs(before[i]),std::abs(after[i]));
                        if (m < 1e-3){ ++tiny; if (d > 1e-15) ++bad; continue; }
                        caseWorst = std::max(caseWorst,d/m);
                        if (d > 1e-12*m) ++bad;
                    }
                    worst = std::max(worst,caseWorst);
                    printf("  D %d %-15s %-9s columns %3d nstep %d: largest relative difference %.2e%s\n",D,s.name,
                           recursive ? "recursive" : "exact",columns,cnf->nstep,caseWorst,bad ? "  FAIL" : "");
                    failed += bad;
                }
            }
        }
        Config_freeAll(cnf);
    }
    // The phases. Both paths keep the phases of the current step of each distinct interval in ws.phases (up, down and,
    // recursive, one step's factors), so after a run it holds those of the last step; last steps around the recursive
    // anchors (0, 16, 32, 64) and the column-block edges. The relative test above cannot see the phase rule (exact and
    // recursive differ by about 1e-16 there), so the rule is held here:
    //  - exact, and the step factors of recursive: the bits of the path before (sin and cos of the same products);
    //  - recursive at an anchor step: the bits of exact (sin and cos there, not a carried product);
    //  - recursive between anchors: carried - some phase differs from exact in a bit, all within 1e-11 (the docs: about 1e-11 between the two) - and within 1e-14
    //    of the path before. Not bit for bit: the carried product a *= b is a complex multiplication whose fused
    //    multiply-adds the compiler places per translation unit (seen: the same bits in the MKL build of wsl-gnu, 1-ulp
    //    differences in its CCEX_BLAS=off build), as for the fixed-size path (A56). The count is printed.
    // An anchor of 15 instead of 16, or recomputing sin and cos at each block's first column, fails it.
    {
        const int D = 300;
        unsigned long long x = 0x2545f4914f6cdd1dULL;
        auto uniform = [&x](double a){ x ^= x << 13; x ^= x >> 7; x ^= x << 17; return a*(2.0*double(x >> 11)/9007199254740992.0-1.0); };
        CceWorkspace ws;
        ws.M = MatrixXcd::Zero(D,D);
        for (int i=0;i<D;i++) ws.M(i,(7*i+3)%D) = 1.0;          // a permutation (7 and 300 coprime): cheap products
        ws.Madj = ws.M.adjoint();
        // cceEvolve binds the eigenvectors (unused on this path): give the solvers something, as the A42 setup does
        const MatrixXcd identity = MatrixXcd::Identity(D,D); const Eigen::VectorXd zeros = Eigen::VectorXd::Zero(D);
        ws.up.setGiven(D,zeros.data(),identity.data()); ws.down.setGiven(D,zeros.data(),identity.data());
        ws.valUp.resize(D); ws.valDown.resize(D);
        for (int i=0;i<D;i++){ ws.valUp(i) = doublec(uniform(300.0),0.0); ws.valDown(i) = doublec(uniform(300.0),0.0); }
        ws.colStart = Eigen::VectorXcd::Zero(D); ws.colStart(5) = 1.0;
        ws.rowStart = Eigen::RowVectorXcd::Zero(D); ws.rowStart(38) = 1.0;
        Config* cnf = Config_init();
        strcpy(cnf->method,"CCE"); strcpy(cnf->quantity,"coherence");
        cnf->nstate = 1; cnf->deltat = 0.0137;
        long compared = 0, bad = 0, recursiveOtherBits = 0;
        double recursiveWorst = 0;
        for (const Sequence& s : sequences){
            Built b(s);
            auto run = [&](bool recursive, bool columns, int nb){
                strcpy(cnf->evolution_update, recursive ? "recursive" : "exact");
                cnf->column_evolution = columns; cnf->column_evolution_min_dim = D; cnf->column_evolution_columns = nb;
                evolve(cnf,&b.pulse,ws);
                return ws.phases;
            };
            for (int nstep : {1,2,16,17,31,32,33,48,63,64,65,77}){
                cnf->nstep = nstep;
                const bool anchor = (nstep-1) % 16 == 0;
                const std::vector<doublec> refExact = run(false,false,1), refRecursive = run(true,false,1);
                for (int columns : {32,1,7,nstep+5}){
                    const std::vector<doublec> exact = run(false,true,columns), recursive = run(true,true,columns);
                    for (int ip=0; ip<=b.pulse.npulse; ip++){
                        if (b.idx[ip] != ip) continue;
                        ++compared;
                        const doublec* re = refExact.data()+size_t(2*D*ip);
                        const doublec* ce = exact.data()+size_t(2*D*ip);
                        const doublec* rr = refRecursive.data()+size_t(4*D*ip);
                        const doublec* cr = recursive.data()+size_t(4*D*ip);
                        const char* why = nullptr;
                        if (std::memcmp(re,ce,2*D*sizeof(doublec))) why = "exact: not the bits of the path before";
                        else if (std::memcmp(rr+2*D,cr+2*D,2*D*sizeof(doublec))) why = "recursive step factors: not the bits of the path before";
                        else if (anchor && std::memcmp(cr,ce,2*D*sizeof(doublec))) why = "recursive at an anchor step: not the bits of exact";
                        // (step 1 is 1 times one step's factor: the bits of exact, legitimately)
                        else if (!anchor && nstep != 2 && !std::memcmp(cr,ce,2*D*sizeof(doublec))) why = "recursive between anchors: the bits of exact (not carried)";
                        else {
                            double toExact = 0, toBefore = 0;
                            for (int k=0;k<2*D;k++){ toExact = std::max(toExact,std::abs(cr[k]-ce[k])); toBefore = std::max(toBefore,std::abs(cr[k]-rr[k])); }
                            if (std::memcmp(rr,cr,2*D*sizeof(doublec))) ++recursiveOtherBits;
                            recursiveWorst = std::max(recursiveWorst,toBefore);
                            if (toExact > 1e-11) why = "recursive: more than 1e-11 from exact";
                            else if (toBefore > 1e-14) why = "recursive: more than 1e-14 from the path before";
                        }
                        if (why){
                            if (!bad) printf("FAIL phases: %s nstep %d columns %d interval %d: %s\n",s.name,nstep,columns,ip,why);
                            ++bad;
                        }
                    }
                }
            }
        }
        Config_freeAll(cnf);
        printf("  phases of the last step: %ld interval comparisons, %ld failed; recursive with other bits than the path before "
               "in %ld (largest difference %.1e)\n",compared,bad,recursiveOtherBits,recursiveWorst);
        failed += bad;
    }
    printf("  %ld outputs, %ld of them with other bits than the path before, %ld below 1e-3 in magnitude\n",outputs,otherBits,tiny);
    if (failed){ printf("FAIL column-block one-state evolution: %ld failures (above)\n",failed); return 1; }
    printf("PASS column-block one-state evolution (A62): within 1e-12 relative of the path before (largest %.2e); "
           "the phase rule held; bit-identical below min_dim\n",worst);
}
