// The fixed-size single-state path of calCoherenceCce (bath dimension 2 and 4, A16) against the dynamic one, bit for bit,
// with evolution_update = exact (P1-04 of review 20261012-233058-452278: the fixed size gave (1, +0) where the dynamic path
// gives (1, -0), and other last bits without -march=native). recursive is not compared: its fixed size may differ in the
// last bits (simulator_cce.cpp, colTimes). The test includes simulator_cce.cpp for evolveSingleState and links the other
// objects of main.out; no MPI runtime is needed.
#include "../../../src/calculation/simulator_cce.cpp"
#include <cstdio>
#include <cstring>
#include <random>

int rank=0,nprocess=1,seed=19;
bool verbosity=false;

static long outputs=0;
// Hahn (npulse 1) or CPMG-2 with an interval that reuses an earlier one's phases
template<int N>
static long compare(const MatrixXcd& ha,const MatrixXcd& hb,const Eigen::VectorXcd& psi,int npulse,int nstep,
                    bool show){
    Eigen::SelfAdjointEigenSolver<MatrixXcd> a(ha),b(hb);
    const MatrixXcd M=a.eigenvectors().adjoint()*b.eigenvectors(),Madj=M.adjoint();
    const Eigen::RowVectorXcd r0=psi.adjoint()*b.eigenvectors();
    const Eigen::VectorXcd c0=a.eigenvectors().adjoint()*psi,valA=a.eigenvalues(),valB=b.eigenvalues();
    double s0[3]={0,0,.25},s1[3]={0,0,.5},s2[3]={0,0,.25};
    double* hahn[]={s1,s1}; double* cpmg[]={s0,s1,s2};
    int hahnIndex[]={0,0},cpmgIndex[]={0,1,0};
    double** sequence=npulse==1 ? hahn : cpmg; int* index=npulse==1 ? hahnIndex : cpmgIndex;
    CceWorkspace wd,wf; std::vector<MatrixXcd> rd(nstep),rf(nstep);
    const SingleStateRun dynamic{valA,valB,sequence,index,npulse,nstep,.1,rd.data(),wd};
    const SingleStateRun fixed{valA,valB,sequence,index,npulse,nstep,.1,rf.data(),wf};
    Eigen::RowVectorXcd r(N),rn(N); Eigen::VectorXcd c(N),cn(N),sc(N);
    evolveSingleState<false>(dynamic,M,Madj,r0,c0,r,rn,c,cn,sc);
    const Eigen::Matrix<doublec,N,N> Mf=M,Madjf=Madj;
    const Eigen::Matrix<doublec,1,N> r0f=r0; Eigen::Matrix<doublec,1,N> rf1,rnf;
    const Eigen::Matrix<doublec,N,1> c0f=c0; Eigen::Matrix<doublec,N,1> cf,cnf,scf;
    evolveSingleState<false>(fixed,Mf,Madjf,r0f,c0f,rf1,rnf,cf,cnf,scf);
    long differ=0;
    for (int k=0;k<nstep;k++){
        ++outputs;
        if (std::memcmp(rd[k].data(),rf[k].data(),sizeof(doublec))==0) continue;
        if (show && differ==0)
            printf("  dimension %d step %d: dynamic (%a, %a), fixed (%a, %a)\n",N,k,rd[k](0,0).real(),rd[k](0,0).imag(),
                   rf[k](0,0).real(),rf[k](0,0).imag());
        ++differ;
    }
    return differ;
}

template<int N>
static MatrixXcd hermitian(std::mt19937_64& g,bool integer){
    std::uniform_real_distribution<double> u(-3,3); MatrixXcd h(N,N);
    for (int i=0;i<N;i++) for (int j=i;j<N;j++){
        double x=u(g),y=i==j ? 0 : u(g);
        if (integer){ x=std::round(x); y=std::round(y); }
        h(i,j)=doublec(x,y); h(j,i)=std::conj(h(i,j));
    }
    return h;
}

int main(){
    long failed=0;
    // the review's case: Ha = [[-2, -1-2i], [-1+2i, 2]], Hb = 0, psi = (1, 0), Hahn; at t = 0 the dynamic path gives (1, -0)
    MatrixXcd ha(2,2); ha<<doublec(-2,0),doublec(-1,-2),doublec(-1,2),doublec(2,0);
    Eigen::VectorXcd psi(2); psi<<1.,0.;
    const long review=compare<2>(ha,MatrixXcd::Zero(2,2),psi,1,33,true);
    if (review) printf("FAIL the review's case: %ld of 33 outputs differ\n",review);
    failed+=review;
    // random Hermitian matrices, integer entries (exact zeros and ties) and not, Hb = 0 every third, Hahn and CPMG-2
    std::mt19937_64 g(1);
    long random=0;
    for (int t=0;t<2000;t++){
        const bool integer=t%2; const int npulse=1+(t/2)%2; const bool zero=t%3==0;
        Eigen::VectorXcd p2=Eigen::VectorXcd::Zero(2),p4=Eigen::VectorXcd::Zero(4); p2(t%2)=1; p4(t%4)=1;
        const MatrixXcd a2=hermitian<2>(g,integer),b2=zero ? MatrixXcd(MatrixXcd::Zero(2,2)) : hermitian<2>(g,integer);
        const MatrixXcd a4=hermitian<4>(g,integer),b4=zero ? MatrixXcd(MatrixXcd::Zero(4,4)) : hermitian<4>(g,integer);
        random+=compare<2>(a2,b2,p2,npulse,40,random==0)+compare<4>(a4,b4,p4,npulse,40,random==0);
    }
    if (random) printf("FAIL random: %ld outputs differ\n",random);
    failed+=random;
    if (failed) return 1;
    printf("PASS fixed-size exact single-state path: %ld outputs bit-identical with the dynamic path (signed zeros included)\n",
           outputs);
}
