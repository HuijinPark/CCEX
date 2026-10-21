// Ownership invariants: C++ subspin members and empty per-rank cluster copies.
#include "physical/defect.h"
#include "common/memory.h"
#include <cassert>
int rank=0, nprocess=1, seed=19;
bool verbosity=false;

int main(int argc,char** argv){
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&nprocess);
    for(int iteration=0;iteration<20;++iteration){
        QubitArray* qubits=QubitArray_init();
        QubitArray_alloc_alphaidx_betaidx(qubits);
        int index=1;QubitArray_set_alphaidx(qubits,&index);QubitArray_set_betaidx(qubits,&index);
        QubitArray_set_alphaidx(qubits,nullptr);QubitArray_set_betaidx(qubits,nullptr);
        assert(qubits->_alphaidx==nullptr && qubits->_betaidx==nullptr);
        QubitArray_freeAll(qubits);
        // A freshly created array owns memory even without any Defect entries.
        DefectArray_freeAll(DefectArray_init());
        BathArray bath{};bath.nspin=3;
        DefectArray* dfa=DefectArray_init();
        DefectArray_allocNaddspins(dfa,3);
        dfa->naddspins[0]=0;dfa->naddspins[1]=1;dfa->naddspins[2]=2;
        DefectArray_allocSubbath(dfa,&bath,2);
        for(int i=1;i<3;++i)for(int j=0;j<dfa->naddspins[i];++j){
            BathSpin* spin=dfa->subbath[i][j];
            assert(spin->quad.size()==0 && spin->hypf_sub.size()==0);
            spin->quad=MatrixXcd::Identity(3,3);
            spin->hypf_sub=MatrixXcd::Constant(3,3,doublec(.25,.5));
            for(int q=0;q<2;++q)spin->hypf[q]=MatrixXcd::Identity(3,3);
        }
        DefectArray_freeSubbath(dfa);assert(dfa->subbath==nullptr);
        DefectArray_freeSubbath(dfa);DefectArray_freeAll(dfa);

        for(int order:{0,2,4}){
            int*** source=(int***)calloc(order+1,sizeof(int**));
            for(int n=0;n<=order;++n){
                const int rows=n==1 ? 2 : 1;
                source[n]=allocInt2d(rows,n+1);source[n][0][0]=rows;
                if(n==1){source[n][1][0]=1;source[n][1][1]=17;}
            }
            int*** local=MPI_getLocalClusters(order,source);
            assert(local!=source && local[0]!=source[0]);
            assert(local[0][0][0]==(rank==0 ? 1 : 0));
            for(int n=2;n<=order;++n)assert(local[n][0][0]==1);
            MPI_freeLocalClusters(order,&local);assert(local==nullptr);
            MPI_freeLocalClusters(order,&local);
            // Freeing the independent copy must not change the source.
            if(order>0)assert(source[1][1][1]==17);
            MPI_freeLocalClusters(order,&source);
        }
    }
    if(rank==0)puts("PASS repeated subspin lifetime and MPI cluster-copy ownership, including empty ranks/orders");
    MPI_Finalize();
}
