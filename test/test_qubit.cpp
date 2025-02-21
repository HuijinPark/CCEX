#include "../include/qubit.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include <iostream>

bool verbosity = true;

int main(){
    QubitArray* qa = QubitArray_init();

    int nqubit = 2;
    int iqubit = 0; // qubit index
    int jqubit = 0; // qubit index

    // Set nqubit #
    QubitArray_setNqubit(qa,nqubit);
    
    // Allocate qubit array memory for nqubit #
    QubitArray_allocQubit(qa);
    
    ////////////////////////////////////////////
    // Set qubit properties
    ////////////////////////////////////////////
    doublec occupied = doublec(1.0,0.0);
    doublec unoccupied = doublec(0.0,0.0);

    printf("\n    >> Qubit properties\n");
    for (iqubit=0; iqubit<nqubit; iqubit++){

        // Qubit spin properties
        char name[10] = "\0"; // name
        snprintf(name,10,"q%d",iqubit);

        // Qubit spin properties
        float spin = 1.0;
        double xyz[3] = {1.0,2.0,3.0}; 
        double gyro = GAMMA_ELECTRON;

        QubitArray_setQubit_i_name(qa,name,iqubit); // name
        QubitArray_setQubit_i_spin(qa,spin,iqubit); // spin
        QubitArray_setQubit_i_xyz(qa,xyz,iqubit); // xyz
        QubitArray_setQubit_i_gyro(qa,gyro,iqubit); // gyro

        // Qubit detuning and overhaus
        double detuning = 0.0;
        double overhaus = 0.0;
        QubitArray_setQubit_i_detuning(qa,detuning,iqubit); // detuning
        QubitArray_setQubit_i_overhaus(qa,overhaus,iqubit); // overhaus

        // Qubit alpha and beta
        int dim = int(2*spin+1);
        MatrixXcd alpha = MatrixXcd::Zero(dim,1);
        MatrixXcd beta = MatrixXcd::Zero(dim,1);
        alpha(0,0) = occupied; alpha(1,0) = unoccupied; alpha(2,0) = unoccupied; //|1>
        beta(0,0) = unoccupied; beta(1,0) = occupied; beta(2,0) = unoccupied; //|0>
        QubitArray_setQubit_i_alpha(qa,alpha,iqubit); //alpha
        QubitArray_setQubit_i_beta(qa,beta,iqubit); //beta

        QubitArray_reportQubit_i(qa,iqubit);
        printf("\n");
    }

    ////////////////////////////////////////////
    // Check dimension
    ////////////////////////////////////////////
    printf("\n");

    int dim = QubitArray_dim(qa);
    
    for (iqubit=0; iqubit<nqubit; iqubit++){
        int dim_i = QubitArray_dimQubit_i(qa,iqubit);
        printf("    QubitArray dimension of qubit[%d] = %d\n",iqubit,dim_i);
    }
    printf("\n    QubitArray total dimension = %d\n",dim);
    
    ////////////////////////////////////////////
    // Set interactions
    ////////////////////////////////////////////

    // alloc
    QubitArray_allocIntmap(qa);

    // Tensors
    MatrixXcd Tzfs(3,3);
    Tzfs <<  100.0, 0.0, 0.0,
             0.0, 200.0, 0.0,
             0.0, 0.0, 300.0;

    MatrixXcd Tint(3,3);
    Tint <<  1.0, 2.0, 3.0,
             4.0, 5.0, 6.0,
             7.0, 8.0, 9.0;

    // zero-field splitting tensor of 0th qubit
    iqubit = 0;
    jqubit = 0;
    QubitArray_setIntmap_i_j(qa,Tzfs,iqubit,jqubit);

    // Spin interaction tensor between 0th and 1st qubit
    iqubit = 0;
    jqubit = 1;
    QubitArray_setIntmap_i_j(qa,Tint,iqubit,jqubit);

    // zero-field splitting tensor of 1st qubit
    iqubit = 1;
    jqubit = 1;
    QubitArray_setIntmap_i_j(qa,Tzfs,iqubit,jqubit);

    ////////////////////////////////////////////
    // Set psia, psib, psi0
    ////////////////////////////////////////////

    //----------------------------------------
    // Case1. Directly set psia, psib, psi0
    //----------------------------------------
    // states of psia , psib
    MatrixXcd psia = MatrixXcd::Zero(dim,1);
    MatrixXcd psib = MatrixXcd::Zero(dim,1);
    
    psia(0,0) = occupied; psia(1,0) = occupied; psia(2,0) = occupied; // |11> + |10> + |1-1>
    psib(dim-3,0) = occupied; psib(dim-2,0) = occupied; psib(dim-1,0) = occupied; // |-11> + |-10> + |-1-1>

    // states of psi0
    MatrixXcd psi0 = psia + psib;

    // set
    QubitArray_setPsia(qa,psia);
    QubitArray_setPsib(qa,psib);
    QubitArray_setPsi0(qa,psi0);

    // report
    printf("\n    >> Directly set psia, psib, psi0\n");
    QubitArray_reportPsiaPsib(qa);
    QubitArray_reportPsi0(qa);
    printf("\n");
    //----------------------------------------
    // Case2. Set psia, psib, psi0 from qubit array
    //----------------------------------------
    //set
    QubitArray_setPsiaPsib_fromQubit(qa); // psia = kron(qubit_i_alpha,qubit_j_alpha, ...)
    QubitArray_setPsi0_fromPsiaPsib(qa); // psi0 = Normalize(psia + psib)

    // report
    printf("\n    >> Set psia, psib from qubit array\n");
    QubitArray_reportPsiaPsib(qa);
    QubitArray_reportPsi0(qa);
    printf("\n");
 
    ////////////////////////////////////////////
    // Report (Again)
    ////////////////////////////////////////////
    QubitArray_report(qa);

    ////////////////////////////////////////////
    // Free
    ////////////////////////////////////////////
    QubitArray_freeAll(qa);
    return 0;
}