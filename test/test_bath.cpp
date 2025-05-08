#include "../include/bath.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include <iostream>
#include <stdlib.h>
#include <float.h>

bool verbosity = false;
int rank = 0;

int main(){
    BathArray* ba = BathArray_init();

    int nspin = 10;
    int nqubit = 3;

    // Set nqubit #
    BathArray_setNspin(ba,nspin);
    
    // Allocate bath memory for nspin #
    BathArray_allocBath(ba,nqubit);

    printf("\n    >> Bath defect related hyperfine tensors\n\n");
    BathArray_reportBath_hypf_sub(ba);
    
    ////////////////////////////////////////////
    // Set bath properties
    ////////////////////////////////////////////

    printf("\n    >> Bath properties\n\n");
    for (int ispin=0; ispin<nspin; ispin++){

        // Bath spin properties
        char name[20] = {"14N"}; // "14N

        // Bath spin properties
        float spin = 1.0;
        double xyz[3] = {rand()%100,rand()%100,rand()%100}; 
        double gyro = GAMMA_ELECTRON;
        float state = -1.0;

        BathArray_setBath_i_name(ba,name,ispin); // name
        BathArray_setBath_i_spin(ba,spin,ispin); // spin
        BathArray_setBath_i_xyz(ba,xyz,ispin); // xyz
        BathArray_setBath_i_gyro(ba,gyro,ispin); // gyro
        BathArray_setBath_i_state(ba,state,ispin); // state

        // Bath detuning and disorder
        double detuning = 0.0;
        double disorder = 0.0;
        BathArray_setBath_i_detuning(ba,detuning,ispin); // detuning
        BathArray_setBath_i_disorder(ba,disorder,ispin); // disorder

        // Hpyerfine tensor between bath spin and qubit spin
        BathArray_allocBath_i_hypf(ba,ispin,nqubit); // hyperfine field
        for (int iqubit=0; iqubit<nqubit; iqubit++){
            MatrixXcd hypf = MatrixXcd::Random(3,3);
            BathArray_setBath_i_hypf_j(ba,hypf,ispin,iqubit); // hyperfine field
        }

        // Quadrupole tensor
        MatrixXcd quad = MatrixXcd::Random(3,3);
        BathArray_setBath_i_quad(ba,quad,ispin); // quadrupole field
        
    }
    BathArray_reportBath(ba);

    printf("\n    >> Bath states\n\n");
    BathArray_reportBath_states(ba);

    printf("\n    >> Bath detuning\n\n");
    BathArray_reportBath_detunings(ba);

    printf("\n    >> Bath disorder\n\n");
    BathArray_reportBath_disorders(ba);

    printf("\n    >> Bath hyperfine tensors\n\n");
    BathArray_reportBath_hypf(ba,nqubit);

    printf("\n    >> Bath quadrupole tensors\n\n");
    BathArray_reportBath_quad(ba);

    printf("\n    >> Bath defect related hyperfine tensors\n\n");
    BathArray_reportBath_hypf_sub(ba);

    ////////////////////////////////////////////
    // Check dimension
    ////////////////////////////////////////////
    printf("\n");

    printf("\n    >> Bath spin dimensions\n\n");
    for (int ispin=0; ispin<nspin; ispin++){
        if (!verbosity && (ispin<3 || ispin>nspin-3)){
            int dim = BathArray_dimBath_i(ba,ispin);
            printf("      BathSpin[%d] = %d\n",ispin,dim);
        }
        if (!verbosity && ispin==3){
            printf("         :\n");
        }
    }

    ////////////////////////////////////////////
    // Bath - Bath interaction tensor
    ////////////////////////////////////////////

    printf("\n    >> Defect : set hyperfine with subspin\n\n");
    int ispinlist[5] = {nspin-5,nspin-4,nspin-3,nspin-2,nspin-1};
    int mainspidxlist[5] = {0,1,2,3,4};

    for (int i=0; i<5; i++){
        int ispin = ispinlist[i];
        int mainspidx = mainspidxlist[i];
        MatrixXcd hypf_sub = MatrixXcd::Random(3,3);
        BathArray_setBath_i_hypf_sub(ba,hypf_sub,ispin); // defect related hyperfine field
        BathArray_setBath_i_mainspidx(ba,mainspidx,ispin); // main spin index
    }

    BathArray_reportBath_hypf_sub(ba);

    printf("\n    >> Bath - Bath interaction tensors\n\n");
    for (int ispin=0; ispin<nspin; ispin++){
        for (int jspin=ispin+1; jspin<nspin; jspin++){
            if (ispin!=jspin){
                MatrixXcd inttensor = BathArray_int_i_j(ba,ispin,jspin);
                double r = BathArray_dist_i_j(ba,ispin,jspin);
                double val = std::abs(inttensor.sum());
                if (val> FLT_EPSILON){
                    char key[20] = "\0";
                    sprintf(key,"(%d,%d), r = %5.2lf",ispin,jspin,r);
                    printInlineMatrixXcd(key,inttensor);
                }
            }
        }
    }

    ////////////////////////////////////////////
    // Connectivity
    ////////////////////////////////////////////
    printf("\n    >> Bath connectivity\n\n");

    int** cmap; float** stmap; 
    float rdip = 40.0; float rdipcut = 10.0;
    BathArray_connectivity(&cmap,&stmap,ba,rdip,rdipcut);
    //////////////////////////////////////////////////
    /* Print strength Map */
    
    printf("    Strength Map (rad kHz)\n");
    printf("\n%10s"," ");
    for (int j=0; j<nspin; j++){
      printf("%10d ",j);
    }   
    printf("\n%10s"," ");
    for (int j=0; j<nspin; j++){
      printf("%10s ","---");
    }   
    printf("\n");
    for (int i=0; i<nspin; i++){
      printf("%10d|",i);
      for (int j=0; j<nspin; j++){
          printf("%10.1f ",stmap[i][j]);
      }   
      printf("\n");
    }
    //////////////////////////////////////////////////
    /* Print connectiyvity Map */
    if (rank==0){
        printf("\n");
        printf("    Connectivity Map\n");
        printf("\n%10s"," ");
        for (int j=0; j<nspin; j++){
        printf("%10d ",j);
        }   
        printf("\n%10s"," ");
        for (int j=0; j<nspin; j++){
        printf("%10s ","---");
        }   
        printf("\n");

       int count = 0;
       for (int i=0; i<nspin; i++){
            printf("%10d|",i);
            for (int j=0; j<nspin; j++){
                printf("%10d ",cmap[i][j]);
                count += cmap[i][j];
            }   
            printf("\n");
       }
       printf("\n");
       printf("     Total pair # : %d\n\n",(count)/2);
    }
    ////////////////////////////////////////////////////
    printf("    ==========================\n");
    printf("    Done\n");
    printf("    ==========================\n");
    ////////////////////////////////////////////
    // Free
    ////////////////////////////////////////////
    BathArray_freeAll(ba);
    return 0;
}