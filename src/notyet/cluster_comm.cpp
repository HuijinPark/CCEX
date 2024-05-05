
#include "../include/utilities.h"
#include "../include/memory.h"
#include "../include/bath.h"

void findConnectivity(int*** cmap, float*** stmap, BathArray* ba, float rdip, float rdipcut){

    // cmap : connectivity map
    // cmap(nspin+1, nspin+1) 2d array
    // if i-th spin and j-th spin are connected,
    //  cmap[i+1][j+1] = 1
    // otherwise, cmap[i+1][j+1] = 0

    // stmap : strength map
    // stmap(nspin+1, nspin+1) 2d array
    // if i-th spin and j-th spin are connected,
    //  stmap[i+1][j+1] = strength
    // otherwise, stmap[i+1][j+1] = 0
    // Definition of strength : interaction_m

    int nspin = BathArray_getNspin(ba);
    int maplength = nspin+1;

    // Connectivity Map and Strength Map
    *cmap = allocInt2d(maplength,maplength);
    *stmap = allocFloat2d(maplength,maplength);

    // Initialize : 0th line is not used but save the legnth of the array
    (*cmap)[0][0] = maplength;
    (*stmap)[0][0] = maplength;

    // Find connectivity Map and strength Map    
    for (int row=1; row<nspin+1; row++){
        for (int col=row+1; col<nspin+1; col++){
            int sp1 = row-1;
            int sp2 = col-1;
            double* xyz1 = ba->bathspin[sp1]->xyz;
            double* xyz2 = ba->bathspin[sp2]->xyz;
            double  gyro1 = ba->bathspin[sp1]->gyro;
            double  gyro2 = ba->bathspin[sp2]->gyro;
            double r = dist(xyz1,xyz2);
            double constant = H_BAR * gyro1 * gyro2 / pow(r,3);
            double strength = constant * (ba->intmapBB[sp1][sp2]).zz;

            if ( r < rdip && r > rdipcut){
                (*cmap)[row][col] = 1;
                (*cmap)[col][row] = 1;
            }

            (*stmap)[row][col] = fabs((float)strength);
            (*stmap)[col][row] = fabs((float)strength);
        }
    }

}