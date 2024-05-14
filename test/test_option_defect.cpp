#include "../include/utilities.h"
#include "../include/option.h"
#include "../include/defect.h"
#include "../include/memory.h"
int rank = 0;
int nprocess = 0;
bool verbosity = false;


int main(int argc, char* argv[]){

    DefectArray* dfa = DefectArray_init();

    char* fccein = allocChar1d(MAX_FILEPATH);
    strcpy(fccein,"../config/cce_defect_v1.json");

    cJSON_readOptionDefectArray(dfa,fccein);

    int ndefect = dfa->ndefect;
    for (int idf=0; idf<ndefect; idf++){
        printf("Defect %d\n",idf);
        printf("  dfname: %s\n",dfa->defect[idf]->dfname);
        printf("  apprx: %d\n",dfa->defect[idf]->apprx);

        int naddspin = dfa->defect[idf]->naddspin;
        for (int i=0; i<naddspin; i++){
            printf("  types[%d]: %s\n",i,dfa->defect[idf]->types[i]);
            printf("  spins[%d]: %f\n",i,dfa->defect[idf]->spins[i]);
            printf("  gyros[%d]: %lf\n",i,dfa->defect[idf]->gyros[i]);
            printf("  eqs[%d]: %lf\n",i,dfa->defect[idf]->eqs[i]);
        }

        int navaax = dfa->defect[idf]->navaax;
        for (int iax=0; iax<navaax; iax++){
            for (int isp=0; isp<naddspin; isp++){
                printf("  rxyzs[%d][%d]: %lf %lf %lf\n",iax,isp,dfa->defect[idf]->rxyzs[iax][isp][0],dfa->defect[idf]->rxyzs[iax][isp][1],dfa->defect[idf]->rxyzs[iax][isp][2]);
                
                char str[100];
                sprintf(str,"  hypf[%d][%d]",iax,isp);
                printInlineMatrixXcd(str,dfa->defect[idf]->hypf[iax][isp]);

                sprintf(str,"  efg[%d][%d]",iax,isp);
                printInlineMatrixXcd(str,dfa->defect[idf]->efg[iax][isp]);
            }

            char str[100];
            sprintf(str,"  zfs[%d]",iax);
            printInlineMatrixXcd(str,dfa->defect[idf]->zfs[iax]);
            printf("  detuning[%d]: %lf\n",iax,dfa->defect[idf]->detuning[iax]);
        }
    }
    DefectArray_freeAll(dfa);

    return 0;
}
