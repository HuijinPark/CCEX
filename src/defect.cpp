#include "../include/defect.h"
#include "../include/memory.h"
#include "../include/hamiltonian.h"



void DefectArray_setPaxesRandom(DefectArray* dfa, BathArray* ba){

    int nbathspin = BathArray_getNspin(ba);
    
    for (int ibs=0; ibs<nbathspin; ibs++){
        char* name = BathArray_getBath_i_name(ba,ibs);
        int idf = DefectArray_findDefectIndex(dfa,name);
        int navaax = DefectArray_getDefect_idf_navaax(dfa,idf);
        int minpax = 1;
        int maxpax = (navaax-1);

        dfa->paxes[ibs] = (rand() % (maxpax))+minpax;
    }

}

void DefectArray_setSubbathStatesRandom(DefectArray* dfa, BathArray* ba){

    int nbathspin = BathArray_getNspin(ba);

    for (int ibs=0; ibs<nbathspin; ibs++){

        char* name = BathArray_getBath_i_name(ba,ibs);
        int idf = DefectArray_findDefectIndex(dfa,name);

        int naddspin = dfa->naddspins[ibs];
        for (int isp=0; isp<naddspin; isp++){
            float S = DefectArray_getDefect_idf_isp_spins(dfa,idf,isp);
            int dim = (int)(2*S+1);
            int r = rand() % dim;
            float ms = S - r;
            dfa->subbath[ibs][isp]->state = ms;
        }
    }

}


/* Low level functions ------------------------------------ */

DefectArray* DefectArray_init(){
    DefectArray* dfa = (DefectArray*)allocArray1d(1,sizeof(DefectArray));
    dfa->ndefect = 0;
    dfa->defect = NULL;
    dfa->paxes = NULL;
    dfa->naddspins = NULL;
    dfa->subbath = NULL;

    return dfa;
}


void DefectArray_setNdefect(DefectArray* dfa, int ndefect){
    dfa->ndefect = ndefect;
}


void DefectArray_allocDefect(DefectArray* dfa){
    dfa->defect = (Defect**)allocArray2d(dfa->ndefect,1,sizeof(Defect));
} 

void DefectArray_allocDefect_idf(DefectArray* dfa, int idf, int navaax, int naddspin){

    Defect* df = dfa->defect[idf];

    df->naddspin = naddspin;

    df->types = allocChar2d(naddspin,MAX_CHARARRAY_LENGTH);
    df->spins = allocFloat1d(naddspin);
    df->gyros = allocDouble1d(naddspin);
    df->eqs = allocDouble1d(naddspin);

    df->navaax = navaax;

    df->rxyzs = allocDouble3d(navaax,naddspin,3);
    df->hypf = allocMatrixXcd2d(navaax,naddspin);
    df->efg = allocMatrixXcd2d(navaax,naddspin);
    df->zfs = allocMatrixXcd1d(navaax);
    df->detuning = allocDouble1d(navaax);
}


void DefectArray_setDefect_idf_dfname(DefectArray* dfa, int idf, char* dfname){
    strcpy(dfa->defect[idf]->dfname,dfname);
}

void DefectArray_setDefect_idf_apprx(DefectArray* dfa, int idf, bool apprx){
    dfa->defect[idf]->apprx = apprx;
}

void DefectArray_setDefect_idf_naddspin(DefectArray* dfa, int idf, int naddspin){
    dfa->defect[idf]->naddspin = naddspin;
}

void DefectArray_setDefect_idf_types(DefectArray* dfa, int idf, char** types){
    for (int i=0; i<dfa->defect[idf]->naddspin; i++){
        strcpy(dfa->defect[idf]->types[i],types[i]);
    }
}

void DefectArray_setDefect_idf_spins(DefectArray* dfa, int idf, float* spins){
    for (int i=0; i<dfa->defect[idf]->naddspin; i++){
        dfa->defect[idf]->spins[i] = spins[i];
    }
}

void DefectArray_setDefect_idf_gyros(DefectArray* dfa, int idf, double* gyros){
    for (int i=0; i<dfa->defect[idf]->naddspin; i++){
        dfa->defect[idf]->gyros[i] = gyros[i];
    }
}

void DefectArray_setDefect_idf_eqs(DefectArray* dfa, int idf, double* eqs){
    for (int i=0; i<dfa->defect[idf]->naddspin; i++){
        dfa->defect[idf]->eqs[i] = eqs[i];
    }
}

void DefectArray_setDefect_idf_navaax(DefectArray* dfa, int idf, int navaax){
    dfa->defect[idf]->navaax = navaax;
}

void DefectArray_setDefect_idf_iax_isp_rxyz(DefectArray* dfa, int idf, int iax, int isp, double* rxyzs){
    dfa->defect[idf]->rxyzs[iax][isp][0] = rxyzs[0];
    dfa->defect[idf]->rxyzs[iax][isp][1] = rxyzs[1];
    dfa->defect[idf]->rxyzs[iax][isp][2] = rxyzs[2];
}
void DefectArray_setDefect_idf_iax_isp_hypf(DefectArray* dfa, int idf, int iax, int isp, MatrixXcd hypf){
    dfa->defect[idf]->hypf[iax][isp] = hypf;
}
void DefectArray_setDefect_idf_iax_isp_efg(DefectArray* dfa, int idf, int iax, int isp, MatrixXcd efg){
    dfa->defect[idf]->efg[iax][isp] = efg;
}
void DefectArray_setDefect_idf_iax_zfs(DefectArray* dfa, int idf, int iax, MatrixXcd zfs){
    dfa->defect[idf]->zfs[iax] = zfs;
}
void DefectArray_setDefect_idf_iax_detuning(DefectArray* dfa, int idf, int iax, double detuning){
    dfa->defect[idf]->detuning[iax] = detuning;
}


void DefectArray_setNaddspins(DefectArray* dfa, BathArray* ba){
    for (int i=0; i<BathArray_getNspin(ba); i++){
        char* dfname = BathArray_getBath_i_name(ba,i);
        int idf = DefectArray_findDefectIndex(dfa,dfname);
        dfa->naddspins[i] = DefectArray_getDefect_idf_naddspin(dfa,idf);
    }
}

int DefectArray_findDefectIndex(DefectArray* dfa, char* dfname){
    for (int idf=0; idf<dfa->ndefect; idf++){
        if (strcasecmp(dfa->defect[idf]->dfname,dfname) == 0){
            return idf;
        }
    }
    return -1;
}

void DefectArray_setPaxes_i(DefectArray* dfa, int ibs, int axis){
    dfa->paxes[ibs] = axis;
}

int DefectArray_getNdefect(DefectArray* dfa){
    return dfa->ndefect;
}

char* DefectArray_getDefect_idf_dfname(DefectArray* dfa, int idf){
    return dfa->defect[idf]->dfname;
}

bool DefectArray_getDefect_idf_apprx(DefectArray* dfa, int idf){
    return dfa->defect[idf]->apprx;
}

int DefectArray_getDefect_idf_naddspin(DefectArray* dfa, int idf){
    return dfa->defect[idf]->naddspin;
}

char* DefectArray_getDefect_idf_isp_types(DefectArray* dfa, int idf, int isp){
    return dfa->defect[idf]->types[isp];
}

float DefectArray_getDefect_idf_isp_spins(DefectArray* dfa, int idf, int isp){
    return dfa->defect[idf]->spins[isp];
}

double DefectArray_getDefect_idf_isp_gyros(DefectArray* dfa, int idf, int isp){
    return dfa->defect[idf]->gyros[isp];

}

double DefectArray_getDefect_idf_isp_eqs(DefectArray* dfa, int idf, int isp){
    return dfa->defect[idf]->eqs[isp];
}

int DefectArray_getDefect_idf_navaax(DefectArray* dfa, int idf){
    return dfa->defect[idf]->navaax;
}

double* DefectArray_getDefect_idf_iax_isp_rxyz(DefectArray* dfa, int idf, int iax, int isp){
    return dfa->defect[idf]->rxyzs[iax][isp];
}

MatrixXcd DefectArray_getDefect_idf_iax_isp_hypf(DefectArray* dfa, int idf, int iax, int isp){
    return dfa->defect[idf]->hypf[iax][isp];
}

MatrixXcd DefectArray_getDefect_idf_iax_isp_efg(DefectArray* dfa, int idf, int iax, int isp){
    return dfa->defect[idf]->efg[iax][isp];
}

MatrixXcd DefectArray_getDefect_idf_iax_isp_zfs(DefectArray* dfa, int idf, int iax){
    return dfa->defect[idf]->zfs[iax];
}

double DefectArray_getDefect_idf_iax_detuning(DefectArray* dfa, int idf, int iax){
    return dfa->defect[idf]->detuning[iax];
}

int DefectArray_getNbathspin(DefectArray* dfa){
    return dfa->nbathspin;
}

int DefectArray_getPaxes_i(DefectArray* dfa, int ibs){
    return dfa->paxes[ibs];
}

int DefectArray_getNaddspins_i(DefectArray* dfa, int ibs){
    return dfa->naddspins[ibs];

}

BathSpin* DefectArray_getSubbath_i_isp(DefectArray* dfa, int ibs, int isp){
    return dfa->subbath[ibs][isp];
}

// free
void DefectArray_freeAll(DefectArray* dfa){

    if (dfa->defect == NULL){
        return;
    }

    for (int idf=0; idf<dfa->ndefect; idf++){

        freeChar2d(&(dfa->defect[idf]->types),dfa->defect[idf]->naddspin);
        freeFloat1d(&(dfa->defect[idf]->spins));
        freeDouble1d(&(dfa->defect[idf]->gyros));
        freeDouble1d(&(dfa->defect[idf]->eqs));

        freeDouble3d(&(dfa->defect[idf]->rxyzs),dfa->defect[idf]->navaax,dfa->defect[idf]->naddspin);
        freeMatrixXcd2d(&(dfa->defect[idf]->hypf),dfa->defect[idf]->navaax);
        freeMatrixXcd2d(&(dfa->defect[idf]->efg),dfa->defect[idf]->navaax);
        freeMatrixXcd1d(&(dfa->defect[idf]->zfs));
        freeDouble1d(&(dfa->defect[idf]->detuning));
    }
    freeArray2d((void***)&(dfa->defect),dfa->ndefect);

    DefectArray_freePaxes(dfa);

    DefectArray_freeSubbath(dfa);
    DefectArray_freeNaddspins(dfa);
    
    freeArray1d((void**)&(dfa));
}

void DefectArray_freePaxes(DefectArray* dfa){

    if (dfa->paxes == NULL){
        return;
    }
    freeInt1d(&(dfa->paxes));
}

void DefectArray_freeNaddspins(DefectArray* dfa){

    if (dfa->naddspins == NULL){
        return;
    }
    freeInt1d(&(dfa->naddspins));
}


void DefectArray_freeSubbath(DefectArray* dfa){

    if (dfa->subbath == NULL){
        return;
    }
    
    for (int i=0; i<dfa->nbathspin; i++){
        for (int j=0; j<dfa->naddspins[i]; j++){
            freeMatrixXcd1d(&(dfa->subbath[i][j]->hypf));
            free(dfa->subbath[i][j]);
        }
        free(dfa->subbath[i]);
    }
    free(dfa->subbath);
    dfa->subbath = NULL;
}


// alloc
void DefectArray_allocPaxes(DefectArray* dfa, int nbathspin){
    dfa->paxes = allocInt1d(nbathspin);
    dfa->nbathspin = nbathspin;
}

void DefectArray_allocNaddspins(DefectArray* dfa, int nbathspin){
    dfa->naddspins = allocInt1d(nbathspin);
    dfa->nbathspin = nbathspin;
}

void DefectArray_allocSubbath(DefectArray* dfa, BathArray* ba, int nqubit){

    int nbathspin = BathArray_getNspin(ba);
    dfa->nbathspin = nbathspin;

    // alloc BathSpin
    dfa->subbath = (BathSpin***)calloc(nbathspin,sizeof(BathSpin**));
    for (int i=0; i<nbathspin; i++){
        int naddspin = dfa->naddspins[i];
        dfa->subbath[i] = (BathSpin**)calloc(naddspin,sizeof(BathSpin*));
        for (int j=0; j<naddspin; j++){
            dfa->subbath[i][j] = (BathSpin*)calloc(1,sizeof(BathSpin));
            dfa->subbath[i][j]->hypf = allocMatrixXcd1d(nqubit);
        }
    }   
    
}

void DefectArray_setSubbath(DefectArray* dfa, BathArray* ba, QubitArray* qa){

    /** What we do here :
     * 1. Set subbath name, xyz, spin, gyro, (eq,) hypf, mainspidx, 
     * 2. Set quad(zfs),hypf_sub (depending on principal axis) 
     * 2. Set detuning = 0.0 (No detuning)
     * 3. Set state = 0.0 and disorder = 0.0 (will be updated later)
     */

    int nbathspin = BathArray_getNspin(ba);

    for (int ibs=0; ibs<nbathspin; ibs++){

        /// BathSpin information (Main spins - electron)
        char* bsname = BathArray_getBath_i_name(ba,ibs);
        double* bsxyz = BathArray_getBath_i_xyz(ba,ibs);

        // Spin principal axis
        int iax = DefectArray_getPaxes_i(dfa,ibs);

        // Defect index
        int naddspin = DefectArray_getNaddspins_i(dfa,ibs);
        int idf = DefectArray_findDefectIndex(dfa,bsname);
        for (int isp=0; isp<naddspin; isp++){
            
            ////////////////////////////////////
            // get BathSpin address (sub)
            ////////////////////////////////////
            BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs,isp);
            ////////////////////////////////////
            // BathSpin : Name
            ////////////////////////////////////
            BathSpin_setName_withType(subspin,bsname,DefectArray_getDefect_idf_isp_types(dfa,idf,isp));

            ////////////////////////////////////
            // BathSpin : Spin properties
            ////////////////////////////////////
            BathSpin_setSpin(subspin,DefectArray_getDefect_idf_isp_spins(dfa,idf,isp));
            BathSpin_setGyro(subspin,DefectArray_getDefect_idf_isp_gyros(dfa,idf,isp));
            // BathSpin_setEq(subspin,DefectArray_getDefect_idf_isp_eqs(dfa,idf,isp));

            ////////////////////////////////////
            // BathSpin : Spin position
            ////////////////////////////////////
            BathSpin_setXyz_fromRxyz(subspin,bsxyz,DefectArray_getDefect_idf_iax_isp_rxyz(dfa,idf,iax,isp));

            ////////////////////////////////////
            // BathSpin : Spin state & Disorder 
            // (This will be changed later)
            ////////////////////////////////////
            BathSpin_setState(subspin,0.0);
            BathSpin_setDisorder(subspin,0.0);
            
            ////////////////////////////////////
            // BathSpin : Detuning (No detuning)
            ////////////////////////////////////
            BathSpin_setDetuning(subspin,0.0);

            ////////////////////////////////////
            // BathSpin : Hyperfine (as much as qubit #)
            ////////////////////////////////////
            // Hyperfine tensor : Qubit-Bath
            int nqubit = QubitArray_getNqubit(qa);
            double* xyz = BathSpin_getXyz(subspin);
            double gyro = BathSpin_getGyro(subspin);
            for (int iq=0; iq<nqubit; iq++){
                double qgyro = QubitArray_getQubit_i_gyro(qa,iq);
                double* qxyz = QubitArray_getQubit_i_xyz(qa,iq);
                BathSpin_setHypf_i(subspin, calPointDipoleTensor(qxyz, xyz, qgyro, gyro),iq);
            }
        
            ////////////////////////////////////
            // BathSpin : Quadrupole & Hyperfine_sub
            // This tensor changes depending on principal axis(pax)
            ////////////////////////////////////
            float spin = BathSpin_getSpin(subspin);
            MatrixXcd efg = DefectArray_getDefect_idf_iax_isp_efg(dfa,idf,iax,isp); // Hartree/Bohr^2
            double eq = DefectArray_getDefect_idf_isp_eqs(dfa,idf,isp); // 10e-30 m^2
            BathSpin_setQuad_fromEFG(subspin,efg,eq,spin); // radkHz
            
            // Hyperfine sub tensor : MainDefect - BathSpin (inter bath)
            BathSpin_setHypfSub(subspin,
            MHZ_TO_RADKHZ(DefectArray_getDefect_idf_iax_isp_hypf(dfa,idf,iax,isp))); //radkHz

            // Maindefect index : MainDefect - BathSpin (inter bath)
            BathSpin_setMainspidx(subspin,ibs); 
        }
    }
}



void updateMainSpins_fromDefectArray(DefectArray* dfa, BathArray* ba){
    int nbathspin = BathArray_getNspin(ba);

    for (int ibs=0; ibs<nbathspin; ibs++){

        char* name = BathArray_getBath_i_name(ba,ibs);
        int idf = DefectArray_findDefectIndex(dfa,name);
        int iax = DefectArray_getPaxes_i(dfa,ibs);

        int naddspin = DefectArray_getNaddspins_i(dfa,ibs);
        for (int isp=0; isp<naddspin; isp++){
        
            ////////////////////////////////////
            // Update mainspin index
            int mainspidx = ibs;
            BathArray_setBath_i_mainspidx(ba,mainspidx,ibs);
            ////////////////////////////////////
            // Update quad(zfs) depending on principal axis
            MatrixXcd zfs = DefectArray_getDefect_idf_iax_isp_zfs(dfa,idf,iax); //MHz
            BathArray_setBath_i_quad(ba,MHZ_TO_RADKHZ(zfs),ibs); //radkHz
            ////////////////////////////////////
            // Update detuning depending on principal axis (add)
            double current_detuning = BathArray_getBath_i_detuning(ba,ibs); //radkHz
            double additional_detuning = DefectArray_getDefect_idf_iax_detuning(dfa,idf,iax); //MHz
            additional_detuning = MHZ_TO_RADKHZ(additional_detuning); //radkHz
            BathArray_setBath_i_detuning(ba,(current_detuning+additional_detuning),ibs); //radkHz
            ////////////////////////////////////
        }
    }
}

void updateDisorder_main_sub(DefectArray* dfa, BathArray* ba){

    if (verbosity && rank==0){
        printSubTitle("Update disorder of BathArray and DefectArray (from defect information)");
        printMessage("The bathspin is actually defect, so the other spins can be exist near the defect.\n\
         We call them as subspins. And we update the disorder originated from the subspins\n\n\
         * mainspin - subspin of this mainspin  : hypf_sub of subspin\n\
         * mainspin - subspin of other mainspin : point-dipole approximation \n");
    }

    int nbathspin = BathArray_getNspin(ba);

    for (int ibs=0; ibs<nbathspin; ibs++){

        BathSpin* mainspin = BathArray_getBath_i(ba,ibs);

        // Current disorder of bathspin(main)
        //
        // Disorder types of main spins:
        // main - other main (already calculated) 
        // main - subspin of the main (hypf_sub) ...(m-ms) (Calculating here)
        // main - subspin of the other main (point-dipole approximation) ...(m-os) (Calculating here)
        double main_current_disorder = BathSpin_getDisorder(mainspin);
        double main_additional_disorder = 0.0;

        // Main spin information
        double* main_xyz = BathSpin_getXyz(mainspin);
        double main_gyro = BathSpin_getGyro(mainspin);
        float main_ms = BathSpin_getState(mainspin);

        // Additional disorder due to sub spins
        for (int ibs_sub = 0; ibs_sub<nbathspin; ibs_sub++){

            int naddspin = DefectArray_getNaddspins_i(dfa,ibs_sub);
            for (int isp=0; isp<naddspin; isp++){

                BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs_sub,isp);

                // Current disorder of subspin
                //
                // Disorder types of sub spins:
                // subspin - mainspin of the subspin (hypf_sub) ...(s-mm) (Calculating here)
                // subspin - mainspin of the other subspin (point-dipole approximation) ...(s-om) (Calculating here)
                // subspin - other subspin of sharing the same mainspn (point-dipole approximation) ...(s-ms) (will be calculated in other function)
                // subspin - other subspin of other main spin (point-dipole approximation) ...(s-os) (will be calculated in other function)
                double sub_current_disorder = BathSpin_getDisorder(subspin);
                double sub_additional_disorder = 0.0;

                // Sub-spin information
                double* sub_xyz = BathSpin_getXyz(subspin);
                double sub_gyro = BathSpin_getGyro(subspin);
                float sub_ms = BathSpin_getState(subspin);
                

                // Get disorders
                if (ibs == ibs_sub){
                    // if main-sub, additional disorder = hypf_sub(2,2)*ms_sub
                    MatrixXcd main_sub_int  = BathSpin_getHypfSub(subspin); // radkHz
                    main_additional_disorder += main_sub_int(2,2).real()*sub_ms; // ...(m-ms)
                    sub_additional_disorder += main_sub_int(2,2).real()*main_ms; // ...(s-mm)
                }
                else{
                    // else, additional disorder = Point-dipole tensor * ms_sub
                    MatrixXcd main_sub_int  = calPointDipoleTensor(main_xyz, sub_xyz, main_gyro, sub_gyro); // radkHz
                    main_additional_disorder += main_sub_int(2,2).real()*sub_ms; // ...(m-os)
                    sub_additional_disorder += main_sub_int(2,2).real()*main_ms; // ...(s-om)
                }

                // Update subspin
                BathSpin_setDisorder(subspin,sub_current_disorder + sub_additional_disorder); //radkHz

                // Print
                if (verbosity && rank==0){
                    char message[500];
                    sprintf(message,"SubSpin[%d][%d].disorder = %+7.2g(Current) %+7.2g(Additional) = %+7.2g (radkHz)"\
                    ,ibs_sub,isp,sub_current_disorder,sub_additional_disorder,BathSpin_getDisorder(subspin));
                    printMessage(message);
                }
            }
        }

        // Update main spin
        BathSpin_setDisorder(mainspin,main_current_disorder+main_additional_disorder); //radkHz

        // Print
        if (verbosity && rank==0){
            printf("\n");
            char message[500];
            sprintf(message,"MainSpin[%d].disorder = %+7.2g(Current) %+7.2g(Additional) = %+7.2g (radkHz)"\
            ,ibs,main_current_disorder,main_additional_disorder,BathSpin_getDisorder(mainspin));
            printMessage(message);
            printf("\n");
        }
    }
}

void updateDisorder_sub_sub(DefectArray* dfa){
    
    // Current disorder of subspin
    //
    // Disorder types of sub spins:
    // subspin - mainspin of the subspin (hypf_sub) ...(s-mm) (will be calculated in other function)
    // subspin - mainspin of the other subspin (point-dipole approximation) ...(s-om) (will be calculated in other function)
    // subspin - other subspin of sharing the same mainspn (point-dipole approximation) ...(s-ms) (Calculating here)
    // subspin - other subspin of other main spin (point-dipole approximation) ...(s-os) (Calculating here)

    if (verbosity && rank==0){
        printSubTitle("Update disorder of Subbath in DefectArray");
        printMessage("The bathspin is actually defect, so the other spins can be exist near the defect.\n\
         We call them as subspins. And we update the disorder originated from the subspins\n\n\
         * subspin - other subspin of sharing the same mainspn : point-dipole approximation\n\
         * subspin - other subspin of other main spin : point-dipole approximation: \n");
    }

    int nbathspin = DefectArray_getNbathspin(dfa);

    for (int ibs_current=0; ibs_current<nbathspin; ibs_current++){

        int naddspin_current = DefectArray_getNaddspins_i(dfa,ibs_current);
        for (int isp_current=0; isp_current<naddspin_current; isp_current++){
            ///////////////////////////////////////////////////
            // Current subspin
            BathSpin* subspin_current = DefectArray_getSubbath_i_isp(dfa,ibs_current,isp_current);
            double sub_current_disorder = BathSpin_getDisorder(subspin_current);
            double sub_current_additional_disorder = 0.0;
            // Current Spin properties
            double* sub_current_xyz = BathSpin_getXyz(subspin_current);
            double sub_current_gyro = BathSpin_getGyro(subspin_current);
            ///////////////////////////////////////////////////

            for (int ibs_walk=0; ibs_walk<nbathspin; ibs_walk++){
                int naddspin_walk = DefectArray_getNaddspins_i(dfa,ibs_walk);
                for (int isp_walk=0; isp_walk<naddspin_walk; isp_walk++){
                    ///////////////////////////////////////////////////
                    // Walk subspin
                    BathSpin* subspin_walk = DefectArray_getSubbath_i_isp(dfa,ibs_walk,isp_walk);
                    if (ibs_current == ibs_walk && isp_current == isp_walk){
                        continue;
                    }else{
                        ///////////////////////////////////////////////////
                        // Calculate the additional disorder
                        ///////////////////////////////////////////////////
                        // Walk Spin properties
                        double* sub_walk_xyz = BathSpin_getXyz(subspin_walk);
                        double sub_walk_gyro = BathSpin_getGyro(subspin_walk);
                        float sub_walk_ms = BathSpin_getState(subspin_walk);
                        // Point-dipole approximation
                        MatrixXcd tensor = calPointDipoleTensor(sub_current_xyz, sub_walk_xyz, sub_current_gyro, sub_walk_gyro);
                        sub_current_additional_disorder += tensor(2,2).real()*sub_walk_ms; // ...(s-os), (s-ms)
                        ///////////////////////////////////////////////////
                    }
                }
            }

            // Update subspin
            BathSpin_setDisorder(subspin_current,\
            sub_current_disorder + sub_current_additional_disorder); //radkHz

            // Print
            if (verbosity && rank==0){
                char message[500];
                sprintf(message,"SubSpin[%d][%d].disorder = %+7.2g(Current) %+7.2g(Additional) = %+7.2g (radkHz)"\
                ,ibs_current,isp_current,sub_current_disorder,sub_current_additional_disorder,BathSpin_getDisorder(subspin_current));
                printMessage(message);
            }
        }
    } 

    if (verbosity && rank==0){
        printf("\n\n");
    }
}

void updateOverhaus_qubit_sub(DefectArray* dfa, QubitArray* qa){
    
    if (verbosity && rank==0){
        printSubTitle("Update the overhauser field of QubitArray from the subspins of DefectArray");
        printMessage(" * qubit - subspin : hypf of subspin (point-dipole approximation) \n");
    }

    int nbathspin = DefectArray_getNbathspin(dfa);
    int nqubit = QubitArray_getNqubit(qa);

    for (int iq=0; iq<nqubit; iq++){

        double current_qubit_overhaus = QubitArray_getQubit_i_overhaus(qa,iq);
        double additional_qubit_overhaus = 0.0;

        for (int ibs=0; ibs<nbathspin; ibs++){
            int naddspin = DefectArray_getNaddspins_i(dfa,ibs);
            for (int isp=0; isp<naddspin; isp++){
                ////////////////////////////////////
                // Current subspin
                BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs,isp);
                MatrixXcd hypf = BathSpin_getHypf_i(subspin,iq); // radkHz
                float state = BathSpin_getState(subspin);
                ////////////////////////////////////
                // Overhauser of qubit-subspin
                double overhaus = hypf(2,2).real()*state; // radkHz
                additional_qubit_overhaus += overhaus; // radkHz
            }
        }

        // Update qubit
        QubitArray_setQubit_i_overhaus(qa,current_qubit_overhaus+additional_qubit_overhaus,iq); // radkHz

        // Print
        if (verbosity && rank==0){
            char message[500];
            sprintf(message,"Qubit[%d].overhaus = %+7.2g(Current) %+7.2g(Additional) = %+7.2g (radkHz)"\
            ,iq,current_qubit_overhaus,additional_qubit_overhaus,QubitArray_getQubit_i_overhaus(qa,iq));
            printMessage(message);
            printf("\n");
        }
    }
}

void DefectArray_reportDefect_idf(DefectArray* dfa, int idf){

    char subtitle[100];
    sprintf(subtitle,"DefectArray->defect[%d]",idf);
    printSubTitle(subtitle);
    printStructElementChar("dfname",DefectArray_getDefect_idf_dfname(dfa,idf));
    printStructElementBool("apprx",DefectArray_getDefect_idf_apprx(dfa,idf));
    printStructElementInt("naddspin",DefectArray_getDefect_idf_naddspin(dfa,idf));
    printf("      %-18s:   %-d - 1 # ( possible axis : 1 ~ %d )\n", "navaax", \
    DefectArray_getDefect_idf_navaax(dfa,idf), DefectArray_getDefect_idf_navaax(dfa,idf)-1);

    printStructElementChar2d("types",dfa->defect[idf]->types,DefectArray_getDefect_idf_naddspin(dfa,idf));
    printStructElementFloat1d("spins",dfa->defect[idf]->spins,DefectArray_getDefect_idf_naddspin(dfa,idf));
    printStructElementDouble1d("gyros",dfa->defect[idf]->gyros,DefectArray_getDefect_idf_naddspin(dfa,idf));
    printStructElementDouble1d("eqs",dfa->defect[idf]->eqs,DefectArray_getDefect_idf_naddspin(dfa,idf));
    printf("\n");

    for (int i=1; i<DefectArray_getDefect_idf_navaax(dfa,idf); i++){
        for (int j=0; j<DefectArray_getDefect_idf_naddspin(dfa,idf); j++){
            char message[100];
            sprintf(message,"rxyz[%d][%d]",i,j);
            printStructElementDouble1d(message,DefectArray_getDefect_idf_iax_isp_rxyz(dfa,idf,i,j),3);            
        }
    }
    printf("\n");

    for (int i=1; i<DefectArray_getDefect_idf_navaax(dfa,idf); i++){
        for (int j=0; j<DefectArray_getDefect_idf_naddspin(dfa,idf); j++){
            char message[100];
            MatrixXcd tensor;
            sprintf(message,"hypf[%d][%d]",i,j);
            tensor = DefectArray_getDefect_idf_iax_isp_hypf(dfa,idf,i,j);
            if (tensor.rows() != 0){
                printInlineMatrixXcd(message,tensor);

                if(j==DefectArray_getDefect_idf_naddspin(dfa,idf)-1 && i==DefectArray_getDefect_idf_navaax(dfa,idf)-1){
                    printf("\n");
                }
            }
        }
    }

    for (int i=1; i<DefectArray_getDefect_idf_navaax(dfa,idf); i++){
        for (int j=0; j<DefectArray_getDefect_idf_naddspin(dfa,idf); j++){
            char message[100];
            MatrixXcd tensor;
            sprintf(message,"efg[%d][%d]",i,j);
            tensor = DefectArray_getDefect_idf_iax_isp_efg(dfa,idf,i,j);
            if (tensor.rows() != 0){
                printInlineMatrixXcd(message,tensor);

                if(j==DefectArray_getDefect_idf_naddspin(dfa,idf)-1 && i==DefectArray_getDefect_idf_navaax(dfa,idf)-1){
                    printf("\n");
                }
            }
        }
    }

    for (int i=1; i<DefectArray_getDefect_idf_navaax(dfa,idf); i++){
        char message[100];
        MatrixXcd tensor;
        sprintf(message,"zfs[%d]",i);
        tensor = DefectArray_getDefect_idf_iax_isp_zfs(dfa,idf,i);
        if (tensor.rows() != 0){
            printInlineMatrixXcd(message,tensor);

            if(i==DefectArray_getDefect_idf_navaax(dfa,idf)-1){
                printf("\n");
            }
        }
    }

    for (int i=1; i<DefectArray_getDefect_idf_navaax(dfa,idf); i++){
        char message[100];
        sprintf(message,"detuning[%d]",i);
        printStructElementDouble(message,DefectArray_getDefect_idf_iax_detuning(dfa,idf,i));
    }
    printf("\n");
    printLine();

}

void DefectArray_reportPaxes(DefectArray* dfa){

    printSubTitle("DefectArray->paxes");

    int nbathspin = DefectArray_getNbathspin(dfa);
    printStructElementInt("nbathspin (#)",nbathspin);
    printLine();
    for (int i=0; i<nbathspin; i++){

        if (verbosity || (i<3 || i>nbathspin-3)){ 
            char message[100];
            sprintf(message,"paxes[%d]",i);
            printStructElementInt(message,DefectArray_getPaxes_i(dfa,i));
        }

        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
    printf("\n");
    printLine();

}

void DefectArray_reportSubbath_states(DefectArray* dfa){

    printSubTitle("DefectArray->subbath->state");

    int nbathspin = DefectArray_getNbathspin(dfa);
    printStructElementInt("nbathspin (#)",nbathspin);
    printLine();
    for (int i=0; i<nbathspin; i++){

        if (verbosity || (i<3 || i>nbathspin-3)){ 
            int naddspin = DefectArray_getNaddspins_i(dfa,i);
            float* states = allocFloat1d(naddspin);

            for (int j=0; j<naddspin; j++){
                BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,i,j);
                states[j] = BathSpin_getState(subspin);
            }

            char message[100];
            sprintf(message,"subbath[%d] : ",i);
            printStructElementFloat1d(message,states,naddspin);

            freeFloat1d(&states);
        }

        if (!verbosity && i==3){
            printf("         :\n");
        }
    }

    printf("\n");
    printLine();

}

void DefectArray_reportNaddspins(DefectArray* dfa){

    printSubTitle("DefectArray->naddspins");

    int nbathspin = DefectArray_getNbathspin(dfa);
    printStructElementInt("nbathspin (#)",nbathspin);
    printLine();
    for (int i=0; i<nbathspin; i++){

        if (verbosity || (i<3 || i>nbathspin-3)){ 
            char message[100];
            sprintf(message,"naddspins[%d]",i);
            printStructElementInt(message,DefectArray_getNaddspins_i(dfa,i));
        }

        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
    printf("\n");
    printLine();

}

void DefectArray_reportSubbath_i_props(DefectArray* dfa, int ibs, int isp){
    BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs,isp);
    char* name = BathSpin_getName(subspin);
    float spin = BathSpin_getSpin(subspin);
    double gyro = BathSpin_getGyro(subspin);
    double* xyz = BathSpin_getXyz(subspin);
    int mainspidx = BathSpin_getMainspidx(subspin);
    printf("      SubBath[%3d][%3d] %7s %7.3lf %7.3lf %7.3lf   ( S = %2.1f, gyro = %7.3lf , mainspidx = %d ) \n"\
    ,ibs,isp,name,xyz[0],xyz[1],xyz[2],spin,gyro, mainspidx);
}

void DefectArray_reportSubbath_hypfs(DefectArray* dfa, int nqubit){

    printSubTitle("SubBath Hyperfine tensors (Qubit-SubBath)");
    
    printMessage("SubBath[ibs][isp].hypf[iq] (radkHz)");
    printf("\n");
    printMessage(" * ibs : BathSpin index (main spin)");
    printMessage(" * isp : SubSpin index (sub spin)");
    printMessage(" * iq  : Qubit index");
    printf("\n");

    int nbathspin = DefectArray_getNbathspin(dfa);
    for (int ibs=0; ibs<nbathspin; ibs++){
        int naddspin = DefectArray_getNaddspins_i(dfa,ibs);

        if (verbosity || (ibs<3 || ibs>nbathspin-3)){ 
            for (int isp=0; isp<naddspin; isp++){
                BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs,isp);

                for (int iq=0; iq<nqubit; iq++){
                    MatrixXcd hypf = BathSpin_getHypf_i(subspin,iq);

                    char message[100];
                    sprintf(message,"SubBath[%d][%d].hypf[%d]",ibs,isp,iq);
                    printInlineMatrixXcd(message,hypf);
                }
            }
        }

        if (!verbosity && ibs==3){
            printf("         :\n");
        }
    }
}

void DefectArray_reportSubbath_quads(DefectArray* dfa){

    printSubTitle("SubBath Quadrupole tensors (Depending on principal axis)");

    int nbathspin = DefectArray_getNbathspin(dfa);
    for (int ibs=0; ibs<nbathspin; ibs++){
        int naddspin = DefectArray_getNaddspins_i(dfa,ibs);

        if (verbosity || (ibs<3 || ibs>nbathspin-3)){ 
            for (int isp=0; isp<naddspin; isp++){
                BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs,isp);
                MatrixXcd quad = BathSpin_getQuad(subspin);

                char message[100];
                sprintf(message,"SubBath[%d][%d].quad",ibs,isp);
                printInlineMatrixXcd(message,quad);
            }
        }

        if (!verbosity && ibs==3){
            printf("         :\n");
        }
    }
}

void DefectArray_reportSubbath_hypf_subs(DefectArray* dfa){

    printSubTitle("SubBath Hyperfine tensors (MainDefect-SubBath)");

    int nbathspin = DefectArray_getNbathspin(dfa);
    for (int ibs=0; ibs<nbathspin; ibs++){
        int naddspin = DefectArray_getNaddspins_i(dfa,ibs);

        if (verbosity || (ibs<3 || ibs>nbathspin-3)){ 
            for (int isp=0; isp<naddspin; isp++){
                BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs,isp);
                MatrixXcd hypf_sub = BathSpin_getHypfSub(subspin);

                char message[100];
                sprintf(message,"SubBath[%d][%d].hypf_sub",ibs,isp);
                printInlineMatrixXcd(message,hypf_sub);
            }
        }

        if (!verbosity && ibs==3){
            printf("         :\n");
        }
    }
}

void DefectArray_reportSubbath_disorders(DefectArray* dfa){

    printSubTitle("SubBath Disorders");

    int nbathspin = DefectArray_getNbathspin(dfa);
    for (int ibs=0; ibs<nbathspin; ibs++){
        int naddspin = DefectArray_getNaddspins_i(dfa,ibs);

        if (verbosity || (ibs<3 || ibs>nbathspin-3)){ 
            for (int isp=0; isp<naddspin; isp++){
                BathSpin* subspin = DefectArray_getSubbath_i_isp(dfa,ibs,isp);
                double disorder = BathSpin_getDisorder(subspin);

                char message[100];
                sprintf(message,"SubBath[%d][%d].disorder",ibs,isp);
                printStructElementDouble(message,disorder);
            }
        }

        if (!verbosity && ibs==3){
            printf("         :\n");
        }
    }
}

void DefectArray_reportSubbath_props(DefectArray* dfa){

    printSubTitle("DefectArray->subbath");

    int nbathspin = DefectArray_getNbathspin(dfa);
    printStructElementInt("nbathspin (#)",nbathspin);
    printLine();
    for (int i=0; i<nbathspin; i++){

        if (verbosity || (i<3 || i>nbathspin-3)){ 
            int naddspin = DefectArray_getNaddspins_i(dfa,i);
            for (int j=0; j<naddspin; j++){
                DefectArray_reportSubbath_i_props(dfa,i,j);
            }
        }

        if (!verbosity && i==3){
            printf("         :\n");
        }
    }

    printf("\n");
    printLine();

}



void DefectArray_reportAll(DefectArray* dfa){
    
}