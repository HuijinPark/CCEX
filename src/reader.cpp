#include "../include/reader.h"
#include "../include/memory.h"
#include "../include/utilities.h"
#include "../include/hamiltonian.h"

void readQubitfile(QubitArray* qa, Config* cnf){

    char* fname = Config_getQubitfile(cnf);
    double* xyz = QubitArray_getQubit_i_xyz(qa,0); // only add the info to 0-th qubit

    // file open
    FILE* data;

    data    =   fopen(fname,"r");

    if (data == NULL && QubitArray_getNqubit(qa) == 0){
        perror("Error(readQubitFile): Cannot open the qubit file\n");
        printf("Current WRONG qubitfile name : %s\n",fname);
        printf("Check if the qubit file exists\n");
        exit(EXIT_FAILURE);
    }else if (data == NULL && QubitArray_getNqubit(qa) != 0){
        if (rank==0){printMessage("QubitArray already set with option ...\n");}
        return;
    }else{
        if (rank==0){
            char message[MAX_FILEPATH];
            sprintf(message,"Read QubitFile : %s ... ",fname);
            printMessage(message);
        }
        
    }

    // read xyz from QubitFile
    size_t max = 1000;
    char lline[max];
    int count = 0;
    while(!feof(data)){
        fgets(lline, max, data);
        count = sscanf(lline, "%lf %lf %lf\n", &(xyz[0]), &(xyz[1]), &(xyz[2]));
        if (count == 3){
            fclose(data);

            if (rank==0){
                printStructElementDouble1d("  Qubit[0].xyz ", xyz, 3);
                printf("\n");
            }
            return;
        }
    }
    fclose(data);

    // if no data
    perror("Error(readQubitFile): Cannot read the qubit file\n");
    printf("Current WRONG qubitfile name : %s\n",fname);
    printf("Check if your qubit file has 'xyz' values in a right form : x y z \n");
    exit(EXIT_FAILURE);
}

void readBathfiles(BathArray* ba, QubitArray* qa, Config* cnf){

    char message[MAX_FILEPATH];

    int     nqubit      = QubitArray_getNqubit(qa);
    int     nbathfiles  = Config_getNbathfiles(cnf);
    double  rbath       = Config_getRbath(cnf);
    double  rbathcut    = Config_getRbathcut(cnf);
    
    // bath spin properties
    int     nspecies    = BathArray_getProp_nspecies(ba);
    char**  names       = BathArray_getProp_names(ba);
    float*  spins       = BathArray_getProp_spins(ba);
    double* gyros       = BathArray_getProp_gyros(ba);

    // Read bath files
    int nspin_prev = 0;
    int nspin = 0;    

    for (int i=0; i<nbathfiles; i++){

        FILE*   data;
        char*   fname       = Config_getBathfiles_i(cnf,i);
        double* bathadjust  = Config_getBathadjust_i(cnf,i);

        // file open        
        data = fopen(fname,"r");
        if(data==NULL){
            perror("Error(readBathFile): Cannot open the bath file\n");
            printf("Current WRONG bathfile name : %s\n",fname);
            printf("Check if the bath file exists\n");
            exit(EXIT_FAILURE);
        }

        // read xyz from BathFile
        double xyz[3] = {0.0,};
        char name[100];

        // Read first line (the number of spins)
        int fline = 1; // the first line is the number of spins
        fscanf(data, "%*lf %*lf %*lf %*s\n");

        while(!feof(data)){
        
            int count = fscanf(data, "%lf %lf %lf %s\n", &xyz[0], &xyz[1], &xyz[2], name);
            xyz[0] += bathadjust[0];
            xyz[1] += bathadjust[1];
            xyz[2] += bathadjust[2];
            fline++; // the line is started from 2

            if (count == 4){
                double r = QubitArray_mindist(xyz, qa);
                if ((r <= rbath) && (r >= rbathcut)){

                    //////////////////////////////
                    // increase the number of spins
                    nspin++;
                    BathArray_setNspin(ba, nspin);                
                    // allocate the BathArray->Bath
                    if (nspin == 0){
                        BathArray_allocBath(ba, nqubit); // alloc bath spins
                        Config_alloc_flines(cnf, nspin); // alloc fline
                    }else{
                        BathArray_reallocBath(ba, nspin-1, nspin, nqubit);
                        Config_realloc_flines(cnf, nspin-1, nspin);
                    }

                    // set the bath
                    BathArray_setBath_i_name(ba, name, nspin-1);
                    BathArray_setBath_i_xyz(ba, xyz, nspin-1);

                    // set fline
                    Config_set_nflines(cnf, nspin);
                    Config_set_flines_i(cnf, fline, nspin-1);

                    // set the bath spin properties
                    int ispeceis = findIndexChar(names,0,nspecies-1,name);

                    if (ispeceis == -1){
                        fprintf(stderr,"Error(readBathFile): Cannot find the species name in the gyro file\n");
                        fprintf(stderr,"The bath spin is : %s %lf %lf %lf \n",name,xyz[0],xyz[1],xyz[2]);
                        exit(EXIT_FAILURE);
                    }else{
                        BathArray_setBath_i_spin(ba, spins[ispeceis], nspin-1);
                        BathArray_setBath_i_gyro(ba, gyros[ispeceis], nspin-1);
                    }


                    // overlap..?
                }
            }
            else{
                perror("Error(readBathFile): BathFile format is NOT x y z name");
                exit(EXIT_FAILURE);
            }
        }
        fclose(data);

        ////////////////////////////////////////////////////////////////////////
        // Print
        if(rank==0){
            printLine();
            sprintf(message,"Read BathFile[%d] : %s\n",i,fname); printMessage(message);
            sprintf(message,"- Number of bath spins : %d (nspin = %d)",nspin-nspin_prev,nspin); printMessage(message);
            for (int isp=nspin_prev; isp<nspin; isp++){
                if (verbosity || (isp<nspin_prev + 3 || isp>nspin-3)){ 
                    // printf("        ( fline %5d )",Config_get_flines_i(cnf,isp));
                    BathArray_reportBath_i_props(ba, isp);
                }
                if (!verbosity && isp==nspin_prev + 3){
                    sprintf(message,"   : \n"); printMessage(message);
                }
            }
            printLine();

            nspin_prev = BathArray_getNspin(ba);
            printf("\n");
        }
        ////////////////////////////////////////////////////////////////////////
    }
    
    if (Config_get_nflines(cnf) != nspin){
        fprintf(stderr,"Error(readBathFile): The number of bath spins is not consistent\n");
        fprintf(stderr,"The number of bath spins from the file : %d\n",nspin);
        fprintf(stderr,"The number of bath spins from the file line : %d\n",Config_get_nflines(cnf));
        exit(EXIT_FAILURE);
    }

}

void setBathStates(BathArray* ba, Config* cnf, int i){

    /**
     * @param ba : BathArray
     * @param cnf : Config
     * @param i : index of the state file ( 1 <= i <= nstate )
     * @details 
     *  If i is zero, this calculation is for the ensemble approach
     *  Therefore, the bath states shouldn't be set
    */

    // Check the number of bath states if i is zero
    if (i == 0){
        // err
        if (rank==0){
            fprintf(stderr,"Error(setBathState): The index of the state file is zero\n");
            fprintf(stderr,"This calculation is for the ensemble approach\n");
            fprintf(stderr,"Therefore, the bath states shouldn't be set\n");
        }
        exit(EXIT_FAILURE);
    }

    // For print
    char message[MAX_FILEPATH];

    // Set filename as fname_i
    char fname[MAX_FILEPATH];
    sprintf(fname,"%s%d",Config_getStatefile(cnf),i);

    // Read the state file
    FILE* data;
    data = fopen(fname,"r");

    // Check the nbathfiles is larger than 1
    int nbathfiles = Config_getNbathfiles(cnf);
    if(data==NULL || nbathfiles > 1){
        
        if (rank==0){
            if (data == NULL){
                sprintf(message,"Warning(readStateFile): Cannot open the state file (%s)",fname); printMessage(message);
            }
            if (nbathfiles > 1){
                sprintf(message,"Warning(readStateFile): The number of bathfiles is larger than 1"); printMessage(message);
            }
            sprintf(message,"The bath state is randomly generated\n"); printMessage(message);
        }
        BathArray_setBathStatesRandom(ba);
        return;

    }else{

        // Read the number of lines
        int idx = 0;
        int count = fscanf(data, "%*s\n");
        int fline = 1; // the first line is the length of the file

        while(!feof(data)){

            //Read the bath state from the file
            //We will get the bath state at the same line of bathfile that we get the bath spins
            fline++;
            // Read the bath state
            float state = 0.0;
            count = fscanf(data, "%f\n", &state);

            if (fline == Config_get_flines_i(cnf,idx)){
                if (count == 1){
                    BathArray_setBath_i_state(ba, state, idx);
                    idx++;
                }else{
                    fprintf(stderr,"Error(readStateFile): StateFile format is NOT `float` \n");
                    exit(EXIT_FAILURE);
                }
            }
        }
        fclose(data);

        int nspin = BathArray_getNspin(ba);
        if (idx == nspin){
            if (rank==0){
                sprintf(message,"Read StateFile : %s\n",fname); printMessage(message);
            }
        }else{
            fprintf(stderr,"Error(readStateFile): The number of bath spins is not consistent\n");
            fprintf(stderr,"The number of bath spins from bathfile : %d\n",idx);
            fprintf(stderr,"The number of bath spins from the file line : %d\n",nspin);
            exit(EXIT_FAILURE);
        }
    }

    return;
}

void setDefectPaxes(DefectArray* dfa, BathArray* ba, Config* cnf){

    // For print
    char message[MAX_FILEPATH];

    // Read the state file
    char* fname = Config_getAvaaxfile(cnf);
    FILE* data;
    data = fopen(fname,"r");

    // Check the nbathfiles is larger than 1
    int nbathfiles = Config_getNbathfiles(cnf);

    if(data==NULL || nbathfiles > 1){
        
        if (rank==0){
            if (data == NULL){
                sprintf(message,"Warning(setPaxes): Cannot open the avaax file (%s)",fname); printMessage(message);
            }
            if (nbathfiles > 1){
                sprintf(message,"Warning(setPaxes): The number of bathfiles is larger than 1\n"); printMessage(message);
            }
            sprintf(message,"The paxes is randomly generated\n"); printMessage(message);
        }
        DefectArray_setPaxesRandom(dfa, ba);
        return;
    }

    // Read the number of lines
    int idx = 0;
    int count = fscanf(data, "%*s\n");
    int fline = 1; // the first line is the length of the file

    while(!feof(data)){

        // Read the bath axes from the file
        // We will get the bath axes at the same line of bathfile that we get the bath spins
        fline++;
        // Read the bath axes
        float paxestmp = 0;
        int paxes = 0;
        count = fscanf(data, "%f\n", &paxestmp);
        paxes = int(paxestmp);
        if (fline == Config_get_flines_i(cnf,idx)){
            
            if (count == 1){
                //////////////////////////////////////////////
                char* dfname = BathArray_getBath_i_name(ba,idx);
                int idf = DefectArray_findDefectIndex(dfa,dfname);
                if (idf == -1){
                    // err if the defect name is not found
                    fprintf(stderr,"Error(setPaxes): Cannot find the defect name in the defect array\n");
                    fprintf(stderr,"The bath spin[%d] is : %s \n",idx,dfname);
                    exit(EXIT_FAILURE);
                }
                int navaax = DefectArray_getDefect_idf_navaax(dfa,idf);
                if (paxes > navaax || paxes < 1){
                    // err if the principal axis is out of range
                    fprintf(stderr,"Error(setPaxes): The principal axis is out of range\n");
                    fprintf(stderr,"The bath spin[%d] is : %s \n",idx,dfname);
                    fprintf(stderr,"The number of principal axis : %d \n",navaax);
                    fprintf(stderr,"The principal axis from the file : %d \n",paxes);
                    exit(EXIT_FAILURE);
                }
                //////////////////////////////////////////////
                DefectArray_setPaxes_i(dfa, idx, paxes);
                idx++;
            }else{
                fprintf(stderr,"Error(setPaxes): AvaaxFile format is NOT `int` \n");
                exit(EXIT_FAILURE);
            }
        }

    }
    fclose(data);

    // Check if the number of bath spins is consistent
    if (Config_get_nflines(cnf) == idx){
        if (rank==0){
            sprintf(message,"Read AvaaxFile : %s\n",fname); printMessage(message);
        }
    }else{
        fprintf(stderr,"Error(setPaxes): The number of bath spins is not consistent\n");
        fprintf(stderr,"The number of bath spins from bathfile : %d\n",idx);
        fprintf(stderr,"The number of bath spins from the file line : %d\n",Config_get_nflines(cnf));
        exit(EXIT_FAILURE);
    }

    return;
}

void setSubbathStates(DefectArray* dfa, BathArray* ba, Config* cnf, int i){

    /**
     * @param ba : BathArray
     * @param dfa : DefectArray
     * @param cnf : Config
     * @param i : index of the state file ( 1 <= i <= nstate )
     * @details 
     *  If i is zero, this calculation is for the ensemble approach
     *  Therefore, the bath states shouldn't be set
    */

    // Check the number of bath states if i is zero
    if (i == 0){
        // err
        if (rank==0){
            fprintf(stderr,"Error(setBathState): The index of the state file is zero\n");
            fprintf(stderr,"This calculation is for the ensemble approach\n");
            fprintf(stderr,"Therefore, the bath states shouldn't be set\n");
        }
        exit(EXIT_FAILURE);
    }

    // For print
    char message[MAX_FILEPATH];

    // File
    char fname[MAX_FILEPATH];
    sprintf(fname,"%s%d",Config_getExstatefile(cnf),i);

    // Read the state file
    FILE* data;
    data = fopen(fname,"r");

    // Check the nbathfiles is larger than 1
    int nbathfiles = Config_getNbathfiles(cnf);
    if (data == NULL || nbathfiles > 1){
        
        if (rank==0){
            if (data == NULL){
                
                sprintf(message,"Warning(setSubbathStates): Cannot open the exstate file (%s)",fname); printMessage(message);
            }
            if (nbathfiles > 1){
                sprintf(message,"Warning(setSubbathStates): The number of bathfiles is larger than 1"); printMessage(message);
            }
            sprintf(message,"The subbath states are randomly generated\n"); printMessage(message);
        }
        DefectArray_setSubbathStatesRandom(dfa, ba);
        return;
    }else{
            
        // Read the number of lines
        int idx = 0;
        int count = fscanf(data, "%*s\n");
        int fline = 1; // the first line is the length of the file

        while(!feof(data)){

            // Read the bath state from the file
            // We will get the bath state at the same line of bathfile that we get the bath spins
            fline++;
            // Read the bath state
            const int max = 1000;
            char lline[max];
            count = fscanf(data, "%s\n", lline);

            if (fline == Config_get_flines_i(cnf,idx)){
                int naddspin = DefectArray_getNaddspins_i(dfa,idx);

                if (naddspin != 0){

                    char* ptr = strtok(lline," ");
                    int isp=0;
                     
                    while (ptr !=NULL){
                        //state add & check
                        float ms = float(atof(ptr));
                        BathSpin* bs = DefectArray_getSubbath_i_isp(dfa,idx,isp);
                        float  S  = BathSpin_getSpin(bs);
                        BathSpin_setState(bs,ms); // Check if the state is valid in this function
                        ptr = strtok(NULL," ");
                        isp++;
                    }

                    if (isp != naddspin){
                        fprintf(stderr,"Error(setSubbathStates): The number of subbath spins is not consistent\n");
                        fprintf(stderr,"The number of subbath spins from the file line : %d\n",isp);
                        fprintf(stderr,"The number of subbath spins from the defect array : %d\n",naddspin);
                        exit(EXIT_FAILURE);
                    }
                }
                idx++;
            }
        }
        fclose(data);

        int nspin = BathArray_getNspin(ba);
        if (idx == nspin){
            if (rank==0){
                sprintf(message,"Read ExstateFile : %s\n",fname); printMessage(message);
            }
        }else{
            fprintf(stderr,"Error(setSubbathStates): The number of bath spins is not consistent\n");
            fprintf(stderr,"The number of bath spins from bathfile : %d\n",idx);
            fprintf(stderr,"The number of bath spins from the file line : %d\n",nspin);
            exit(EXIT_FAILURE);
        }
    }


    return;
}

void readGyrofile(BathArray* ba, Config* cnf){
    
    char* fname = Config_getGyrofile(cnf);

    FILE* data;
    data = fopen(fname,"r");

    if(data==NULL){
        fprintf(stderr,"Error(readGyroFile): Cannot open the gyro file\n"); 
        fprintf(stderr,"Current WRONG gyrofile name : %s\n",fname);
        fprintf(stderr,"Check if the gyro file exists\n");
        exit(EXIT_FAILURE);
    }

    // read gyrofile set spin properties
    int nspecies = 0;
    while(!feof(data)){
        char    name[MAX_CHARARRAY_LENGTH] = "\0";
        float   spin = 0.0;
        double  gyro = 0.0;
        int     count = fscanf(data, "%s %f %lf\n", name, &spin, &gyro);
        if (count == 3){
            nspecies++;
            BathArray_setProp_nspecies(ba, nspecies);
            if (nspecies == 0){ 
                BathArray_allocProp(ba); }
            else{
                BathArray_reallocProp(ba, nspecies-1, nspecies);
            }
            BathArray_setProp_names_i(ba, name, nspecies-1);
            BathArray_setProp_spins_i(ba, spin, nspecies-1);
            BathArray_setProp_gyros_i(ba, gyro, nspecies-1);
        }
        else{
            fprintf(stderr,"Error(readGyroFile): GyroFile format is NOT name spin gyro\n");
            exit(EXIT_FAILURE);
        }
    }
    fclose(data);
    if(rank==0){printf("        Read GyroFile : %s\n",fname);}
    if(rank==0){BathArray_reportSpinProperties(ba);}
    if(rank==0){printf("\n");}
}

void readHftensorfile(BathArray* ba, QubitArray* qa, Config* cnf){

    // Required bath spin properties 
    int nspecies = BathArray_getProp_nspecies(ba);
    char** names = BathArray_getProp_names(ba);

    // hf-related parameters
    int hf_readmode = Config_getHf_readmode(cnf);
    char* hf_tensorfile = Config_getHf_tensorfile(cnf);
    double hf_cutoff = Config_getHf_cutoff(cnf);
    int hf_ignore_oor = Config_getHf_ignore_oor(cnf);

    double DefectTotSpin = Config_getDefectTotSpin(cnf);
    double CorrTotSpin = Config_getCorrTotSpin(cnf);
    double SpinFactor = 0.0;

    if (rank==0){
        printSubTitle("Read the Hyperfine interaction from DFT inputfile...");
    }

    //the information in Atensor file
    if (hf_readmode==0){
        printf("      %-18s:   %4d ( Point-dipole tensor )\n\n", "HF Readmode ", hf_readmode);
        BathArray_setBathHypfs(ba,qa); 
    }
    else if (hf_readmode==1 || hf_readmode==2 || hf_readmode==3){

        if (rank==0){

            if (hf_readmode==1){
                printf("      %-18s:   %d ( Fermi-contact term + point-dipole tensor )\n", "HF Readmode ", hf_readmode);
            }
            else if (hf_readmode==2){
                printf("      %-18s:   %d ( DFT dipolar tensor )\n", "HF Readmode ", hf_readmode);
            }
            else if (hf_readmode==3){
                printf("      %-18s:   %d ( Fermi-contact term + DFT dipole tensor )\n", "HF Readmode ", hf_readmode);
            }
            printStructElementChar("HF Tensor file ", hf_tensorfile);
            printStructElementDouble("HF Cutoff ", hf_cutoff);
            printStructElementInt("HF ignoring range ", hf_ignore_oor);
            printf("\n");

        }

        //////////////////////////////////////////////////////////////////
        // Tensor information
        //////////////////////////////////////////////////////////////////

        // main information
        double*  A_Gfactor;     //g-factor___ ;
        double** A_Etc;         //etc____ ; Axx, Axy, Axz,..., Ayz, Azz
        double** AtensorArray;  //x,y,z,Axx,Axy,Axz,...,Ayz,Azz in wiDefect file

        // Boundary condition (vertex version)
        double** A_vertex;      //vertex about A-tensor
        double** A_center;      //center of plane about A-tensor
        double** A_normal;      //normal vector of plane about A-tensor
        char* vertex_condi[8] = {"v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"};

        // Boundary condition (xyz version)
        double  A_MinDif[3] = {0.0,0.0,0.0};    // MinDif[A]
        double  A_MaxDif[3] = {0.0,0.0,0.0};    // MaxDif[A]

        // Check if the BD is used with vertex or xyz version
        bool isVertex = false;

        ////////////////////////////////////////////////////////////
        // Read the tensor files
        ////////////////////////////////////////////////////////////
    
        int version = READ_Tensor_ver(hf_tensorfile,&SpinFactor ,DefectTotSpin, &CorrTotSpin);

        if (rank==0){
            printHfInfo_version(version, DefectTotSpin, CorrTotSpin, SpinFactor);
        }

        //the information to checking BD (in A-&Q-tensor)
        isVertex = \
        READ_BD_vertex(hf_tensorfile, &(A_vertex), &(A_center), &(A_normal), (vertex_condi));

        if (!isVertex){
            bool Check_min = READ_Tensor_array(hf_tensorfile,&(A_MinDif),"MinDif[A]",3);
            bool Check_max = READ_Tensor_array(hf_tensorfile,&(A_MaxDif),"MaxDif[A]",3);
            if ((Check_min==false) || (Check_max==false)){
                printf("\tThe Boundary Condition (vertex and Range) can not read from tensor file (%s)\n",hf_tensorfile);
                printf("\tPlz, Check the Tensor file !!!\n");
                exit(1);
            }
        }

        if (rank==0 && verbosity){
            printHfInfo_BD(A_vertex, A_center, A_normal, A_MinDif, A_MaxDif, isVertex);
        }

        // Read gfactor and etc
        READ_Tensor_const(hf_tensorfile,names,nspecies,&(A_Gfactor),"g-factor___"); // length = nspecies+1
        READ_Tensor_etc(hf_tensorfile,names,nspecies,&(A_Etc),"etc____",1); // length = nspecies+1

        if (rank==0 && verbosity){
            printHfInfo_etc(A_Etc, A_Gfactor, names, nspecies, 0);
        }

        // Read the tensors
        READ_Tensor(hf_tensorfile,&(AtensorArray),13); // length = AtensorArray[0][0]

        if (rank==0 && verbosity){
            printHfInfo_tensor(AtensorArray, 0);
        }
        
        ////////////////////////////////////////////////////////////
        // Set Hyperfine tensor to bath array
        ////////////////////////////////////////////////////////////
        int nspin = BathArray_getNspin(ba);
        int nqubit = QubitArray_getNqubit(qa);

        if (nqubit > 1){
            fprintf(stderr,"Error(readHftensorfile): The number of qubit is more than 1\n");
            fprintf(stderr,"In this case, DONOT use Atensorfile\n");
            fprintf(stderr,"Instead of that, put the tensor directly.\n");
            fprintf(stderr,"Recommanded keyword is Bath : {}, use cce.in as json file\n");
            exit(EXIT_FAILURE);
        }

        for (int iqubit=0; iqubit<nqubit; iqubit++){
            for (int ispin=0; ispin<nspin; ispin++){

                // resultant paramter
                double fc = 0.0;
                MatrixXcd Adip = MatrixXcd::Zero(3,3);
                MatrixXcd Atot = MatrixXcd::Zero(3,3);

                // qubit information
                float qspin = QubitArray_getQubit_i_spin(qa,iqubit);
                double qgyro = QubitArray_getQubit_i_gyro(qa,iqubit);
                double* qxyz = QubitArray_getQubit_i_xyz(qa,iqubit);

                if (DefectTotSpin != qspin){
                    fprintf(stderr,"Error(readHftensorfile): The defect spin is not matched with the qubit spin\n");
                    fprintf(stderr,"DefectTotSpin : %lf\n",DefectTotSpin);
                    fprintf(stderr,"QubitSpin : %lf\n",qspin);
                    exit(EXIT_FAILURE);
                }

                // bath spin information
                double spgyro = BathArray_getBath_i_gyro(ba,ispin);
                double* spxyz = BathArray_getBath_i_xyz(ba,ispin);
                char* spname = BathArray_getBath_i_name(ba,ispin);

                // g-factor of the bath spin
                int speciesIdx = findIndexChar(names,0,nspecies-1,spname);
                double gfactor = A_Gfactor[speciesIdx+1]; // length = nspecies+1

                // relative distance between qubit and bath spin
                double rxyz[3];
                rxyz[0] = spxyz[0] - qxyz[0];
                rxyz[1] = spxyz[1] - qxyz[1];
                rxyz[2] = spxyz[2] - qxyz[2];
                
                // Check if the bath spin is in the maximum, minimum range
                bool isExist = false;
                double err = 0.001; // tolerance for the range

                if (isVertex){
                    isExist = CheckBD_vertex(rxyz, A_vertex, A_center, A_normal, err);
                }else{
                    isExist = CheckBD_Range(rxyz, A_MinDif, A_MaxDif, err);
                }

                // Check if the bath spin is matched with the relative distance of the file
                bool isFound = false;
                int tensorIdx = 0;

                if (isExist){
                    //AtesnorArray[][0~2]  : atomic position
                    isFound = FIND_AtomPosi(rxyz, AtensorArray, &tensorIdx);
                }

                // read the tensor information
                if (isExist && isFound){
                    //AtesnorArray[][3]    : fermi contact term
                    //AtesnorArray[][4~12] : Atensor value
                    fc = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][3] * gfactor * SpinFactor);
                    Adip(0,0) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][4] * gfactor * SpinFactor);
                    Adip(0,1) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][5] * gfactor * SpinFactor);
                    Adip(0,2) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][6] * gfactor * SpinFactor);
                    Adip(1,0) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][7] * gfactor * SpinFactor);
                    Adip(1,1) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][8] * gfactor * SpinFactor);
                    Adip(1,2) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][9] * gfactor * SpinFactor);
                    Adip(2,0) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][10] * gfactor * SpinFactor);
                    Adip(2,1) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][11] * gfactor * SpinFactor);
                    Adip(2,2) = MHZ_TO_RADKHZ(AtensorArray[tensorIdx][12] * gfactor * SpinFactor);

                    if (hf_readmode == 1){
                        // fermi contact term + Point-dipole approximation
                        Adip = calPointDipoleTensor(qxyz, spxyz, qgyro, spgyro);
                        Atot = Adip + fc * MatrixXcd::Identity(3,3);
                    }else if (hf_readmode == 2){
                        // DFT dipolar tensor
                        Atot = Adip;
                    }else if (hf_readmode == 3){
                        // fermi contact term + DFT dipolar tensor
                        Atot = Adip + fc * MatrixXcd::Identity(3,3);
                    }
                    BathArray_setBath_i_hypf_j(ba, Atot, ispin, iqubit);

                }else if (isExist && !isFound){
                    if (hf_ignore_oor==0){
                        //not Ignore the mis-match
                        fprintf(stderr,"Error(readHftensorfile): The nuclear spin is within the range, but doesn't exist in A-file!!\n");
                        fprintf(stderr,"A_MinDif : %lf %lf %lf\n",A_MinDif[0],A_MinDif[1],A_MinDif[2]);
                        fprintf(stderr,"A_MaxDif : %lf %lf %lf\n",A_MaxDif[0],A_MaxDif[1],A_MaxDif[2]);
                        fprintf(stderr,"difXYZ : %lf %lf %lf\n", rxyz[0],rxyz[1],rxyz[2]);
                        fprintf(stderr,"Nuclear spin : %lf %lf %lf\n",spxyz[0],spxyz[1],spxyz[2]);
                        fprintf(stderr,"In the A-file --> %d-line is the same!!", tensorIdx);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }

        ////////////////////////////////////////////////////////////
        // AInfo memory free
        ////////////////////////////////////////////////////////////

        int length_AtensorArray = (int)AtensorArray[0][0];
        int length_others = nspecies+1;
        
        freeDouble1d(&A_Gfactor);
        freeDouble2d(&A_Etc,length_others);
        freeDouble2d(&AtensorArray, length_AtensorArray);

        if (isVertex){
            freeDouble2d(&A_vertex, 8);
            freeDouble2d(&A_center, 8);
            freeDouble2d(&A_normal, 8);
        }

    }

    if (rank==0){
        int nqubit = QubitArray_getNqubit(qa);
        BathArray_reportBath_hypf(ba, nqubit); //report updated hyperfine tensors
    }
}

void readQdtensorfile(BathArray* ba, QubitArray* qa, Config* cnf){

    // Required bath spin properties 
    int nspecies = BathArray_getProp_nspecies(ba);
    char** names = BathArray_getProp_names(ba);

    // qd-related parameters
    int qd_readmode = Config_getQd_readmode(cnf);
    char* qd_tensorfile = Config_getQd_tensorfile(cnf);
    char* qd_tensorfile_woqubit = Config_getQd_tensorfile_woqubit(cnf);

    //the information in Atensor file
    if (qd_readmode == 0){
        ; // ignore the quadrupole interaction
    }
    else if (qd_readmode == 1 || qd_readmode == 2 || qd_readmode == 3 || qd_readmode == 4){ 
        ; // to be added
    }
}


// void Hamil_quad(MatrixXcd** H_nQ , MatrixXcd** Pmat, double** spPos, double spInfo[][3], Gyro* Gy, Hamil_data* Hdata, INFOBATH* Info){
//     //unit : [rad/msec]

// //    printf("Quadrupole calculation\n");
//     int line = 0;

//     //1) input the Hamiltonian in each array 0 to n-1
//     for (int i = 1; i < (int)spInfo[0][0]; i++){
//         for (int j = 1; j < (int)spInfo[0][0]; j++){
//             int    mSize = (int)spInfo[j][2];
//             double gamma = spInfo[j][0];
//             double numI  = spInfo[j][1]; 
//             MatrixXcd tmp(mSize,mSize);
//             if (j == i){
//                 if (numI > 0.50){
//                     switch (Hdata->Hamil_Quad){
//                     case 0:
//                         //ignore the quadrupole interaction
//                         break;
//                     case 1:
//                         Hamil_quad_mkexp(&tmp,Pmat[j-1],spInfo[j],Hdata,spPos[j],Gy);
//                         //Hamil_quad_mkexp(&tmp,Pmat[j-1],spInfo[j],Hdata,spPos[j];
//                         break;
//                     case 2: case 3: case 4:
//                         {
//                                MatrixXcd EFG(3,3);
//                                EFG = MatrixXcd::Zero(3,3);    
//                                MatrixXcd tmp_xx(mSize,mSize),tmp_xy(mSize,mSize),tmp_xz(mSize,mSize);
//                                MatrixXcd tmp_yx(mSize,mSize),tmp_yy(mSize,mSize),tmp_yz(mSize,mSize);
//                                MatrixXcd tmp_zx(mSize,mSize),tmp_zy(mSize,mSize),tmp_zz(mSize,mSize);

//                             double eQ   = 0.0; //Value Would be added;   //eQ = 20.44[Q/millibarn = Q * 10^-31 m^2] , in DFT eQ = 2.044 * 10^-30 [m^2]
//                             Hamil_quad_mkdft(&EFG,&eQ,Gy,qbInfo,spPos[j],Hdata,Info);
//                             //double quadConst = eQ * pow(10,-30)/((6*numI) * (2*numI-1)) * 27.211386 /pow((0.5291772 * pow(10,-10)),2)/(4.135667 * pow(10,-15)); //Hz
//                             //quadConst = quadConst * pow(10,-6);  //Hz to MHz
//                             //quadConst = quadConst * 2 * M_PI * pow(10,3); // MHz to rad/msec
//                             //
//                             //double quadConst = eQ * 1e-30 / ((6*numI) * (2*numI-1)) * 27.211386 /pow((0.5291772 * 1e-10),2)/(4.135667 * 1e-15); //Hz
//                             double quadConst = eQ * 1e+5 / ((6*numI) * (2*numI-1)) * 27.211386 /pow(0.5291772,2)/(4.135667); //Hz
//                             quadConst = quadConst * 1e-6;  //Hz to MHz
//                             quadConst = quadConst * 2 * M_PI * 1e+3; // MHz to rad/msec



//                             ////Checking for quadrupole parameter
//                             //double quadConst_test = eQ * pow(10,-30) * 27.211386 /pow((0.5291772 * pow(10,-10)),2)/(4.135667 * pow(10,-15)); //Hz
//                             ////quadConst_test = quadConst_test * pow(10,-6); //Hz to MHz;
//                             //quadConst_test = quadConst_test * 1e-6; //Hz to MHz;
//                             //std::cout << "Quad : Vzz [Hartree/Bohr_radius^2]  = "<< EFG(2,2) << std::endl;
//                             //std::cout << "Quad : Zeta [MHz] = "<<quadConst_test*(EFG(2,2) - (EFG(0,0)+EFG(1,1))/2) << std::endl;
//                             //std::cout << "Quad : Cq [MHz]   = "<< quadConst_test*EFG(2,2) << std::endl;
//                             //printf("\n");


//                             tmp_xx = EFG(0,0)*(1.5*(Pmat[j-1][1]*Pmat[j-1][1]*2) - numI*(numI+1)*Pmat[j-1][0]); //
//                             tmp_yy = EFG(1,1)*(1.5*(Pmat[j-1][2]*Pmat[j-1][2]*2) - numI*(numI+1)*Pmat[j-1][0]);
//                             tmp_zz = EFG(2,2)*(1.5*(Pmat[j-1][3]*Pmat[j-1][3]*2) - numI*(numI+1)*Pmat[j-1][0]);

//                             tmp_xy = EFG(0,1)*(1.5*((Pmat[j-1][1]*Pmat[j-1][2])+(Pmat[j-1][2]*Pmat[j-1][1])));
//                             tmp_yx = EFG(1,0)*(1.5*((Pmat[j-1][2]*Pmat[j-1][1])+(Pmat[j-1][1]*Pmat[j-1][2])));
                            
//                             tmp_xz = EFG(0,2)*(1.5*((Pmat[j-1][1]*Pmat[j-1][3])+(Pmat[j-1][3]*Pmat[j-1][1])));
//                             tmp_zx = EFG(2,0)*(1.5*((Pmat[j-1][3]*Pmat[j-1][1])+(Pmat[j-1][1]*Pmat[j-1][3])));

//                             tmp_yz = EFG(1,2)*(1.5*((Pmat[j-1][2]*Pmat[j-1][3])+(Pmat[j-1][3]*Pmat[j-1][2])));
//                             tmp_zy = EFG(2,1)*(1.5*((Pmat[j-1][3]*Pmat[j-1][2])+(Pmat[j-1][2]*Pmat[j-1][3])));

//                             tmp = quadConst * ( tmp_xx + tmp_xy + tmp_xz 
//                                               + tmp_yx + tmp_yy + tmp_yz
//                                               + tmp_zx + tmp_zy + tmp_zz);

// //                            printf("spinNum = %lf, quad_mkdft, eQ = %lf, quadConst = %lf\n",numI, eQ, quadConst);
//                         }
//                         break;
//                     default:
//                         printf("error, Hdata->Hamil_Quad(%d) is not neither 0,1,2,3,4!!\n",Hdata->Hamil_Quad);exit(1);
//                         //printf("error, Hdata->Hamil_Quad(%d) is not neither 0,1,2\n",Hdata->Hamil_Quad);exit(1);
//                         break;
//                     }

//                 }else if (numI == 0.50){
//                     tmp = MatrixXcd::Zero(mSize,mSize);
// //                    printf("spinNum = %lf, no quadrupole term\n",numI);

//                 }else{
//                     printf("Error(Hamil_quad) spin is something else not 1/2 or more than 1/2\n");exit(1);
//                 }

//             }else{
//                 tmp = Pmat[j-1][0];
//             }

//             if (j == 1){
//                 (*H_nQ)[line] = tmp;
//             }else if(j > 1){
//                 (*H_nQ)[line] = kron((*H_nQ)[line],tmp);
//             }
// //            printf("(%dx%d) = \n",tmp.rows(),tmp.cols());
// //            std::cout<<tmp<<std::endl;
//         }
// //        printf("H_nQ[%d] =(%dx%d)\n",line,((*H_nQ)[line]).rows(),((*H_nQ)[line]).cols());
// //        std::cout<<(*H_nQ)[line]<<std::endl;
//         line++;
//     }

//     if (line != (int)spInfo[0][0]-1){printf("Error(Hamil_quad)\n");exit(1);}

// }

// void Hamil_quad_mkexp(MatrixXcd* tmp, MatrixXcd* Pmat, double spInfo[3], Hamil_data* Hdata, double spin2[],Gyro* Gy){
//     //unit : [rad/msec]

//     int i = (int)spin2[3];
//     double numI  = spInfo[1]; 

//     if (i >= Hdata->Quad_exp[0][0] ){
//         printf(" i & Hdata->Quad_exp[0][0] : %d %d\n", i, Hdata->Quad_exp[0][0]);
//         printf("There are some wrong experiment Quadrupole moment value\n");
//         printf("The Atom name : %s\n",Gy->Atom_name[int(spin2[3])]);
//         for (int j=1; j<Hdata->Quad_exp[0][0]; j++){
//             printf("Quad_exp(%d) : %lf(P||), %lf, %lf\n",j, Hdata->Quad_exp[j][0], Hdata->Quad_exp[j][1], Hdata->Quad_exp[j][2]);
//         }
//         printf("if the atom has no Quadrupole moment value,\n");
//         printf("insert Quad_exp(xx) = 0!!\n");
//         exit(1);
//     }

//     //1.only Vzz
//     ////double Q_zz  =   -1 * Hdata->Quad_exp[i] * pow(10,3) * 2 * M_PI; // MHz to rad/msec
//     //double Q_zz  =   -1 * Hdata->Quad_exp[i] * 1e+3 * 2 * M_PI; // MHz to rad/msec
//     //(*tmp) = Q_zz * Pmat[3] * Pmat[3];
//     //
//     ////double Pzz  = Hdata->Quad_exp[i][0] * pow(10,3)*2*M_PI;
//     //double Pzz  = Hdata->Quad_exp[i][0] * 1e+3 *2*M_PI;
//     //(*tmp) = (Pzz/3) * (3 *Pmat[3] * Pmat[3]);
//     //
//     //////////////////////////////////////////////

//     //2.diagonal term; Vxx & Vyy & Vzz
//     ////double Q_xx  =   -1 * Hdata->Quad_exp[i][0] * pow(10,3) * 2 * M_PI;
//     ////double Q_yy  =   -1 * Hdata->Quad_exp[i][1] * pow(10,3) * 2 * M_PI;
//     ////double Q_zz  =   -1 * Hdata->Quad_exp[i][2] * pow(10,3) * 2 * M_PI;
//     //double Q_xx  =   -1 * Hdata->Quad_exp[i][0] * 1e+3 * 2 * M_PI;
//     //double Q_yy  =   -1 * Hdata->Quad_exp[i][1] * 1e+3 * 2 * M_PI;
//     //double Q_zz  =   -1 * Hdata->Quad_exp[i][2] * 1e+3 * 2 * M_PI;
//     //(*tmp) = Q_xx * Pmat[1] * Pmat[1] + Q_yy * Pmat[2] * Pmat[2] + Q_zz * Pmat[3] * Pmat[3];
//     //
//     //////////////////////////////////////////////

    
//     //3. using P|| & asymmetry parameter
//     ////double Pzz = 1 * Hdata->Quad_exp[i][0] * pow(10,3) * 2 * M_PI; //unit: MHz -> rad/msec
//     //double Pzz = 1 * Hdata->Quad_exp[i][0] * 1e+3 * 2 * M_PI; //unit: MHz -> rad/msec
//     //double asym = Hdata->Quad_exp[i][1]; //no unit
//     //(*tmp) = (Pzz/3) * ( 3*Pmat[3]*Pmat[3] - numI*(numI+1)*Pmat[0] + asym * (Pmat[1]*Pmat[1] - Pmat[2]*Pmat[2]) );
//     //
//     //eq : eQ/(4I*(2I-1))*[Vzz(3I_z^2-I^2) + (V_xx - V_yy)(I_x^2 - I_y^2)]
//     //////////////////////////////////////////////
    
//     //4. using P|| & asymmetry parameter
//     //double Pzz = 1 * Hdata->Quad_exp[i][2] * pow(10,3) * 2 * M_PI; //unit: MHz -> rad/msec
//     double Pzz = 1 * Hdata->Quad_exp[i][2] * 1e+3 * 2 * M_PI; //unit: MHz -> rad/msec
//     double VxxPVzz = Hdata->Quad_exp[i][0] ; //no unit
//     double VyyPVzz = Hdata->Quad_exp[i][1] ; //no unit
    
//     (*tmp) = (Pzz/3) * ( (3*Pmat[3]*Pmat[3] - numI*(numI+1)*Pmat[0]) + VxxPVzz * (3*Pmat[1]*Pmat[1] -numI*(numI+1)*Pmat[0]) +VyyPVzz * (3*Pmat[1]*Pmat[1] -numI*(numI+1)*Pmat[0]));
//     //
//     //eq : eQ/(4I*(2I-1))*[Vzz(3I_z^2-I^2) + Vxx(3I_z^2-I^2)+ Vyy(3I_z^2-I^2)]
//     //////////////////////////////////////////////
// }

// void Hamil_quad_mkdft(MatrixXcd* H_nEFG,double* eQ, Gyro* Gy, double spin1[], double spin2[], Hamil_data* Hdata, INFOBATH* Info){
//     //unit : [rad/msec]
//     //spin1 : qubit , spin2 : nuclear spin bath, The difference is spin2 - spin1 

//     //Calculation difference 
//     double difXYZ[3];
//     difXYZ[0] = spin2[0]-spin1[0];
//     difXYZ[1] = spin2[1]-spin1[1];
//     difXYZ[2] = spin2[2]-spin1[2];
//     //printf("Quad : difXYZ(%s) = (%lf, %lf, %lf)\n",isoAtom,difXYZ[0],difXYZ[1],difXYZ[2]);

//     //about isoAtom
//     //char isoAtom[strlen(Gy->Atom_name[int(spin2[3])])+1];
//     //strcpy(isoAtom, Gy->Atom_name[int(spin2[3])]);
//     //printf("isoAtom : %s\n",isoAtom);
    
//     //using below four function!!
//     //1. CheckRange     :   Check the difXYZ is in the Range
//     //2. FIND_AtomPosi  :   Find and Check the index in Tensor using difXYZ
    
//     (*eQ) = Hdata->Q_eQ[int(spin2[3])];
//     int num=0; //index num of tensor value
//     switch (Hdata->Hamil_Quad){
//     case 2:
//         //only consider the defect cell Q-tensor
//         //about Quadrupole file with Defect 
//         //if(CheckRange(difXYZ,Hdata->Q_MinDifwiDef,Hdata->Q_MaxDifwiDef))
//         if(CheckRange(difXYZ,Hdata,1)){
//             if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWiDef,&num)){
//                 (*H_nEFG)(0,0) = Hdata->QtensorArrayWiDef[num][3];   (*H_nEFG)(0,1) = Hdata->QtensorArrayWiDef[num][4];   (*H_nEFG)(0,2) = Hdata->QtensorArrayWiDef[num][5];
//                 (*H_nEFG)(1,0) = Hdata->QtensorArrayWiDef[num][6];   (*H_nEFG)(1,1) = Hdata->QtensorArrayWiDef[num][7];   (*H_nEFG)(1,2) = Hdata->QtensorArrayWiDef[num][8];
//                 (*H_nEFG)(2,0) = Hdata->QtensorArrayWiDef[num][9];   (*H_nEFG)(2,1) = Hdata->QtensorArrayWiDef[num][10];  (*H_nEFG)(2,2) = Hdata->QtensorArrayWiDef[num][11];

//                 //(*H_nEFG)(0,0) = Mat[0][0]; (*H_nEFG)(0,1) = Mat[0][1]; (*H_nEFG)(0,2) = Mat[0][2];
//                 //(*H_nEFG)(1,0) = Mat[1][0]; (*H_nEFG)(1,1) = Mat[1][1]; (*H_nEFG)(1,2) = Mat[1][2];
//                 //(*H_nEFG)(2,0) = Mat[2][0]; (*H_nEFG)(2,1) = Mat[2][1]; (*H_nEFG)(2,2) = Mat[2][2];
//             }else{
//                 printf("error in EFG(Opt:%d), the atom is within the range but doesn't exist in defect EFG file \n",Hdata->Hamil_Quad);
//                 exit(1);
//             }
//         //out of range atom
//         //about etc value in quadrupole file
//         }else{
//             (*H_nEFG)(0,0) = Hdata->Q_Etc[int(spin2[3])][0];     (*H_nEFG)(0,1) = Hdata->Q_Etc[int(spin2[3])][1];     (*H_nEFG)(0,2) = Hdata->Q_Etc[int(spin2[3])][2];
//             (*H_nEFG)(1,0) = Hdata->Q_Etc[int(spin2[3])][3];     (*H_nEFG)(1,1) = Hdata->Q_Etc[int(spin2[3])][4];     (*H_nEFG)(1,2) = Hdata->Q_Etc[int(spin2[3])][5];
//             (*H_nEFG)(2,0) = Hdata->Q_Etc[int(spin2[3])][6];     (*H_nEFG)(2,1) = Hdata->Q_Etc[int(spin2[3])][7];     (*H_nEFG)(2,2) = Hdata->Q_Etc[int(spin2[3])][8];
//         }
//         break;

//     case 3:
//         //consider the defect cell Q-tensor
//         //about Quadrupole file with Defect 
//         //if(CheckRange(difXYZ,Hdata->Q_MinDifwiDef,Hdata->Q_MaxDifwiDef))
//         if(CheckRange(difXYZ,Hdata,1)){
//             if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWiDef,&num)){
//                 (*H_nEFG)(0,0) = Hdata->QtensorArrayWiDef[num][3];       (*H_nEFG)(0,1) = Hdata->QtensorArrayWiDef[num][4];       (*H_nEFG)(0,2) = Hdata->QtensorArrayWiDef[num][5];
//                 (*H_nEFG)(1,0) = Hdata->QtensorArrayWiDef[num][6];       (*H_nEFG)(1,1) = Hdata->QtensorArrayWiDef[num][7];       (*H_nEFG)(1,2) = Hdata->QtensorArrayWiDef[num][8];
//                 (*H_nEFG)(2,0) = Hdata->QtensorArrayWiDef[num][9];       (*H_nEFG)(2,1) = Hdata->QtensorArrayWiDef[num][10];      (*H_nEFG)(2,2) = Hdata->QtensorArrayWiDef[num][11];
//             }else{
//                         printf("error1 in EFG(Opt:%d), the atom is withIn the range of defect cell(curve or bubble),\n",Hdata->Hamil_Quad);
//                         printf("                           but is not the same atoms in the defect cell(curve or bubble)!!\n");
//                 exit(1);
//             }
//         //consider the wrinkle cell Q-tensor along x,y,z-axis (2D & 3D)
//         //about Quadrupole file without Defect position
//         }else{
//             if ((Hdata->cell_para[0] == 0.0) && (Hdata->cell_para[1]==0.0) && (Hdata->cell_para[2]==0.0)){
//                 printf("Something is wrong to Hdata structure!!\n");
//                 printf("Hdata->cell_para : %lf %lf %lf\n",Hdata->cell_para[0],Hdata->cell_para[1],Hdata->cell_para[2]);
//             }
//             if(Hdata->cell_para[0]!=0){difXYZ[0]=ReDefinediff(difXYZ[0],Hdata->Q_MinDifwoDef[0],Hdata->Q_MaxDifwoDef[0],Hdata->cell_para[0]);}
//             if(Hdata->cell_para[1]!=0){difXYZ[1]=ReDefinediff(difXYZ[1],Hdata->Q_MinDifwoDef[1],Hdata->Q_MaxDifwoDef[1],Hdata->cell_para[1]);}
//             if(Hdata->cell_para[2]!=0){difXYZ[2]=ReDefinediff(difXYZ[2],Hdata->Q_MinDifwoDef[2],Hdata->Q_MaxDifwoDef[2],Hdata->cell_para[2]);}

//             //if(CheckRange(difXYZ,Hdata->Q_MinDifwoDef, Hdata->Q_MaxDifwoDef ))
//             if(CheckRange(difXYZ,Hdata,2)){
//                 if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWoDef,&num)){
//                     (*H_nEFG)(0,0) = Hdata->QtensorArrayWoDef[num][3];       (*H_nEFG)(0,1) = Hdata->QtensorArrayWoDef[num][4];      (*H_nEFG)(0,2) = Hdata->QtensorArrayWoDef[num][5];
//                     (*H_nEFG)(1,0) = Hdata->QtensorArrayWoDef[num][6];       (*H_nEFG)(1,1) = Hdata->QtensorArrayWoDef[num][7];      (*H_nEFG)(1,2) = Hdata->QtensorArrayWoDef[num][8];
//                     (*H_nEFG)(2,0) = Hdata->QtensorArrayWoDef[num][9];       (*H_nEFG)(2,1) = Hdata->QtensorArrayWoDef[num][10];     (*H_nEFG)(2,2) = Hdata->QtensorArrayWoDef[num][11];
//                 }else{
//                         printf("error2 in EFG(Opt:%d), the atom is withOut the range of defect cell,\n",Hdata->Hamil_Quad);
//                         printf("                           and is withIn the range of No defect cell(curve or bubble),\n");
//                         printf("                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
//                     exit(1);
//                 }
        
//             //out of range atom
//             //about etc value in quadrupole file
//             }else{
//                 (*H_nEFG)(0,0) = Hdata->Q_Etc[int(spin2[3])][0];     (*H_nEFG)(0,1) = Hdata->Q_Etc[int(spin2[3])][1];     (*H_nEFG)(0,2) = Hdata->Q_Etc[int(spin2[3])][2];
//                 (*H_nEFG)(1,0) = Hdata->Q_Etc[int(spin2[3])][3];     (*H_nEFG)(1,1) = Hdata->Q_Etc[int(spin2[3])][4];     (*H_nEFG)(1,2) = Hdata->Q_Etc[int(spin2[3])][5];
//                 (*H_nEFG)(2,0) = Hdata->Q_Etc[int(spin2[3])][6];     (*H_nEFG)(2,1) = Hdata->Q_Etc[int(spin2[3])][7];     (*H_nEFG)(2,2) = Hdata->Q_Etc[int(spin2[3])][8];
//             }
//         }
//         break;

//     case 4:
//         //consider the special defect cell Q-tensor
//         //about Quadrupole file with Defect 
//         //if(CheckRange(difXYZ,Hdata->Q_MinDifwiDef,Hdata->Q_MaxDifwiDef))
//         if(CheckRange(difXYZ,Hdata,1)){
//             if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWiDef,&num)){
//                 //printf("case 1.1 : w/i defect cell in defect cell range\n");
//                 (*H_nEFG)(0,0) = Hdata->QtensorArrayWiDef[num][3];   (*H_nEFG)(0,1) = Hdata->QtensorArrayWiDef[num][4];  (*H_nEFG)(0,2) = Hdata->QtensorArrayWiDef[num][5];
//                 (*H_nEFG)(1,0) = Hdata->QtensorArrayWiDef[num][6];   (*H_nEFG)(1,1) = Hdata->QtensorArrayWiDef[num][7];  (*H_nEFG)(1,2) = Hdata->QtensorArrayWiDef[num][8];
//                 (*H_nEFG)(2,0) = Hdata->QtensorArrayWiDef[num][9];   (*H_nEFG)(2,1) = Hdata->QtensorArrayWiDef[num][10]; (*H_nEFG)(2,2) = Hdata->QtensorArrayWiDef[num][11];
//             }else{
//                 //consider the defect cell Q-tensor in the special defect cell
//                 //about Quadrupole file without Defect position
//                 //if(CheckRange(difXYZ,Hdata->Q_MinDifwoDef, Hdata->Q_MaxDifwoDef ))
//                 if(CheckRange(difXYZ,Hdata,2)){
//                     if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWoDef,&num)){
//                         //printf("case 1.2 : w/o defect cell in defect cell range\n");
//                         (*H_nEFG)(0,0) = Hdata->QtensorArrayWoDef[num][3];   (*H_nEFG)(0,1) = Hdata->QtensorArrayWoDef[num][4];  (*H_nEFG)(0,2) = Hdata->QtensorArrayWoDef[num][5];
//                         (*H_nEFG)(1,0) = Hdata->QtensorArrayWoDef[num][6];   (*H_nEFG)(1,1) = Hdata->QtensorArrayWoDef[num][7];  (*H_nEFG)(1,2) = Hdata->QtensorArrayWoDef[num][8];
//                         (*H_nEFG)(2,0) = Hdata->QtensorArrayWoDef[num][9];   (*H_nEFG)(2,1) = Hdata->QtensorArrayWoDef[num][10]; (*H_nEFG)(2,2) = Hdata->QtensorArrayWoDef[num][11];
//                     }else{
//                         printf("error1 in EFG(Opt:%d), the atom is withIn the range of special defect cell,\n",Hdata->Hamil_Quad);
//                         printf("                           and is withIn the range of No defect cell(curve or bubble),\n");
//                         printf("                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
//                         exit(1);
//                     }
//                 }else{
//                     if ((Hdata->cell_para[0] == 0.0) && (Hdata->cell_para[1]==0.0) && (Hdata->cell_para[2]==0.0)){
//                         printf("Something is wrong to Hdata structure!!\n");
//                         printf("Hdata->cell_para : %lf %lf %lf\n",Hdata->cell_para[0],Hdata->cell_para[1],Hdata->cell_para[2]);
//                     }
//                     if(Hdata->cell_para[0]!=0){difXYZ[0]=ReDefinediff(difXYZ[0],Hdata->Q_MinDifwoDef[0],Hdata->Q_MaxDifwoDef[0],Hdata->cell_para[0]);}
//                     if(Hdata->cell_para[1]!=0){difXYZ[1]=ReDefinediff(difXYZ[1],Hdata->Q_MinDifwoDef[1],Hdata->Q_MaxDifwoDef[1],Hdata->cell_para[1]);}
//                     if(Hdata->cell_para[2]!=0){difXYZ[2]=ReDefinediff(difXYZ[2],Hdata->Q_MinDifwoDef[2],Hdata->Q_MaxDifwoDef[2],Hdata->cell_para[2]);}

//                     //if(CheckRange(difXYZ,Hdata->Q_MinDifwoDef, Hdata->Q_MaxDifwoDef ))
//                     if(CheckRange(difXYZ,Hdata,2)){
//                         if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWoDef,&num)){
//                             //printf("case 1.3 : w/o defect cell(curve) in defect cell range\n");
//                             (*H_nEFG)(0,0) = Hdata->QtensorArrayWoDef[num][3];   (*H_nEFG)(0,1) = Hdata->QtensorArrayWoDef[num][4];  (*H_nEFG)(0,2) = Hdata->QtensorArrayWoDef[num][5];
//                             (*H_nEFG)(1,0) = Hdata->QtensorArrayWoDef[num][6];   (*H_nEFG)(1,1) = Hdata->QtensorArrayWoDef[num][7];  (*H_nEFG)(1,2) = Hdata->QtensorArrayWoDef[num][8];
//                             (*H_nEFG)(2,0) = Hdata->QtensorArrayWoDef[num][9];   (*H_nEFG)(2,1) = Hdata->QtensorArrayWoDef[num][10]; (*H_nEFG)(2,2) = Hdata->QtensorArrayWoDef[num][11];
//                         }else{
//                             printf("error2 in EFG(Opt:%d), the atom is withIn the range of special defect cell,\n",Hdata->Hamil_Quad);
//                             printf("                           and is withOut the range of No defect cell(curve or bubble),\n");
//                             printf("                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
//                             exit(1);
//                         }
        
//                     //out of range atom
//                     }else{
//                         //printf("case 1.4 : out range atom in defect cell range\n");
//                         (*H_nEFG)(0,0) = Hdata->Q_Etc[int(spin2[3])][0]; (*H_nEFG)(0,1) = Hdata->Q_Etc[int(spin2[3])][1];  (*H_nEFG)(0,2) = Hdata->Q_Etc[int(spin2[3])][2];
//                         (*H_nEFG)(1,0) = Hdata->Q_Etc[int(spin2[3])][3]; (*H_nEFG)(1,1) = Hdata->Q_Etc[int(spin2[3])][4];  (*H_nEFG)(1,2) = Hdata->Q_Etc[int(spin2[3])][5];
//                         (*H_nEFG)(2,0) = Hdata->Q_Etc[int(spin2[3])][6]; (*H_nEFG)(2,1) = Hdata->Q_Etc[int(spin2[3])][7];  (*H_nEFG)(2,2) = Hdata->Q_Etc[int(spin2[3])][8];
//                     }
//                 }
//             }

//         //consider the bubble cell Q-tensor along x,y,z-axis (2D & 3D)
//         }else{
//             //consider the no defect cell Q-tensor in out of the special defect cell
//             //if(CheckRange(difXYZ,Hdata->Q_MinDifwoDef, Hdata->Q_MaxDifwoDef ))
//             if(CheckRange(difXYZ,Hdata,2)){
//                 if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWoDef,&num)){
//                     //printf("case 2.1 : w/o defect cell in No defect cell range\n");
//                     (*H_nEFG)(0,0) = Hdata->QtensorArrayWoDef[num][3];   (*H_nEFG)(0,1) = Hdata->QtensorArrayWoDef[num][4];   (*H_nEFG)(0,2) = Hdata->QtensorArrayWoDef[num][5];
//                     (*H_nEFG)(1,0) = Hdata->QtensorArrayWoDef[num][6];   (*H_nEFG)(1,1) = Hdata->QtensorArrayWoDef[num][7];   (*H_nEFG)(1,2) = Hdata->QtensorArrayWoDef[num][8];
//                     (*H_nEFG)(2,0) = Hdata->QtensorArrayWoDef[num][9];   (*H_nEFG)(2,1) = Hdata->QtensorArrayWoDef[num][10];  (*H_nEFG)(2,2) = Hdata->QtensorArrayWoDef[num][11];
//                 }else{
//                     printf("error3 in EFG(Opt:%d), the atom is withOut the range of special defect cell,\n",Hdata->Hamil_Quad);
//                     printf("                           and is withIn the range of No defect cell(curve or bubble),\n");
//                     printf("                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
//                     exit(1);
//                 }
//             }else{
//                 if ((Hdata->cell_para[0] == 0.0) && (Hdata->cell_para[1]==0.0) && (Hdata->cell_para[2]==0.0)){
//                     printf("Something is wrong to Hdata structure!!\n");
//                     printf("Hdata->cell_para : %lf %lf %lf\n",Hdata->cell_para[0],Hdata->cell_para[1],Hdata->cell_para[2]);
//                 }
//                 if(Hdata->cell_para[0]!=0){difXYZ[0]=ReDefinediff(difXYZ[0],Hdata->Q_MinDifwoDef[0],Hdata->Q_MaxDifwoDef[0],Hdata->cell_para[0]);}
//                 if(Hdata->cell_para[1]!=0){difXYZ[1]=ReDefinediff(difXYZ[1],Hdata->Q_MinDifwoDef[1],Hdata->Q_MaxDifwoDef[1],Hdata->cell_para[1]);}
//                 if(Hdata->cell_para[2]!=0){difXYZ[2]=ReDefinediff(difXYZ[2],Hdata->Q_MinDifwoDef[2],Hdata->Q_MaxDifwoDef[2],Hdata->cell_para[2]);}

//                 //if(CheckRange(difXYZ,Hdata->Q_MinDifwoDef, Hdata->Q_MaxDifwoDef ))
//                 if(CheckRange(difXYZ,Hdata,2)){
//                     if(FIND_AtomPosi(difXYZ, Hdata->QtensorArrayWoDef,&num)){
//                         //printf("case 2.2 : w/o defect cell(curve) in No defect cell range\n");
//                         (*H_nEFG)(0,0) = Hdata->QtensorArrayWoDef[num][3];   (*H_nEFG)(0,1) = Hdata->QtensorArrayWoDef[num][4];  (*H_nEFG)(0,2) = Hdata->QtensorArrayWoDef[num][5];
//                         (*H_nEFG)(1,0) = Hdata->QtensorArrayWoDef[num][6];   (*H_nEFG)(1,1) = Hdata->QtensorArrayWoDef[num][7];  (*H_nEFG)(1,2) = Hdata->QtensorArrayWoDef[num][8];
//                         (*H_nEFG)(2,0) = Hdata->QtensorArrayWoDef[num][9];   (*H_nEFG)(2,1) = Hdata->QtensorArrayWoDef[num][10]; (*H_nEFG)(2,2) = Hdata->QtensorArrayWoDef[num][11];
//                     }else{
//                         printf("error4 in EFG(Opt:%d), the atom is withOut the range of special defect cell,\n",Hdata->Hamil_Quad);
//                         printf("                           and is withOut the range of No defect cell(curve or bubble),\n");
//                         printf("                           but is not the same atoms in the No defect cell(curve or bubble)!!\n");
                        
//                         exit(1);
//                     }
        
//                 //out of range atom
//                 }else{
//                     //printf("case 2.3 : out fo range atom in No defect cell range\n");
//                     (*H_nEFG)(0,0) = Hdata->Q_Etc[int(spin2[3])][0]; (*H_nEFG)(0,1) = Hdata->Q_Etc[int(spin2[3])][1];  (*H_nEFG)(0,2) = Hdata->Q_Etc[int(spin2[3])][2];
//                     (*H_nEFG)(1,0) = Hdata->Q_Etc[int(spin2[3])][3]; (*H_nEFG)(1,1) = Hdata->Q_Etc[int(spin2[3])][4];  (*H_nEFG)(1,2) = Hdata->Q_Etc[int(spin2[3])][5];
//                     (*H_nEFG)(2,0) = Hdata->Q_Etc[int(spin2[3])][6]; (*H_nEFG)(2,1) = Hdata->Q_Etc[int(spin2[3])][7];  (*H_nEFG)(2,2) = Hdata->Q_Etc[int(spin2[3])][8];
//                 }
//             }
//         }
//         break;

//     default:
//         printf("error, Hamil_Quad option is weird\n");
//         exit(1);
//         break;
//     }
// }


///////////////////////////////////////////////////////////////////
//This script is related with the judging the range and atom
//
//In A-&Q-tensor, we should judge the atom is in the range or not.
//And also we match the nuclear spin and tensor data!!
///////////////////////////////////////////////////////////////////


bool READ_BD_vertex(const char* inputfile, double*** vertex, double*** center, double*** normal, char** vertex_condi){

    //char* vertex_condi = {"v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8"};
    
    //Allocate the vertex
    double tmp[3] = {0.0,}; bool Vertex_OK = false;
    (*vertex) = allocDouble2d(8,3);

    //Find the vertex from inputfile
    for (int i = 0; i <8; i++){
        Vertex_OK = READ_Tensor_array(inputfile, &(tmp), vertex_condi[i],3);

        //Fail to find the vertex
        if (Vertex_OK == false){
            for (int j=0; j<8; j++){free((*vertex)[j]);}
            free(*vertex); vertex = NULL;
            return false;
        }

        //copy the data into vertex
        for (int j=0; j<3; j++){
            (*vertex)[i][j] = tmp[j];
        }
        tmp[0] = 0.0; tmp[1] = 0.0; tmp[2] = 0.0;
    }
    
    //Sort the vertex along value
    // vertex is like as,
    //  0 : (-x,-y,-z)    1 : (+x,-y,-z)
    //  2 : (-x,+y,-z)    3 : (+x,+y,-z)
    //  4 : (-x,-y,+z)    5 : (+x,-y,+z)
    //  6 : (-x,+y,+z)    7 : (+x,+y,+z)
    //
    //  QuickSort_2d(&(*vertex),2,0,7);
    //  QuickSort_2d(&(*vertex),1,0,7);
    //  QuickSort_2d(&(*vertex),0,0,7);


    //Create the Plane information
    (*center) = allocDouble2d(8,3);
    (*normal) = allocDouble2d(8,3);
    CreatePlaneInfo(*vertex, center, normal);
    
    return true;
}



////////////////////
//related Hyperfine
////////////////////

//Read the tensor(3x3) value in tensor file
void READ_Tensor(const char* inputfile, double*** Tensor,int numCol){
    FILE* data;
    
    //check file exit
    data = fopen(inputfile, "r");
    if(data==NULL){printf("\tThere is no interaction-file, %s!!\n",inputfile);exit(1);}

    //allocation Tensor
    (*Tensor) = allocDouble2d(1,numCol);
    (*Tensor)[0][0]=1;

    double test[numCol];
    memset(test, 0, numCol*sizeof(test[0]));
    //read and save the tensor value
    const int max=1000;
    char lline[max];
    while(!feof(data)){
        fgets(lline,max,data);

        char *ptr = strtok(lline," ");
        int count=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            double tempNum=atof(ptr);
            if (tempNum != 0){
                test[count]=tempNum;
                count++;
            }else{
                if (isStringDouble(ptr) != 0){test[count]=tempNum;count++;}
            }
            ptr = strtok(NULL, " ");
            if (count==numCol){break;}
        }

        //if the number of value is same, write the value in Tensor array
        if (count==numCol){
            int NowRow=(int)(*Tensor)[0][0];
            reallocDouble2d(Tensor, NowRow, NowRow+1,numCol);
            for (int i=0; i<numCol;i++){
                (*Tensor)[NowRow][i]=test[i];
            }
            (*Tensor)[0][0]=NowRow+1;
            
            ////for test the Tensor
            //for (int i=0;i<numCol;i++){printf("%lf ",(*Tensor)[NowRow-1][i]);}
            //printf("\n");
        }
        //initalize the test array
        memset(test, 0, numCol*sizeof(test[0])); 

    }
    rewind(data);
    fclose(data);

    if ((*Tensor)[0][0]==1){
        printf("\tThe Tensor tensor can not read from tensor file (%s)\n",inputfile);
        printf("\tPlz, Check the Tensor file !!!\n");
    }
}

//Read the etc value in tensor file
void READ_Tensor_etc(const char* inputfile, char** names, int nspecies, double*** Tensor,char* condition, int numCol){
    FILE* data;
    
    //check file exit
    data = fopen(inputfile, "r");
    if(data==NULL){printf("\tThere is no interaction-file, %s!!\n",inputfile);exit(1);}

    //allocate the memory of Array
    (*Tensor) = allocDouble2d(nspecies+1,numCol);
    (*Tensor)[0][0]=nspecies+1;

    //read and save the tensor value
    const int max=1000;
    char lline[max];
    bool CheckIs=true;
    while(!feof(data)){
        fgets(lline,max,data);
        bool CheckCondi=false;
        bool CheckAtom=false;
        char *ptr = strtok(lline," ");
        int count=0;
        int col=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            //find the conditions
            if (strncmp(ptr,condition,strlen(condition))==0){
                CheckCondi=true;
            }
            if (CheckAtom == false){
                for (int i=1; i<nspecies+1;i++){
                    char *Atom = names[i-1];
                    if (strncmp(ptr, Atom, strlen(Atom))==0){
                        CheckAtom=true;
                        count=i;
                        break;
                    }
                }
            }
            // write the value in array
            if (CheckCondi && CheckAtom){
                double tempNum=atof(ptr);
                if (tempNum != 0){
                    (*Tensor)[count][col]=tempNum;col++;;
                }else{
                    if (isStringDouble(ptr) !=0){
                        (*Tensor)[count][col]=tempNum;col++;
                    }
                }
            }
            ptr = strtok(NULL, " ");
            if (numCol==col){CheckIs=false;break;}
        }
    }
    ////print array for test
    //for (int i=1; i<Gy->Gyro_size;i++){
    //    printf("Tensor[%d] : ",i);
    //    for (int j=0; j<numCol; j++){
    //        printf(" %lf ",(*Tensor)[i][j]);
    //    }
    //    printf("\n");
    //}

    rewind(data);
    fclose(data);

    //Check the conditions is in the input file
    if (CheckIs){
        printf("\tThe Tensor Condition (%s) can not read from tensor file (%s)\n",condition,inputfile);
        printf("\tPlz, Check the Tensor file !!!\n");
        exit(0);
    }
}

//Read the constant value in tensor file (g-factor, eQ)
void READ_Tensor_const(const char* inputfile, char** names, int nspecies, double** Array,char* condition){
    FILE* data;
    
    //check file exit
    data = fopen(inputfile, "r");
    if(data==NULL){printf("\tThere is no interaction-file, %s!!\n",inputfile);exit(1);}

    //allocate the memory of Array
    *Array=(double*)calloc((nspecies+1),sizeof(double));
    (*Array)[0] = nspecies+1;
    //read and save the tensor value
    const int max=1000;
    char lline[max];
    bool CheckIs=true;
    while(!feof(data)){
        fgets(lline,max,data);
        bool CheckCondi=false;
        bool CheckAtom=false;
        char *ptr = strtok(lline," ");
        int count=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            //find the conditions
            if (strncmp(ptr,condition,strlen(condition))==0){
                CheckCondi=true;
            }
            if (CheckAtom==false){
                for (int i=1; i<nspecies+1;i++){
                    char *Atom = names[i-1];
                    if (strncmp(ptr, Atom, strlen(Atom))==0){
                        CheckAtom=true;
                        count=i;
                        break;
                    }
                }
            }
            // write the value in array
            if (CheckCondi && CheckAtom){
                double tempNum=atof(ptr);
                CheckIs=false;
                if (tempNum != 0){
                    (*Array)[count]=tempNum;break;
                }else{
                    if (isStringDouble(ptr) !=0){(*Array)[count]=tempNum;break;}
                }
            }
            ptr = strtok(NULL, " ");
        }
    }
    //print array for test
    //for (int i=1; i<Gy->Gyro_size;i++){printf("Array : %lf\n",(*Array)[i]);}

    rewind(data);
    fclose(data);

    //Check the conditions is in the input file
    if (CheckIs){
        printf("\tThe Tensor Condition (%s) can not read from tensor file (%s)\n",condition,inputfile);
        printf("\tPlz, Check the Tensor file !!!\n");
        exit(0);
    }
}

//Read the condition's value in tensor file (MinDif, MaxDif)
bool READ_Tensor_array(const char* inputfile, double(*Array)[3],char* condition,int numCol){
    
    //check file exit
    FILE* data;
    data = fopen(inputfile, "r");
    if(data==NULL){printf("\tThere is no interaction-file, %s!!\n",inputfile);exit(1);}

    //read and save the tensor value
    const int max=1000;
    char lline[max];
    while(!feof(data)){
        fgets(lline,max,data);
        bool CheckOption=false;
        char *ptr = strtok(lline," ");
        int count=0;
        //count the value in tensor file, and tempory save the value in test
        while(ptr !=NULL){
            if (strncmp(ptr,condition,strlen(condition))==0){
                CheckOption=true;
            }
            if (CheckOption){
                double tempNum=atof(ptr);
                if (tempNum != 0){
                    (*Array)[count]=tempNum;
                    count++;
                }else{
                    if (isStringDouble(ptr) !=0){(*Array)[count]=tempNum;count++;}
                }
            }
            ptr = strtok(NULL, " ");
            if (count==numCol){rewind(data);fclose(data);return true;} //success to get data
        }
    }
    rewind(data);
    fclose(data);
    
    //fail to get data
    return false; 

    //printf("\tThe Tensor Condition (%s) can not read from tensor file (%s)\n",condition,inputfile);
    //printf("\tPlz, Check the Tensor file !!!\n");
    //exit(0);
}

//Read the version tag of A-tensor of Q-tensor file (Q.E or VASP)
int READ_Tensor_ver(const char* inputfile,double* SpinFactor, double DefectTotSpin, double* CorrTotSpin){
    
    //check file exit
    FILE* data;
    data = fopen(inputfile, "r");
    if(data==NULL){printf("\tThere is no interaction-file, %s!!\n",inputfile);exit(1);}

    //version tag name
    char QE[] = "QUANTUMESPRESSO";  // will return 0; (with no-find)
    char VASP[] = "VASP";           // will return 1;
    
    //read and save the tensor value
    const int max=1000;     char lline[max];
    bool CheckOption=false;
    while(!feof(data)){
        fgets(lline,max,data);
        char *ptr = strtok(lline," ");
        int count=0;

        while(ptr !=NULL){
            //QE check
            if (strncmp(ptr,QE,strlen(QE))==0){
                CheckOption=true;
                if (rank==0){
                    printf("\t\t  Successfully, The version tag is read from Tensorfile (Q.E)\n");
                }
                break; // Q.E. version
            }
            
            //VASP check
            if (strncmp(ptr,VASP,strlen(VASP))==0){
                CheckOption=true;
                if (rank==0){
                    printf("\t\t  Successfully, The version tag is read from Tensorfile (VASP)\n");
                }
                //calculate the spin factor with total spin & correlation
                if (*CorrTotSpin != 0.0){
                    //consider other central spin value
                    (*SpinFactor) = (*CorrTotSpin) / DefectTotSpin;
                    if (rank==0){
                        printf("\t\t  Are you sure that the Central Spin is different from A-tensor file?\n");
                        printf("\t\t  To check the detail option, set '-v' options\n");
                    }
                }else{
                    //just use total spin written in A-tensor file
                    //(*SpinFactor) = DefectTotSpin; 
                    (*SpinFactor) = 1;
                    if (rank==0){
                        printf("\t\t  You will use the raw data of A-tensor file\n");
                    }
                }
                rewind(data);fclose(data);
                return 1;   // VASP version
            }
            
            ptr = strtok(NULL, " ");
        }
    }
    rewind(data);
    fclose(data);

    //if you can search the version tag, 
    //The 'Version Tag' will be recongize as 'Q.E.' version
    if (CheckOption == false && rank == 0){
        printf("\t\t  The version tage is not finded from Tensorfile, Now we will use 'Q.E.'\n");
    }
    
    //calculate the spin factor with total spin & correlation
    (*CorrTotSpin) = 0.5; //automatically, set S=1/2 in Q.E.
    (*SpinFactor) = (*CorrTotSpin) / DefectTotSpin;
    if (DefectTotSpin != 1 && rank == 0){
        printf("\t\t  DefectTotSpin is different from default value(1)!!\n");
        printf("\t\t  To check the detail option, set '-v' options\n");
    }
    return 0;   // Q.E.

}

//print information from inputfile

//Atensor
void printHfInfo_version(int version, double DefectTotSpin, double CorrTotSpin, double SpinFactor){

    printf("\n");

    //print A-tensor version (Q.E. or VASP)
    if (version == 0){ //Q.E.
        printf("\t\tYour Hyperfine file version is 'Q.E.'\n");
        printf("\t\tCentral Tot spin & Correlation spin = %lf & %lf\n",DefectTotSpin, CorrTotSpin);
        printf("\t\t-->SpinFactor = (Correlation Spin) / (Central Tot spin) = %lf\n",SpinFactor);

    }else if (version == 1){    //VASP
        printf("\t\tYour Hyperfine file version is 'VASP'\n");
        if (CorrTotSpin != 0.0){
            printf("\t\tCentral Tot spin & Correlation spin = %lf & %lf\n",DefectTotSpin, CorrTotSpin);
            printf("\t\t-->SpinFactor = (Correlation Spin) / (Central Tot spin) = %lf\n",SpinFactor);
        }else{
            printf("\t\tNow doesn't consider the SpinFactor(%lf) due to using raw data of A-file!!\n",SpinFactor);
        }

    }else{
        printf("\tThe version option(%d) is weird!!\n",version);
        printf("\tNeed to check the reading version from A-tensor File");
        exit(1);
    }
    printf("\n");
}

void printHfInfo_BD(double** vertex, double** center, double** normal, double MinDif[3], double MaxDif[3], bool Usingvertex){
    if (Usingvertex == true){
        for (int i =0; i<8; i++){
            printf("\tv%d : %+10lf  %+10lf  %+10lf (angstrom)\n",(i+1),vertex[i][0], vertex[i][1], vertex[i][2]);
        }
        printf("\n");

        for (int i =0; i<8; i++){
            printf("\tv%d : %+10lf  %+10lf  %+10lf (angstrom)\n",(i+1),center[i][0], center[i][1], center[i][2]);
        }
        printf("\n");

        for (int i =0; i<8; i++){
            printf("\tv%d : %+10lf  %+10lf  %+10lf (angstrom)\n",(i+1),normal[i][0], normal[i][1], normal[i][2]);
        }
        printf("\n");


    }else{
        printf("\t[Min Range] (angstrom) : %+10lf  %+10lf  %+10lf\n",MinDif[0], MinDif[1],MinDif[2]);
        printf("\t[Max Range] (angstrom) : %+10lf  %+10lf  %+10lf\n",MaxDif[0], MaxDif[1],MaxDif[2]);
        printf("\n");
    }
}

void printHfInfo_etc(double** Etc, double* Gfactor, char** names, int nspecies, int mode){

    char string1[200];
    char string2[200];

    switch(mode){
        case 0:
            strcpy(string1,"\t[g-factor] (no-unit) :\n");
            strcpy(string2,"\t[A etc] (MHz) :\n\tIso  \t\tA-tensor\n");
            break;
        case 1:
            strcpy(string1,"\t[eQ] (Q/milibarn*10^-1) :\n");
            strcpy(string2,"\t[Q etc] (MHz) :\n\tIso   \t\tEFG-tensor\n");
            break;
    }
    
    printf("%s",string1);
    for (int i =1; i<nspecies+1; i++){
        printf("\t%s  =  %12lf\n",names[i-1],Gfactor[i]);
    }
    printf("\n");

    //print etc
    printf("%s",string2);
    for (int i =1; i<nspecies+1; i++){
        printf("\t%s  : ",names[i-1]);
        for (int j=0; j<1; j++){
            printf(" %12lf ",Etc[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printHfInfo_tensor(double** tensorArray, int mode){

    char string1[200];
    
    switch(mode){
        case 0:
            printf("\t[A-tensor] (MHz) :\n");
            printf("\tion  \t\t\t\t\tposition\t\t\t\t\t\tFermi\t\t\t\tA-tensor\n");
            break;
        case 1:
            printf("\t[Q-tensor wi defect] (angstrom) & (Hatree/Bohr_radius^2) :\n");
            printf("\tion  \t\t\t\t\tposition\t\t\t\t\t\t\tEFG-tensor\n");
            break;

        case 2:
            printf("\t[Q-tensor wo defect] (angstrom) & (Hatree/Bohr_radius^2) :\n");
            printf("\tion  \t\t\t\t\tposition\t\t\t\t\t\t\tEFG-tensor\n");
    }   

    for (int i=1; i< int(tensorArray[0][0]); i++){
        //num of ion
        printf("\t%3d  \t",i);
        //ion position
        for (int j =0; j<3; j++){
            printf("%+12lf ",tensorArray[i][j]);
        }
        printf("    ");
        //fermi contact
        printf("%12lf ",tensorArray[i][3]);
        printf("    ");
        //A-tensor
        for (int j=4; j<13; j++){
            printf("%12lf ",tensorArray[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

////
//For checking the atom and tensor's posi
bool CheckPosition(const double Posi[3],double* refPosi){

    bool Is_there=false;
    double err=0.05;

    if(fabs(Posi[0] - refPosi[0]) <= err){
        if(fabs(Posi[1] - refPosi[1]) <= err){
            if(fabs(Posi[2] - refPosi[2]) <= err){
                Is_there=true;   
            }
        }
    }
    return Is_there;
}
bool FIND_AtomPosi(const double refPosi[3], double** TensorValue, int* num){

    //TensorValue ==> " x, y, z, ~~~~~"
    
    int totcol=(int)TensorValue[0][0];
    for(int i=1;i<totcol;i++){
        if(CheckPosition(refPosi,TensorValue[i])){
            *num=i;
            return true;
        }
    }
    return false;
}

//Check the BD using maxRange & minRange
bool CheckBD_Range(const double Posi[3],const double minRange[3],const double maxRange[3],double err){

    bool Is_there=false;

    if(Posi[0] >= minRange[0]-err && Posi[0] <= maxRange[0]+err){
        if(Posi[1] >= minRange[1]-err && Posi[1] <= maxRange[1]+err){
            if(Posi[2] >= minRange[2]-err && Posi[2] <= maxRange[2]+err){
                Is_there=true;
            }
        }
    }
         
    return Is_there;
}

//Check the BD using 8-vertex
bool CheckBD_vertex(const double Posi[3],double** vertex, double** center, double** normal, double err){

    bool Is_there = false;

    //Set the BD to use cutting the atoms
    //
    //Max positive & negative value
    double OuterPosiBD[3] = {0,};    double OuterNegaBD[3] = {0,}; 
    
    //find the max & min value of vertex
    for (int i=0; i<8; i++){
        for (int j=0; j<3; j++){
            if (vertex[i][j] > OuterPosiBD[j]) {OuterPosiBD[j] = vertex[i][j];}
            if (vertex[i][j] < OuterNegaBD[j]) {OuterNegaBD[j] = vertex[i][j];}
        }
    }

    //1. Check that the atom is out of Boundary
    bool tmp_there = false;
    if(Posi[0] >= OuterNegaBD[0]-err && Posi[0] <= OuterPosiBD[0]+err){
        if(Posi[1] >= OuterNegaBD[1]-err && Posi[1] <= OuterPosiBD[1]+err){
            if(Posi[2] >= OuterNegaBD[2]-err && Posi[2] <= OuterPosiBD[2]+err){
                tmp_there=true;
            }
        }
    }
    if (tmp_there == false){return false;} //the atom is out of Outer BD

    //2. Check that the atom is in the unique BD (related with Hexagonal shape, ...)
    Is_there = CheckingUniqeBD(Posi, center, normal);

    //  if (Is_there == false){
    //      printf("\tCaution!, The Atom is between OuterBD and InnerBD!!\n");
    //      printf("\tIf the shape of A- & Q-tensor is orthgonal, Something is werid!!\n");
    //      printf("\tCheck plz!\n\n");
    //  }else{
    //      printf("\tSucess the checking! Posi : %lf %lf %lf\n",Posi[0],Posi[1],Posi[2]);
    //  }

    return Is_there;
}

//
bool CheckingUniqeBD(const double Posi[3], double** center, double** normal){
    //Ref : https://github.com/fgomez03/hexacheck/blob/main/hexacheck.py
    //
    //This function find the atom is in the BD or not
    //To find that, we use the normal vector and center of plane

    //It checks if all dot products are non-positive
    //return 'true' if they are and -1 otherwise.
    double c =0; double vect[3] = {0.0,};
    for (int i=0; i<6; i++){
        for (int j=0; j<3; j++){
            vect[j] = Posi[j] - center[i][j];
        }

        c = vector_dot(vect, normal[i]);
        
        if (c>0){return false;} //point is outside
    }
    return true; //the point is inside
}

//Create the normal vector and center position
void CreatePlaneInfo(double** vertex, double*** center, double*** normal){
    //Ref : https://github.com/fgomez03/hexacheck/blob/main/hexacheck.py
    //
    //The vertex of plane is like as,
    //  0 : (-x,-y,-z)    1 : (+x,-y,-z)
    //  2 : (-x,+y,-z)    3 : (+x,+y,-z)
    //  4 : (-x,-y,+z)    5 : (+x,-y,+z)
    //  6 : (-x,+y,+z)    7 : (+x,+y,+z)
    //    
    //            z
    //            |
    //       .6------7  y
    //     .' |    .'|.'
    //    4---+--5'  |    
    //    |   |  | --|----->x 
    //    |  .2--+---3
    //    |.'    | .'
    //    0------1'

    //Create the center
    //  0 : (0,1,2,3), 1 : (0,2,4,6), 2 : (0,1,4,5)
    //  3 : (1,3,5,7), 4 : (2,3,6,7), 5 : (4,5,6,7)
    CalCenter(&((*center)[0]), vertex[0], vertex[1], vertex[2], vertex[3]);
    CalCenter(&((*center)[1]), vertex[0], vertex[2], vertex[4], vertex[6]);
    CalCenter(&((*center)[2]), vertex[0], vertex[1], vertex[4], vertex[5]);
    CalCenter(&((*center)[3]), vertex[1], vertex[3], vertex[5], vertex[3]);
    CalCenter(&((*center)[4]), vertex[2], vertex[3], vertex[6], vertex[7]);
    CalCenter(&((*center)[5]), vertex[4], vertex[5], vertex[6], vertex[7]);

    //Create the normal
    // 0 : (2-0) x (1-0), 1 : (4-0) x (2-0), 2 : (1-0) x (4-0)
    // 3 : (3-1) x (5-1), 4 : (2-3) x (7-3), 5 : (5-4) x (6-4)
    double u[3] = {0.0,}; double v[3] = {0.0,};

    vector_diff(&(u),vertex[2],vertex[0]); vector_diff(&v,vertex[1],vertex[0]); 
    CalNormal(&((*normal)[0]), u, v);
    vector_diff(&u,vertex[4],vertex[0]); vector_diff(&v,vertex[2],vertex[0]); 
    CalNormal(&((*normal)[1]), u, v);
    vector_diff(&u,vertex[1],vertex[0]); vector_diff(&v,vertex[4],vertex[0]); 
    CalNormal(&((*normal)[2]), u, v);
    vector_diff(&u,vertex[3],vertex[1]); vector_diff(&v,vertex[5],vertex[1]); 
    CalNormal(&((*normal)[3]), u, v);
    vector_diff(&u,vertex[2],vertex[3]); vector_diff(&v,vertex[7],vertex[3]); 
    CalNormal(&((*normal)[4]), u, v);
    vector_diff(&u,vertex[5],vertex[4]); vector_diff(&v,vertex[6],vertex[4]); 
    CalNormal(&((*normal)[5]), u, v);

}

////
//For Q-tensor without defect, re-define the difXYZ
double ReDefinediff(double difXYZ,const  double minRange,const double maxRange,const double Copy_Length){

    if(maxRange-minRange > Copy_Length){
        printf("Some problem in Copy_Length\n");
        printf("maxRange-minRange : %lf\nCopy_Length : %lf\n",maxRange-minRange,Copy_Length);
        exit(1);
    }

    double err=0.05;
    int check_exit=0;

    while(difXYZ > maxRange+err || difXYZ < minRange-err){
        int pre_exit=check_exit;

        if(difXYZ > maxRange+err){
            difXYZ=difXYZ-Copy_Length;
            check_exit++;
            continue;
        }
        if(difXYZ < minRange-err){
            difXYZ=difXYZ+Copy_Length;
            check_exit--;
            continue;
        }
        
        if(pre_exit == check_exit ){
            printf("error, spin standard position is weird\n");
            printf("Need to check the Copy_Length & curve atom position\n");
            exit(1);
        }

    }
    double newdifXYZ=difXYZ;

    return newdifXYZ;
}

////
//calculate the center position using 4-vertex
void CalCenter(double** center, double p1[], double p2[], double p3[], double p4[]){
    for (int i =0; i<3; i++){
        (*center)[i] = (p1[i] + p2[i] + p3[i] + p4[i])/4;
    }
}
//calculate the normal vector using 2-vector
void CalNormal(double** normal, double u[], double v[]){
    vector_cross(&(*normal), u, v, true);
}

////
//Vector cross & dot
void vector_diff(double (*result)[3], double p1[], double p2[]){
    //making vector 
    //result = p1 - p2
    for (int i=0; i<3; i++){
        (*result)[i] = p1[i] - p2[i];
    }
}

double vector_dot(double u[], double v[]){
    double result = 0.0;
    for (int i=0; i<3; i++){
        result += u[i]*v[i];
    }
    return result;
}
void vector_cross(double** cross_P, double u[], double v[], bool norm){
    //vector cross product
    //cross_P = u x v
    (*cross_P)[0] = u[1] * v[2] - u[2] * v[1];
    (*cross_P)[1] = u[2] * v[0] - u[0] * v[2];
    (*cross_P)[2] = u[0] * v[1] - u[1] * v[0];

    //nomalization
    if (norm == true){
        double val=pow((*cross_P)[0],2) + pow((*cross_P)[1],2) + pow((*cross_P)[2],2);
        (*cross_P)[0] /= val;
        (*cross_P)[1] /= val;
        (*cross_P)[2] /= val;
    }
}

