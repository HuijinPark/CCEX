#include "../include/option.h"
#include "../include/json.h"
#include "../include/memory.h"

#include <fstream>
#include <unistd.h>

/**
 * @brief Read the option from the input file
 * @details Read &General tag  options
 * @note
 *      free the 1d values from cJSON should be free if it is not string
*/
void cJSON_readOptionConfig(Config* cnf, char* fccein){

    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            printf("Error before: %s\n", cJSON_GetErrorPtr());
        }
        exit(EXIT_FAILURE);
        freeChar1d(&data);
    }else{
		;//printf("%s",cJSON_Print(root));
	}

    ////////////////////////////////////////////////////////////////////////
    // General Options
    ////////////////////////////////////////////////////////////////////////
    if (rank==0){ 
        printMessage("Read Config Options ...\n");
    }

    if (rank==0){
        printMessage("  - General option-related keys : ");
        printMessage("    [ method, quantity, order, bfield, rbath, rdip, deltat, nstep, rbathcut, rdipcut, nstate, seed ] \n");
    }
    char* method = cJSON_ReadString(root,"method",true,"cce");
    Config_setMethod(cnf,method); // Current possible options : cce, gcce, pcce, dsj, dsjitb, itb, mecce, megcce, pmecce

    char* quantity = cJSON_ReadString(root,"quantity",true,"coherence");
    Config_setQuantity(cnf,quantity); // Current possible options : coherence, noise, dm

    int order = cJSON_ReadInt(root,"order",false,-1);
    Config_setOrder(cnf,order);    

    // Magnetic field : This value can be set by [bx,by,bz] or bz
    float* bfield = cJSON_ReadFloat1d(root,"bfield",true,NULL,3);

    if (bfield != NULL){
        // set bfield as [bx,by,bz]
        Config_setBfield(cnf,bfield);
        freeFloat1d(&bfield);
    }else{
        // set bfield as [0,0,bz]
        float bfieldtry[3] = {0.0,0.0,0.0};
        float bfield_z = cJSON_ReadFloat(root,"bfield",false,-1);
        bfieldtry[2] = bfield_z;
        Config_setBfield(cnf,bfieldtry);
    }
    ////////////////////////////////////////////////////////////////////////
    
    float rbath = cJSON_ReadFloat(root,"rbath",false,-1);
    Config_setRbath(cnf,rbath);

    float rdip = cJSON_ReadFloat(root,"rdip",false,-1);
    Config_setRdip(cnf,rdip);

    float deltat = cJSON_ReadFloat(root,"deltat",false,-1);
    Config_setDeltat(cnf,deltat);

    int nstep = cJSON_ReadInt(root,"nstep",false,-1);
    Config_setNstep(cnf,nstep);

    float rbathcut = cJSON_ReadFloat(root,"rbathcut",true,0.0);
    Config_setRbathcut(cnf,rbathcut);

    float rdipcut = cJSON_ReadFloat(root,"rdipcut",true,0.0);
    Config_setRdipcut(cnf,rdipcut);

    int nstate = cJSON_ReadInt(root,"nstate",true,0);
    Config_setNstate(cnf,nstate);

    int seed = cJSON_ReadInt(root,"seed",true,-1);
    if (seed == -1) {
        seed = time(NULL);
    }
    Config_setSeed(cnf,seed);
    srand(seed);

    ////////////////////////////////////////////////////////////////////////
    // Qubit and BathFiles
    ////////////////////////////////////////////////////////////////////////
    if (rank==0){
        printf("\n");
        printMessage("  - File-related keys :");
        printMessage("    [ qubitfile, gyrofile, bathfile, bathadjust, avaaxfile, statefile, exstatefile ] \n");    
    }

    // qubitfile
    char* qubitfile = cJSON_ReadFilePath(root,"qubitfile",true,NULL);
    if (qubitfile != NULL) {
        Config_allocQubitfile(cnf);
        Config_setQubitfile(cnf,qubitfile);
    }

    // gyrofile
    char* gyrofile = cJSON_ReadFilePath(root,"gyrofile",true,NULL);
    if (gyrofile != NULL) {
        Config_allocGyrofile(cnf);
        Config_setGyrofile(cnf,gyrofile);
    }

    // bathfiles
    int length = 0;
    char** bathfiles = cJSON_ReadFilePath1d(&length, root, "bathfile", true, NULL); // get length
    if (bathfiles != NULL) {
        Config_setNbathfiles(cnf,length);
        Config_allocBathfiles(cnf);
        Config_allocBathadjust(cnf);
        for (int i = 0; i < length; i++){
            Config_setBathfiles_i(cnf,bathfiles[i],i);
        }        
        freeChar2d(&bathfiles,length);
    }
    
    // bathadjust
    double** bathadjustDefault = allocDouble2d(length,3); // all elements are 0.0
    double** bathadjust = cJSON_ReadDouble2d(root,"bathadjust",true,bathadjustDefault,length,3);
    
    for (int i = 0; i < length; i++){
        Config_setBathadjust_i(cnf,bathadjust[i],i);
    }

    freeDouble2d(&bathadjustDefault,length);
    freeDouble2d(&bathadjust,length);

    char* avaaxfile = cJSON_ReadFilePath(root,"avaaxfile",true,NULL);
    if (avaaxfile != NULL) {
        Config_allocAvaaxfile(cnf);
        Config_setAvaaxfile(cnf,avaaxfile);
    }

    char* statefile = cJSON_ReadFilePath(root,"statefile",true,NULL);
    if (statefile != NULL) {
        Config_allocStatefile(cnf);
        Config_setStatefile(cnf,statefile);
    }

    char* exstatefile = cJSON_ReadFilePath(root,"exstatefile",true,NULL);
    if (exstatefile != NULL) {
        Config_allocExstatefile(cnf);
        Config_setExstatefile(cnf,exstatefile);
    }

    ////////////////////////////////////////////////////////////////////////
    // Tensor file central spin option
    ////////////////////////////////////////////////////////////////////////
    if (rank==0){
        printf("\n");
        printMessage("  - Tensorfile-related keys :");
        printMessage("    [ DefectTotSpin, CorrTotSpin, ");
        printMessage("      hf_readmode, hf_tensorfile, hf_cutoff, hf_ignore_oor, ");
        printMessage("      qd_readmode, qd_tensorfile, qd_tensorfile_woqubit, qd_cellpara ] \n");
    }

    double DefectTotSpin = cJSON_ReadDouble(root,"DefectTotSpin",true,1.0);
    Config_setDefectTotSpin(cnf,DefectTotSpin);

    double CorrTotSpin = cJSON_ReadDouble(root,"CorrTotSpin",true,0.0);
    Config_setCorrTotSpin(cnf,CorrTotSpin);
    
    ////////////////////////////////////////////////////////////////////////
    // DFT Hyperfine tensor
    ////////////////////////////////////////////////////////////////////////
    int hf_readmode = cJSON_ReadInt(root,"hf_readmode",true,0);
    Config_setHf_readmode(cnf,hf_readmode);

    if (hf_readmode != 0){
        char* hf_tensorfile = cJSON_ReadFilePath(root,"hf_tensorfile",false,NULL);
        Config_allocHf_tensorfile(cnf);
        Config_setHf_tensorfile(cnf,hf_tensorfile);

        double hf_cutoff = cJSON_ReadDouble(root,"hf_cutoff",true,0.0);
        Config_setHf_cutoff(cnf,hf_cutoff);

        int hf_ignore_oor = cJSON_ReadInt(root,"hf_ignore_oor",true,0.0);
        Config_setHf_ignore_oor(cnf,hf_ignore_oor);
    }    

    ////////////////////////////////////////////////////////////////////////
    // DFT Quadrupole tensor
    ////////////////////////////////////////////////////////////////////////
    int qd_readmode = cJSON_ReadInt(root,"qd_readmode",true,0);
    Config_setQd_readmode(cnf,qd_readmode);

    char* qd_tensorfile = NULL;
    char* qd_tensorfile_woqubit = NULL;
    double* qd_cellpara = NULL;

    switch (qd_readmode){
        case 1:
            if (rank==0){
                printf("Error! qd readmode = 1 is removed option. use 0,2-4\n");
            }
            exit(1);
            break;
        case 2:
            qd_tensorfile = cJSON_ReadFilePath(root,"qd_tensorfile",false,NULL);
            Config_allocQd_tensorfile(cnf);
            Config_setQd_tensorfile(cnf,qd_tensorfile);

            break;
        case 3: case 4:
            qd_tensorfile = cJSON_ReadFilePath(root,"qd_tensorfile",false,NULL);
            Config_allocQd_tensorfile(cnf);
            Config_setQd_tensorfile(cnf,qd_tensorfile);

            qd_tensorfile_woqubit = cJSON_ReadFilePath(root,"qd_tensorfile_woqubit",true,NULL);
            Config_allocQd_tensorfile_woqubit(cnf);
            Config_setQd_tensorfile_woqubit(cnf,qd_tensorfile_woqubit);

            double* qd_cellpara = cJSON_ReadDouble1d(root, "qd_cellpara", "false", NULL,3);
            Config_setQd_cellpara(cnf,qd_cellpara);
            Config_allocQd_tensorfile(cnf);
            freeDouble1d(&qd_cellpara);

            break;

    }

    ////////////////////////////////////////////////////////////////////////
    if (rank==0){
        printf("\n");
        printMessage("  - Additional Hamiltonian-related keys :");
        printMessage("    [ hfmedi, knight ] \n");
    }

    bool hfmediDefault = false;
    bool hfmedi = cJSON_ReadBool(root,"hfmedi",true,hfmediDefault);
    Config_setHfmedi(cnf, hfmedi);

    bool knightDefault = false;
    bool knight = cJSON_ReadBool(root,"knight",true,knightDefault);
    Config_setKnight(cnf, knight);

    cJSON_Delete(root);    
    freeChar1d(&data);
}

/**
 * @brief Read the option from the input file
 * @details Read &Qubit tag  options
*/
void cJSON_readOptionQubitArray(QubitArray* qa, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Qubit Options ...\n");
        printMessage("  - Read values of main-key 'Qubit'");
        printMessage("    sub-key of 'Qubit' : [ nqubit, qubit, intmap, psia, psib, psi0, overhaus, alphaidx, betaidx ] \n");
        printMessage("  - Read values of sub-key 'qubit'");
        printMessage("    sub-sub-key : [ name, spin, gyro, xyz, detuning, alphams, betams ] \n");
        printMessage("  - Read values of sub-key 'intmap'");
        printMessage("    sub-sub-key : [ between, tensor ] \n");
        printMessage("  - Read values of main-key 'qubitfile', 'qzfs', 'qspin', 'qalphams', 'qbetams'");
    }

    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            printf("Error before: %s\n", cJSON_GetErrorPtr());
        }
        exit(EXIT_FAILURE);
        freeChar1d(&data);
    }

    int nqubitDefault = 1;

    char nameDefault[MAX_CHARARRAY_LENGTH] = "q";
    float spinDefault = 1.0;
    double gyroDefault = GAMMA_ELECTRON;
    double detuningDefault = 0.0;
    float alphaMsDefault = 1.0;
    float betaMsDefault = 0.0;
    
    MatrixXcd intmapDefault = MatrixXcd::Zero(3,3);
    MatrixXcd psiDefault;
    int* alphaidxDefault = NULL;
    int* betaidxDefault = NULL;
    bool overhausDefault = false;
    
    // read qubitfile (priority)
    char* _qubitfile = cJSON_ReadFilePath(root,"qubitfile",true,NULL);
    
    if (_qubitfile != NULL) {
        // when the qubitfile exist : doesn't read "Qubit section"
        // It is possible to read only one qubit

        ////////////////////////////////////////////////////////////////////////
        // set nqubit
        ////////////////////////////////////////////////////////////////////////
        // set one qubit
        QubitArray_setNqubit(qa,nqubitDefault);
        QubitArray_allocQubit(qa);
        
        ////////////////////////////////////////////////////////////////////////
        // set qubit properties by nqubit
        ////////////////////////////////////////////////////////////////////////
        // set qubit properties
        snprintf(nameDefault,MAX_CHARARRAY_LENGTH,"q%d",0);
        QubitArray_setQubit_i_name(qa,nameDefault,0);
        QubitArray_setQubit_i_gyro(qa,gyroDefault,0);
        QubitArray_setQubit_i_detuning(qa,detuningDefault,0);
        float spin = cJSON_ReadFloat(root,"qspin",true,spinDefault);
        QubitArray_setQubit_i_spin(qa,spin,0);

        // set the qubit state (alpha, beta)
        float alphams = cJSON_ReadFloat(root,"qalphams",true,alphaMsDefault);
        float betams = cJSON_ReadFloat(root,"qbetams",true,betaMsDefault);
        QubitArray_setQubit_i_alpha_fromMs(qa,alphams,0);
        QubitArray_setQubit_i_beta_fromMs(qa,betams,0);

        ////////////////////////////////////////////////////////////////////////
        //  Intmap
        ////////////////////////////////////////////////////////////////////////
        // set the options of qubit-related Hamiltonian tensor
        //
        // alloc Interaction map
        QubitArray_allocIntmap(qa);

        // read interaction map
        MatrixXcd tensor = cJSON_ReadTensor(root,"qzfs",true,intmapDefault);
        tensor = KHZ_TO_RADKHZ(tensor);

        // Set ZFS in qubit intmap 
        QubitArray_setIntmap_i_j(qa,tensor,0,0); // qubit index:0,0

        // mediatedTerm IO
        cJSON_Delete(root);
        return;
    }else{
        // when the qubitfile doesn't exist : read "Qubit section"
        // It is possible to read multiple qubits
        cJSON* QubitSection = cJSON_GetObjectItem(root,"Qubit");

        int iqubit = 0;
    
        ////////////////////////////////////////////////////////////////////////
        // set nqubit
        ////////////////////////////////////////////////////////////////////////
        int nqubit = cJSON_ReadInt(QubitSection,"nqubit",true,nqubitDefault);
        QubitArray_setNqubit(qa,nqubit);
        QubitArray_allocQubit(qa);
        
        ////////////////////////////////////////////////////////////////////////
        // set qubit properties by nqubit
        ////////////////////////////////////////////////////////////////////////
        cJSON* qubitArray = cJSON_GetObjectItem(QubitSection,"qubit");
        cJSON* qubit;

        cJSON_ArrayForEach(qubit, qubitArray){
            if (qubit == NULL) {
                fprintf(stderr, "Error: The length of \"qubit\" is different from nqubit(%d) \n", nqubit);
                exit(EXIT_FAILURE);
            }
            // read qubit properties from the input file
            snprintf(nameDefault,MAX_CHARARRAY_LENGTH,"q%d",iqubit);
            char* name = cJSON_ReadString(qubit,"name",true,nameDefault);
            float spin = cJSON_ReadFloat(qubit,"spin",true,spinDefault);
            double gyro = cJSON_ReadDouble(qubit,"gyro",true,gyroDefault);
            double* xyz = cJSON_ReadDouble1d(qubit,"xyz",false,NULL,3);
            double detuning = cJSON_ReadDouble(qubit,"detuning",true,detuningDefault); //kHz
            float alphams = cJSON_ReadFloat(qubit,"alphams",true,alphaMsDefault); 
            float betams = cJSON_ReadFloat(qubit,"betams",true,betaMsDefault); 

            // set qubit properties
            QubitArray_setQubit_i_name(qa,name,iqubit);
            QubitArray_setQubit_i_spin(qa,spin,iqubit);
            QubitArray_setQubit_i_xyz(qa,xyz,iqubit);
            QubitArray_setQubit_i_gyro(qa,gyro,iqubit);
            QubitArray_setQubit_i_detuning(qa,KHZ_TO_RADKHZ(detuning),iqubit); // radkHz
            QubitArray_setQubit_i_alpha_fromMs(qa,alphams,iqubit);
            printf("[DEBUG] a2 :\n");
            QubitArray_setQubit_i_beta_fromMs(qa,betams,iqubit);
            printf("[DEBUG] std::cout << qa->qubit[iqubit]->beta << std::endl;\n");
            std::cout << qa->qubit[iqubit]->beta << std::endl;
            std::cout << betams << std::endl;
            iqubit++;
            freeDouble1d(&xyz);
        }
        
        ////////////////////////////////////////////////////////////////////////
        //  Intmap
        ////////////////////////////////////////////////////////////////////////
        // alloc Interaction map
        QubitArray_allocIntmap(qa);
        QubitArray_setIntmap_dipAuto(qa);

        // Set the initial interaction map
        // as dipolar interaction tensor in radkHz
        QubitArray_setIntmap_dipAuto(qa);

        // read interaction map
        cJSON* intmapArray = cJSON_GetObjectItem(QubitSection,"intmap");
        cJSON* intmap;
        cJSON_ArrayForEach(intmap, intmapArray){

            int qubit1_idx; char* qubit1_name;
            int qubit2_idx; char* qubit2_name;
            cJSON* between = cJSON_GetObjectItem(intmap,"between");
            
            // read which qubits are interacting
            if (cJSON_IsArray(between)) {
                qubit1_name = cJSON_GetStringValue(cJSON_GetArrayItem(between,0));
                qubit2_name = cJSON_GetStringValue(cJSON_GetArrayItem(between,1));
                qubit1_idx = QubitArray_getQubitIdx_fromName(qa,qubit1_name);
                qubit2_idx = QubitArray_getQubitIdx_fromName(qa,qubit2_name);
                // check if the qubit name exists
                if (qubit1_idx == -1 || qubit2_idx == -1 ) {
                    fprintf(stderr, "Error: The qubit name is not found in the input file\n");
                    exit(EXIT_FAILURE);
                }
                if (qubit1_idx > qubit2_idx){
                    fprintf(stderr, "Error: First qubit index is greater than second qubit index\n");
                    fprintf(stderr, "       qubit1_name = %s, qubit2_name = %s\n", qubit1_name, qubit2_name);
                    fprintf(stderr, "       qubit1_idx = %d, qubit2_idx = %d\n", qubit1_idx, qubit2_idx);
                    fprintf(stderr, "       qubit1_idx should be less than qubit2_idx\n");
                    exit(EXIT_FAILURE);
                }
            }else{
                fprintf(stderr, "Error: \"between\" is not found in the input file\n");
                exit(EXIT_FAILURE);
            }

            // read tensor properties from the input file
            MatrixXcd tensor = cJSON_ReadTensor(intmap,"tensor",true,intmapDefault);
            tensor = KHZ_TO_RADKHZ(tensor);

            // set Interaction map
            QubitArray_setIntmap_i_j(qa,tensor,qubit1_idx,qubit2_idx);
        }

        ////////////////////////////////////////////////////////////////////////
        //  psia, psib, psi0
        ////////////////////////////////////////////////////////////////////////
        // read psia, psib, psi0
        int dim = QubitArray_dim(qa);
  
        double* psiamat = cJSON_ReadDouble1d(QubitSection,"psia",true,NULL,dim);
        double* psibmat = cJSON_ReadDouble1d(QubitSection,"psib",true,NULL,dim);
        double* psi0mat = cJSON_ReadDouble1d(QubitSection,"psi0",true,NULL,dim);

        // psiDefault is the 0 dimension matrix
        if (psiamat != NULL && psibmat != NULL) {
            QubitArray_setPsia(qa,Double1dToMatrixXcd(psiamat,dim));
            QubitArray_setPsib(qa,Double1dToMatrixXcd(psibmat,dim));
            freeDouble1d(&psiamat);
            freeDouble1d(&psibmat);
        }else{
            QubitArray_setPsia(qa,psiDefault);
            QubitArray_setPsib(qa,psiDefault);
        }

        if (psi0mat != NULL) {
            QubitArray_setPsi0(qa,Double1dToMatrixXcd(psi0mat,dim));
            freeDouble1d(&psi0mat);
        }else{
            QubitArray_setPsi0(qa,psiDefault);
        }

        ////////////////////////////////////////////////////////////////////////
        //  etc options
        //////////////////////////////////////////////////////////////////////// 

        // overhaus
        bool overhaus = cJSON_ReadBool(QubitSection,"overhaus",true,overhausDefault);
        QubitArray_setOverhaus(qa,overhaus);

        // alphaidx, betaidx
        int alphaidx = cJSON_ReadInt(QubitSection,"alphaidx",true,-1);
        int betaidx = cJSON_ReadInt(QubitSection,"betaidx",true,-1);
        if (alphaidx == -1 || betaidx == -1) {
            if (alphaidx != betaidx) {
                fprintf(stderr, "Error: alphaidx and betaidx should be set together\n");
                exit(EXIT_FAILURE);
            }
            QubitArray_alloc_alphaidx_betaidx(qa);
            QubitArray_set_alphaidx(qa,alphaidxDefault);
            QubitArray_set_betaidx(qa,betaidxDefault);
        }else{
            QubitArray_alloc_alphaidx_betaidx(qa);
            QubitArray_set_alphaidx(qa,&alphaidx);
            QubitArray_set_betaidx(qa,&betaidx);
        }
    }    
    cJSON_Delete(root);
    freeChar1d(&data);
}

/**
 * @brief Read the option from the input file
 * @details Read &Cluster tag  options

*/
void cJSON_readOptionCluster(Cluster* clus, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Cluster Options ...");
        printMessage("  [ order, method, addsubclus, nk ] \n");
    }
    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            printf("Error before: %s\n", cJSON_GetErrorPtr());
        }
        exit(EXIT_FAILURE);
        freeChar1d(&data);
    }

    int order = cJSON_ReadInt(root,"order",false, -1); // read twice
    Cluster_setOrder(clus,order);

    char* method = cJSON_ReadString(root,"method",true,"cce"); // read twice
    Cluster_setMethod(clus,method);

    bool addsubclus = cJSON_ReadBool(root,"addsubclus",true,true);
    Cluster_setAddsubclus(clus,addsubclus);

    // nk
    Cluster_allocNk(clus);
    int* nk = cJSON_ReadInt1d(root,"nk",true,NULL,order+1);
    if (nk!=NULL){
        Cluster_setNk(clus,nk);
    }else{;} 
    freeInt1d(&nk);
    // all element would be "0" in which we consider all pairs within rdip


    if (strcasecmp(Cluster_getMethod(clus),"pcce") == 0 || strcasecmp(Cluster_getMethod(clus),"pmecce") == 0 || strcasecmp(Cluster_getMethod(clus),"pgcce") == 0 || strcasecmp(Cluster_getMethod(clus),"pmegcce") == 0){

        if (rank==0){
            printMessage("  [ sK , max_trial, max_iter, kmeans_pp, iter_detail ] \n");
        }

        int* defaultSk        = NULL;
        int default_max_trial = 30000;
        int default_max_iter  = 30000;

        int sK           = cJSON_ReadInt(root,"sK"          , false, *defaultSk ); 
        int max_trial    = cJSON_ReadInt(root,"max_trial"   , true , default_max_trial); 
        int max_iter     = cJSON_ReadInt(root,"max_iter"    , true , default_max_iter); 
        bool kmeans_pp   = cJSON_ReadBool(root,"kmeans_pp"  , true , true);
        bool iter_detail = cJSON_ReadBool(root,"iter_detail", true , false);

        Cluster_setSk(clus,sK);
        Cluster_setMax_trial(clus,max_trial);
        Cluster_setMax_iter(clus,max_iter);
        Cluster_setKmeans_pp(clus,kmeans_pp);
        Cluster_setIter_detail(clus,iter_detail);
    }

    cJSON_Delete(root);
    freeChar1d(&data);
}

////////////////////////// 
// Pulse HS version  //
////////////////////////// 
double parse_fraction(const std::string& frac_str) {
    int num = 0, denom = 1;
    if (sscanf(frac_str.c_str(), "%d/%d", &num, &denom) == 2) {
        return static_cast<double>(num) / denom;
    }
    return std::stod(frac_str);
}

char parse_axis(cJSON* axis_json, int row_index) {
//char parse_axis(cJSON* axis_json) {
    if (!cJSON_IsString(axis_json)) {
        fprintf(stderr, "[Error] Axis at row %d is not a string.\n", row_index);
        exit(EXIT_FAILURE);
    }

    const char* axis_str = axis_json->valuestring;
    if (!axis_str || strlen(axis_str) == 0) {
        fprintf(stderr, "[Error] Axis at row %d is empty.\n", row_index);
        exit(EXIT_FAILURE);
    }

    // Check if the string is one of the allowed axes: X, Y, Z, I
    if (strlen(axis_str) > 1) {
        fprintf(stderr, "[Error] Axis at row %d is too long. Only single characters are allowed.\n", row_index);
        exit(EXIT_FAILURE);
    }

    char c = toupper(axis_str[0]);

    if (c == 'X' || c == 'Y' || c == 'Z' || c == 'I') {
        return c;
    }

    fprintf(stderr, "[Error] Invalid axis '%s'. Only 'X', 'Y', 'Z', 'I' are allowed.\n", axis_str);
    exit(EXIT_FAILURE);
}

std::string read_json_file(const std::string& fccein) {
    std::ifstream file(fccein, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open ccein file: " + fccein);
    file.seekg(0, std::ios::end);
    std::string contents(file.tellg(), '\0');
    file.seekg(0, std::ios::beg);
    file.read(&contents[0], contents.size());
    return contents;
}

void cJSON_readOptionPulse(Pulse* pulse, char* fccein) {

    if (rank == 0) {
        printf("\n");
        printMessage("Read Pulse Options ...");
        printMessage("  [ npulse, pulsename, sequence ] \n");
    }

    char*  data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank == 0) {
            printf("Error before: %s\n", cJSON_GetErrorPtr());
        }
        exit(EXIT_FAILURE);
        freeChar1d(&data);
    }

    bool made = false;

    // Read "pulse number" (Default: -1, i.e. not set)
    int npulse = cJSON_ReadInt(root,"npulse",false, -1);
    Pulse_setNpulse(pulse,npulse);

    // ================= //
    // New Parameters !! //
    // Input pulse duration time (unit : ns)
    double pulse_time     = cJSON_ReadDouble(root, "pulse_time", true, 0.0);
    Pulse_setPulseTime(pulse,pulse_time);
    // Input pulse detuning offset (default : 0)
    double detuning_factor = cJSON_ReadDouble(root, "detuning_factor", true, 1.0);
    Pulse_setPulseDetuningFactor(pulse,detuning_factor);
    // ================= //

    char* pulsename = cJSON_ReadString(root,"pulsename",true,"None");
    Pulse_setPulsename(pulse,pulsename);

    Pulse_allocSequence(pulse);
    Pulse_allocFracs(pulse);
    Pulse_allocAxes(pulse);
    Pulse_allocAngles(pulse);
    Pulse_allocSequenceIndices(pulse);

    cJSON* sequence_array = cJSON_GetObjectItem(root, "sequence");
    if (sequence_array == NULL){
        allocateDefaultSequence(pulse);
        made = true;
    }
    if (strcasecmp(pulse->pulsename, "Manual") == 0) {
        parse_manual_sequence(pulse, sequence_array);
    }

    // New tag~! 2025.05.09 
    // Read Bath pulse sequence !! 
    char* bpulse_defect = cJSON_ReadString(root,"bpulse_defect",true,"NoneNoNe");
    double bpulse_energy_shift = cJSON_ReadDouble(root, "bpulse_energy_shift", true, 0.0);
    double bspin          = cJSON_ReadDouble(root, "bspin", true, 0.5);
    double balphams       = cJSON_ReadDouble(root, "balphams", true,  0.5);
    double bbetams        = cJSON_ReadDouble(root, "bbetams",  true, -0.5);

    int bnpulse = cJSON_ReadInt(root,"bnpulse",true, 0);
    Pulse_setBNpulse(pulse,bnpulse);
    Pulse_setBPulse_Defect(pulse,bpulse_defect);
    Pulse_setBPulse_EnergyShift(pulse, bpulse_energy_shift);
    Pulse_setBPulse_bspin(pulse, bspin);
    Pulse_setBPulse_alphams(pulse, balphams);
    Pulse_setBPulse_betams(pulse, bbetams);
    Pulse_setBPulse_bpsia_fromMs(pulse);
    Pulse_setBPulse_bpsib_fromMs(pulse);

    ////////////////////////////////
    // Generate Bath pulse sequence
    Pulse_allocBSequence(pulse);
    Pulse_allocBFracs(pulse);
    Pulse_allocBAxes(pulse);
    Pulse_allocBAngles(pulse);
    if ((pulse->bnpulse != 0)){
        //Pulse_allocBSequenceIndices(pulse);
        cJSON* bsequence_array = cJSON_GetObjectItem(root, "bsequence");
        parse_manual_bsequence(pulse, bsequence_array);
    }else{
        pulse->bsequence[0][0]  = 0.0;
        pulse->bsequence[0][1]  = 1.0;
        pulse->bsequence[0][2]  = 1.0;
        pulse->bpulse_fracs[0]  = 1.0;
        pulse->bpulse_axes[0]   = 'I';
        pulse->bpulse_angles[0] = 0.0;
    }

    //////////////////////////
    // Generate total sequence 
    int      total_npulse;
    double** total_sequence;     double* total_pulse_fracs;
    double*  total_pulse_angles; char*   total_pulse_axes;
    double* total_bpulse_angles; char*   total_bpulse_axes;

    build_tot_sequence(pulse, &total_sequence, &total_pulse_fracs, 
                        &total_pulse_axes, &total_pulse_angles, 
                        &total_bpulse_axes, &total_bpulse_angles, &total_npulse);

    Pulse_setTotNPulse(pulse, total_npulse);
    Pulse_setTotPulseFracs(pulse, total_pulse_fracs);
    Pulse_setTotPulseAxes(pulse, total_pulse_axes);
    Pulse_setTotPulseAngles(pulse, total_pulse_angles);
    Pulse_setTotBPulseAxes(pulse, total_bpulse_axes);
    Pulse_setTotBPulseAngles(pulse, total_bpulse_angles);
    Pulse_setTotSequence(pulse, total_sequence);

    Pulse_allocTotalSequenceIndices(pulse);
    assign_total_sequence_indices(pulse);

    //////////////////////////
    cJSON_Delete(root);
    freeChar1d(&data);

}

//////////////////////////////////////////////////////
void cJSON_readOptionOutput(Output* op, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Output Options ...");
        printMessage("  [ savemode, outfile ] ");
    }
    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            printf("Error before: %s\n", cJSON_GetErrorPtr());
        }
        exit(EXIT_FAILURE);
        freeChar1d(&data);
    }

    char savemodeDefault[MAX_CHARARRAY_LENGTH] = "normal";
    char* savemode = cJSON_ReadString(root,"savemode",true,savemodeDefault);
    Output_setSavemode(op,savemode);
;
    char* outfileDefault = NULL;
    char* outfile = cJSON_ReadString(root,"outfile",true,outfileDefault);
	if (outfile != NULL){
	    Output_allocOutfile(op);
    	Output_setOutfile(op,outfile);
	}

    cJSON_Delete(root);
    freeChar1d(&data);
}


/**
 * @brief Read the option from the input file
 * @details Read &Defect tag  options
 * @note Unit of the input file :
 * - rxyzs      : [Angstrom]
 * - hypf       : [MHz]
 * - efg        : [Hartree/Bohr^2]
 * - zfs        : [MHz]
 * - gyros      : [radkHz/G]
 * - eqs        : [10^-30 m^2]
 * - detuning   : [MHz]
*/
void cJSON_readOptionDefectArray(DefectArray* dfa, char* fccein){

    if (rank==0){
        printf("\n");
        printMessage("Read Defect Options ...\n");
        printMessage("  - Read values of main-key 'Defect' (Array format)");
        printMessage("    sub-key : [ dfname, naddspin, navaax, apprx,   ( type : single value) ");
        printMessage("                types, spins, gyros, eqs,          ( type : array ) ");
        printMessage("                rxyzs, hypf, efg, zfs, detuning ]  ( type : [ axis, spname, array ] ) \n");
    }

    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        if (rank==0){
            printf("Error before: %s\n", cJSON_GetErrorPtr());
        }
        exit(EXIT_FAILURE);
        freeChar1d(&data);
    }

    int ndefect = 0;

    cJSON* defectArray = cJSON_GetObjectItem(root,"Defect");
    cJSON* defect;

    // Read the number of defect information
    cJSON_ArrayForEach(defect, defectArray){ndefect++;}

    // // Allocate the defect array
    DefectArray_setNdefect(dfa,ndefect);
    DefectArray_allocDefect(dfa);

    // Default values
    bool apprxDefault = true; // do approximation (do not clusterize)
    int naddspinDefault = 0;
    int navaaxDefault = 0;

    // Set the each defect information
    int idf = 0;
    int length = 0;
    cJSON_ArrayForEach(defect,defectArray){

        ////////////////////////////////////////////////////////////////////////
        // General information && Allocation
        ////////////////////////////////////////////////////////////////////////
        int naddspin = cJSON_ReadInt(defect,"naddspin",true,naddspinDefault);
        int navaax = cJSON_ReadInt(defect,"navaax",true,navaaxDefault);
        
        DefectArray_setDefect_idf_naddspin(dfa,idf,naddspin);
        DefectArray_setDefect_idf_navaax(dfa,idf,navaax+1); // 0th : main spin
        DefectArray_allocDefect_idf(dfa,idf,navaax+1,naddspin);

        char* dfname = cJSON_ReadString(defect,"dfname",false,NULL);
        DefectArray_setDefect_idf_dfname(dfa,idf,dfname);

        bool apprx = cJSON_ReadBool(defect,"apprx",true,apprxDefault);
        DefectArray_setDefect_idf_apprx(dfa,idf,apprx);

        ////////////////////////////////////////////////////////////////////////
        // Spin information
        ////////////////////////////////////////////////////////////////////////
        bool HaveDefault = true;
        if (naddspin >= 1){HaveDefault=false;}

        // If there is nuclear spin on defect site, 
        // then the following values should be included in your cce.json:
        // types, spins, gyros (So the "HaveDefault" value would become false)
        char** types = cJSON_ReadString1d(defect,"types",HaveDefault,NULL, naddspin);
        DefectArray_setDefect_idf_types(dfa,idf,types);
        freeChar2d(&types,naddspin);
        
        float* spins = cJSON_ReadFloat1d(defect,"spins",HaveDefault,NULL,naddspin);    
        DefectArray_setDefect_idf_spins(dfa,idf,spins);
        freeFloat1d(&spins);


        double* gyros = cJSON_ReadDouble1d(defect,"gyros",HaveDefault,NULL,naddspin);        
        DefectArray_setDefect_idf_gyros(dfa,idf,gyros);
        freeDouble1d(&gyros);

        double* eqs = cJSON_ReadDouble1d(defect,"eqs",true,NULL,naddspin);
        if (eqs != NULL){
            DefectArray_setDefect_idf_eqs(dfa,idf,eqs);
            freeDouble1d(&eqs);
        }else{
            // set default value
            double* eqsDefault = allocDouble1d(naddspin);
            DefectArray_setDefect_idf_eqs(dfa,idf,eqsDefault);
            freeDouble1d(&eqsDefault);
        }

        ////////////////////////////////////////////////////////////////////////
        // Relative position && tensor information for each axis/additional spins
        ////////////////////////////////////////////////////////////////////////

        // Read rxyzs other wise set 0.0
        cJSON_ReadDefectInfo_IntCharDoubleArray(defect,"rxyzs",3,&(dfa->defect[idf]->rxyzs),dfa->defect[idf]->types,navaax+1,naddspin);

        // Read hypf, otherwise set 0.0
        cJSON_ReadDefectInfo_IntCharMatrixXcd2d(defect,"hypf",9,&(dfa->defect[idf]->hypf),dfa->defect[idf]->types,navaax+1,naddspin);

        // Read quad, otherwise set 0.0
        cJSON_ReadDefectInfo_IntCharMatrixXcd2d(defect,"efg",9,&(dfa->defect[idf]->efg),dfa->defect[idf]->types,navaax+1,naddspin);

        // Read zfs, otherwise set 0.0
        cJSON_ReadDefectInfo_IntCharMatrixXcd1d(defect,"zfs",9,&(dfa->defect[idf]->zfs),navaax+1);

        // Read detuning, otherwise set 0.0
        cJSON_ReadDefectInfo_IntCharDouble(defect,"detuning",&(dfa->defect[idf]->detuning),navaax+1);

        // ////////////////////////////////////////////////////////////////////////
        idf++;
    }
}

////////////////////////////////////////////////////////////////
char* cJSON_ReadFccein(char* fccein){

    FILE* inputFile = fopen(fccein, "r");
    
    if (inputFile == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    
    char line[1024]; 
    char* data = allocChar1d(1); 
    data[0] = '\0'; 
    size_t currentLength = 0; 

    while (fgets(line, sizeof(line), inputFile) != NULL) {
        char* commentStart = strpbrk(line, "!#");
        if (commentStart != NULL) {
            *commentStart = '\0'; 
        }

        size_t lineLength = strlen(line);
        if (lineLength > 0) {
            reallocChar1d(&data, currentLength + lineLength + 1);            
            strcat(data, line);
            currentLength += lineLength;
        }
    }

    fclose(inputFile);
    return data;
}

char* cJSON_ReadFilePath(cJSON* root, char* key, bool _default, char* default_value){
    cJSON* item = cJSON_GetObjectItem(root, key);

    if (cJSON_IsString(item) && (item->valuestring != NULL)) {
        if (access(item->valuestring, R_OK) == 0) {
            return item->valuestring;
        }else{
            // fprintf(stderr, "Warning: %s cannot open/read (%d)\n",key,access(item->valuestring, R_OK));
            // fprintf(stderr, "Current path: %s\n",item->valuestring);
            return item->valuestring;
        }
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

char** cJSON_ReadFilePath1d(int* length, cJSON* root, char* key, bool _default, char** default_value){
    cJSON* item = cJSON_GetObjectItem(root, key);

    if (cJSON_IsArray(item)) {
        char** array = NULL;
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){

            ////access 
            //if (access(itemElement->valuestring, R_OK) != 0) {
            //    fprintf(stderr, "Error: %s is not found in the input file\n",key);
            //    fprintf(stderr, "Current path: %s\n",itemElement->valuestring);
            //    exit(EXIT_FAILURE);
            //}

            if (i==0){
                array = allocChar2d(1,MAX_FILEPATH);
            }else{
                reallocChar2d(&array,i,i+1,MAX_FILEPATH);
            }
            strcpy(array[i],itemElement->valuestring);
            i++;
        }
        *length = i;
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;
            
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

char* cJSON_ReadString(cJSON* root, char* key, bool _default, char* default_value){
    cJSON* item = cJSON_GetObjectItem(root, key);
    if (cJSON_IsString(item) && (item->valuestring != NULL)) {
        return item->valuestring;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

char** cJSON_ReadString1d(cJSON* root, char* key, bool _default, char** default_value, int size){

    cJSON* item = cJSON_GetObjectItem(root, key);

    if (cJSON_IsArray(item)) {
        char** array = allocChar2d(size,MAX_CHARARRAY_LENGTH);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if (i>=size){
                fprintf(stderr, "Error: %s is too long in the input file\n",key);
                exit(EXIT_FAILURE);
            }else{
                strcpy(array[i],itemElement->valuestring);
            }
            i++;
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;
            
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

int cJSON_ReadInt(cJSON* root, char* key, bool _default, int default_value){

    cJSON* item = cJSON_GetObjectItem(root, key);
    if (cJSON_IsNumber(item)) {
        return item->valueint;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

int* cJSON_ReadInt1d(cJSON* root, char* key, bool _default, int* default_value, int size){
    cJSON* item = cJSON_GetObjectItem(root, key);
    
    if (cJSON_IsArray(item)) {
        int* array = allocInt1d(size);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= size){
                fprintf(stderr, "Error: %s is too long in the input file\n",key);
                exit(EXIT_FAILURE);
            }else{
                array[i] = itemElement->valueint;
            }
            i++;
        }
        if (i != size) {
            fprintf(stderr, "Error: %s is too short(%d) in the input file (expected size = %d)\n",key,i,size);
            exit(EXIT_FAILURE);
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;
            
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

double cJSON_ReadDouble(cJSON* root, char* key, bool _default, double default_value){
    cJSON* item = cJSON_GetObjectItem(root, key);
    if (cJSON_IsNumber(item)) {
        return item->valuedouble;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;    
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

double* cJSON_ReadDouble1d(cJSON* root, char* key, bool _default, double* default_value, int size){
    cJSON* item = cJSON_GetObjectItem(root, key);
    
    if (cJSON_IsArray(item)) {
        double* array = allocDouble1d(size);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= size){
                fprintf(stderr, "Error: %s is too long in the input file\n",key);
                exit(EXIT_FAILURE);
            }else{
                array[i] = itemElement->valuedouble;
            }
            i++;
        }
        if (i != size) {
            fprintf(stderr, "Error: %s is too short(%d) in the input file (expected size = %d)\n",key,i,size);
            exit(EXIT_FAILURE);
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;
            
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}


double** cJSON_ReadDouble2d(cJSON* root, char* key, bool _default, double** default_value, int row, int col){

    cJSON* item = cJSON_GetObjectItem(root, key);
    double** array = allocDouble2d(row,col);

    if (cJSON_IsArray(item)) {
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= row){
                fprintf(stderr, "Error: %s is too long in the input file\n",key);
                exit(EXIT_FAILURE);
            }else{
                int j = 0;
                cJSON* itemElement_j;
                cJSON_ArrayForEach(itemElement_j, itemElement){
                    if(j >= col){
                        fprintf(stderr, "Error: %s is too long in the input file\n",key);
                        exit(EXIT_FAILURE);
                    }else{
                        array[i][j] = itemElement_j->valuedouble;
                    }
                    j++;
                }
            }
            i++;
        }
        return array;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            copyDouble2d(array,(const double**)default_value,row,col);
            return array;
            
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

float cJSON_ReadFloat(cJSON* root, char* key, bool _default, float default_value){
    cJSON* item = cJSON_GetObjectItem(root, key);
    if (cJSON_IsNumber(item)) {
        return item->valuedouble;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;    
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

float* cJSON_ReadFloat1d(cJSON* root, char* key, bool _default, float* default_value, int size){
    cJSON* item = cJSON_GetObjectItem(root, key);
    
    if (cJSON_IsArray(item)) {
        float* array = allocFloat1d(size);
        int i = 0;
        cJSON* itemElement;
        cJSON_ArrayForEach(itemElement, item){
            if(i >= size){
                fprintf(stderr, "Error: %s is too long in the input file\n",key);
                exit(EXIT_FAILURE);
            }else{
                array[i] = itemElement->valuedouble;
            }
            i++;
        }
        return array;
    }else{
        if(_default){
            if (rank==0 && strcasecmp(key,"bfield")!=0){printf("           %s : use default value .. \n",key);}
            return default_value;
            
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

bool cJSON_ReadBool(cJSON* root, char* key, bool _default, bool default_value){
    cJSON* item = cJSON_GetObjectItem(root, key);
    if (cJSON_IsBool(item)) {
        return item->valueint;
    }else{
        if(_default){
            if (rank==0){printf("           %s : use default value .. \n",key);}
            return default_value;    
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

MatrixXcd cJSON_ReadTensor(cJSON* root, char* key, bool _default, MatrixXcd default_value){

    int row = 3;
    int col = 3;
    MatrixXcd mat = MatrixXcd::Zero(row,col);

    cJSON* tensor = cJSON_GetObjectItem(root,key);

    if (cJSON_IsArray(tensor)) {

        for (int i = 0; i < cJSON_GetArraySize(tensor); ++i) {

            if (i >= row){
                fprintf(stderr, "Error: %s is too long in the input file (>3)\n",key);
                exit(EXIT_FAILURE);
            }

            cJSON* tensor_i = cJSON_GetArrayItem(tensor, i);
            cJSON_Print(tensor_i);
            if (cJSON_IsArray(tensor_i)) {
                for (int j = 0; j < cJSON_GetArraySize(tensor_i); ++j) {

                    if (j >= col){
                        fprintf(stderr, "Error: %s is too long in the input file (>3)\n",key);
                        exit(EXIT_FAILURE);
                    }

                    cJSON* tensor_i_j = cJSON_GetArrayItem(tensor_i, j);
                    
                    cJSON_Print(tensor_i_j);
                    if (cJSON_IsNumber(tensor_i_j)) {
                        mat(i,j) = tensor_i_j->valuedouble;
                    }else{
                        fprintf(stderr, "Error: %s type error, it should be tensor\n",key);
                        exit(EXIT_FAILURE);
                    }   
                }    
            }else{
                fprintf(stderr, "Error: %s type error, it should be tensor\n",key);
                exit(EXIT_FAILURE);
            }
        }
        return mat;
    }else{
        if(_default){
            return default_value;
        }else{
            fprintf(stderr, "Error: %s is not found in the input file\n",key);
            exit(EXIT_FAILURE);
        }
    }
}

void cJSON_ReadDefectInfo_IntCharDoubleArray(cJSON* root, char* key, int valuecount, double**** array, char** types, int navaax, int naddspin){

    // itemArray2d : (int, char, doubleArray) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);

    if (itemArray2d == NULL && (navaax == 0 && naddspin == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        for (int isp=0; isp<naddspin; isp++){
            for (int j=0; j<valuecount; j++){
                (*array)[iax][isp][j] = 0.0;
            }
        }
    }

    // set the array from the input file
    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, doubleArray
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            fprintf(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            fprintf(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            exit(EXIT_FAILURE);
        }

        // Find axis index
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax > navaax) {
            fprintf(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            fprintf(stderr, "Error: %s[%d] is out of range\n", key, i);
            exit(EXIT_FAILURE);
        }
        
        // Find spin index
        char* spname = cJSON_GetArrayItem(itemArray1d, 1)->valuestring;
        int isp = findIndexChar(types,0,naddspin-1,spname);
        if (isp == -1) {
            fprintf(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            fprintf(stderr, "Error: %s is not found in the input file\n",spname);
            exit(EXIT_FAILURE);
        }
        
        // Find double array values
        cJSON* values = cJSON_GetArrayItem(itemArray1d, 2);

        if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != valuecount) {
            fprintf(stderr, "Error : cJSON_ReadDefectInfo_IntCharDoubleArray");
            fprintf(stderr, "The third element of '%s[%d]' should be an array or the length is not %d\n", key, i,valuecount);
            exit(EXIT_FAILURE);
        }
    
        for (int j = 0; j < valuecount; j++) {
            (*array)[iax][isp][j] = cJSON_GetArrayItem(values, j)->valuedouble;
        }
    }
}

void cJSON_ReadDefectInfo_IntCharMatrixXcd2d(cJSON* root, char* key, int valuecount, MatrixXcd*** array, char** types, int navaax, int naddspin){

    if (valuecount != 9){
        fprintf(stderr, "Error: valuecount should be 9\n");
        exit(EXIT_FAILURE);
    }

    // itemArray2d : (int, char, MatrixXcd) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);

    if (itemArray2d == NULL && (navaax == 0 && naddspin == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        for (int isp=0; isp<naddspin; isp++){
            (*array)[iax][isp] = MatrixXcd::Zero(3,3);
        }
    }



    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, MatrixXcd
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            fprintf(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            exit(EXIT_FAILURE);
        }

        // Find axis index
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax > navaax) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            fprintf(stderr, "Error: %s[%d] is out of range\n", key, i);
            exit(EXIT_FAILURE);
        }

        // Find spin index
        char* spname = cJSON_GetArrayItem(itemArray1d, 1)->valuestring;
        int isp = findIndexChar(types,0,naddspin-1,spname);
        if (isp == -1) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            fprintf(stderr, "Error: %s is not found in the input file\n",spname);
            exit(EXIT_FAILURE);
        }

        // Find double array values
        cJSON* values = cJSON_GetArrayItem(itemArray1d, 2);

        if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != valuecount) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd2d");
            fprintf(stderr, "The third element of '%s[%d]' should be an array or the length is not %d\n", key, i,valuecount);
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < valuecount; j++) {
            int row = j / 3;
            int col = j % 3;
            double value = cJSON_GetArrayItem(values, j)->valuedouble;
            (*array)[iax][isp](row,col) = doublec(value,0.0);
        }
    }
}

void cJSON_ReadDefectInfo_IntCharMatrixXcd1d(cJSON* root, char* key, int valuecount, MatrixXcd** array, int navaax){

    if (valuecount != 9){
        fprintf(stderr, "Error: valuecount should be 9\n");
        exit(EXIT_FAILURE);
    }

    // itemArray2d : (int, char, MatrixXcd) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);

    if (itemArray2d == NULL && (navaax == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        (*array)[iax] = MatrixXcd::Zero(3,3);
    }

    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, MatrixXcd
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            fprintf(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            exit(EXIT_FAILURE);
        }

        // Find axis index
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax > navaax) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            fprintf(stderr, "Error: %s[%d] is out of range\n", key, i);
            exit(EXIT_FAILURE);
        }

        // Find spin index
        char* spname = cJSON_GetArrayItem(itemArray1d, 1)->valuestring;
        if (strcasecmp(spname,"e") != 0) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            fprintf(stderr, "Error: key : %s, %s is not found in the input file\n",key,spname);
            exit(EXIT_FAILURE);
        }

        // Find double array values
        cJSON* values = cJSON_GetArrayItem(itemArray1d, 2);

        if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != valuecount) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharMatrixXcd1d");
            fprintf(stderr, "The third element of '%s[%d]' should be an array or the length is not %d\n", key, i,valuecount);
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < valuecount; j++) {
            int row = j / 3;
            int col = j % 3;
            double value = cJSON_GetArrayItem(values, j)->valuedouble;
            (*array)[iax](row,col) = doublec(value,0.0);
        }
    }
}

void cJSON_ReadDefectInfo_IntCharDouble(cJSON* root, char* key, double** array, int navaax){

    // itemArray2d : (int, char, MatrixXcd) * n
    cJSON* itemArray2d = cJSON_GetObjectItem(root, key);

    if (itemArray2d == NULL && (navaax == 0)){
        return;
    }

    // Initialize the array
    for (int iax=0; iax<navaax; iax++){ 
        (*array)[iax] = 0.0;
    }

    int itemArray2dCount = cJSON_GetArraySize(itemArray2d);

    for (int i = 0; i < itemArray2dCount; i++) {

        // itemArray1d : int, char, MatrixXcd
        cJSON* itemArray1d = cJSON_GetArrayItem(itemArray2d, i);

        if (!cJSON_IsArray(itemArray1d) || cJSON_GetArraySize(itemArray1d) != 3) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharDouble");
            fprintf(stderr, "Each '%s[%d]' entry should be an array of three elements\n", key, i);
            exit(EXIT_FAILURE);
        }

        // Find axis index
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax > navaax) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharDouble");
            fprintf(stderr, "Error: %s[%d] is out of range\n", key, i);
            exit(EXIT_FAILURE);
        }

        // Find spin index
        char* spname = cJSON_GetArrayItem(itemArray1d, 1)->valuestring;
        if (strcasecmp(spname,"e") != 0) {
            fprintf(stderr, "cJSON_ReadDefectInfo_IntCharDouble");
            fprintf(stderr, "Error: key : %s, %s is not found in the input file\n",key,spname);
            exit(EXIT_FAILURE);
        }

        // Find double array values
        (*array)[iax] = cJSON_GetArrayItem(itemArray1d, 2)->valuedouble;
    }
}

void parse_manual_sequence(Pulse* pulse, cJSON* sequence_array) {
    if (!(sequence_array && cJSON_IsArray(sequence_array))) {
        fprintf(stderr, "[Error] sequence_array is missing or not a valid JSON array. [case pulse name = manual]\n");
        exit(EXIT_FAILURE);
    }

    int count = cJSON_GetArraySize(sequence_array);
    if (count != (pulse->npulse) + 1) {
        throw std::runtime_error("\n\n !! Mismatch between sequence length and npulse !! \n");
    }

    double total_frac = 0.0;
    for (int i = 0; i < (pulse->npulse) + 1; ++i) {
        cJSON* item = cJSON_GetArrayItem(sequence_array, i);
        if (!cJSON_IsArray(item)) continue;

        cJSON* frac  = cJSON_GetArrayItem(item, 0);
        cJSON* axis  = cJSON_GetArrayItem(item, 1);
        cJSON* angle = cJSON_GetArrayItem(item, 2);

        if (!frac || !axis || !angle) continue;

        pulse->sequence[i][0] = total_frac;
        double fraction       = parse_fraction(frac->valuestring);
        total_frac           += fraction;
        pulse->sequence[i][1] = total_frac;
        pulse->sequence[i][2] = fraction;

        pulse->pulse_fracs[i]  = fraction;
        std::cout << "pulse->pulse_fracs[" << i << "] = " << pulse->pulse_fracs[i] << std::endl;
        pulse->pulse_axes[i]   = parse_axis(axis, i);
        pulse->pulse_angles[i] = static_cast<double>(angle->valuedouble);
    }

    assign_sequence_indices(pulse);

    if (total_frac - 1.0 > 1e-9) {
        throw std::runtime_error("Total Fraction is not 1");
    }
}

void parse_manual_bsequence(Pulse* pulse, cJSON* bsequence_array) {
    if (!(bsequence_array && cJSON_IsArray(bsequence_array))) {
        fprintf(stderr, "[Error] bsequence_array is missing or not a valid JSON array. [case pulse name = manual]\n");
        exit(EXIT_FAILURE);
    }

    int count = cJSON_GetArraySize(bsequence_array);
    if (count != (pulse->bnpulse) + 1) {
        throw std::runtime_error("\n\n  !! Mismatch between bsequence length and bnpulse ||\n");
    }

    double btotal_frac = 0.0;
    for (int i = 0; i < (pulse->bnpulse) + 1; ++i) {
        cJSON* item = cJSON_GetArrayItem(bsequence_array, i);
        if (!cJSON_IsArray(item)) continue;

        cJSON* bfrac  = cJSON_GetArrayItem(item, 0);
        cJSON* baxis  = cJSON_GetArrayItem(item, 1);
        cJSON* bangle = cJSON_GetArrayItem(item, 2);

        if (!bfrac || !baxis || !bangle) continue;

        pulse->bsequence[i][0] = btotal_frac;
        double bfraction       = parse_fraction(bfrac->valuestring);
        btotal_frac           += bfraction;
        pulse->bsequence[i][1] = btotal_frac;
        pulse->bsequence[i][2] = bfraction;

        pulse->bpulse_fracs[i]  = bfraction;
        pulse->bpulse_axes[i]   = parse_axis(baxis, i);
        pulse->bpulse_angles[i] = static_cast<double>(bangle->valuedouble);
    }

    //assign_bsequence_indices(pulse);

    if (btotal_frac - 1.0 > 1e-9) {
        throw std::runtime_error("Total Fraction is not 1");
    }
}

void build_tot_sequence(Pulse* pulse, double*** total_sequence, double** total_fracs, char** total_pulse_axes, double** total_pulse_angles, char** total_bpulse_axes, double** total_bpulse_angles, int* total_npulse){

    const double*  pulse_fracs  =   pulse->pulse_fracs;
    const double*  bpulse_fracs =  pulse->bpulse_fracs;
    const int            npulse =        pulse->npulse;
    const int           bnpulse =       pulse->bnpulse;
    const char*      pulse_axes =    pulse->pulse_axes;
    const double*  pulse_angles =  pulse->pulse_angles;
    const char*     bpulse_axes =   pulse->bpulse_axes;
    const double* bpulse_angles = pulse->bpulse_angles;


    double pulse_times[512], bpulse_times[512];
    double edges[1024];
    int edge_count = 0;

    // build absolute pulse times
    double acc = 0.0;
    for (int i = 0; i < npulse; ++i){
        acc += pulse_fracs[i];
        pulse_times[i] = acc;
        edges[edge_count++] = acc;
    }

    acc = 0.0;
    for (int i = 0; i < bnpulse; ++i){
        acc += bpulse_fracs[i];
        bpulse_times[i] = acc;
        edges[edge_count++] = acc;
    }

    // sort edges
    for (int i = 0; i < edge_count - 1; ++i){
        for (int j = i + 1; j < edge_count; ++j){
            if (edges[i] > edges[j]){
                double tmp = edges[i]; edges[i] = edges[j]; edges[j] = tmp;
            }
        }
    }

    // remove duplicates
    double unique_edges[1024];
    int unique_count = 0;
    if (edge_count == 0 || edges[0] > 1e-12)
        unique_edges[unique_count++] = 0.0; // 시작점 보장

    for (int i = 0; i < edge_count; ++i){
        if (unique_count == 0 || fabs(edges[i] - unique_edges[unique_count - 1]) > 1e-12) {
            unique_edges[unique_count++] = edges[i];
        }
    }

    if (fabs(unique_edges[unique_count - 1] - 1.0) > 1e-12)
        unique_edges[unique_count++] = 1.0;  // 종료점 보장

    int seg_count = unique_count - 1;

    *total_fracs          = allocDouble1d(seg_count);
    *total_pulse_axes     = allocChar1d(seg_count);
    *total_pulse_angles   = allocDouble1d(seg_count);
    *total_bpulse_axes    = allocChar1d(seg_count);
    *total_bpulse_angles  = allocDouble1d(seg_count);
    *total_sequence       = allocDouble2d(seg_count, 3);

    int valid_seg_count = 0;

    for (int i = 0; i < seg_count; ++i){
        double t1 = unique_edges[i];
        double t2 = unique_edges[i + 1];
        double duration = t2 - t1;

        if (duration < 1e-12) continue;

        (*total_fracs)[valid_seg_count] = duration;

        // Qubit pulse detection at t2
        int qidx = -1;
        double qubit_t = 0.0;
        for (int q = 0; q < npulse; ++q){
            qubit_t += pulse_fracs[q];
            if (fabs(qubit_t - t2) < 1e-12){
                if (!(pulse_axes[q] == 'I' || fabs(pulse_angles[q]) < 1e-6)) {
                    qidx = q;
                }
                break;
            }
        }

        // Bath pulse detection at t2
        int bidx = -1;
        double bath_t = 0.0;
        for (int b = 0; b < bnpulse; ++b){
            bath_t += bpulse_fracs[b];
            if (fabs(bath_t - t2) < 0.01){
                bidx = b;
                break;
            } 
        }

        // assign qubit pulse info
        if (qidx == -1) {
            (*total_pulse_axes)[valid_seg_count] = 'I';
            (*total_pulse_angles)[valid_seg_count] = 0.0;
        } else {
            (*total_pulse_axes)[valid_seg_count] = pulse_axes[qidx];
            (*total_pulse_angles)[valid_seg_count] = pulse_angles[qidx];
        }
        // assign bath pulse info
        if (bidx == -1 || bidx >= bnpulse){
            (*total_bpulse_axes)[valid_seg_count] = 'I';
            (*total_bpulse_angles)[valid_seg_count] = 0.0;
        } else{
            (*total_bpulse_axes)[valid_seg_count] = bpulse_axes[bidx];
            (*total_bpulse_angles)[valid_seg_count] = bpulse_angles[bidx];
        }

        ++valid_seg_count;
    }
    *total_npulse   = valid_seg_count - 1;
    *total_sequence = allocDouble2d(*total_npulse+1, 3);
    double tempfrac = 0.0;
    for (int i=0; i<=*total_npulse; i++){
        (*total_sequence)[i][0] = tempfrac;
        tempfrac += (*total_fracs)[i];
        (*total_sequence)[i][1] = tempfrac;
        (*total_sequence)[i][2] = (*total_fracs)[i];
    }
}


/**
 * @brief Read jump operator options for ME-CCE from the input JSON file
 * @details Reads the "jump_operators" array from the JSON input.
 *          Each element should have: "bath_name", "operator", "rate"
 *
 * Example JSON:
 * {
 *     "jump_operators": [
 *         {"bath_name": "E_up",   "operator": "+-", "rate": 0.001},
 *         {"bath_name": "E_up",   "operator": "z",  "rate": 0.002},
 *         {"bath_name": "E_down", "operator": "+-", "rate": 0.001}
 *     ]
 * }
 *
 * @note Rate unit in input file: MHz (converted to radkHz internally)
 */
void cJSON_readOptionJumpOperators(JumpOperatorArray* joa, char* fccein){

    char* data = cJSON_ReadFccein(fccein);
    cJSON* root = cJSON_Parse(data);

    if (root == NULL){
        freeChar1d(&data);
        return; // No JSON or parse error - just return with empty joa
    }

    cJSON* jumpArray = cJSON_GetObjectItem(root, "jump_operators");
    if (jumpArray == NULL || !cJSON_IsArray(jumpArray)){
        // No jump operators defined - this is fine for non-ME-CCE methods
        cJSON_Delete(root);
        freeChar1d(&data);
        return;
    }

    // Count the number of jump operators
    int njump = 0;
    cJSON* item;
    cJSON_ArrayForEach(item, jumpArray){ njump++; }

    if (njump == 0){
        cJSON_Delete(root);
        freeChar1d(&data);
        return;
    }

    if (rank == 0){
        printf("\n");
        printMessage("Read Jump Operator Options (ME-CCE) ...");
        char msg[200];
        snprintf(msg, 200, "  Found %d jump operator(s)\n", njump);
        printMessage(msg);
    }

    // Allocate
    JumpOperatorArray_alloc(joa, njump);

    // Read each jump operator
    int idx = 0;
    cJSON_ArrayForEach(item, jumpArray){
        char* bath_name = cJSON_ReadString(item, "bath_name", false, NULL);
        char* op_type   = cJSON_ReadString(item, "operator",  false, NULL);
        double rate     = cJSON_ReadDouble(item, "rate",      false, 0.0);

        if (bath_name == NULL || op_type == NULL){
            fprintf(stderr, "Error: jump_operators[%d] must have 'bath_name' and 'operator' fields\n", idx);
            exit(EXIT_FAILURE);
        }

        // Convert rate from MHz to radkHz
        double rate_radkhz = MHZ_TO_RADKHZ(rate);

        JumpOperatorArray_setOp(joa, idx, bath_name, op_type, rate_radkhz);

        if (rank == 0){
            printf("    [%d] bath_name: %s, operator: %s, rate: %.6e MHz (%.6e radkHz)\n",
                   idx, bath_name, op_type, rate, rate_radkhz);
        }

        idx++;
    }

    cJSON_Delete(root);
    freeChar1d(&data);
}
