#include "../include/option.h"
#include "../include/json.h"
#include "../include/memory.h"

#include <cmath>
#include <fstream>
#include <unistd.h>

/**
 * @brief Read a 3-vector axis, which must be exactly three JSON numbers
 * @details Deliberately NOT cJSON_ReadDouble1d. That helper checks the array length but
 *          then takes valuedouble unconditionally, and valuedouble is 0 for a string, a
 *          null or a boolean -- so "qubit_axis": [1,"1",1] would quietly become [1,0,1],
 *          a different but perfectly valid-looking axis, and the whole run would be
 *          silently wrong. An axis is not a place to guess.
 *
 *          File-local on purpose: this is stricter than the shared cJSON_Read* family,
 *          and tightening those would change how every other key is parsed.
 * @param section_name name of the enclosing section, for the error message
 * @param out [in,out] left at its default value if the key is absent and not required
*/
static void readAxis3(cJSON* section, const char* section_name, const char* key,
                     double out[3], bool required){

    cJSON* item = cJSON_GetObjectItem(section,key);
    if (item == NULL){
        if (required){
            fprintf(stderr,"Error: %s.%s is required\n",section_name,key);
            exit(EXIT_FAILURE);
        }
        return; // absent and optional : keep the default
    }

    if (!cJSON_IsArray(item) || cJSON_GetArraySize(item) != 3){
        fprintf(stderr,"Error: %s.%s must be an array of exactly 3 numbers\n",section_name,key);
        exit(EXIT_FAILURE);
    }

    int i = 0;
    cJSON* element;
    cJSON_ArrayForEach(element,item){
        if (!cJSON_IsNumber(element)){
            fprintf(stderr,"Error: %s.%s[%d] is not a number\n",section_name,key,i);
            fprintf(stderr,"Write the axis as plain numbers, e.g. [0, 0, 1] -- not strings, null or booleans.\n");
            exit(EXIT_FAILURE);
        }
        out[i] = element->valuedouble;
        i++;
    }
}

// A JSON boolean and nothing else. cJSON_ReadBool falls back to its default for any
// other type, which would turn "enabled": "true" into a silent no-op.
static bool readStrictBool(cJSON* section, const char* section_name, const char* key, bool default_value){

    cJSON* item = cJSON_GetObjectItem(section,key);
    if (item == NULL){ return default_value; }
    if (!cJSON_IsBool(item)){
        fprintf(stderr,"Error: %s.%s must be a boolean (true or false)\n",section_name,key);
        exit(EXIT_FAILURE);
    }
    return (bool)item->valueint;
}

static DefectCoordinateFrame readDefectCoordinateFrame(cJSON* defect, const char* dfname,
                                                       bool* was_explicit){

    cJSON* item = cJSON_GetObjectItem(defect,"coordinate_frame");
    *was_explicit = (item != NULL);

    // Backward compatibility: legacy Defect tensors were consumed as already being in
    // the computational frame, because there was no Defect frame conversion at all.
    if (item == NULL){ return DEFECT_COORDINATE_FRAME_QUBIT; }

    if (!cJSON_IsString(item) || item->valuestring == NULL){
        fprintf(stderr,"Error: Defect[%s].coordinate_frame must be \"bath\" or \"qubit\"\n",dfname);
        exit(EXIT_FAILURE);
    }
    if (strcasecmp(item->valuestring,"bath") == 0){ return DEFECT_COORDINATE_FRAME_BATH; }
    if (strcasecmp(item->valuestring,"qubit") == 0){ return DEFECT_COORDINATE_FRAME_QUBIT; }

    fprintf(stderr,"Error: Defect[%s].coordinate_frame = \"%s\" is not supported\n",
            dfname,item->valuestring);
    fprintf(stderr,"Use \"bath\" for crystal/bath-file components or \"qubit\" for computational-frame components.\n");
    exit(EXIT_FAILURE);
}

static double readDefectFrequencyScaleToMHz(cJSON* object, const char* dfname,
                                            const char* field_name){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; } // documented default
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        fprintf(stderr,"Error: Defect[%s].%s.unit must be a string\n",dfname,field_name);
        exit(EXIT_FAILURE);
    }

    if (strcasecmp(item->valuestring,"Hz")  == 0){ return 1.0e-6; }
    if (strcasecmp(item->valuestring,"kHz") == 0){ return 1.0e-3; }
    if (strcasecmp(item->valuestring,"MHz") == 0){ return 1.0; }
    if (strcasecmp(item->valuestring,"GHz") == 0){ return 1.0e3; }

    fprintf(stderr,"Error: Defect[%s].%s.unit = \"%s\" is not supported\n",
            dfname,field_name,item->valuestring);
    fprintf(stderr,"Supported frequency units are Hz, kHz, MHz and GHz.\n");
    exit(EXIT_FAILURE);
}

static double readFiniteNumber(cJSON* item, const char* what){
    if (!cJSON_IsNumber(item) || !std::isfinite(item->valuedouble)){
        fprintf(stderr,"Error: %s must be a finite JSON number\n",what);
        exit(EXIT_FAILURE);
    }
    return item->valuedouble;
}

static MatrixXcd readDefectTensor3x3(cJSON* item, const char* what, double scale){

    if (!cJSON_IsArray(item)){
        fprintf(stderr,"Error: %s must be a flat 9-number or nested 3x3 array\n",what);
        exit(EXIT_FAILURE);
    }

    MatrixXcd tensor = MatrixXcd::Zero(3,3);
    int outer_size = cJSON_GetArraySize(item);

    if (outer_size == 9){
        for (int k=0; k<9; k++){
            char element_label[220];
            snprintf(element_label,sizeof(element_label),"%s[%d]",what,k);
            tensor(k/3,k%3) = doublec(scale * readFiniteNumber(cJSON_GetArrayItem(item,k),element_label),0.0);
        }
        return tensor;
    }

    if (outer_size == 3){
        for (int i=0; i<3; i++){
            cJSON* row = cJSON_GetArrayItem(item,i);
            if (!cJSON_IsArray(row) || cJSON_GetArraySize(row) != 3){
                fprintf(stderr,"Error: %s row %d must contain exactly 3 numbers\n",what,i);
                exit(EXIT_FAILURE);
            }
            for (int j=0; j<3; j++){
                char element_label[220];
                snprintf(element_label,sizeof(element_label),"%s[%d][%d]",what,i,j);
                tensor(i,j) = doublec(scale * readFiniteNumber(cJSON_GetArrayItem(row,j),element_label),0.0);
            }
        }
        return tensor;
    }

    fprintf(stderr,"Error: %s must be a flat 9-number or nested 3x3 array\n",what);
    exit(EXIT_FAILURE);
}

static void readDefectZfs(cJSON* defect, DefectArray* dfa, int idf, int nconfig,
                          bool coordinate_frame_explicit){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"zfs");

    if (item == NULL){
        cJSON_ReadDefectInfo_IntCharMatrixXcd1d(defect,(char*)"zfs",9,&df->zfs,nconfig+1);
        DefectArray_setDefect_idf_zfs_input_mode(dfa,idf,DEFECT_ZFS_NONE);
        return;
    }

    if (cJSON_IsArray(item)){
        cJSON_ReadDefectInfo_IntCharMatrixXcd1d(defect,(char*)"zfs",9,&df->zfs,nconfig+1);
        DefectArray_setDefect_idf_zfs_input_mode(dfa,idf,DEFECT_ZFS_INDEXED_LEGACY);
        return;
    }

    if (!cJSON_IsObject(item)){
        fprintf(stderr,"Error: Defect[%s].zfs must be an object or the legacy indexed array\n",df->dfname);
        exit(EXIT_FAILURE);
    }
    if (!coordinate_frame_explicit){
        fprintf(stderr,"Error: Defect[%s] uses an object-form zfs but has no coordinate_frame\n",df->dfname);
        fprintf(stderr,"Set coordinate_frame to \"bath\" or \"qubit\" so the tensor basis is explicit.\n");
        exit(EXIT_FAILURE);
    }
    if (nconfig < 1){
        fprintf(stderr,"Error: Defect[%s] has an object-form zfs but navaax = %d\n",df->dfname,nconfig);
        fprintf(stderr,"At least one indexed configuration is required.\n");
        exit(EXIT_FAILURE);
    }

    cJSON* D_item = cJSON_GetObjectItem(item,"D");
    cJSON* E_item = cJSON_GetObjectItem(item,"E");
    cJSON* tensor_item = cJSON_GetObjectItem(item,"tensor");
    if ((D_item != NULL) == (tensor_item != NULL)){
        fprintf(stderr,"Error: Defect[%s].zfs must contain exactly one of D or tensor\n",df->dfname);
        exit(EXIT_FAILURE);
    }
    if (tensor_item != NULL && E_item != NULL){
        fprintf(stderr,"Error: Defect[%s].zfs.E cannot be combined with zfs.tensor\n",df->dfname);
        exit(EXIT_FAILURE);
    }

    double scale_to_mhz = readDefectFrequencyScaleToMHz(item,df->dfname,"zfs");
    MatrixXcd zfs = MatrixXcd::Zero(3,3);
    DefectZfsInputMode mode;

    if (D_item != NULL){
        if (!df->axis_set){
            fprintf(stderr,"Error: Defect[%s].zfs uses D/E form but Defect[%s].axis is missing\n",
                    df->dfname,df->dfname);
            fprintf(stderr,"The symmetry axis is required to construct the Cartesian ZFS tensor.\n");
            exit(EXIT_FAILURE);
        }

        char what[160];
        snprintf(what,sizeof(what),"Defect[%s].zfs.D",df->dfname);
        double D_mhz = scale_to_mhz * readFiniteNumber(D_item,what);
        double E_mhz = 0.0;
        if (E_item != NULL){
            snprintf(what,sizeof(what),"Defect[%s].zfs.E",df->dfname);
            E_mhz = scale_to_mhz * readFiniteNumber(E_item,what);
        }
        if (fabs(E_mhz) > 1.0e-12){
            fprintf(stderr,"Error: Defect[%s].zfs.E is non-zero (%g MHz)\n",df->dfname,E_mhz);
            fprintf(stderr,"An axis fixes z but not the transverse x/y directions needed by E. Use zfs.tensor instead.\n");
            exit(EXIT_FAILURE);
        }

        // H_ZFS = S . [D (n n^T - I/3)] . S. For n=[1,1,1]/sqrt(3),
        // D=2870 MHz gives the familiar +/-956.6667 MHz off-diagonal tensor.
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                double value = D_mhz * (df->axis[i]*df->axis[j] - ((i == j) ? 1.0/3.0 : 0.0));
                zfs(i,j) = doublec(value,0.0);
            }
        }
        mode = DEFECT_ZFS_SHARED_DE;
    }else{
        char tensor_label[160];
        snprintf(tensor_label,sizeof(tensor_label),"Defect[%s].zfs.tensor",df->dfname);
        zfs = readDefectTensor3x3(tensor_item,tensor_label,scale_to_mhz);
        mode = DEFECT_ZFS_SHARED_TENSOR;
    }

    df->zfs[0] = MatrixXcd::Zero(3,3);
    for (int iax=1; iax<=nconfig; iax++){ df->zfs[iax] = zfs; }
    DefectArray_setDefect_idf_zfs_input_mode(dfa,idf,mode);
}

static void readDefectDetuning(cJSON* defect, DefectArray* dfa, int idf, int nconfig){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"detuning");

    // The absent and legacy-array paths preserve the original MHz behavior exactly.
    if (item == NULL || cJSON_IsArray(item)){
        cJSON_ReadDefectInfo_IntCharDouble(defect,(char*)"detuning",&df->detuning,nconfig+1);
        return;
    }

    if (!cJSON_IsObject(item)){
        fprintf(stderr,"Error: Defect[%s].detuning must be an object or the legacy indexed array\n",
                df->dfname);
        exit(EXIT_FAILURE);
    }

    cJSON* values = cJSON_GetObjectItem(item,"values");
    if (!cJSON_IsArray(values)){
        fprintf(stderr,"Error: Defect[%s].detuning.values must be an indexed array\n",df->dfname);
        fprintf(stderr,"Example: \"detuning\": {\"values\": [[1,\"e\",-2.14]], \"unit\": \"MHz\"}\n");
        exit(EXIT_FAILURE);
    }

    for (int iax=0; iax<=nconfig; iax++){ df->detuning[iax] = 0.0; }
    double scale_to_mhz = readDefectFrequencyScaleToMHz(item,df->dfname,"detuning");

    int nentry = cJSON_GetArraySize(values);
    for (int i=0; i<nentry; i++){
        cJSON* entry = cJSON_GetArrayItem(values,i);
        if (!cJSON_IsArray(entry) || cJSON_GetArraySize(entry) != 3){
            fprintf(stderr,"Error: Defect[%s].detuning.values[%d] must be [configuration, \"e\", value]\n",
                    df->dfname,i);
            exit(EXIT_FAILURE);
        }

        cJSON* iax_item = cJSON_GetArrayItem(entry,0);
        if (!cJSON_IsNumber(iax_item)
            || iax_item->valuedouble != (double)iax_item->valueint
            || iax_item->valueint < 1 || iax_item->valueint > nconfig){
            fprintf(stderr,"Error: Defect[%s].detuning.values[%d] configuration must be an integer in 1..%d\n",
                    df->dfname,i,nconfig);
            exit(EXIT_FAILURE);
        }

        cJSON* spin_item = cJSON_GetArrayItem(entry,1);
        if (!cJSON_IsString(spin_item) || spin_item->valuestring == NULL
            || strcasecmp(spin_item->valuestring,"e") != 0){
            fprintf(stderr,"Error: Defect[%s].detuning.values[%d] spin label must be \"e\"\n",
                    df->dfname,i);
            exit(EXIT_FAILURE);
        }

        char what[180];
        snprintf(what,sizeof(what),"Defect[%s].detuning.values[%d][2]",df->dfname,i);
        df->detuning[iax_item->valueint] =
            scale_to_mhz * readFiniteNumber(cJSON_GetArrayItem(entry,2),what);
    }
}

static double readDefectGyroScaleToRadkHzPerG(cJSON* object, const char* dfname){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        fprintf(stderr,"Error: Defect[%s].gyros.unit must be a string\n",dfname);
        exit(EXIT_FAILURE);
    }

    const char* unit = item->valuestring;
    if (strcasecmp(unit,"radkHz/G") == 0 || strcasecmp(unit,"rad/ms/G") == 0){ return 1.0; }
    if (strcasecmp(unit,"kHz/G") == 0){ return 2.0*M_PI; }
    if (strcasecmp(unit,"MHz/G") == 0){ return 2.0*M_PI*1.0e3; }
    if (strcasecmp(unit,"MHz/T") == 0){ return 2.0*M_PI*1.0e-1; }
    if (strcasecmp(unit,"Hz/T") == 0){ return 2.0*M_PI*1.0e-7; }
    if (strcasecmp(unit,"radMHz/T") == 0 || strcasecmp(unit,"rad/us/T") == 0){ return 1.0e-1; }
    if (strcasecmp(unit,"rad/s/T") == 0){ return 1.0e-7; }

    fprintf(stderr,"Error: Defect[%s].gyros.unit = \"%s\" is not supported\n",dfname,unit);
    fprintf(stderr,"Supported gyro units are radkHz/G, rad/ms/G, kHz/G, MHz/G, MHz/T, Hz/T, radMHz/T and rad/s/T.\n");
    exit(EXIT_FAILURE);
}

static double readDefectEqScaleTo1eMinus30m2(cJSON* object, const char* dfname){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        fprintf(stderr,"Error: Defect[%s].eqs.unit must be a string\n",dfname);
        exit(EXIT_FAILURE);
    }

    const char* unit = item->valuestring;
    if (strcasecmp(unit,"1e-30 m^2") == 0 || strcasecmp(unit,"1e-30m^2") == 0
        || strcasecmp(unit,"10^-30 m^2") == 0 || strcasecmp(unit,"10^-30m^2") == 0
        || strcasecmp(unit,"fm^2") == 0 || strcasecmp(unit,"fm2") == 0
        || strcasecmp(unit,"10 mbarn") == 0){ return 1.0; }
    if (strcasecmp(unit,"m^2") == 0 || strcasecmp(unit,"m2") == 0){ return 1.0e30; }
    if (strcasecmp(unit,"cm^2") == 0 || strcasecmp(unit,"cm2") == 0){ return 1.0e26; }
    if (strcasecmp(unit,"barn") == 0 || strcasecmp(unit,"b") == 0){ return 1.0e2; }
    if (strcasecmp(unit,"mbarn") == 0 || strcasecmp(unit,"millibarn") == 0
        || strcasecmp(unit,"mb") == 0){ return 1.0e-1; }

    fprintf(stderr,"Error: Defect[%s].eqs.unit = \"%s\" is not supported\n",dfname,unit);
    fprintf(stderr,"Supported quadrupole-moment units are 1e-30 m^2, fm^2, m^2, cm^2, barn and mbarn.\n");
    exit(EXIT_FAILURE);
}

static double readDefectEfgScaleToHartreePerBohr2(cJSON* object, const char* dfname){

    cJSON* item = cJSON_GetObjectItem(object,"unit");
    if (item == NULL){ return 1.0; }
    if (!cJSON_IsString(item) || item->valuestring == NULL){
        fprintf(stderr,"Error: Defect[%s].efg.unit must be a string\n",dfname);
        exit(EXIT_FAILURE);
    }

    const char* unit = item->valuestring;
    if (strcasecmp(unit,"Hartree/Bohr^2") == 0 || strcasecmp(unit,"Ha/Bohr^2") == 0
        || strcasecmp(unit,"Ha/a0^2") == 0 || strcasecmp(unit,"a.u.") == 0
        || strcasecmp(unit,"au") == 0 || strcasecmp(unit,"atomic_unit") == 0){ return 1.0; }

    // Keep these constants identical to BathSpin_setQuad_fromEFG. Numerically, one
    // eV of electron potential energy per m^2 corresponds to one V/m^2; the electron
    // charge is implicit in the stored eQ convention used by that routine.
    const double hartree_to_ev = 27.211386;
    const double bohr_to_m = 0.5291772e-10;
    const double volt_per_m2_to_au = bohr_to_m*bohr_to_m/hartree_to_ev;

    if (strcasecmp(unit,"V/m^2") == 0 || strcasecmp(unit,"V/m2") == 0){
        return volt_per_m2_to_au;
    }
    if (strcasecmp(unit,"V/angstrom^2") == 0 || strcasecmp(unit,"V/angstrom2") == 0
        || strcasecmp(unit,"V/A^2") == 0 || strcasecmp(unit,"V/A2") == 0){
        return 1.0e20*volt_per_m2_to_au;
    }

    fprintf(stderr,"Error: Defect[%s].efg.unit = \"%s\" is not supported\n",dfname,unit);
    fprintf(stderr,"Supported EFG units are Hartree/Bohr^2 (a.u.), V/m^2 and V/angstrom^2.\n");
    exit(EXIT_FAILURE);
}

static cJSON* readDefectObjectValues(cJSON* object, const char* dfname,
                                     const char* field_name, int expected_count){

    cJSON* values = cJSON_GetObjectItem(object,"values");
    if (!cJSON_IsArray(values) || cJSON_GetArraySize(values) != expected_count){
        fprintf(stderr,"Error: Defect[%s].%s.values must be an array of exactly %d numbers\n",
                dfname,field_name,expected_count);
        exit(EXIT_FAILURE);
    }
    return values;
}

static void readDefectGyros(cJSON* defect, DefectArray* dfa, int idf,
                            int naddspin, bool have_default){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"gyros");

    if (item == NULL || cJSON_IsArray(item)){
        double* values = cJSON_ReadDouble1d(defect,(char*)"gyros",have_default,NULL,naddspin);
        DefectArray_setDefect_idf_gyros(dfa,idf,values);
        freeDouble1d(&values);
        return;
    }
    if (!cJSON_IsObject(item)){
        fprintf(stderr,"Error: Defect[%s].gyros must be an object or the legacy value array\n",df->dfname);
        exit(EXIT_FAILURE);
    }

    cJSON* values = readDefectObjectValues(item,df->dfname,"gyros",naddspin);
    double scale = readDefectGyroScaleToRadkHzPerG(item,df->dfname);
    for (int isp=0; isp<naddspin; isp++){
        char what[180];
        snprintf(what,sizeof(what),"Defect[%s].gyros.values[%d]",df->dfname,isp);
        df->gyros[isp] = scale * readFiniteNumber(cJSON_GetArrayItem(values,isp),what);
    }
}

static void readDefectEqs(cJSON* defect, DefectArray* dfa, int idf, int naddspin){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,"eqs");

    if (item == NULL || cJSON_IsArray(item)){
        double* values = cJSON_ReadDouble1d(defect,(char*)"eqs",true,NULL,naddspin);
        if (values != NULL){
            DefectArray_setDefect_idf_eqs(dfa,idf,values);
            freeDouble1d(&values);
        }else{
            double* defaults = allocDouble1d(naddspin);
            DefectArray_setDefect_idf_eqs(dfa,idf,defaults);
            freeDouble1d(&defaults);
        }
        return;
    }
    if (!cJSON_IsObject(item)){
        fprintf(stderr,"Error: Defect[%s].eqs must be an object or the legacy value array\n",df->dfname);
        exit(EXIT_FAILURE);
    }

    cJSON* values = readDefectObjectValues(item,df->dfname,"eqs",naddspin);
    double scale = readDefectEqScaleTo1eMinus30m2(item,df->dfname);
    for (int isp=0; isp<naddspin; isp++){
        char what[180];
        snprintf(what,sizeof(what),"Defect[%s].eqs.values[%d]",df->dfname,isp);
        df->eqs[isp] = scale * readFiniteNumber(cJSON_GetArrayItem(values,isp),what);
    }
}

static int findDefectSpinType(Defect* df, const char* spin_name){
    for (int isp=0; isp<df->naddspin; isp++){
        if (strcasecmp(df->types[isp],spin_name) == 0){ return isp; }
    }
    return -1;
}

static void readDefectIndexedTensorWithUnit(cJSON* defect, DefectArray* dfa, int idf,
                                            const char* field_name, MatrixXcd*** target,
                                            int nconfig, int naddspin,
                                            bool coordinate_frame_explicit){

    Defect* df = dfa->defect[idf];
    cJSON* item = cJSON_GetObjectItem(defect,field_name);

    if (item == NULL || cJSON_IsArray(item)){
        cJSON_ReadDefectInfo_IntCharMatrixXcd2d(defect,(char*)field_name,9,target,
                                                df->types,nconfig+1,naddspin);
        return;
    }
    if (!cJSON_IsObject(item)){
        fprintf(stderr,"Error: Defect[%s].%s must be an object or the legacy indexed array\n",
                df->dfname,field_name);
        exit(EXIT_FAILURE);
    }
    if (!coordinate_frame_explicit){
        fprintf(stderr,"Error: Defect[%s] uses object-form %s but has no coordinate_frame\n",
                df->dfname,field_name);
        fprintf(stderr,"Set coordinate_frame to \"bath\" or \"qubit\" so the tensor basis is explicit.\n");
        exit(EXIT_FAILURE);
    }

    cJSON* values = cJSON_GetObjectItem(item,"values");
    if (!cJSON_IsArray(values)){
        fprintf(stderr,"Error: Defect[%s].%s.values must be an indexed tensor array\n",
                df->dfname,field_name);
        exit(EXIT_FAILURE);
    }

    double scale = (strcasecmp(field_name,"hypf") == 0)
                 ? readDefectFrequencyScaleToMHz(item,df->dfname,"hypf")
                 : readDefectEfgScaleToHartreePerBohr2(item,df->dfname);

    for (int iax=0; iax<=nconfig; iax++){
        for (int isp=0; isp<naddspin; isp++){ (*target)[iax][isp] = MatrixXcd::Zero(3,3); }
    }

    int nentry = cJSON_GetArraySize(values);
    for (int i=0; i<nentry; i++){
        cJSON* entry = cJSON_GetArrayItem(values,i);
        if (!cJSON_IsArray(entry) || cJSON_GetArraySize(entry) != 3){
            fprintf(stderr,"Error: Defect[%s].%s.values[%d] must be [configuration, spin, tensor]\n",
                    df->dfname,field_name,i);
            exit(EXIT_FAILURE);
        }

        cJSON* iax_item = cJSON_GetArrayItem(entry,0);
        if (!cJSON_IsNumber(iax_item)
            || iax_item->valuedouble != (double)iax_item->valueint
            || iax_item->valueint < 1 || iax_item->valueint > nconfig){
            fprintf(stderr,"Error: Defect[%s].%s.values[%d] configuration must be an integer in 1..%d\n",
                    df->dfname,field_name,i,nconfig);
            exit(EXIT_FAILURE);
        }

        cJSON* spin_item = cJSON_GetArrayItem(entry,1);
        if (!cJSON_IsString(spin_item) || spin_item->valuestring == NULL){
            fprintf(stderr,"Error: Defect[%s].%s.values[%d] spin label must be a string\n",
                    df->dfname,field_name,i);
            exit(EXIT_FAILURE);
        }
        int isp = findDefectSpinType(df,spin_item->valuestring);
        if (isp < 0){
            fprintf(stderr,"Error: Defect[%s].%s.values[%d] spin type \"%s\" is not in types\n",
                    df->dfname,field_name,i,spin_item->valuestring);
            exit(EXIT_FAILURE);
        }

        char tensor_label[220];
        snprintf(tensor_label,sizeof(tensor_label),"Defect[%s].%s.values[%d][2]",
                 df->dfname,field_name,i);
        (*target)[iax_item->valueint][isp] =
            readDefectTensor3x3(cJSON_GetArrayItem(entry,2),tensor_label,scale);
    }
}

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
        printMessage("    [ method, quantity, propagator, evolution, order, bfield, ");
        printMessage("      rbath, rdip, deltat, nstep, rbathcut, rdipcut, nstate, seed ] \n");
    }
    char* method = cJSON_ReadString(root,"method",true,"cce");
    Config_setMethod(cnf,method); // Current possible options : cce, gcce, pcce, dsj, dsjitb, itb

    char* quantity = cJSON_ReadString(root,"quantity",true,"coherence");
    Config_setQuantity(cnf,quantity); // Current possible options : coherence, noise, dm

    char* propagator = cJSON_ReadString(root,"propagator",true,"eigen");
    Config_setPropagator(cnf,propagator); // gCCE only : eigen (default) | expm (legacy)

    // vector is the DEFAULT, but it only exists where its preconditions hold: gCCE, a pure
    // rho0 (nstate > 0) and propagator=eigen. Do NOT decide that here -- this function
    // runs from inside main.cpp's getopt loop, so -m and -N have not been applied yet.
    // Only record whether the key was written; Config_resolveEvolution settles it once
    // every source is in.
    char* evolution = cJSON_ReadString(root,"evolution",true,"vector");
    Config_setEvolution(cnf,evolution); // gCCE only : vector (default) | matrix (legacy)
    cnf->evolution_isdefault = (cJSON_GetObjectItem(root,"evolution") == NULL);

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
    // magnetic_field_axis : removed
    ////////////////////////////////////////////////////////////////////////
    // It only ever existed in unreleased development work. Silently ignoring it would
    // leave an old input running with a field direction it no longer describes, so the
    // key is refused outright rather than skipped.
    if (cJSON_GetObjectItem(root,"magnetic_field_axis") != NULL){
        fprintf(stderr,"Error: \"magnetic_field_axis\" was removed.\n");
        fprintf(stderr,"Use \"bfield\": [0,0,Bz] with coordinate_frame_rotation.\n");
        exit(EXIT_FAILURE);
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
        // Without a seed the bath states come from a clock-seeded srand, so the same
        // input run twice gives different numbers. Print the value that was used --
        // reproducing the run afterwards needs it.
        if (rank==0){
            fprintf(stderr,"Warning: no \"seed\" in the input file -- seeding from time(NULL) (%d).\n"
                           "         This run is NOT reproducible. Set \"seed\" to make it so.\n", seed);
        }
    }
    Config_setSeed(cnf,seed);
    srand(seed);

    ////////////////////////////////////////////////////////////////////////
    // Qubit and BathFiles
    ////////////////////////////////////////////////////////////////////////
    if (rank==0){
        printf("\n");
        printMessage("  - File-related keys :");
        printMessage("    [ qubitfile, gyrofile, bathfile, bathadjust, coordinate_frame_rotation, ");
        printMessage("      defect_axis_reference, avaaxfile, statefile, exstatefile ] \n");
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

    ////////////////////////////////////////////////////////////////////////
    // Coordinate frame rotation
    ////////////////////////////////////////////////////////////////////////
    // Absent section, or "enabled": false, leaves Config at the values set by
    // Config_init -- rotation off, both axes on +z -- so the run is bit-for-bit
    // the legacy one.
    //
    // "bath_coordinate_rotation" was the name while the feature only moved bath
    // POSITIONS. It now also decides the magnetic-field direction, so the name is
    // wrong; it is still accepted, with a warning, so that inputs written against
    // the earlier version keep working.
    const char* ROTKEY = "coordinate_frame_rotation";
    cJSON* RotSection    = cJSON_GetObjectItem(root,ROTKEY);
    cJSON* RotSectionOld = cJSON_GetObjectItem(root,"bath_coordinate_rotation");

    if (RotSection != NULL && RotSectionOld != NULL){
        fprintf(stderr,"Error: both \"coordinate_frame_rotation\" and \"bath_coordinate_rotation\" are present\n");
        fprintf(stderr,"They are the same option under two names. Keep \"coordinate_frame_rotation\".\n");
        exit(EXIT_FAILURE);
    }

    if (RotSection == NULL && RotSectionOld != NULL){
        RotSection = RotSectionOld;
        ROTKEY = "bath_coordinate_rotation";
        if (rank==0){
            fprintf(stderr,"Warning: \"bath_coordinate_rotation\" is deprecated -- rename it to\n"
                           "         \"coordinate_frame_rotation\". The behaviour is unchanged.\n");
        }
    }

    if (RotSection != NULL){

        // A non-object -- "coordinate_frame_rotation": true, or a path string -- has no
        // children, so every sub-key below would read as missing and the run would go on
        // silently WITHOUT the rotation the input clearly asked for.
        if (!cJSON_IsObject(RotSection)){
            fprintf(stderr,"Error: %s must be an object, e.g.\n",ROTKEY);
            fprintf(stderr,"  \"%s\" : { \"enabled\" : true, \"bath_axis\" : [0,0,1], \"qubit_axis\" : [1,1,1] }\n",ROTKEY);
            exit(EXIT_FAILURE);
        }

        bool rot_enabled = readStrictBool(RotSection,ROTKEY,"enabled",false);
        Config_setRot_enabled(cnf,rot_enabled);

        if (rot_enabled){
            double bath_axis[3]  = {0.0,0.0,1.0};
            double qubit_axis[3] = {0.0,0.0,1.0};

            readAxis3(RotSection,ROTKEY,"bath_axis" ,bath_axis ,false);
            readAxis3(RotSection,ROTKEY,"qubit_axis",qubit_axis,false);

            Config_setRot_bath_axis(cnf,bath_axis);
            Config_setRot_qubit_axis(cnf,qubit_axis);

            // Which qubit the rotation is taken about. Optional while nqubit is
            // restricted to 1, where there is nothing to disambiguate; when given it
            // must name an existing qubit exactly (resolved in readBathfiles, where
            // the QubitArray exists).
            // Which basis the external tensor COMPONENTS are written in. No silent
            // default: reading a bath-frame tensor as if it were already rotated (or
            // the reverse) is a double / missing rotation that nothing downstream can
            // detect, so it has to be stated. Required only where it has an effect --
            // enforced below, once hf_readmode and qd_readmode have been read too.
            cJSON* hfFrame = cJSON_GetObjectItem(RotSection,"hf_tensor_frame");
            if (hfFrame != NULL){
                if (!cJSON_IsString(hfFrame) || hfFrame->valuestring == NULL){
                    fprintf(stderr,"Error: %s.hf_tensor_frame must be \"bath\" or \"qubit\"\n",ROTKEY);
                    exit(EXIT_FAILURE);
                }
                Config_setRot_hf_tensor_frame(cnf,hfFrame->valuestring);
            }

            cJSON* qdFrame = cJSON_GetObjectItem(RotSection,"qd_tensor_frame");
            if (qdFrame != NULL){
                if (!cJSON_IsString(qdFrame) || qdFrame->valuestring == NULL){
                    fprintf(stderr,"Error: %s.qd_tensor_frame must be \"bath\" or \"qubit\"\n",ROTKEY);
                    exit(EXIT_FAILURE);
                }
                Config_setRot_qd_tensor_frame(cnf,qdFrame->valuestring);
            }

            // Which frame every Qubit.xyz is written in. One value for the whole
            // QubitArray -- per-qubit frames are not supported, and would not make
            // sense: a single R and r0 move every qubit and every bath spin together.
            cJSON* posFrame = cJSON_GetObjectItem(RotSection,"qubit_position_frame");
            if (posFrame != NULL){
                if (!cJSON_IsString(posFrame) || posFrame->valuestring == NULL){
                    fprintf(stderr,"Error: %s.qubit_position_frame must be the string \"bath\"\n",ROTKEY);
                    exit(EXIT_FAILURE);
                }
                Config_setRot_qubit_position_frame(cnf,posFrame->valuestring);
            }

            cJSON* refItem = cJSON_GetObjectItem(RotSection,"reference_qubit");
            if (refItem != NULL){
                if (!cJSON_IsString(refItem) || refItem->valuestring == NULL){
                    fprintf(stderr,"Error: %s.reference_qubit must be a qubit name (a string)\n",ROTKEY);
                    exit(EXIT_FAILURE);
                }
                Config_setRot_reference_qubit(cnf,refItem->valuestring);
            }
        }
    }

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

            // Both keys are REQUIRED here, and were meant to be. Passing them through as
            // optional handed Config_setQd_tensorfile_woqubit and Config_setQd_cellpara a
            // NULL, which reach strcpy and memcpy -- so an input that omitted either one
            // segfaulted instead of being told what was missing. "false" as the _default
            // argument was a string literal, which is a non-null pointer and therefore
            // true, so qd_cellpara asked for the opposite of what it meant.
            qd_tensorfile_woqubit = cJSON_ReadFilePath(root,"qd_tensorfile_woqubit",false,NULL);
            Config_allocQd_tensorfile_woqubit(cnf);
            Config_setQd_tensorfile_woqubit(cnf,qd_tensorfile_woqubit);

            double* qd_cellpara = cJSON_ReadDouble1d(root, "qd_cellpara", false, NULL,3);
            Config_setQd_cellpara(cnf,qd_cellpara);
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

    ////////////////////////////////////////////////////////////////////////
    // Defect axis reference
    ////////////////////////////////////////////////////////////////////////
    // Read here, after the rotation section, because the only value it can take names
    // an axis that only that section defines.
    cJSON* DefAxisItem = cJSON_GetObjectItem(root,"defect_axis_reference");

    if (DefAxisItem != NULL){

        if (!cJSON_IsString(DefAxisItem) || DefAxisItem->valuestring == NULL){
            fprintf(stderr,"Error: defect_axis_reference must be a string\n");
            fprintf(stderr,"The only accepted value is \"qubit_axis\".\n");
            exit(EXIT_FAILURE);
        }
        if (strcasecmp(DefAxisItem->valuestring,"qubit_axis") != 0){
            fprintf(stderr,"Error: defect_axis_reference = \"%s\" is not available\n",
                    DefAxisItem->valuestring);
            fprintf(stderr,"The only accepted value is \"qubit_axis\".\n");
            exit(EXIT_FAILURE);
        }
        if (!Config_getRot_enabled(cnf)){
            fprintf(stderr,"Error: defect_axis_reference = \"qubit_axis\" needs a qubit axis,\n");
            fprintf(stderr,"and coordinate_frame_rotation is what defines one in the source frame.\n");
            fprintf(stderr,"Enable coordinate_frame_rotation, or drop defect_axis_reference.\n");
            exit(EXIT_FAILURE);
        }
        Config_setDefect_axis_reference(cnf,DEFECT_AXIS_REFERENCE_QUBIT_AXIS);
    }

    ////////////////////////////////////////////////////////////////////////
    // Cross-checks that need both the rotation section and the readmodes
    ////////////////////////////////////////////////////////////////////////
    if (Config_getRot_enabled(cnf)){

        // readmode 1 discards the file's anisotropic tensor and rebuilds it as a
        // point-dipole tensor from the source geometry, keeping only the isotropic
        // Fermi-contact term, which has no frame. So hf_tensor_frame decides nothing
        // there -- not required, and said out loud if it is given anyway.
        if (Config_getHf_readmode(cnf) == 1 && Config_getRot_hf_tensor_frame(cnf)[0] != '\0' && rank==0){
            printMessage("  Note: hf_tensor_frame is not applicable at hf_readmode = 1.");
            printMessage("        Only the isotropic Fermi-contact term is taken from the file, and the");
            printMessage("        point-dipole tensor is transformed from the bath frame regardless.\n");
        }

        if ((Config_getHf_readmode(cnf) == 2 || Config_getHf_readmode(cnf) == 3)
            && Config_getRot_hf_tensor_frame(cnf)[0] == '\0'){
            fprintf(stderr,"Error: hf_readmode = %d with the coordinate frame rotation on, but\n",
                    Config_getHf_readmode(cnf));
            fprintf(stderr,"coordinate_frame_rotation.hf_tensor_frame is not set.\n");
            fprintf(stderr,"State the basis the tensor COMPONENTS in \"hf_tensorfile\" are written in:\n");
            fprintf(stderr,"  \"bath\"  : the original bath Cartesian basis -- they get transformed\n");
            fprintf(stderr,"  \"qubit\" : already the qubit-aligned basis -- they are left alone\n");
            fprintf(stderr,"(The POSITION columns of that file are always read in the original frame.)\n");
            exit(EXIT_FAILURE);
        }
        if (Config_getQd_readmode(cnf) != 0 && Config_getRot_qd_tensor_frame(cnf)[0] == '\0'){
            fprintf(stderr,"Error: qd_readmode = %d with the coordinate frame rotation on, but\n",
                    Config_getQd_readmode(cnf));
            fprintf(stderr,"coordinate_frame_rotation.qd_tensor_frame is not set.\n");
            fprintf(stderr,"State the basis the tensor COMPONENTS in \"qd_tensorfile\" are written in:\n");
            fprintf(stderr,"  \"bath\"  : the original bath Cartesian basis -- they get transformed\n");
            fprintf(stderr,"  \"qubit\" : already the qubit-aligned basis -- they are left alone\n");
            exit(EXIT_FAILURE);
        }
    }

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
        printMessage("    sub-sub-key : [ between, tensor, tensor_frame ] \n");
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
        // "qzfs" has no tensor_frame sub-key, and this path is single-qubit only, so it
        // stays where the legacy code left it: read as written, never rotated.
        QubitArray_setIntmap_i_j_frame(qa,INTMAP_EXPLICIT_UNSPECIFIED,0,0);

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
            QubitArray_setQubit_i_beta_fromMs(qa,betams,iqubit);
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

            // Which basis the components are written in. Recorded, never guessed: a ZFS
            // tensor and a pair tensor look the same, and so does one that has already
            // been rotated. Whether the omission is fatal depends on nqubit and on
            // whether the rotation is on, which validateCoordinateFrameRotationInputs
            // decides once the whole input is in.
            IntmapProvenance frame = INTMAP_EXPLICIT_UNSPECIFIED;
            cJSON* tframe = cJSON_GetObjectItem(intmap,"tensor_frame");
            if (tframe != NULL){
                if (!cJSON_IsString(tframe) || tframe->valuestring == NULL){
                    fprintf(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].tensor_frame must be \"bath\" or \"qubit\"\n",
                            qubit1_name,qubit2_name);
                    exit(EXIT_FAILURE);
                }
                if (strcasecmp(tframe->valuestring,"bath") == 0){
                    frame = INTMAP_EXPLICIT_BATH;
                }else if (strcasecmp(tframe->valuestring,"qubit") == 0){
                    frame = INTMAP_EXPLICIT_QUBIT;
                }else{
                    fprintf(stderr,"Error: Qubit.intmap[\"%s\",\"%s\"].tensor_frame = \"%s\" is not available\n",
                            qubit1_name,qubit2_name,tframe->valuestring);
                    fprintf(stderr,"  \"bath\"  : components are in the original bath Cartesian basis\n");
                    fprintf(stderr,"  \"qubit\" : components are already in the common computational basis\n");
                    fprintf(stderr,"            (reference-qubit-aligned; NOT an individual local frame)\n");
                    exit(EXIT_FAILURE);
                }
            }

            // set Interaction map
            QubitArray_setIntmap_i_j(qa,tensor,qubit1_idx,qubit2_idx);
            QubitArray_setIntmap_i_j_frame(qa,frame,qubit1_idx,qubit2_idx);
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


    if (strcasecmp(Cluster_getMethod(clus),"pcce") == 0){

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

    bool pulseiter = cJSON_ReadBool(root,"pulseiter",true,false);
    Pulse_setPulseiter(pulse, pulseiter);

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
    Pulse_allocAxes(pulse);
    Pulse_allocAngles(pulse);
    Pulse_allocSequenceIndices(pulse);

    cJSON* sequence_array = cJSON_GetObjectItem(root, "sequence");
    if (sequence_array == NULL){
        allocateDefaultSequence(pulse);
        made = true;
    }
    if (strcasecmp(pulse->pulsename, "Manual") == 0) {
        if (!(sequence_array && cJSON_IsArray(sequence_array))) {
            fprintf(stderr, "[Error] sequence_array is missing or not a valid JSON array. [case pulse name = manual]\n");
            exit(EXIT_FAILURE);
        }

        int count = cJSON_GetArraySize(sequence_array);
        if (count != (pulse->npulse) + 1) {
            throw std::runtime_error("Mismatch between sequence length and npulse");
        }

        double total_frac = 0.0;
        for (int i = 0; i < (pulse->npulse)+1; ++i) {
            cJSON* item = cJSON_GetArrayItem(sequence_array, i);
            if (!cJSON_IsArray(item)) continue;

            cJSON* frac  = cJSON_GetArrayItem(item, 0);
            cJSON* axis  = cJSON_GetArrayItem(item, 1);
            cJSON* angle = cJSON_GetArrayItem(item, 2);

            if (!frac || !axis || !angle) continue;

            pulse->sequence[i][0] = total_frac;               
            double fraction       = parse_fraction(frac->valuestring);
            total_frac            = total_frac + fraction;
            pulse->sequence[i][1] = total_frac ;               
            pulse->sequence[i][2] = fraction;                               // fraction
            pulse->pulse_axes[i]    = parse_axis(axis, i);                  // axis
            pulse->pulse_angles[i]  = static_cast<double>(angle->valueint); // angle
        }
        assign_sequence_indices(pulse);
        if (total_frac - 1.0 > 1e-9){
            throw std::runtime_error("Total Fraction is not 1: ");
        }
    }

    //////////////////////////
    cJSON_Delete(root);
    freeChar1d(&data);
}

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
        printMessage("    sub-key : [ dfname, naddspin, navaax, apprx, coordinate_frame, ( type : single value) ");
        printMessage("                axis, types, spins, gyros, eqs,    ( type : array ) ");
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

        bool coordinate_frame_explicit = false;
        DefectCoordinateFrame coordinate_frame =
            readDefectCoordinateFrame(defect,dfname,&coordinate_frame_explicit);
        DefectArray_setDefect_idf_coordinate_frame(dfa,idf,coordinate_frame);

        cJSON* axis_item = cJSON_GetObjectItem(defect,"axis");
        if (axis_item != NULL){
            if (!coordinate_frame_explicit){
                fprintf(stderr,"Error: Defect[%s].axis is present but coordinate_frame is missing\n",dfname);
                fprintf(stderr,"Set coordinate_frame to \"bath\" or \"qubit\" so the axis basis is explicit.\n");
                exit(EXIT_FAILURE);
            }
            double axis_raw[3] = {0.0,0.0,0.0};
            double axis_normalized[3];
            char section_label[160];
            char axis_label[160];
            snprintf(section_label,sizeof(section_label),"Defect[%s]",dfname);
            snprintf(axis_label,sizeof(axis_label),"Defect[%s].axis",dfname);
            readAxis3(defect,section_label,"axis",axis_raw,true);
            normalizeAxisOrDie(axis_raw,axis_normalized,axis_label);
            DefectArray_setDefect_idf_axis(dfa,idf,axis_normalized);
        }

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


        readDefectGyros(defect,dfa,idf,naddspin,HaveDefault);
        readDefectEqs(defect,dfa,idf,naddspin);

        ////////////////////////////////////////////////////////////////////////
        // Relative position && tensor information for each axis/additional spins
        ////////////////////////////////////////////////////////////////////////

        // Read rxyzs other wise set 0.0
        cJSON_ReadDefectInfo_IntCharDoubleArray(defect,"rxyzs",3,&(dfa->defect[idf]->rxyzs),dfa->defect[idf]->types,navaax+1,naddspin);

        // Unit-aware objects are converted to the legacy internal units while parsing.
        // Their tensor basis is explicit for the same reason as object-form zfs.
        readDefectIndexedTensorWithUnit(defect,dfa,idf,"hypf",&(dfa->defect[idf]->hypf),
                                        navaax,naddspin,coordinate_frame_explicit);
        readDefectIndexedTensorWithUnit(defect,dfa,idf,"efg",&(dfa->defect[idf]->efg),
                                        navaax,naddspin,coordinate_frame_explicit);

        // Object forms provide one physical tensor shared by every indexed
        // configuration. The legacy list remains per-index and is unchanged.
        readDefectZfs(defect,dfa,idf,navaax,coordinate_frame_explicit);

        // The object form adds an explicit frequency unit; the legacy array stays MHz.
        readDefectDetuning(defect,dfa,idf,navaax);

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
            // The path is returned anyway, so an unreadable file does not stop the run --
            // whatever the caller meant to fill from it just keeps its default. Silence
            // here is how a qubit ends up at the cell origin and the run still "works".
            if (rank==0){
                fprintf(stderr,"Warning: %s : cannot read \"%s\" -- continuing with defaults\n",
                        key, item->valuestring);
            }
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

        // Find axis index.
        //
        // NOTE the parameter called navaax here is the caller's navaax + 1: index 0 is
        // reserved for the main spin, so DefectArray_allocDefect_idf allocates navaax+1
        // entries and the valid range is 0 .. navaax.
        //
        // The bound was `iax > navaax` before, which accepted one index past the end of
        // every array this reads into -- rxyzs, hypf, efg, zfs and detuning alike. An
        // input addressing configuration navaax+1 therefore ran, and corrupted memory
        // while it did. It is refused now, which is the one way this branch can stop a
        // v1.1.0 input that used to complete.
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax >= navaax) {
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

        // Find axis index. navaax here is the caller's navaax + 1; valid range is
        // 0 .. navaax. See cJSON_ReadDefectInfo_IntCharDoubleArray for why this is
        // `>=` and not `>`.
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax >= navaax) {
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

        // Find axis index. navaax here is the caller's navaax + 1; valid range is
        // 0 .. navaax. See cJSON_ReadDefectInfo_IntCharDoubleArray for why this is
        // `>=` and not `>`.
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax >= navaax) {
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

        // Find axis index. navaax here is the caller's navaax + 1; valid range is
        // 0 .. navaax. See cJSON_ReadDefectInfo_IntCharDoubleArray for why this is
        // `>=` and not `>`.
        int iax = cJSON_GetArrayItem(itemArray1d, 0)->valueint;
        if (iax < 0 || iax >= navaax) {
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
