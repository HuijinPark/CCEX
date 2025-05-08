// ReadOption
#include "../include/option.h"

// Structures
#include "../include/qubit.h"
#include "../include/bath.h"
#include "../include/partition.h"
#include "../include/defect.h"
#include "../include/cluster.h"
#include "../include/pulse.h"
#include "../include/output.h"
#include "../include/general.h"

int seed; // for random genearator
int rank,nprocess;
bool verbosity = true;

int main(int argc, char* argv[]){

    //Set up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get process id
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess); // Get number of processes

    BathArray* ba  = BathArray_init();
    Cluster*   cls = Cluster_init();
    cJSON_readOptionCluster(cls, "../config/cce_general.json");

    int nspin =  30;
    int nqubit = 2;
    int sK = 6;

    // Set nqubit #
    BathArray_setNspin(ba,nspin);
    
    // Allocate bath memory for nspin #
    BathArray_allocBath(ba,nqubit);

    //printf("\n    >> Bath defect related hyperfine tensors\n\n");
    //BathArray_reportBath_hypf_sub(ba);
    
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
    ////////////////////////////////////////////
    // Set qubit properties
    ////////////////////////////////////////////
    QubitArray* qa = QubitArray_init();

    //int nqubit = 1;
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

        //QubitArray_reportQubit_i(qa,iqubit);
        printf("\n");
    }
    //BathArray_reportBath(ba);


    // Change "bath" ==>  "Point" type. //
    printf("NSPIN: %d\n", ba->nspin);
    Partition_info pinfo;
    set_spinfinite(ba, qa, sK, &pinfo);
    BathArray_reportBath(ba);
    Cluster_report(cls);
    printf("NSPIN: %d\n", ba->nspin);

    ////// Run constrained k-means algorithm! //
    if (strcmp(cls->method, "pcce")==0){
        int ClusterNum   = pinfo.npartition;
        int TotalSpinNum = pinfo.partition_nspin;
        printf("ClusterNum: %d\n", ClusterNum);
        printf("TotalSpinNum: %d\n", TotalSpinNum);
        Point* best_centers      = (Point*) malloc(sizeof(Point) * ClusterNum);
        int*   best_assigned_idx = (int*)   malloc(sizeof(int)   * TotalSpinNum);
       
        ba->centers = best_centers;
        ba->assigned_idx = best_assigned_idx;

        simulator_cluster_partition(ba, &pinfo, cls);
    
        // Find the "clusinfo" data. 
        //Cluster_clusterize(cls, ba, cnf);
    }


    MPI_Finalize();
}
