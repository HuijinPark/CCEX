/**
 * Unit test for the bath-coordinate rotation geometry.
 *
 * Covers everything that can be decided without running a calculation: the matrix
 * itself, the fixed point, and what the transformation preserves. The parts that
 * need a real run -- legacy-identical output, bath order and labels, multiple bath
 * files, the fatal-error paths -- are in run.sh next to this file.
 *
 * Build and run:  bash test/rotation/run.sh
 */

#include "../../include/utilities.h"
#include "../../include/general.h"
#include "../../include/hamiltonian.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

// Defined in main.cpp for the real binary; the test provides its own.
bool verbosity = false;
int  rank = 0;
int  nprocess = 1;

static int nrun  = 0;
static int nfail = 0;

static void check(bool ok, const char* what){
    nrun++;
    if (ok){
        printf("  ok    %s\n",what);
    }else{
        nfail++;
        printf("  FAIL  %s\n",what);
    }
}

static bool matClose(const double A[3][3], const double B[3][3], double tol){
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            if (fabs(A[i][j]-B[i][j]) > tol){ return false; }
        }
    }
    return true;
}

static void printMat(const char* name, const double R[3][3]){
    printf("        %s =\n",name);
    for (int i=0; i<3; i++){
        printf("          [ %18.15lf  %18.15lf  %18.15lf ]\n",R[i][0],R[i][1],R[i][2]);
    }
}

static double norm3(const double v[3]){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static double det3(const double R[3][3]){
    return R[0][0]*(R[1][1]*R[2][2] - R[1][2]*R[2][1])
         - R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0])
         + R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0]);
}

static double orthoErr(const double R[3][3]){
    double e = 0.0;
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            double rtr = R[0][i]*R[0][j] + R[1][i]*R[1][j] + R[2][i]*R[2][j];
            double err = fabs(rtr - ((i==j) ? 1.0 : 0.0));
            if (err > e){ e = err; }
        }
    }
    return e;
}

// max|R*normalize(q) - [0,0,1]|
static double zAxisErr(const double R[3][3], const double q[3]){
    double n = norm3(q);
    double e = 0.0;
    for (int i=0; i<3; i++){
        double val = (R[i][0]*q[0] + R[i][1]*q[1] + R[i][2]*q[2]) / n;
        double err = fabs(val - ((i==2) ? 1.0 : 0.0));
        if (err > e){ e = err; }
    }
    return e;
}

int main(){

    printf("\n=== bath-coordinate rotation : geometry ===\n\n");

    ////////////////////////////////////////////////////////////////////////
    // [3] bath_axis=[001], qubit_axis=[001] : identity
    ////////////////////////////////////////////////////////////////////////
    printf("[3] identity transformation, bath_axis=[001] qubit_axis=[001]\n");
    {
        double b[3] = {0.0,0.0,1.0};
        double q[3] = {0.0,0.0,1.0};
        double R[3][3];
        double I[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

        buildQubitAlignedRotation(b,q,R);
        check(matClose(R,I,0.0),"R is exactly the identity matrix");
    }

    ////////////////////////////////////////////////////////////////////////
    // [4] bath_axis=[001], qubit_axis=[111] : the required sign convention
    // [5] R * normalize([111]) = [001]
    ////////////////////////////////////////////////////////////////////////
    printf("\n[4] bath_axis=[001] qubit_axis=[111] : expected matrix\n");
    {
        double b[3] = {0.0,0.0,1.0};
        double q[3] = {1.0,1.0,1.0};
        double R[3][3];

        double s2 = sqrt(2.0), s3 = sqrt(3.0), s6 = sqrt(6.0);
        double E[3][3] = {
            { 1.0/s2, -1.0/s2,  0.0    },   // ex = [ 1,-1, 0]/sqrt(2)
            { 1.0/s6,  1.0/s6, -2.0/s6 },   // ey = [ 1, 1,-2]/sqrt(6)
            { 1.0/s3,  1.0/s3,  1.0/s3 }    // ez = [ 1, 1, 1]/sqrt(3)
        };

        buildQubitAlignedRotation(b,q,R);
        printMat("R",R);
        check(matClose(R,E,1e-14),"R matches ex=[1,-1,0]/sqrt2, ey=[1,1,-2]/sqrt6, ez=[1,1,1]/sqrt3");

        printf("\n[5] R * normalize([111]) = [001]\n");
        check(zAxisErr(R,q) < 1e-14,"the qubit axis lands on +z");

        // The axes are directions, so a rescaled input must give the same frame.
        double q2[3] = {5.0,5.0,5.0};
        double b2[3] = {0.0,0.0,7.0};
        double R2[3][3];
        buildQubitAlignedRotation(b2,q2,R2);
        check(matClose(R,R2,1e-15),"unnormalized input [0,0,7] / [5,5,5] gives the same R");
    }

    ////////////////////////////////////////////////////////////////////////
    // Proper rotation for a range of axis pairs (R^T R = I, det = +1, R q = z)
    ////////////////////////////////////////////////////////////////////////
    printf("\n[*] R is a proper rotation for several axis pairs\n");
    {
        double pairs[5][2][3] = {
            {{0.0,0.0,1.0},{1.0,1.0,1.0}},
            {{0.0,0.0,1.0},{1.0,1.0,0.0}},
            {{0.0,1.0,0.0},{1.0,1.0,1.0}},
            {{1.0,2.0,3.0},{3.0,-1.0,2.0}},
            {{0.0,0.0,1.0},{0.0,0.0,1.0}}
        };

        for (int k=0; k<5; k++){
            double R[3][3];
            buildQubitAlignedRotation(pairs[k][0],pairs[k][1],R);

            char what[160];
            snprintf(what,sizeof(what),
                     "b=[%g,%g,%g] q=[%g,%g,%g] : |R^T R - I| = %.1e, det = %.15g, |Rq-z| = %.1e",
                     pairs[k][0][0],pairs[k][0][1],pairs[k][0][2],
                     pairs[k][1][0],pairs[k][1][1],pairs[k][1][2],
                     orthoErr(R),det3(R),zAxisErr(R,pairs[k][1]));

            check(orthoErr(R) < 1e-14 && fabs(det3(R)-1.0) < 1e-14 && zAxisErr(R,pairs[k][1]) < 1e-14, what);
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // [6] a non-zero rotation origin stays fixed
    // [7] qubit-bath distances are preserved
    // [8] bath-bath distances are preserved
    ////////////////////////////////////////////////////////////////////////
    printf("\n[6,7,8] fixed point and preserved distances (r0 = [10,20,30])\n");
    {
        double b[3]  = {0.0,0.0,1.0};
        double q[3]  = {1.0,1.0,1.0};
        double r0[3] = {10.0,20.0,30.0};
        double R[3][3];
        buildQubitAlignedRotation(b,q,R);

        const int n = 6;
        double p0[n][3] = {   // bath positions, as r0 + offset
            {12.0, 20.0, 30.0},
            {10.0, 23.0, 30.0},
            {10.0, 20.0, 34.0},
            {11.0, 21.0, 31.0},
            { 8.0, 21.5, 33.0},
            {15.0, 19.0, 32.0}
        };
        double p[n][3];
        for (int i=0; i<n; i++){
            for (int k=0; k<3; k++){ p[i][k] = p0[i][k]; }
            rotateAboutPoint(p[i],R,r0);
        }

        // [6]
        double r0moved[3] = {r0[0],r0[1],r0[2]};
        rotateAboutPoint(r0moved,R,r0);
        check(fabs(r0moved[0]-r0[0]) < 1e-15 &&
              fabs(r0moved[1]-r0[1]) < 1e-15 &&
              fabs(r0moved[2]-r0[2]) < 1e-15,"the rotation origin r0 maps to itself");

        // [7]
        double worst_qb = 0.0;
        for (int i=0; i<n; i++){
            double e = fabs(dist(p0[i],r0) - dist(p[i],r0));
            if (e > worst_qb){ worst_qb = e; }
        }
        {
            char what[120];
            snprintf(what,sizeof(what),"every qubit-bath distance preserved (worst = %.1e)",worst_qb);
            check(worst_qb < 1e-13,what);
        }

        // [8]
        double worst_bb = 0.0;
        for (int i=0; i<n; i++){
            for (int j=i+1; j<n; j++){
                double e = fabs(dist(p0[i],p0[j]) - dist(p[i],p[j]));
                if (e > worst_bb){ worst_bb = e; }
            }
        }
        {
            char what[120];
            snprintf(what,sizeof(what),"every bath-bath distance preserved (worst = %.1e)",worst_bb);
            check(worst_bb < 1e-13,what);
        }

        // A spin sitting exactly on the qubit axis must land exactly on the new z axis:
        // r0 + (1,1,1) -> r0 + (0,0,sqrt(3)). This is the whole point of the transform.
        double onaxis[3] = {11.0,21.0,31.0};
        rotateAboutPoint(onaxis,R,r0);
        check(fabs(onaxis[0]-r0[0]) < 1e-14 &&
              fabs(onaxis[1]-r0[1]) < 1e-14 &&
              fabs(onaxis[2]-r0[2]-sqrt(3.0)) < 1e-14,
              "a spin on the [111] qubit axis lands on the new z axis, at z = r0z + sqrt(3)");
    }

    ////////////////////////////////////////////////////////////////////////
    // [13,14] a rank-2 tensor transforms as R * T * R^T
    ////////////////////////////////////////////////////////////////////////
    printf("\n[13,14] non-diagonal tensor : T_new = R * T * R^T\n");
    {
        double b[3] = {0.0,0.0,1.0};
        double q[3] = {1.0,1.0,1.0};
        double R[3][3];
        buildQubitAlignedRotation(b,q,R);

        // A real quadrupole tensor, lifted from the h-BN V_B run in
        // test/code_verification/ -- symmetric, and genuinely non-diagonal.
        double T0[3][3] = {
            { -70.76,  0.00, -1.21 },
            {   0.00,-70.31,  1.66 },
            {  -1.21,  1.66, 141.07 }
        };
        MatrixXcd T = MatrixXcd::Zero(3,3);
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){ T(i,j) = doublec(T0[i][j],0.0); }
        }

        // Expected value computed here, by hand, from the definition -- NOT by calling
        // the same helper the test is checking.
        double E[3][3];
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                double sum = 0.0;
                for (int k=0; k<3; k++){
                    for (int l=0; l<3; l++){ sum += R[i][k]*T0[k][l]*R[j][l]; }
                }
                E[i][j] = sum;
            }
        }

        MatrixXcd Trot = rotateTensor(T,R);
        double worst = 0.0;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                worst = std::max(worst,fabs(Trot(i,j).real() - E[i][j]));
                worst = std::max(worst,fabs(Trot(i,j).imag()));
            }
        }
        {
            char what[140];
            snprintf(what,sizeof(what),"rotateTensor matches the hand-computed R*T*R^T (worst = %.1e)",worst);
            check(worst < 1e-12,what);
        }

        // A similarity transform by an orthogonal matrix keeps the trace and the
        // Frobenius norm. Those two are what the integration test can still compare
        // after the log has rounded the components to two decimals.
        double tr0 = T0[0][0]+T0[1][1]+T0[2][2];
        double tr1 = Trot(0,0).real()+Trot(1,1).real()+Trot(2,2).real();
        double f0 = 0.0, f1 = 0.0;
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){ f0 += T0[i][j]*T0[i][j]; f1 += Trot(i,j).real()*Trot(i,j).real(); }
        }
        check(fabs(tr0-tr1) < 1e-10 && fabs(sqrt(f0)-sqrt(f1)) < 1e-10,
              "trace and Frobenius norm are unchanged by the transform");

        // Off-diagonal elements really do move -- otherwise the check above is vacuous.
        check(fabs(Trot(0,1).real() - T0[0][1]) > 1.0,
              "the tensor is genuinely re-expressed, not left alone");

        // An unset tensor (BathArray leaves those 0x0) must pass through untouched.
        MatrixXcd empty;
        check(rotateTensor(empty,R).rows() == 0,"an unset 0x0 tensor is returned untouched");
    }

    ////////////////////////////////////////////////////////////////////////
    // [16] the point-dipole tensor is covariant :
    //      A(R*r1, R*r2) = R * A(r1,r2) * R^T
    ////////////////////////////////////////////////////////////////////////
    printf("\n[16] point-dipole tensor : recomputing after the rotation == rotating the tensor\n");
    {
        double b[3] = {0.0,0.0,1.0};
        double q[3] = {1.0,1.0,1.0};
        double r0[3] = {10.0,20.0,30.0};
        double R[3][3];
        buildQubitAlignedRotation(b,q,R);

        double gq = -17608.597050;  // electron
        double gb =   6.728284;     // 13C

        double pts[4][3] = {
            {12.0, 20.0, 30.0},
            {11.0, 21.0, 31.0},
            { 8.0, 21.5, 33.0},
            {15.0, 19.0, 32.0}
        };

        double worst = 0.0;
        for (int k=0; k<4; k++){
            double p[3] = {pts[k][0],pts[k][1],pts[k][2]};

            // (a) tensor from the source geometry, then transformed
            MatrixXcd A_src = calPointDipoleTensor(r0,p,gq,gb);
            MatrixXcd A_then_rot = rotateTensor(A_src,R);

            // (b) geometry transformed first, then the tensor recomputed. r0 is the
            //     fixed point, so the qubit itself does not move.
            double p_rot[3] = {p[0],p[1],p[2]};
            rotateAboutPoint(p_rot,R,r0);
            MatrixXcd A_rot_geom = calPointDipoleTensor(r0,p_rot,gq,gb);

            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    worst = std::max(worst,std::abs(A_then_rot(i,j) - A_rot_geom(i,j)));
                }
            }
        }
        {
            // Relative to tensor components of order 1e2-1e4 radkHz, this is round-off.
            char what[150];
            snprintf(what,sizeof(what),
                     "R*A(r)*R^T equals A(R*r) for every off-diagonal element too (worst = %.1e radkHz)",worst);
            check(worst < 1e-9,what);
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // [1] Config defaults : rotation off, both axes on +z
    ////////////////////////////////////////////////////////////////////////
    printf("\n[1] Config defaults\n");
    {
        Config* cnf = Config_init();
        double* ba  = Config_getRot_bath_axis(cnf);
        double* qax = Config_getRot_qubit_axis(cnf);

        check(Config_getRot_enabled(cnf) == false,"rotation is off by default");
        check(ba[0]==0.0  && ba[1]==0.0  && ba[2]==1.0 ,"default bath_axis  = [0,0,1]");
        check(qax[0]==0.0 && qax[1]==0.0 && qax[2]==1.0,"default qubit_axis = [0,0,1]");
        // Config_freeAll would free members this test never allocated; the leak of one
        // Config in a test binary that is about to exit is not worth guarding against.
    }

    ////////////////////////////////////////////////////////////////////////
    printf("\n------------------------------------------------------------\n");
    printf("  %d checks, %d failed\n\n",nrun,nfail);
    return (nfail == 0) ? 0 : 1;
}
