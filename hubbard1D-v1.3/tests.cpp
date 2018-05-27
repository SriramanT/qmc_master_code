#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <typeinfo>
#include <cmath>
#include <random>
#include <time.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <fstream>

#include "matrixgen.h"
#include "prints.h"
#include "green.h"
#include "UDV.h"
#include "VDU.h"

int main()
{
    //  SET PARAMETERS.
    const int N = 4;  //  # sites
    const double dt = 0.125;  //  time subinterval width. error scales as dt^2
    const double beta = 1;  //  inverse temperature
    const int L = beta / dt;  //  # slices
    const double t = 1.;  //  hopping
    const double U = 4.;  //  interaction
    const double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;  //  HS transformation parameter
    const double mu = 0.;  //  chemical potential
    
    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 12345;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    //  to accept or not to accept, hence the question.
    const int totalMCSweeps = 10;
    //  number of measurements will be totalMCSweeps * L ,
    //  i.e. we make a measurement every lattice sweep (N steps), then find correlations in post-processing
    const int totalMCSteps = totalMCSweeps * N * L;
//    const int greenAfreshFreq = L + 1 ;  //   how often to calculate Green's functions afresh
    const int k = L / 2;  //  k = # intervals. Must be commensurate with L and k < L
    
    
    // -- INITIALIZATION ---
    
    
    //  HOPPING + CHEM. POTENTIAL MATRIX FOR 1D CHAIN W/ PBCs.
    const Eigen::MatrixXd K = genHoppingMatrix(N);
    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Eigen::MatrixXd h = genHsMatrix(L, N);  //    HS field h_{l = 1,...,L ; i = 1,...,N}
    //  COMPUTE MATRIX 'PREFACTOR' OF THE B-MATRICES e^(t dt K).
    const Eigen::MatrixXd BpreFactor = (t * dt * K).exp();
    
    //  GENERATE THE B-MATRICES.
    Eigen::MatrixXd Bup[L];
    genBmatrix(Bup, true, nu, N, L, h, BpreFactor, dt, mu);
    
    //  ALLOCATE MEMORY TO STORE THE PARTIAL PRODUCTS
    
    Eigen::MatrixXd Us[L];
    Eigen::MatrixXd Vs[L];
    Eigen::MatrixXd Ds[L];
    
    Green GreenUp(N, L);
    GreenUp.computeStableGreenNaiveR(Bup, L -1, k);
    std::cout << GreenUp.getG() << std::endl;
    GreenUp.computeStableGreenNaiveL(Bup, L -1, k);
    std::cout << GreenUp.getG() << std::endl;
    GreenUp.computeGreenNaive(Bup, L -1);
    std::cout << GreenUp.getG() << std::endl;
    GreenUp.storeVDU(Bup, k, Us, Ds, Vs);
    for (int la = 0; la < k ; la++)
    {
        std::cout << Ds[la] << std::endl;
    }
//    Eigen::MatrixXd Bdown[L];
//    genBmatrix(Bdown, false, nu, N, L, h, BpreFactor, dt, mu);
//    Eigen::MatrixXd B = Bup[0];
//    std::cout << B << std::endl;
//    VDU Bmat(B);
//    Bmat.printMatrixToDecompose();
//    Eigen::MatrixXd Q = Bmat.QR_and_getU();
//    std::cout << Q << std::endl;
//    Eigen::MatrixXd D = Bmat.getD();
//    std::cout << D << std::endl;
//    Eigen::MatrixXd V = Bmat.getV();
//    std::cout << V << std::endl;
//    std::cout << V * D * Q << std::endl;
    
    std::cout << 2 % 4 << std::endl;
    
    return 0;
    
}
