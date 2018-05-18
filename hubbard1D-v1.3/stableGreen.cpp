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

int main()
{
    //  TOOGLE DEBUGGING MODE.
    const bool debug = false;
    
    //  SET PARAMETERS.
    const int N = 4;  //  # sites
    const double dt = 0.125;  //  time subinterval width. error scales as dt^2
    const double beta = 9.5;  //  inverse temperature
    const int L = beta / dt;  //  # slices
    const double t = 1.;  //  hopping
    const double U = 5.;  //  interaction
    const double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;  //  HS transformation parameter
    const double mu = 0.;  //  chemical potential
    
    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 12345;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double decisionMaker;
    //  to accept or not to accept, hence the question.
    const int totalMCSweeps = 10;
    //  number of measurements will be totalMCSweeps * L ,
    //  i.e. we make a measurement every lattice sweep (N steps), then find correlations in post-processing
    const int totalMCSteps = totalMCSweeps * N * L;
    const int greenAfreshFreq = L + 1 ;  //   how often to calculate Green's functions afresh
    const int k = L / 2;  //  k = # intervals. Must be commensurate with L and k < L
    //NOTE:
    //  measured in SPATIAL lattice sweeps,
    //  i.e. greenAfreshFreq = L corresponds to 1 sweep.
    //  Comment the default definition and comment the following lines to test.
    //    const int greenAfreshFreq = 1;  //  corresponds to computing always computing the Green's function from scratch after N steps
    //    const int greenAfreshFreq >= totalMCSweeps * N;  //  greenAfreshFreq = totalMCSweeps * N corresponds to the updates + wrapping version
    
    
    // -- INITIALIZATION ---
    
    
    //  HOPPING + CHEM. POTENTIAL MATRIX FOR 1D CHAIN W/ PBCs.
    const Eigen::MatrixXd K = genHoppingMatrix(N, mu);
    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Eigen::MatrixXd h = genHsMatrix(L, N);  //    HS field h_{l = 1,...,L ; i = 1,...,N}
    //  COMPUTE MATRIX 'PREFACTOR' OF THE B-MATRICES e^(t dt K).
    const Eigen::MatrixXd BpreFactor = (t * dt * K).exp();
    
    initialPrint(N, dt, beta, L, t, U, mu, totalMCSteps, debug, K, BpreFactor, h);
    
    //  GENERATE THE B-MATRICES.
    Eigen::MatrixXd Bup[L];
    genBmatrix(Bup, true, nu, N, L, h, BpreFactor);
    Eigen::MatrixXd Bdown[L];
    genBmatrix(Bdown, false, nu, N, L, h, BpreFactor);
    
    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green GreenUp(N, L);
    GreenUp.computeGreenNaive(Bup, 1);  //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    Eigen::MatrixXd Gup = GreenUp.getG();
    Green GreenDown(N, L);
    GreenDown.computeGreenNaive(Bdown, 1);  //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    Eigen::MatrixXd Gdown = GreenDown.getG();
    
    std::cout << "\nGup\n" << Gup << std::endl;
    std::cout << "\nGdown\n" << Gdown << std::endl;
    
    GreenUp.computeStableGreenNaive(Bup, 1, k);
    Gup = GreenUp.getG();
    std::cout << "\nGup\n" << Gup << std::endl;
    
    return 0;
    
}
