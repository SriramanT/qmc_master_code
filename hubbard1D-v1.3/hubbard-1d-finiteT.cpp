//
//  hubbard-1d-finiteT.cpp
//  
//
//  Created by Francisco Brito on 09/05/2018.
//
//  This program simulates the Hubbard model on a 1D chain
//  using auxiliary field (or determinant) Quantum Monte Carlo.
//  The used notation is based on the lecture notes "Numerical Methods for Quantum Monte Carlo
//  Simulations of the Hubbard Model by Zhaojun Bai, Wenbin Chen, Richard Scalettar, and
//  Ichitaro Yamazaki (2006)
//

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <typeinfo>
#include <cmath>
#include <random>
#include <time.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>
#include <fstream>

#include "matrixgen.h"
#include "prints.h"
#include "green.h"

int main()
{
    //  TOOGLE DEBUGGING MODE.
    const bool debug = false;
    
    //  SET PARAMETERS.
    const int N = 50;  //  # sites
    const double dt = 0.1;  //  time subinterval width. error scales as dt^2
    const double beta = 1.;  //  inverse temperature
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
    const int k = 2;  //  k = # intervals. Must be commensurate with L and k < L
        //NOTE:
        //  measured in SPATIAL lattice sweeps,
        //  i.e. greenAfreshFreq = L corresponds to 1 sweep.
        //  Comment the default definition and comment the following lines to test.
//    const int greenAfreshFreq = 1;  //  corresponds to computing always computing the Green's function from scratch
//    const int greenAfreshFreq = totalMCSweeps * N;  //  greenAfreshFreq = totalMCSweeps * N corresponds to the updates + wrapping version
    
    
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
    GreenUp.computeGreenNaive(Bup, L - 1);  //  start at l = L - 1
    Eigen::MatrixXd Gup = GreenUp.getG();
    Green GreenDown(N, L);
    GreenDown.computeGreenNaive(Bdown, L - 1);  //  start at l = L - 1
    Eigen::MatrixXd Gdown = GreenDown.getG();
    
    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    Eigen::VectorXd uUp;
    Eigen::VectorXd uDown;
    Eigen::VectorXd wUp;
    Eigen::VectorXd wDown;
    double alphaUp;
    double alphaDown;
    double dUp;
    double dDown;
    double accRatio;

    //  INITIALIZE ARRAY TO STORE THE WEIGHT OF THE ACCEPTED CONFIGURATIONS.
    double weights[totalMCSteps];
    weights[0] = GreenUp.getM().determinant() * GreenDown.getM().determinant();

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE, THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0;
    int i = 0;
    int latticeSweepUntilAfresh = 0;
    
    
    // --- MC LOOP ---

    
    std::cout << "\n\nMC loop started. Progress:\n\n";
    for (int step = 0; step < totalMCSteps; step++)
    {
        //  DISPLAY PROGRESS OF THE RUN.
        if ( (step + 1)  % (totalMCSteps/8) == 0 )
        {
            std::cout << (step + 1)*1. / totalMCSteps * 100 << " %" << std::endl;
        }

        //  COMPUTE THE ACCEPTANCE RATIO.
        alphaUp = ( exp( -2 * h(l, i) * nu ) - 1 );
        alphaDown = ( exp( 2 * h(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup(i, i) ) );
        dDown = ( 1 + alphaDown  * ( 1 - Gdown(i, i) ) );
        accRatio = dUp * dDown;

        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
        decisionMaker = dis(gen);

        if (decisionMaker <= accRatio)
        {
            //  STORE WEIGHT
            if (step < totalMCSteps - 1)
            {
                weights[step + 1] = accRatio * weights[step];
            }
            //  FLIP A SPIN
            h(l, i) *= -1;

            //  RANK-ONE UPDATE -> O(N^2)
            uUp = uSigma(N, Gup, i);
            uDown = uSigma(N, Gdown, i);
            wUp = wSigma(N, Gup, i);
            wDown = wSigma(N, Gdown, i);
            //NOTE:
            //  kroneckerProduct is the tensor product of two vectors. .eval() is needed because of an Eigen quirk. Check documentation.
            Gup -= alphaUp/( 1 + alphaUp  * ( 1 - Gup(i, i) ) ) * kroneckerProduct(uUp, wUp.transpose() ).eval();
            Gdown -= alphaDown/( 1 + alphaDown * ( 1 - Gdown(i, i) ) ) * kroneckerProduct(uDown, wDown.transpose() ).eval();
        }
        else
        {
            //  STORE WEIGHT
            if (step < totalMCSteps - 1)
            {
                weights[step + 1] = weights[step];
            }
        }
        
        
        // --- COMPUTE WRAPPED GREEN'S FUNCTIONS. MEASUREMENTS. ---
        
        
        if (i < N - 1)
        {   //  CONTINUE LOOPING THROUGH THE SPATIAL LATTICE
            i += 1;
        }
        else
        {   //  EITHER WRAP OR COMPUTE GREEN'S FUNCTIONS FROM SCRATCH. MAKE MEASUREMENTS
            latticeSweepUntilAfresh += 1;
            
            
            //  MEASUREMENTS
            
            
            //  make a measurement
            
            //  DEAL WITH THE GREEN'S FUNCTIONS.
            
                //  REBUILD B-MATRICES.
            Bup[l] = regenB(true, nu, N, h.row(l), BpreFactor);
            Bdown[l] = regenB(false, nu, N, h.row(l), BpreFactor);

                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == greenAfreshFreq)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                Green GreenUp(N, L);
                Green GreenDown(N, L);
                if (l == 0)
                {   //  THIS PARTICULAR CASE IS WHEN WE GO BACK TO THE ORIGINAL ORDER M = B_{L} B_{L-1} ... B_{0}.
                    GreenUp.computeGreenNaive(Bup, l);
                    GreenDown.computeGreenNaive(Bdown, l);
                    
                }
                else
                {   //  THIS ONE USES M = B_l B_{l-1} ... B_{0} B_{L} B_{L-1} TO COMPUTE THE GREEN'S FUNCTIONS.
                    GreenUp.computeGreenNaive(Bup, l);
                    GreenDown.computeGreenNaive(Bdown, l);
                }
                Gup = GreenUp.getG();
                Gdown = GreenDown.getG();
                latticeSweepUntilAfresh = 0;
            }
            else
            {   //  WRAPPING.
                Gup = Bup[l] * Gup * Bup[l].inverse();
                Gdown = Bdown[l] * Gdown * Bdown[l].inverse();
            }
            if (l < L - 1)
            {
                l += 1;
                i = 0;
            }
            else
            {
                l = 0;
                i = 0;
            }
        }
    }   //  END OF MC LOOP.

    //  SAVE OUTPUT.
    std::ofstream file1("plots/simulationParameters.txt");
    if (file1.is_open())
    {
        file1 << N << '\n';
        file1 << dt << '\n';
        file1 << beta << '\n';
        file1 << L << '\n';
        file1 << t << '\n';
        file1 << U << '\n';
        file1 << mu << '\n';
        file1 << totalMCSweeps << '\n';
        file1 << greenAfreshFreq << '\n' ;
    }
    std::ofstream file2("plots/weights.txt");
    if (file2.is_open())
    {
        for (int s = 0; s < totalMCSteps; s++)
        {
            file2 << weights[s] << '\n';
        }
    }
    return 0;
}
