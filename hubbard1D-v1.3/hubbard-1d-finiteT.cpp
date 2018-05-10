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
    //  TOOGLE DEBUGGING MODE
    const bool debug = false;
    
    //  SET PARAMETERS
    const int N = 4;  //  number of sites
    const double dt = 0.125;  //  time subinterval width. error scales as dt^2
    const double beta = 1.;  //  inverse temperature (imaginary time of equivalent maximum temperature of the (d+1) classical system)
    const int L = beta/dt;  //  number of imaginary time subintervals
    const double t = 1.;  //  hopping
    const double U = 4.;  //  interaction
    const double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;  //  HS parameter
    const double mu = 0.;  //  chemical potential
    
    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES
    const int seed = 12345;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double decisionMaker;  //  to accept or not to accept, hence the question.
    const int totalMCSweeps = 1000;  //  same as the number of measurements, i.e. we make a measurement every sweep (NL sweeps), then find correlations in post-processing
    const int totalMCSteps = totalMCSweeps * N * L;
    const int greenAfreshFreq = 1;  //   how often to calculate Green functions afresh (measured in SPATIAL lattice sweeps)
    int latticeSweep = 0;
    
    
    // --- INITIALIZATION ---
    
    
    //  HOPPING FOR 1D CHAIN W/ PBCs
    const Eigen::MatrixXd K = genHoppingMatrix(N);
    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY
    Eigen::MatrixXd h = genHsMatrix(L, N);  //    HS field h_{l = 1,...,L ; i = 1,...,N}
    //  COMPUTE MATRIX 'PREFACTOR' OF THE B-MATRICES e^{t\delta\tau K}
    const Eigen::MatrixXd BpreFactor = (t * dt * K).exp();
    
    initialPrint(N, dt, beta, L, t, U, mu, totalMCSteps, debug, K, BpreFactor, h);
    
    //  GENERATE THE B-MATRICES
    Eigen::MatrixXd Bup[L];
    genBmatrix(Bup, true, nu, N, L, h, BpreFactor);
    Eigen::MatrixXd Bdown[L];
    genBmatrix(Bdown, false, nu, N, L, h, BpreFactor);

    //  SPIN-UP GREEN FUNCTION
    Green GreenUp(N);
    GreenUp.computeGreenNaive(Bup, L);
    Eigen::MatrixXd Mup = GreenUp.getM();
    Eigen::MatrixXd Gup = GreenUp.getG();
    //  SPIN-DOWN GREEN FUNCTION
    Green GreenDown(N);
    GreenDown.computeGreenNaive(Bdown, L);
    Eigen::MatrixXd Mdown = GreenDown.getM();
    Eigen::MatrixXd Gdown = GreenDown.getG();
    
    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO
    Eigen::VectorXd uUp;
    Eigen::VectorXd uDown;
    Eigen::VectorXd wUp;
    Eigen::VectorXd wDown;
    double alphaUp;
    double alphaDown;
    double dUp;
    double dDown;
    double accRatio;

    //  INITIALIZE ARRAY TO STORE THE WEIGHT OF THE ACCEPTED CONFIGURATIONS
    double weights[totalMCSteps];
    weights[0] = Mup.determinant() * Mdown.determinant();
    
    std::cout << weights[0] << std::endl;

    //  INITIALIZE (l, i) <- (0, 0).
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE, THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0;
    int i = 0;


    // --- MC LOOP ---


    std::cout << "\n\nMC loop started. Progress:\n\n";
    for (int step = 0; step < totalMCSteps; step++)
    {
//        //  DISPLAY PROGRESS OF THE RUN
//        if ( (step + 1)  % (totalMCSteps/8) == 0 )
//        {
//            std::cout << (step + 1)*1. / totalMCSteps * 100 << " %" << std::endl;
//        }

        //  COMPUTE THE ACCEPTANCE RATIO
        alphaUp = ( exp( -2 * h(l, i) * nu ) - 1 );
        alphaDown = ( exp( 2 * h(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup(i, i) ) );
        dDown = ( 1 + alphaDown  * ( 1 - Gdown(i, i) ) );
        accRatio = dUp * dDown;

        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP
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

            //  RANK-ONE UPDATE O(N^2).
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

        if (i < N - 1)
        {
            i += 1;
        }
        else
        {
            latticeSweep += 1;

            //  REBUILD B-MATRICES
            Bup[l] = regenB(true, nu, N, h.row(l), BpreFactor);
            Bdown[l] = regenB(false, nu, N, h.row(l), BpreFactor);

            if (l < L - 1)
            {
                //  DECIDE WHETHER TO COMPUTE GREEN MATRIX AFRESH OR TO WRAP
                if (latticeSweep == greenAfreshFreq)
                {
                    //  COMPUTE GREEN'S MATRIX AFRESH
                    
                    //  SPIN-UP GREEN FUNCTION
                    Green GreenUp(N);
                    GreenUp.computeWrappedGreenNaive(Bup, L, l);
                    Gup = GreenUp.getG();
                    //  SPIN-DOWN GREEN FUNCTION
                    Green GreenDown(N);
                    GreenDown.computeWrappedGreenNaive(Bdown, L, l);
                    Gdown = GreenDown.getG();
                    
//                    std::cout << "\nMup\n" << GreenUp.getM() * Gup;
//                    std::cout << "\nMdown\n" << GreenDown.getM() * Gdown;
                    
                    latticeSweep = 0;
                }
                else
                {
                    //  WRAPPING
                    Gup = Bup[l] * Gup * Bup[l].inverse();
                    Gdown = Bdown[l] * Gdown * Bdown[l].inverse();
                }
                l += 1;
                i = 0;
            }
            else
            {
                //  DECIDE WHETHER TO COMPUTE GREEN MATRIX AFRESH OR TO WRAP
                if (latticeSweep == greenAfreshFreq)
                {
                    //  COMPUTE GREEN'S MATRIX AFRESH
                    
                    //  SPIN-UP GREEN FUNCTION
                    Green GreenUp(N);
                    GreenUp.computeGreenNaive(Bup, L);
                    Gup = GreenUp.getG();
                    //  SPIN-DOWN GREEN FUNCTION
                    Green GreenDown(N);
                    GreenDown.computeGreenNaive(Bdown, L);
                    Gdown = GreenDown.getG();
                    
//                    std::cout << "\nMup\n" << GreenUp.getM() * Gup;
//                    std::cout << "\nMdown\n" << GreenDown.getM() * Gdown;
                    
                    latticeSweep = 0;
                }
                else
                {
                    //  WRAPPING
                    Gup = Bup[l] * Gup * Bup[l].inverse();
                    Gdown = Bdown[l] * Gdown * Bdown[l].inverse();
                }
                l = 0;
                i = 0;
            }
        }

    }

    //  SAVE OUTPUT
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
