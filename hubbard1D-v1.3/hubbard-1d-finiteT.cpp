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
    const double beta = 8;  //  inverse temperature
    const int L = beta / dt;  //  # slices
    const double t = 1.;  //  hopping
    const double U = 4.;  //  interaction
    const double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;  //  HS transformation parameter
    const double mu = 0.4 * U;  //  chemical potential
    const int greenAfreshFreq = 4 ;  //   how often to calculate Green's functions afresh
    //NOTE:
    //  measured in SPATIAL lattice sweeps,
    //  i.e. greenAfreshFreq = L corresponds to 1 sweep.
    //    const int greenAfreshFreq = 1;  //  always compute the Green's function from scratch after N steps
    //    const int greenAfreshFreq = totalMCSweeps * N;  //  updates + wrapping version
    const int newL = L / 4;  //  newL = # intervals in which the product of B's is divided to stabilize.
    //  Must be commensurate with L and k < L
    
    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 12345;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double decisionMaker;
        //  to accept or not to accept, hence the question.
    const int totalMCSweeps = 65536; //   65536 to reproduce
        //  number of measurements will be totalMCSweeps * L ,
        //  i.e. we make a measurement every lattice sweep (N steps), then find correlations in post-processing
    const int totalMCSteps = totalMCSweeps * N * L;

    
    // -- INITIALIZATION ---
    
    
    //  HOPPING + CHEM. POTENTIAL MATRIX FOR 1D CHAIN W/ PBCs.
    const Eigen::MatrixXd K = genHoppingMatrix(N);
    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Eigen::MatrixXd h = genHsMatrix(L, N);  //    HS field h_{l = 1,...,L ; i = 1,...,N}
    //  COMPUTE MATRIX 'PREFACTOR' OF THE B-MATRICES e^(t dt K).
    const Eigen::MatrixXd BpreFactor = (t * dt * K).exp();
    
    initialPrint(N, dt, beta, L, t, U, mu, totalMCSteps, debug, K, BpreFactor, h);
    
    //  GENERATE THE B-MATRICES.
    Eigen::MatrixXd Bup[L];
    genBmatrix(Bup, true, nu, N, L, h, BpreFactor, dt, mu);
    Eigen::MatrixXd Bdown[L];
    genBmatrix(Bdown, false, nu, N, L, h, BpreFactor, dt, mu);

    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green GreenUp(N, L);
    GreenUp.computeGreenNaive(Bup, L - 1);  //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    Eigen::MatrixXd Gup = GreenUp.getG();
    Green GreenDown(N, L);
    GreenDown.computeGreenNaive(Bdown, L - 1);  //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
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
    double* weights = new double[totalMCSweeps * L];
    double weight = GreenUp.getM().determinant() * GreenDown.getM().determinant();
    weights[0] = weight;

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE, THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0;
    int i = 0;
    int latticeSweepUntilAfresh = 0;
    int sweep = 0;
    
    
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

        if (decisionMaker <= abs(accRatio) )
        {
            //  KEEP TRACK OF WEIGHT
            weight = accRatio * weight;
            //  FLIP A SPIN
            h(l, i) *= -1;

            //  RANK-ONE UPDATE -> O(N^2)
            uUp = uSigma(N, Gup, i) * alphaUp / ( 1 + alphaUp  * ( 1 - Gup(i, i) ) );
            uDown = uSigma(N, Gdown, i) * alphaDown /( 1 + alphaDown * ( 1 - Gdown(i, i) ) ) ;
            wUp = wSigma(N, Gup, i);
            wDown = wSigma(N, Gdown, i);
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < N; y++)
                {
                    Gup(x, y) -= uUp(x) * wUp(y);
                    Gdown(x, y) -= uDown(x) * wDown(y);
                }
            }
        }
        
        
        // --- COMPUTE WRAPPED GREEN'S FUNCTIONS. MEASUREMENTS. ---
        
        
        if (i < N - 1)
        {   //  CONTINUE LOOPING THROUGH THE SPATIAL LATTICE
            i += 1;
        }
        else
        {
            //  MEASUREMENTS
            sweep += 1;
            
                //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
            weights[sweep] = weight;
            
            
            //  EITHER WRAP OR COMPUTE GREEN'S FUNCTIONS FROM SCRATCH.
            latticeSweepUntilAfresh += 1;
            //  DEAL WITH THE GREEN'S FUNCTIONS.
            
                //  REBUILD B-MATRICES.
            Bup[l] = regenB(true, nu, N, h.row(l), BpreFactor, dt, mu);
            Bdown[l] = regenB(false, nu, N, h.row(l), BpreFactor, dt, mu);

                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == greenAfreshFreq)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                Green GreenUp(N, L);
                Green GreenDown(N, L);
                GreenUp.computeStableGreenNaive(Bup, l, newL);
                GreenDown.computeStableGreenNaive(Bdown, l, newL);
                //  Uncomment to compute the product in the naive, unstable manner
//                GreenUp.computeGreenNaive(Bup, l);
//                GreenDown.computeGreenNaive(Bdown, l);
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
        file1 << greenAfreshFreq << '\n';
        file1 << newL << '\n';
    }
    std::ofstream file2("plots/weights.txt");
    if (file2.is_open())
    {
        for (int s = 0; s < totalMCSweeps; s++)
        {
            file2 << weights[s] << '\t' << std::copysign( 1. , weights[s] ) << '\n';
        }
    }
    return 0;
}

//int main()
//{
//    int nSamples = 10
//    for (int sim = 0; sim < nSamples; sim++)
//    {
//
//    }
//
//    return 0
//}
