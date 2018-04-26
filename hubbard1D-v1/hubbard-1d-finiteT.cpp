//
//  hubbard-1d-finiteT.cpp
//  
//
//  Created by Francisco Brito on 09/03/2018.
//
//  This code constitutes a first attempt at simulating the Hubbard model on a 1D chain
//  using auxiliary field quantum Monte Carlo.
//  The notation in use is based on the lecture notes "Numerical Methods for Quantum Monte Carlo
//  Simulations of the Hubbard Model by Zhaojun Bai, Wenbin Chen, Richard Scalettar, and
//  Ichitaro Yamazaki
//

//  Standard header files
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

//  Header files created for the program
#include "auxFunctions.h"
#include "parameters.h"
#include "matrices.h"

// --- DEFINITIONS ---

#define NSITES 50

//  Physical parameters
double dt;                                   //  time subinterval width. error scales as dt^2
double beta;                                 //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
int L;                                       //  number of imaginary time subintervals
double t;                                    //  hopping parameters
double U;                                    //  interaction energy
double nu;                                   //  Hubbard Stratonovich transformation parameter

//  Some general use variables
int i;
int j;
int l;
int m;
int step;
int i_chosen;
int l_chosen;

//  Toggle prints
bool printsOn = false;                                               // 0 - NO PRINTS, 1 - PRINTS

//  Things to do with random numbers throughout the code
const static int seed = 12345;                              //  set a seed
std::mt19937 gen(seed);                                     //  mt19937 algorithm to generate random numbers
std::uniform_real_distribution<> dis(0.0, 1.0);             //  set the distribution from which to draw numbers
double decisionMaker;                                       //  to accept or not to accept, hence the question.

//  Set Monte Carlo-specific variables
const static int totalMCSteps = 10;
const static int W = 3000;                                  //  warm-up steps
const static int autoCorrTime = 500;                        //  auto-correlation time
const static int M = (totalMCSteps - W)/autoCorrTime;       //  number of measurements


// --- MAIN ---


int main()
{
    std::cout << "\n\nAuxiliary field QMC for 1D Hubbard chain\n" << std::endl;
    std::cout << "\nNumber of sites: " << NSITES << std::endl;
    
    //  Set physical parameters and Trotter parameter (error scales as dt^2)
    dt = 0.1;                                                   //  time subinterval width. error scales as dt^2
    beta = 1.;                                                  //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
    U = 5;                                                      //  interaction energy
    
    
    L = beta/dt;                                                //  number of imaginary time subintervals
    t = 1.;                                                     //  hopping parameter (default to 1, and measure energy in units of hopping parameter)
    nu = acosh( exp( U * dt / 2 ) );                            //  Hubbard Stratonovich transformation parameter
    printParameters();
    std::cout << "Number of MC sweeps: " << totalMCSteps / (NSITES*L) << " (" << totalMCSteps << " steps)" << "\n\n";
    
    // --- INITIALIZATION ---
    Eigen::MatrixXd h;
    Eigen::MatrixXd K;
    Eigen::MatrixXd BpOld[L]; // B plus (spin-up)
    Eigen::MatrixXd BmOld[L]; // B minus (spin-down)
    Eigen::MatrixXd BpNew[L]; // B plus (spin-up)
    Eigen::MatrixXd BmNew[L]; // B minus (spin-up)
    Eigen::MatrixXd MPlusOld = Eigen::MatrixXd::Identity(NSITES,NSITES);
    Eigen::MatrixXd MMinusOld = Eigen::MatrixXd::Identity(NSITES,NSITES);
    Eigen::MatrixXd MPlusNew;
    Eigen::MatrixXd MMinusNew;
    Eigen::MatrixXd GreenPlus;
    Eigen::MatrixXd GreenMinus;
    double detNewPlus;
    double detNewMinus;
    double acceptanceRatio;
    //  Acceptance ratio computed by update of Green's functions
    double r;
    double alphaPlus;
    double alphaMinus;
    double dPlus;
    double dMinus;
    Eigen::VectorXd uPlus;
    Eigen::VectorXd uMinus;
    Eigen::VectorXd wPlus;
    Eigen::VectorXd wMinus;
    
    // Initialize the HS field at all time slices in [0, beta] with +1 and -1 randomly
    h = generateHSmatrix(L, NSITES); // HS field h_{l = 1,... L, i = 1, ..., NSITES}
    
    // Hopping matrix 1D chain w/ PBC
    K = createHoppingMatrix(NSITES);
    
    // Compute matrix 'prefactor' of the B-matrices e^{t\delta\tau K}
    const static Eigen::MatrixXd B_preFactor = (t * dt * K).exp();
  
    // Some prints for debugging. Toggle them if you want at the top.
    if (printsOn == true)
    {
        std::cout << "\n K = \n\n" << K << "\n\n";
        std::cout << "exp (t * dt * K) = \n\n" << B_preFactor << "\n\n";
        std::cout << "Hubbard Stratonovich Field\n\n h = \n\n" << h << "\n\n";
    }
    
    // Build the B-matrices. We need a copy of each to perform the updates (naively)
    for (l = 0; l < L; l++)
    {
        BpOld[l] = build_Bmatrix(true, nu, NSITES, h.row(l), B_preFactor);      //  true for spin up
        BmOld[l] = build_Bmatrix(false, nu, NSITES, h.row(l), B_preFactor);     //  false for spin down
        BpNew[l] = BpOld[l];
        BmNew[l] = BmOld[l];
    }
    
    //  Build the M-matrices and Green's function (matrix)
    for (l = 0; l < L; l++)
    {
        MPlusOld *= BpOld[L - 1 - l];
        MMinusOld *= BmOld[L - 1 - l];
    }
    MPlusOld += Eigen::MatrixXd::Identity(NSITES,NSITES);
    MMinusOld += Eigen::MatrixXd::Identity(NSITES,NSITES);

    GreenPlus = MPlusOld.inverse();
    GreenMinus = MMinusOld.inverse();
    
    //  Initialize determinants and acceptance ratio
    double detOldPlus = MPlusOld.determinant();
    double detOldMinus = MMinusOld.determinant();
    double detsProdOld = detOldPlus * detOldMinus;
    double detsProdNew = detsProdOld;
    
    //  Initialize arrays to store measurements
    Eigen::VectorXd weightsNaive(totalMCSteps);
    Eigen::VectorXd weightsUpdate(totalMCSteps);
    weightsUpdate(0) = detsProdNew;
    
    //  Inititialize entry of HS field matrix to (0, 0)
    l_chosen = 0;
    i_chosen = 0;
    
    //  Initialize number of measurements
    m = 0;
    
    
    // --- MC LOOP ---
    

    for (step = 0; step < totalMCSteps; step++)
    {
        //  To check progress of the run
        if (printsOn == true)
        {
            if ((step+1) % (totalMCSteps/10) == 0)
            {
                std::cout << "\nMC loop: " << (step + 1)*1. / totalMCSteps * 100 << " %" << std::endl;
            }
        }
        
        //  For the updates
        uPlus = uSigma(NSITES, GreenPlus, i_chosen);
        uMinus = uSigma(NSITES, GreenMinus, i_chosen);
        wPlus = wSigma(NSITES, GreenPlus, i_chosen);
        wMinus = wSigma(NSITES, GreenMinus, i_chosen);
        
        //  Save weight of the configuration to check convergence
        weightsNaive(step) = detsProdNew;
        
        //  Print current Green's matrices
        if (printsOn == true)
        {
            std::cout << "\n\nGreenPlus\n\n" << GreenPlus << std::endl;
            std::cout << "\n\nGreenMinus\n\n" << GreenMinus << std::endl;
        }
        
        //  Flip a 'spin'
        h(l_chosen, i_chosen) *= -1;
        
        //  Rebuild B-matrices
        BpNew[l_chosen] = build_Bmatrix(true, nu, NSITES, h.row(l_chosen), B_preFactor);
        BmNew[l_chosen] = build_Bmatrix(false, nu, NSITES, h.row(l_chosen), B_preFactor);
        
        //  Rebuild M-matrices
        MPlusNew = Eigen::MatrixXd::Identity(NSITES,NSITES);
        MMinusNew = Eigen::MatrixXd::Identity(NSITES,NSITES);
        for (l = 0; l < L; l++)
        {
            MPlusNew *= BpNew[L - 1 - l];
            MMinusNew *= BmNew[L - 1 - l];
        }
        MPlusNew += Eigen::MatrixXd::Identity(NSITES,NSITES);
        MMinusNew += Eigen::MatrixXd::Identity(NSITES,NSITES);

        //  Compute the acceptance ratio (brute force)
        detNewPlus = MPlusNew.determinant();
        detNewMinus = MMinusNew.determinant();
        detsProdNew = detNewPlus * detNewMinus;
        acceptanceRatio = detsProdNew / detsProdOld;
        if (printsOn == true)
        {
            std::cout << "\n\nacceptance ratio: " << acceptanceRatio << " (naive)" << std::endl;
        }
        
        //  Compute the acceptance ratio (via updates)
        alphaPlus = ( exp( -2 * (-1) * h(l_chosen,i_chosen) * nu ) - 1 );
        alphaMinus = ( exp( 2 * (-1) * h(l_chosen,i_chosen) * nu ) - 1 );
        dPlus = ( 1 + alphaPlus  * ( 1 - GreenPlus(i_chosen, i_chosen) ) );
        dMinus = ( 1 + alphaMinus  * ( 1 - GreenMinus(i_chosen, i_chosen) ) );
        r = dPlus * dMinus;
        if (printsOn == true)
        {
            std::cout << "\n\nacceptance ratio: " << r << " (via Green's function update)" << std::endl;
        }

        //  Draw random number to decide whether or not to accept the move
        decisionMaker = dis(gen);
        
        if (decisionMaker <= acceptanceRatio)
        {
            //  Save weight of configuration
            if (step < totalMCSteps - 1)
            {
                weightsUpdate(step + 1) = r * weightsUpdate(step) ;
            }
            //update everything
            BpOld[l_chosen] = BpNew[l_chosen];
            BmOld[l_chosen] = BmNew[l_chosen];
            MPlusOld = MPlusNew;
            MMinusOld = MMinusNew;
            detsProdOld = detsProdNew;

            // Rank-one update --> O ( 2 * N^2 )
            GreenPlus -= alphaPlus/( 1 + alphaPlus  * ( 1 - GreenPlus(i_chosen, i_chosen) ) ) * kroneckerProduct( uPlus , wPlus.transpose() ).eval();
            GreenMinus -= alphaMinus/( 1 + alphaMinus  * ( 1 - GreenMinus(i_chosen, i_chosen) ) ) * kroneckerProduct( uMinus , wMinus.transpose() ).eval();

            if (printsOn == true)
            {
                std::cout << "\n\nO(N^2) update of the Green's function" << std::endl;
                std::cout << "\n\nGreenPlus\n\n" << GreenPlus << std::endl;
                std::cout << "\n\nGreenMinus\n\n" << GreenMinus << std::endl;
                std::cout << "\n\nBrute force update\n\n" << std::endl;
                std::cout << MPlusNew.inverse() << std::endl << std::endl;
                std::cout << MMinusNew.inverse() << std::endl << std::endl;
            }
            
        }
        else
        {
            //  Save weight of configuration
            if (step < totalMCSteps - 1)
            {
                weightsUpdate(step + 1) = weightsUpdate(step) ;
            }
            
            // revert changes (switch new -> old wrt the above statements)
            BpNew[l_chosen] = BpOld[l_chosen];
            BmNew[l_chosen] = BmOld[l_chosen];
            MPlusNew = MPlusOld;
            MMinusNew = MMinusOld;
            h(l_chosen, i_chosen) *= -1;
            detsProdNew = detsProdOld;
        }
        
        //  loop through (l, i)
        if (i_chosen < NSITES - 1)
        {
            i_chosen += 1;
            //l_chosen = l_chosen;
        }
        else
        {
            if (l_chosen < L - 1)
            {
                // Wrapping
                GreenPlus = BpNew[l_chosen] * GreenPlus * BpNew[l_chosen].inverse();
                GreenMinus = BmNew[l_chosen] * GreenMinus * BmNew[l_chosen].inverse();

                l_chosen += 1;
                i_chosen = 0;
            }
            else
            {
                // Back to original order
                GreenPlus = BpNew[l_chosen] * GreenPlus * BpNew[l_chosen].inverse();
                GreenMinus = BmNew[l_chosen] * GreenMinus * BmNew[l_chosen].inverse();
                
                l_chosen = 0;
                i_chosen = 0;
            }
        }
        
        
    }
    
    // Save weights of accepted configurations to file

    std::ofstream file1("plots/weightsNaive.txt");
    if (file1.is_open())
    {
        file1 << weightsNaive << '\n';
    }
    
    std::ofstream file2("plots/weightsUpdate.txt");
    if (file2.is_open())
    {
        file2 << weightsUpdate << '\n';
    }

    return 0;
}

