//
//  hubbard-1d-finiteT.cpp -> main!
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
#include "prints.h"

// --- DEFINITIONS ---

#define NSITES 50

//  Physical parameters
double dt;                                   //  time subinterval width. error scales as dt^2
double beta;                                 //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
int L;                                       //  number of imaginary time subintervals
double t;                                    //  hopping parameters
double U;                                    //  interaction energy
double nu;                                   //  Hubbard Stratonovich transformation parameter
double mu;                                   //  chemical potential

//  Some general use variables
int i;
int j;
int l;
int m;
int step;
int sweep = 0;
int i_chosen;
int l_chosen;

//  Toggle prints for debugging
bool printsOn = false;                                      // false - NO PRINTS, true - PRINTS

//  Things to do with random numbers throughout the code
const static int seed = 12345;                              //  set a seed
std::mt19937 gen(seed);                                     //  mt19937 algorithm to generate random numbers
std::uniform_real_distribution<> dis(0.0, 1.0);             //  set the distribution from which to draw numbers
double decisionMaker;                                       //  to accept or not to accept, hence the question.


// --- MAIN ---


int main()
{
    //  Set Monte Carlo-specific variables
    const static int totalMCSweeps = 10;
    
    //  How often to calculate Green functions afresh
    int freq = 1;
    
    //  Set physical parameters and Trotter parameter (error scales as dt^2)
    dt = 0.1;                                                   //  time subinterval width. error scales as dt^2
    beta = 1.;                                                  //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
    U = 2.;                                                     //  interaction energy
    t = 1.;                                                     //  hopping parameter (default to 1, and measure energy in units of hopping parameter)
    mu = 0;                                                     //  chemical potential
    
    //  Variables that depend on the parameters
    L = beta/dt;                                                //  number of imaginary time subintervals
    nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12 ;       //  Hubbard Stratonovich transformation parameter : nu = acosh( exp( U * dt / 2 ) );
    const static int totalMCSteps = totalMCSweeps * NSITES * L;
    
    printWelcome(NSITES);
    printParameters();
    printMC(totalMCSteps, NSITES, L);
    
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
    Eigen::VectorXd uPlus;
    Eigen::VectorXd uMinus;
    Eigen::VectorXd wPlus;
    Eigen::VectorXd wMinus;
    double detNewPlus;
    double detNewMinus;
    double acceptanceRatio;
    double r;
    double alphaPlus;
    double alphaMinus;
    double dPlus;
    double dMinus;

    // Initialize the HS field at all time slices in [0, beta] with +1 and -1 randomly
    h = generateHSmatrix(L, NSITES); // HS field h_{l = 1,... L, i = 1, ..., NSITES}
    
    // Hopping matrix 1D chain w/ PBC
    K = createHoppingMatrix(NSITES);
    
    // Compute matrix 'prefactor' of the B-matrices e^{t\delta\tau K} * e^{d\tau \mu}
    const static Eigen::MatrixXd B_preFactor = (t * dt * K).exp() * exp(dt * mu);
  
    // Some prints for debugging. Toggle them if you want at the top.
    printStartingMatrices(K, B_preFactor, h, printsOn);

    // Build the B-matrices. We need a copy of each to perform the updates naively
    for (l = 0; l < L; l++)
    {
        BpOld[l] = build_Bmatrix(true, nu, NSITES, h.row(l), B_preFactor);      //  true for spin up
        BmOld[l] = build_Bmatrix(false, nu, NSITES, h.row(l), B_preFactor);     //  false for spin down
        BpNew[l] = BpOld[l];
        BmNew[l] = BmOld[l];
    }

    //  Build the M-matrices and Green's matrix
    for (l = 0; l < L; l++)
    {
        MPlusOld *= BpOld[L - 1 - l];
        MMinusOld *= BmOld[L - 1 - l];
    }
    MPlusOld += Eigen::MatrixXd::Identity(NSITES,NSITES);
    MMinusOld += Eigen::MatrixXd::Identity(NSITES,NSITES);

    //  Initialize the inverses
    GreenPlus = MPlusOld.inverse();
    GreenMinus = MMinusOld.inverse();

    //  Initialize determinants
    double detOldPlus = MPlusOld.determinant();
    double detOldMinus = MMinusOld.determinant();
    double detsProdOld = detOldPlus * detOldMinus;
    double detsProdNew = detsProdOld;

    //  Initialize arrays to store measurements
    Eigen::VectorXd weightsNaive(totalMCSteps);
    Eigen::VectorXd weightsUpdate(totalMCSteps);
    weightsUpdate(0) = detsProdNew;

    //  Inititialize entry of HS field matrix to (0, 0).
    //  For each imaginary time slice l, we loop over all i's, then we change slice and do the same. And so on.
    l_chosen = 0;
    i_chosen = 0;


    // --- MC LOOP ---

    std::cout << "\n\nMC loop started. Please wait.\n\n";

    for (step = 0; step < totalMCSteps; step++)
    {
        if ( (step + 1)  % (totalMCSteps/10) == 0 )
        {
            std::cout << "\nProgress: " << (step + 1)*1. / totalMCSteps * 100 << " %" << std::endl;
        }

        //  Vectors that are used to construct the rank-one update
        uPlus = uSigma(NSITES, GreenPlus, i_chosen);
        uMinus = uSigma(NSITES, GreenMinus, i_chosen);
        wPlus = wSigma(NSITES, GreenPlus, i_chosen);
        wMinus = wSigma(NSITES, GreenMinus, i_chosen);

        //  Save weight of the configuration to check convergence
        weightsNaive(step) = detsProdNew;

        //  Print current Green's matrices (obtained by updates)
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
            //  Save weight of configuration (obtained by updates)
            if (step < totalMCSteps - 1)
            {
                weightsUpdate(step + 1) = r * weightsUpdate(step) ;
            }
            //  Update everything
            BpOld[l_chosen] = BpNew[l_chosen];
            BmOld[l_chosen] = BmNew[l_chosen];
            MPlusOld = MPlusNew;
            MMinusOld = MMinusNew;
            detsProdOld = detsProdNew;

            //  Rank-one update --> O ( 2 * N^2 ). kroneckerProduct is the tensor product of two vectors. .eval() is needed because of an Eigen quirk. Check documentation.
            GreenPlus -= alphaPlus/( 1 + alphaPlus  * ( 1 - GreenPlus(i_chosen, i_chosen) ) ) * kroneckerProduct( uPlus , wPlus.transpose() ).eval();
            GreenMinus -= alphaMinus/( 1 + alphaMinus  * ( 1 - GreenMinus(i_chosen, i_chosen) ) ) * kroneckerProduct( uMinus , wMinus.transpose() ).eval();

            //  Compare the updated Green matrices and the ones computed by brute force. They should match
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
            //  Save weight of configuration (obtained by updates)
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
            i_chosen += 1;  // l_chosen = l_chosen (it doesn't change)
        }
        else
        {
            if (l_chosen < L - 1)
            {

                //  Wrapping
                GreenPlus = BpNew[l_chosen] * GreenPlus * BpNew[l_chosen].inverse();
                GreenMinus = BmNew[l_chosen] * GreenMinus * BmNew[l_chosen].inverse();

                //  New im-time slice
                l_chosen += 1;
                i_chosen = 0;
            }
            else
            {
                sweep += 1;

                if (sweep % freq == 0)
                {
                    //  Compute Green's functions afresh
                    GreenPlus = MPlusNew.inverse();
                    GreenMinus = MMinusNew.inverse();
                }
                else
                {
                    //  This particular wrapping takes us back to the original order
                    GreenPlus = BpNew[l_chosen] * GreenPlus * BpNew[l_chosen].inverse();
                    GreenMinus = BmNew[l_chosen] * GreenMinus * BmNew[l_chosen].inverse();
                }


                //  New sweep
                l_chosen = 0;
                i_chosen = 0;
            }
        }



    }

    //  Save parameters of the simulation to file

    std::ofstream file("plots/simulationParameters.txt");
    if (file.is_open())
    {
        file << NSITES << '\n';
        file << dt << '\n';
        file << beta << '\n';
        file << L << '\n';
        file << t << '\n';
        file << U << '\n';
        file << mu << '\n';
        file << totalMCSteps*1./(NSITES*L) << '\n'; //sweeps
        file << freq << '\n' ; // frequency of naive computation
    }

    //  Save weights of accepted configurations to file

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

