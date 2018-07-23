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

#define NSITES 4

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
    const static int totalMCSweeps = 1;                   //  same as the number of measurements,
                                                            //  i.e. we make a measurement every sweep, then find correlations in post-process
    
    //  How often to calculate Green functions afresh
    int freq = 1;
    
    //  Set physical parameters and Trotter parameter (error scales as dt^2)
    dt = 0.125;                                                   //  time subinterval width. error scales as dt^2
    beta = 1.;                                                  //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
    U = 4.;                                                     //  interaction energy
    t = 1.;                                                     //  hopping parameter (default to 1, and measure energy in units of hopping parameter)
    mu = 0.;                                                   //  chemical potential
    
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
    Eigen::MatrixXd BPlus[L]; // B plus (spin-up)
    Eigen::MatrixXd BMinus[L]; // B minus (spin-down)
    Eigen::MatrixXd MPlus = Eigen::MatrixXd::Identity(NSITES,NSITES);
    Eigen::MatrixXd MMinus = Eigen::MatrixXd::Identity(NSITES,NSITES);
    Eigen::MatrixXd GreenPlus;
    Eigen::MatrixXd GreenMinus;
    Eigen::VectorXd uPlus;
    Eigen::VectorXd uMinus;
    Eigen::VectorXd wPlus;
    Eigen::VectorXd wMinus;
    double r;
    double alphaPlus;
    double alphaMinus;
    double dPlus;
    double dMinus;
    Eigen::VectorXd electronDensity = Eigen::VectorXd::Zero(totalMCSweeps);
    
    // Initialize the HS field at all time slices in [0, beta] with +1 and -1 randomly
    h = generateHSmatrix(L, NSITES); // HS field h_{l = 1,... L, i = 1, ..., NSITES}
    
    // Hopping matrix 1D chain w/ PBC
    K = createHoppingMatrix(NSITES);
    
    // Compute matrix 'prefactor' of the B-matrices e^{t\delta\tau K}
    const static Eigen::MatrixXd B_preFactor = (t * dt * K).exp();
  
    // Some prints for debugging. Toggle them if you want at the top.
    printStartingMatrices(K, B_preFactor, h, printsOn);

    // Build the B-matrices. We need a copy of each to perform the updates naively
    for (l = 0; l < L; l++)
    {
        BPlus[l] = build_Bmatrix(true, nu, NSITES, h.row(l), B_preFactor, dt, mu);      //  true for spin up
        BMinus[l] = build_Bmatrix(false, nu, NSITES, h.row(l), B_preFactor, dt, mu);     //  false for spin down
    }

    //  Build the M-matrices and Green's matrix
    for (l = 0; l < L; l++)
    {
        MPlus *= BPlus[L - 1 - l];
        MMinus *= BMinus[L - 1 - l];
    }
    MPlus += Eigen::MatrixXd::Identity(NSITES,NSITES);
    MMinus += Eigen::MatrixXd::Identity(NSITES,NSITES);

    //  Initialize the inverses
    GreenPlus = MPlus.inverse();
    GreenMinus = MMinus.inverse();

    //  Initialize determinants
    double detPlus = MPlus.determinant();
    double detMinus = MMinus.determinant();
    double detsProd = detPlus * detMinus;
    
    std::cout << GreenPlus << std::endl;

    //  Initialize arrays to store measurements
    Eigen::VectorXd weights(totalMCSteps);
    weights(0) = detsProd;
    std::cout << weights(0) << std::endl;

    //  Inititialize entry of HS field matrix to (0, 0).
    //  For each imaginary time slice l, we loop over all i's, then we change slice and do the same. And so on.
    l_chosen = 0;
    i_chosen = 0;
//
//
//    // --- MC LOOP ---
//
//    std::cout << "\n\nMC loop started. Please wait.\n\n";
//
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
//
//        //  Print current Green's matrices (obtained by updates)
//        if (printsOn == true)
//        {
//            std::cout << "\n\nGreenPlus\n\n" << GreenPlus << std::endl;
//            std::cout << "\n\nGreenMinus\n\n" << GreenMinus << std::endl;
//        }
//
//        //  Compute the acceptance ratio
        alphaPlus = ( exp( -2 * h(l_chosen,i_chosen) * nu ) - 1 );
        alphaMinus = ( exp( 2 * h(l_chosen,i_chosen) * nu ) - 1 );
        dPlus = ( 1 + alphaPlus  * ( 1 - GreenPlus(i_chosen, i_chosen) ) );
        dMinus = ( 1 + alphaMinus  * ( 1 - GreenMinus(i_chosen, i_chosen) ) );
        r = dPlus * dMinus;
//        if (printsOn == true)
//        {
            std::cout << "\n\nacceptance ratio: " << r << " (via Green's function update)" << std::endl;
//        }
//
        //  Draw random number to decide whether or not to accept the move
        decisionMaker = dis(gen);

        if (decisionMaker <= r)
        {
            //  Save weight of configuration (obtained by updates)
            if (step < totalMCSteps - 1)
            {
                weights(step + 1) = r * weights(step) ;
            }
            //  Flip a 'spin'
            h(l_chosen, i_chosen) *= -1;

            //  Rank-one update --> O ( 2 * N^2 ). kroneckerProduct is the tensor product of two vectors. .eval() is needed because of an Eigen quirk. Check documentation.
            GreenPlus -= alphaPlus/( 1 + alphaPlus  * ( 1 - GreenPlus(i_chosen, i_chosen) ) ) * kroneckerProduct( uPlus , wPlus.transpose() ).eval();
            GreenMinus -= alphaMinus/( 1 + alphaMinus  * ( 1 - GreenMinus(i_chosen, i_chosen) ) ) * kroneckerProduct( uMinus , wMinus.transpose() ).eval();

        }
        else
        {
            //  Save weight of configuration (obtained by updates)
            if (step < totalMCSteps - 1)
            {
                weights(step + 1) = weights(step) ;
            }

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
                //  Rebuild B-matrices
                BPlus[l_chosen] = build_Bmatrix(true, nu, NSITES, h.row(l_chosen), B_preFactor, dt, mu);
                BMinus[l_chosen] = build_Bmatrix(false, nu, NSITES, h.row(l_chosen), B_preFactor, dt, mu);

                //  Wrapping
                GreenPlus = BPlus[l_chosen] * GreenPlus * BPlus[l_chosen].inverse();
                GreenMinus = BMinus[l_chosen] * GreenMinus * BMinus[l_chosen].inverse();

                //  New im-time slice
                l_chosen += 1;
                i_chosen = 0;
            }
            else
            {
//                //  Measurements
//                electronDensity(sweep) = 1 - ( GreenPlus.trace() + GreenMinus.trace() ) / 2 / NSITES;

                sweep += 1;

                //  Rebuild B-matrices
                BPlus[l_chosen] = build_Bmatrix(true, nu, NSITES, h.row(l_chosen), B_preFactor, dt, mu);
                BMinus[l_chosen] = build_Bmatrix(false, nu, NSITES, h.row(l_chosen), B_preFactor, dt, mu);

                if (sweep % freq == 0)
                {
                    //  Rebuild M-matrices
                    MPlus = Eigen::MatrixXd::Identity(NSITES,NSITES);
                    MMinus = Eigen::MatrixXd::Identity(NSITES,NSITES);
                    for (l = 0; l < L; l++)
                    {
                        MPlus *= BPlus[L - 1 - l];
                        MMinus *= BMinus[L - 1 - l];
                    }
                    MPlus += Eigen::MatrixXd::Identity(NSITES,NSITES);
                    MMinus += Eigen::MatrixXd::Identity(NSITES,NSITES);

                    //  Compute Green's functions afresh
                    GreenPlus = MPlus.inverse();
                    GreenMinus = MMinus.inverse();
                }
                else
                {
                    //  This particular wrapping takes us back to the original order
                    GreenPlus = BPlus[l_chosen] * GreenPlus * BPlus[l_chosen].inverse();
                    GreenMinus = BMinus[l_chosen] * GreenMinus * BMinus[l_chosen].inverse();
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

    std::ofstream file1("plots/weights.txt");
    if (file1.is_open())
    {
        file1 << weights << '\n';
    }

//    //  Save measurements
//
//    std::ofstream file2("plots/electronDensity.txt");
//    if (file2.is_open())
//    {
//        file2 << electronDensity << '\n';
//    }

    return 0;
}

