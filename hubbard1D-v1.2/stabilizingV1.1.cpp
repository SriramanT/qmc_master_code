//
//  stabilizing_experiments.cpp
//  
//
//  Created by Francisco Brito on 06/05/2018.
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
#include "UDV.h"

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

//  Toggle prints for debugging
bool printsOn = false;                                      // false - NO PRINTS, true - PRINTS

int l;
int i;
int j;

int main()
{
    //  Set physical parameters and Trotter parameter (error scales as dt^2)
    dt = 0.125;                                                   //  time subinterval width. error scales as dt^2
    beta = 10.;                                                  //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
    U = 4.;                                                     //  interaction energy
    t = 1.;                                                     //  hopping parameter (default to 1, and measure energy in units of hopping parameter)
    mu = 0.;                                                   //  chemical potential
    
    //  Variables that depend on the parameters
    L = beta/dt;                                                //  number of imaginary time subintervals
    std::cout << "\nL\n" << L << std::endl;
    int k = 40;                                                 //  k = # intervals. Must be commensurate with L and k < L
    std::cout << "\nk\n" << k << std::endl;
    nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12 ;       //  Hubbard Stratonovich transformation parameter : nu = acosh( exp( U * dt / 2 ) );
    
    // --- INITIALIZATION ---
    Eigen::MatrixXd h;
    Eigen::MatrixXd K;
    Eigen::MatrixXd BPlus[L]; // B plus (spin-up)
    Eigen::MatrixXd GreenPlus;
    Eigen::MatrixXd GreenPlusHouseholder;
    
    // Initialize the HS field at all time slices in [0, beta] with +1 and -1 randomly
    h = generateHSmatrix(L, NSITES); // HS field h_{l = 1,... L, i = 1, ..., NSITES}
    
    // Hopping matrix 1D chain w/ PBC
    K = createHoppingMatrix(NSITES);
    
    // Compute matrix 'prefactor' of the B-matrices e^{t\delta\tau K}
    const static Eigen::MatrixXd B_preFactor = (t * dt * K).exp();
    
    // Some prints for debugging. Toggle them if you want at the top.
    printStartingMatrices(K, B_preFactor, h, printsOn);
    
    // Build the B-matrices
    for (l = L - 1; l >= 0; l--)
    {
        BPlus[l] = build_Bmatrix(true, nu, NSITES, h.row(l), B_preFactor, dt, mu);      //  true for spin up
    }
    
    GreenPlus = computeGreenNaive(BPlus, L, NSITES);
    
    std::cout << "\nGreenPlus\n" << GreenPlus << std::endl;
    
    GreenPlusHouseholder = computeGreen(BPlus, L, k, NSITES);
    
    std::cout << "\nGreenPlusSmart\n" << GreenPlusHouseholder << std::endl;
    std::cout << "\nDifference\n" << GreenPlus - GreenPlusHouseholder << std::endl;
    
    return 0;
}
