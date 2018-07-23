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
int p;

class UDV
{
    Eigen::MatrixXd m;
    Eigen::MatrixXd R;
public:
    UDV(Eigen::MatrixXd toDecompose);
    void printMatrixToDecompose();
    Eigen::MatrixXd QR_and_getU();
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getV();
};

UDV::UDV(Eigen::MatrixXd toDecompose)
{
    m = toDecompose;
}

void UDV::printMatrixToDecompose()
{
    std::cout << "\nMatrixToDecompose\n" << m << std::endl;
}

Eigen::MatrixXd UDV::QR_and_getU()
{
    Eigen::HouseholderQR<Eigen::MatrixXd> qrPartial(m);
    R = qrPartial.matrixQR().template triangularView<Eigen::Upper>();
    
    return qrPartial.householderQ();
}

Eigen::MatrixXd UDV::getD()
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(NSITES, NSITES);
    D.diagonal() = R.diagonal();
    return D;
}

Eigen::MatrixXd UDV::getV()
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(NSITES, NSITES);
    for (i = 0; i < NSITES; i++)
    {
        for (p = i + 1; p < NSITES; p++) // unit upper triangular -> loop starts in i + 1
        {
            V(i, p) = R(i, p) / R(i, i);
        }
    }
    return V;
}

class UDVpartialPivoting
{
    Eigen::MatrixXd m;
    Eigen::MatrixXd R;
    Eigen::MatrixXd P;
public:
    UDVpartialPivoting(Eigen::MatrixXd toDecompose);
    void printMatrixToDecompose();
    Eigen::MatrixXd QR_and_getU();
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getV();
    Eigen::MatrixXd getP();
};

UDVpartialPivoting::UDVpartialPivoting(Eigen::MatrixXd toDecompose)
{
    m = toDecompose;
}

void UDVpartialPivoting::printMatrixToDecompose()
{
    std::cout << "\nMatrixToDecompose\n" << m << std::endl;
}

Eigen::MatrixXd UDVpartialPivoting::QR_and_getU()
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrPartial(m);
    R = qrPartial.matrixQR().template triangularView<Eigen::Upper>();
    P = qrPartial.colsPermutation();
    return qrPartial.householderQ();
}

Eigen::MatrixXd UDVpartialPivoting::getD()
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(NSITES, NSITES);
    D.diagonal() = R.diagonal();
    return D;
}

Eigen::MatrixXd UDVpartialPivoting::getV()
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(NSITES, NSITES);
    for (i = 0; i < NSITES; i++)
    {
        for (p = i + 1; p < NSITES; p++) // unit upper triangular -> loop starts in i + 1
        {
            V(i, p) = R(i, p) / R(i, i);
        }
    }
    return V;
}

Eigen::MatrixXd UDVpartialPivoting::getP()
{
    return P;
}



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
    Eigen::MatrixXd MPlus = Eigen::MatrixXd::Identity(NSITES,NSITES);
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
    
    // Build the B-matrices. We need a copy of each to perform the updates naively
    for (l = 0; l < L; l++)
    {
        BPlus[l] = build_Bmatrix(true, nu, NSITES, h.row(l), B_preFactor, dt, mu);      //  true for spin up
//        std::cout << "\nBPlus_" << l << std::endl << BPlus[l] << std::endl;
    }
    
    Eigen::MatrixXd partialProdBs;
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity(NSITES, NSITES);
    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(NSITES, NSITES);
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(NSITES, NSITES);
    Eigen::MatrixXd U_ = Eigen::MatrixXd::Identity(NSITES, NSITES);
    Eigen::MatrixXd D_ = Eigen::MatrixXd::Identity(NSITES, NSITES);
    Eigen::MatrixXd V_ = Eigen::MatrixXd::Identity(NSITES, NSITES);
    
    //  Build the M-matrices and Green's matrix
    
    for (j = 0; j < k; j++)
    {
        partialProdBs = Eigen::MatrixXd::Identity(NSITES,NSITES);
        for (l = 0; l < L/k; l++)
        {
            partialProdBs = BPlus[j * (L/k) + l] * partialProdBs;
        }
        
        //  Householder
        
        UDV m1(partialProdBs * U * D);
        m1.printMatrixToDecompose();
        U = m1.QR_and_getU();
        std::cout << "\nU'\n" << U << std::endl;
        D = m1.getD();
        std::cout << "\nD'\n" << D << std::endl;
        V = m1.getV() * V;
        std::cout << "\nV'\n" << V << std::endl;
        
//        //  Column pivoting
//
//        UDVpartialPivoting m2(partialProdBs * U_ * D_);
//        m2.printMatrixToDecompose();
//        U_ = m2.QR_and_getU();
//        std::cout << "\nU'\n" << U_ << std::endl;
//        D_ = m2.getD();
//        std::cout << "\nD'\n" << D_ << std::endl;
//        V_ = m2.getV() * V_;
//        std::cout << "\nV'\n" << V_ << std::endl;

        MPlus = partialProdBs * MPlus;
    }
    
    std::cout << "Compute Green" << std::endl;

    UDV middleMatrix(U.inverse() * V.inverse() + D);
    middleMatrix.printMatrixToDecompose();
    U = U * middleMatrix.QR_and_getU();
    std::cout << "\nU U'\n" << U << std::endl;
    D = middleMatrix.getD();
    //  Compute inverse of diagonal
    for (i = 0; i < NSITES; i++)
    {
        D(i, i) = 1 / D(i, i);
    }
    std::cout << "\nD'^-1\n" << D << std::endl;
    V = middleMatrix.getV() * V;
    std::cout << "\nV' V\n" << V << std::endl;
    
    

    GreenPlusHouseholder = V.inverse() * D * U.inverse();

    MPlus += Eigen::MatrixXd::Identity(NSITES,NSITES);

    GreenPlus = MPlus.inverse();

    std::cout << "\nMPlus\n" << MPlus << std::endl;
    std::cout << "\nGreenPlus\n" << GreenPlus << std::endl;
    std::cout << "\nGreenPlusSmart\n" << GreenPlusHouseholder << std::endl;
    std::cout << "\nDifference\n" << GreenPlus - GreenPlusHouseholder << std::endl;
    
    return 0;
}
