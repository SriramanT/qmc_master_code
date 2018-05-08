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
bool printsOn = true;                                      // false - NO PRINTS, true - PRINTS

int main()
{
    int l;
    int i;
    int j;

    
    //  Set physical parameters and Trotter parameter (error scales as dt^2)
    dt = 0.125;                                                   //  time subinterval width. error scales as dt^2
    beta = 8.;                                                  //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
    U = 4.;                                                     //  interaction energy
    t = 1.;                                                     //  hopping parameter (default to 1, and measure energy in units of hopping parameter)
    mu = 0.;                                                   //  chemical potential
    
    //  Variables that depend on the parameters
    L = beta/dt;                                                //  number of imaginary time subintervals
    int k = L / 10 ;                                            // number of intervals
    nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12 ;       //  Hubbard Stratonovich transformation parameter : nu = acosh( exp( U * dt / 2 ) );
    
    // --- INITIALIZATION ---
    Eigen::MatrixXd h;
    Eigen::MatrixXd K;
    Eigen::MatrixXd BPlus[L]; // B plus (spin-up)
    Eigen::MatrixXd BMinus[L]; // B minus (spin-down)
    Eigen::MatrixXd MPlus = Eigen::MatrixXd::Identity(NSITES,NSITES);
    Eigen::MatrixXd MPlusLeft = Eigen::MatrixXd::Identity(NSITES,NSITES);
    Eigen::MatrixXd MPlusRight = Eigen::MatrixXd::Identity(NSITES,NSITES);
    Eigen::MatrixXd GreenPlus;
    Eigen::MatrixXd GreenPlusSmart;
    Eigen::MatrixXd GreenMinus;
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
    printStartingMatrices(K, B_preFactor, h, printsOn);
    
    // Build the B-matrices. We need a copy of each to perform the updates naively
    for (l = 0; l < L; l++)
    {
        BPlus[l] = build_Bmatrix(true, nu, NSITES, h.row(l), B_preFactor, dt, mu);      //  true for spin up
        std::cout << "\nBPlus_" << l << std::endl << BPlus[l] << std::endl;
    }
    
    //  Build the M-matrices and Green's matrix
    for (l = 0; l < k; l++)
    {
        MPlusLeft *= BPlus[L - 1 - l];
    }
    
    for (l = k; l < L; l++)
    {
        MPlusRight *= BPlus[L - 1 - l];
    }
    
    MPlus = (MPlusLeft*MPlusRight);
    
    std::cout << "\nMPlus" << std::endl << MPlus << std::endl;
    
//    Eigen::MatrixXd Q(NSITES, NSITES), R(NSITES, NSITES);
//    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(NSITES, NSITES);
//    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(NSITES, NSITES);
//    Eigen::HouseholderQR<Eigen::MatrixXd> qrMPlus(MPlus);
//    Q = qrMPlus.householderQ();
//    R= qrMPlus.matrixQR().template triangularView<Eigen::Upper>();
//    D.diagonal() = R.diagonal();
//    for (i = 0; i < NSITES; i++)
//    {
//        for (j = i + 1; j < NSITES; j++) // unit upper triangular -> loop starts in i + 1
//        {
//            V(i, j) = R(i, j) / R(i, i);
//        }
//    }
//    std::cout << "\nThe unitary matrix Q is:\n" << Q<< "\n\n";
//    std::cout << "\nThe diagonal matrix D is:\n" << D << "\n\n";
//    std::cout << "\nThe unit upper triangular matrix V is:\n" << V << "\n\n";
//    std::cout << "Q * R - MPlus \n" << Q * D * V - MPlus  << "\n\n";
//
//    Eigen::MatrixXd toDecompose = Q.inverse() * V.inverse() + D;
//
//    std::cout << "toDecompose\n" << toDecompose << "\n\n";
//
//    Eigen::MatrixXd Q1(NSITES, NSITES), R1(NSITES, NSITES);
//    Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(NSITES, NSITES);
//    Eigen::MatrixXd V1 = Eigen::MatrixXd::Identity(NSITES, NSITES);
//    Eigen::HouseholderQR<Eigen::MatrixXd> qr_toDecompose(toDecompose);
//    Q1 = qr_toDecompose.householderQ();
//    R1= qr_toDecompose.matrixQR().template triangularView<Eigen::Upper>();
//    D1.diagonal() = R1.diagonal();
//    for (i = 0; i < NSITES; i++)
//    {
//        for (j = i + 1; j < NSITES; j++) // unit upper triangular -> loop starts in i + 1
//        {
//            V1(i, j) = R1(i, j) / R1(i, i);
//        }
//    }
//    std::cout << "\nThe unitary matrix Q1 is:\n" << Q1<< "\n\n";
//    std::cout << "\nThe diagonal matrix D1 is:\n" << D1 << "\n\n";
//    std::cout << "\nThe unit upper triangular matrix V is:\n" << V1 << "\n\n";
//    std::cout << "Q * R - toDecompose \n" << Q1 * D1 * V1 - toDecompose << "\n\n";
    
    Eigen::MatrixXd Q(NSITES, NSITES), R(NSITES, NSITES), P(NSITES, NSITES);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(NSITES, NSITES);
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(NSITES, NSITES);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrMPlus(MPlus);
    Q = qrMPlus.householderQ();
//    R = qrMPlus.matrixQR().template triangularView<Eigen::Upper>();
    R = qrMPlus.matrixR().template triangularView<Eigen::Upper>();
    P = qrMPlus.colsPermutation();
    D.diagonal() = R.diagonal();
    for (i = 0; i < NSITES; i++)
    {
        for (j = i + 1; j < NSITES; j++) // unit upper triangular -> loop starts in i + 1
        {
            V(i, j) = R(i, j) / R(i, i);
        }
    }
    std::cout << "\nThe unitary matrix Q is:\n" << Q<< "\n\n";
    std::cout << "\nThe diagonal matrix D is:\n" << D << "\n\n";
    std::cout << "\nThe unit upper triangular matrix V is:\n" << V << "\n\n";
    std::cout << "Q * R - MPlus \n" << Q * D * V * P.transpose() - MPlus  << "\n\n";

    Eigen::MatrixXd toDecompose = Q.inverse() * V.inverse() + D;

    std::cout << "toDecompose\n" << toDecompose << "\n\n";

    Eigen::MatrixXd Q1(NSITES, NSITES), R1(NSITES, NSITES), P1(NSITES, NSITES);
    Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(NSITES, NSITES);
    Eigen::MatrixXd V1 = Eigen::MatrixXd::Identity(NSITES, NSITES);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_toDecompose(toDecompose);
    Q1 = qr_toDecompose.householderQ();
//    R1 = qr_toDecompose.matrixQR().template triangularView<Eigen::Upper>();
    R1 = qr_toDecompose.matrixR().template triangularView<Eigen::Upper>();
    P1 = qr_toDecompose.colsPermutation();
    D1.diagonal() = R1.diagonal();
    for (i = 0; i < NSITES; i++)
    {
        for (j = i + 1; j < NSITES; j++) // unit upper triangular -> loop starts in i + 1
        {
            V1(i, j) = R1(i, j) / R1(i, i);
        }
    }
    std::cout << "\nThe unitary matrix Q1 is:\n" << Q1<< "\n\n";
    std::cout << "\nThe diagonal matrix D1 is:\n" << D1 << "\n\n";
    std::cout << "\nThe unit upper triangular matrix V is:\n" << V1 << "\n\n";
    std::cout << "Q * R - toDecompose \n" << Q1 * D1 * V1 * P1.transpose() - toDecompose << "\n\n";
    
//    std::cout << "\nMPlusLeft\n" << std::endl << MPlusLeft << "\n\n";
//
//    Eigen::MatrixXd QLeft(NSITES, NSITES), RLeft(NSITES, NSITES);
//    Eigen::MatrixXd DLeft = Eigen::MatrixXd::Zero(NSITES, NSITES);
//    Eigen::MatrixXd VLeft = Eigen::MatrixXd::Identity(NSITES, NSITES);
//    Eigen::HouseholderQR<Eigen::MatrixXd> qrMPlusLeft(MPlusLeft);
//    QLeft = qrMPlusLeft.householderQ();
//    RLeft = qrMPlusLeft.matrixQR().template triangularView<Eigen::Upper>();
//    DLeft.diagonal() = RLeft.diagonal();
//    for (i = 0; i < NSITES; i++)
//    {
//        for (j = i + 1; j < NSITES; j++) // unit upper triangular -> loop starts in i + 1
//        {
//            VLeft(i, j) = RLeft(i, j) / RLeft(i, i);
//        }
//    }
//    std::cout << "\nThe unitary matrix QLeft is:\n" << QLeft << "\n\n";
//    std::cout << "\nThe diagonal matrix DLeft is:\n" << DLeft << "\n\n";
//    std::cout << "\nThe unit upper triangular matrix VLeft is:\n" << VLeft << "\n\n";
//    std::cout << "Q * R - MPlusLeft \n" << QLeft * DLeft * VLeft - MPlusLeft  << "\n\n";
//
//    std::cout << "\nMPlusRight\n" << std::endl << MPlusLeft << "\n\n";
//
//    Eigen::MatrixXd QRight(NSITES, NSITES), RRight(NSITES, NSITES);
//    Eigen::MatrixXd DRight = Eigen::MatrixXd::Zero(NSITES, NSITES);
//    Eigen::MatrixXd VRight = Eigen::MatrixXd::Identity(NSITES, NSITES);
//    Eigen::HouseholderQR<Eigen::MatrixXd> qrMPlusRight(MPlusRight);
//    QRight = qrMPlusRight.householderQ();
//    RRight= qrMPlusRight.matrixQR().template triangularView<Eigen::Upper>();
//    DRight.diagonal() = RRight.diagonal();
//    for (i = 0; i < NSITES; i++)
//    {
//        for (j = i + 1; j < NSITES; j++) // unit upper triangular -> loop starts in i + 1
//        {
//            VRight(i, j) = RRight(i, j) / RRight(i, i);
//        }
//    }
//
//    std::cout << "\nThe unitary matrix Q is:\n" << QRight << "\n\n";
//    std::cout << "\nThe diagonal matrix D is:\n" << DRight << "\n\n";
//    std::cout << "\nThe unit upper triangular matrix R is:\n" << VRight << "\n\n";
//    std::cout << "Q * R - MPlusRight \n" << QRight * DRight * VRight - MPlusRight  << "\n\n";
//
    //  Initialize the inverses
    
    MPlus += Eigen::MatrixXd::Identity(NSITES,NSITES);
    
    GreenPlus = MPlus.inverse();
    
    GreenPlusSmart = ( V1 * P1.transpose() * V * P.transpose() ).inverse() * D1.inverse() * ( Q * Q1 ).inverse();
    
//    GreenPlusSmart = QLeft.inverse() * ( QLeft.inverse() * QRight.inverse() + DLeft * VLeft * VRight * DRight ).inverse() * QRight.inverse();
    
    std::cout << "\nGreenPlus" << std::endl << GreenPlus << std::endl;
    
    std::cout << "\nGreenPlusSmart" << std::endl << GreenPlusSmart << std::endl;
    
    std::cout << "\nGreenPlus * MPlus\n" << GreenPlus * MPlus << std::endl;

    std::cout << "\nGreenPlusSmart * MPlus\n" << GreenPlusSmart * MPlus << std::endl;
    
    return 0;
}
