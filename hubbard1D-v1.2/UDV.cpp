//
//  UDV.cpp
//  
//
//  Created by Francisco Brito on 08/05/2018.
//

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "UDV.h"

int a;
int b;
int c;
int d;

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
    R = qrPartial.matrixQR().triangularView<Eigen::Upper>();
    
    return qrPartial.householderQ();
}

Eigen::MatrixXd UDV::getD()
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(m.rows(), m.rows());
    D.diagonal() = R.diagonal();
    return D;
}

Eigen::MatrixXd UDV::getV()
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(m.rows(), m.rows());
    for (a = 0; a < m.rows(); a++)
    {
        for (b = a + 1; b < m.rows(); b++) // unit upper triangular -> loop starts in i + 1
        {
            V(a, b) = R(a, b) / R(a, a);
        }
    }
    return V;
}

Eigen::MatrixXd computeGreen(Eigen::MatrixXd Bs[], int L, int k, int nSites)
{
    Eigen::MatrixXd partialProdBs;
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity(nSites, nSites);
    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(nSites, nSites);
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(nSites, nSites);
    
    //  Compute partial products, UDV's and multiply
    for (c = 0; c < k; c++)
    {
        partialProdBs = Eigen::MatrixXd::Identity(nSites, nSites);
        for (d = 0; d < L/k; d++)
        {
            partialProdBs = Bs[c * (L/k) + d] * partialProdBs;
        }
        
        //  Householder QR
        
        UDV m1(partialProdBs * U * D);
        m1.printMatrixToDecompose();
        U = m1.QR_and_getU();
        std::cout << "\nU'\n" << U << std::endl;
        D = m1.getD();
        std::cout << "\nD'\n" << D << std::endl;
        V = m1.getV() * V;
        std::cout << "\nV'\n" << V << std::endl;
        
    }
    
    std::cout << "\nCompute Green\n" << std::endl;
    
    UDV middleMatrix(U.inverse() * V.inverse() + D);
    middleMatrix.printMatrixToDecompose();
    U = U * middleMatrix.QR_and_getU();
    std::cout << "\nU U'\n" << U << std::endl;
    D = middleMatrix.getD();
    //  Compute inverse of diagonal
    for (c = 0; c < nSites; c++)
    {
        D(c, c) = 1 / D(c, c);
    }
    std::cout << "\nD'^-1\n" << D << std::endl;
    V = middleMatrix.getV() * V;
    std::cout << "\nV' V\n" << V << std::endl;
    
    return V.inverse() * D * U.inverse(); //    Green's function
    
    
}

Eigen::MatrixXd computeGreenNaive(Eigen::MatrixXd Bs[], int L, int nSites)
{
    int slice;
    Eigen::MatrixXd M = Eigen::MatrixXd::Identity(nSites, nSites);
    for (slice = L - 1; slice >= 0; slice--)
    {
        M *= Bs[slice];
    }
    
    M += Eigen::MatrixXd::Identity(nSites, nSites);
    
    return M.inverse();
}



