//
//  green.cpp
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include "green.h"
#include "UDV.h"

Green::Green(int nSites, int nSlices)
{
    N = nSites;
    L = nSlices;
    M = Eigen::MatrixXd::Identity(nSites, nSites);
}

void Green::computeGreenNaive(Eigen::MatrixXd* Bs, int l)
{
    int slice;
    for (slice = l; slice >= 0; slice--)
    {
        M *= Bs[slice];
    }
    
    for (slice = L - 1; slice > l; slice--)
    {
        M *= Bs[slice];
    }
    M += Eigen::MatrixXd::Identity(N, N);
    G = M.inverse();
}

void Green::computeStableGreenNaive(Eigen::MatrixXd* Bs, int k)
{
    int c, d;
    Eigen::MatrixXd partialProdBs;
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(N, N);
    
    //  Compute partial products, UDV's and multiply
    for (c = 0; c < k; c++)
    {
        partialProdBs = Eigen::MatrixXd::Identity(N, N);
        for (d = 0; d < L/k; d++)
        {
            partialProdBs = Bs[c * (L/k) + d] * partialProdBs;
        }
        //  Householder QR
        UDV m1(partialProdBs * U * D);
        U = m1.QR_and_getU();
        D = m1.getD();
        V = m1.getV() * V;
    }
    UDV middleMatrix(U.inverse() * V.inverse() + D);
    U = U * middleMatrix.QR_and_getU();
    D = middleMatrix.getD();
    //  Compute inverse of diagonal
    for (c = 0; c < N; c++)
    {
        D(c, c) = 1 / D(c, c);
    }
    V = middleMatrix.getV() * V;
    
    G = V.inverse() * D * U.inverse(); //    Green's function
}

Eigen::MatrixXd Green::getG()
{
    return G;
}

void Green::resetM()
{
    M = Eigen::MatrixXd::Identity(N, N);
}

Eigen::MatrixXd Green::getM()
{
    return M;
}

void Green::printGreen(bool spin)
{
    if (spin == true)
    {
        std::cout << "\nGreenUp\n\n" << G << std::endl;
    }
    else
    {
        std::cout << "\nGreenDown\n\n" << G << std::endl;
    }
}

void Green::printM(bool spin)
{
    if (spin == true)
    {
        std::cout << "\nMUp\n\n" << M << std::endl;
    }
    else
    {
        std::cout << "\nMDown\n\n" << M << std::endl;
    }
}
