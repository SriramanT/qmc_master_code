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
}

void Green::computeGreenNaive(Eigen::MatrixXd* Bs, int l)
{
    M = Eigen::MatrixXd::Identity(N, N);
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

void Green::computeStableGreenNaive(Eigen::MatrixXd* Bs, int l, int newL)
{
    int site;
    int interval;
    int slice = l;
    int countSlice = 0;
    Eigen::MatrixXd partialProdBs;
    Eigen::MatrixXd U = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd D = Eigen::MatrixXd::Identity(N, N);
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(N, N);
    
    //  Compute partial products, UDV's and multiply
    for (interval = 0; interval < newL; interval++)
    {
        //  Initialiaze partial product
        partialProdBs = Eigen::MatrixXd::Identity(N, N);
        for (countSlice = 0; countSlice < L / newL; countSlice++)
        {
            if ( slice == ( L - 1 ) )
            {
                slice = 0;
            }
            else
            {
                slice += 1;
            }
            //  Compute partial product of (L / newL) matrices
            partialProdBs = Bs[slice] * partialProdBs;
        }
        //  Householder QR applied to (previous partial product) * U * D
        //  where U, and D are from the current partial product's decomposition
        UDV m1(partialProdBs * U * D);
        U = m1.QR_and_getU();
        D = m1.getD();
        V = m1.getV() * V;
    }
    //  Compute Green's function.
    //  See eq. (4.3) of "Stable simulations of models of interacting electrons"
    //  by Loh Jr and Gubernatis
    UDV middleMatrix(U.inverse() * V.inverse() + D);
    U = U * middleMatrix.QR_and_getU(); //  U U'
    D = middleMatrix.getD();
    //  Compute inverse of the diagonal matrix in the decomposition
    for (site = 0; site < N; site++)
    {
        D(site, site) = 1 / D(site, site);
    }
    V = middleMatrix.getV() * V;    //  V' V
    //  Final form of eq. (4.3)
    G = V.inverse() * D * U.inverse();
}

Eigen::MatrixXd Green::getG()
{
    return G;
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
