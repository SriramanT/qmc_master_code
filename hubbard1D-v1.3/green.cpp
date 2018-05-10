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

Green::Green(int nSites)
{
    M = Eigen::MatrixXd::Identity(nSites, nSites);
}
    
void Green::computeGreenNaive(Eigen::MatrixXd* Bs, int L)
{
    int slice;
    int nSites = Bs[0].rows();
    for (slice = 0; slice < L; slice++)
    {
        M *= Bs[L - 1 - slice];
    }
    M += Eigen::MatrixXd::Identity(nSites, nSites);
    G = M.inverse();
}

void Green::computeWrappedGreenNaive(Eigen::MatrixXd* Bs, int L, int l)
{
    int slice;
    int nSites = Bs[0].rows();
    for (slice = l; slice >= 0; slice--)
    {
        M *= Bs[slice];
    }
    
    for (slice = L - 1; slice > l; slice--)
    {
        M *= Bs[slice];
    }
    M += Eigen::MatrixXd::Identity(nSites, nSites);
    G = M.inverse();
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
