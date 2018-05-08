//
//  parameters.cpp
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#include <iostream>
#include "parameters.h"
#include <Eigen/Dense>

void printParameters()
{
    std::cout << "dt: " << dt << std::endl;
    std::cout << "beta: " << beta << std::endl;
    std::cout << "L: " << L << std::endl;
    std::cout << "t: " << t << std::endl;
    std::cout << "U: " << U << std::endl;
    std::cout << "mu: " << mu << std::endl;
}

void printWelcome(int nSites)
{
    std::cout << "\n\nDQMC for 1D Hubbard chain\n" << std::endl;
    std::cout << "\nNumber of sites: " << nSites << std::endl;
}

void printMC(int totalMCSteps, int nSites, int L)
{
    std::cout << "Number of MC sweeps: " << totalMCSteps / (nSites*L) << " (" << totalMCSteps << " steps)" << "\n\n";
}

void printStartingMatrices(Eigen::MatrixXd K, Eigen::MatrixXd B_preFactor, Eigen::MatrixXd h, bool yesOrNo)
{
    if (yesOrNo == true)
    {
        std::cout << "\n K = \n\n" << K << "\n\n";
        std::cout << "exp (t * dt * K) = \n\n" << B_preFactor << "\n\n";
        std::cout << "Hubbard Stratonovich Field\n\n h = \n\n" << h << "\n\n";
    }
}
