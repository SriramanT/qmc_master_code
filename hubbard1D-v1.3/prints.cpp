//
//  parameters.cpp
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#include <iostream>
#include <Eigen/Dense>

void initialPrint(int N, double dt, double beta, int L, double t, double U, double mu, int totalMCSteps, bool debug, Eigen::MatrixXd K, Eigen::MatrixXd B_preFactor, Eigen::MatrixXd h)
{
    std::cout << "\nNumber of sites N: " << N << std::endl;
    std::cout << "Trotter breakup parameter dt: " << dt << std::endl;
    std::cout << "Inverse temperature beta: " << beta << std::endl;
    std::cout << "Number of imaginary time slices L: " << L << std::endl;
    std::cout << "Hopping parameter t: " << t << std::endl;
    std::cout << "On-site interaction U: " << U << std::endl;
    std::cout << "Chemical potential mu: " << mu << std::endl;
    std::cout << "Number of (lattice + imaginary time) sweeps: " << totalMCSteps / N / L << " (" << totalMCSteps << " steps)" << std::endl;
    if (debug == true)
    {
        std::cout << "\n K = \n\n" << K << "\n\n";
        std::cout << "exp (t * dt * K) = \n\n" << B_preFactor << "\n\n";
        std::cout << "Hubbard Stratonovich Matrix\n\n h = \n\n" << h << "\n\n";
    }
}

