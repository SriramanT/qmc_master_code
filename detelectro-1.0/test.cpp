//
//  main.cpp
//  
//
//  Created by Francisco Brito on 08/06/2018.
//
//  This program simulates the Hubbard model for an arbitrary geometry lattice
//  using auxiliary field (or determinant) Quantum Monte Carlo: in particular, the BSS algorithm.
//  The used notation is based on the lecture notes "Numerical Methods for Quantum Monte Carlo
//  Simulations of the Hubbard Model by Zhaojun Bai, Wenbin Chen, Richard Scalettar, and
//  Ichitaro Yamazaki (2009)
//

#include <iostream>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <fstream>
#include <iomanip>
#include <string>

#include "matrixgen.h"


int main(int argc, char **argv)
{
//    //  SET SIMULATION PARAMETERS.
//    const int N = 2;  //  # sites
//
    try
    {
        const int N = std::stoi(argv[1]);
        // -- INITIALIZATION ---
        
        //  HOPPING MATRIX FOR ARBITRARY GEOMETRY
        Geometry< N > K;
        //  1D HUBBARD CHAIN
        K.oneDimensionalChainPBC();
        
        std::cout << K.matrix() << std::endl ;
    }
    catch (std::invalid_argument const &ex)
    {
        std::cerr << "Invalid number " << argv[1] << '\n';
    }
    


    
    return 0;
}

