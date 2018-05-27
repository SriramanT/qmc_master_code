//
//  checkerboard.cpp
//  
//
//  Created by Francisco Brito on 23/05/2018.
//

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "matrixgen.h"

int main()
{
    int N = 2;  //  # sites
    const Eigen::MatrixXd K = genHoppingMatrix(N);
    std::cout << K.exp() << std::endl;

    return 0;
}

