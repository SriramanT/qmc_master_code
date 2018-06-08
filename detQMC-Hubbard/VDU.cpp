//
//  VDU.cpp
//  
//
//  Created by Francisco Brito on 21/05/2018.
//

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "VDU.h"

VDU::VDU(Eigen::MatrixXd toDecompose)
{
    m = toDecompose.colwise().reverse().transpose();
}

void VDU::printMatrixToDecompose()
{
    std::cout << "\nMatrixToDecompose\n" << m << std::endl;
}

Eigen::MatrixXd VDU::QR_and_getU()
{
    Eigen::MatrixXd Q;
    Eigen::HouseholderQR<Eigen::MatrixXd> qrHH(m);
    R = qrHH.matrixQR().triangularView<Eigen::Upper>();
    R = R.transpose().colwise().reverse().rowwise().reverse().eval();
    Q = qrHH.householderQ();
    return Q.transpose().colwise().reverse();
}

Eigen::MatrixXd VDU::getD()
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(m.rows(), m.rows());
    D.diagonal() = R.diagonal();
    return D;
}

Eigen::MatrixXd VDU::getV()
{
    int a;
    int b;
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(m.rows(), m.rows());
    for (a = 0; a < m.rows(); a++)
    {
        for (b = a + 1; b < m.rows(); b++) // unit upper triangular -> loop starts in i + 1
        {
            V(a, b) = R(a, b) / R(b, b);
        }
    }
    return V;
}


