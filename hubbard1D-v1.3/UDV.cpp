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


