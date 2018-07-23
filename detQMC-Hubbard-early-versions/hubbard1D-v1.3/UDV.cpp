//
//  UDV.cpp
//  
//
//  Created by Francisco Brito on 08/05/2018.
//
//  BEWARE! You must run the method .QR_and_getU() to actually
//  do the QR decomposition. The method returns U. To get
//  D, and V, you must run .getD() and .getV().

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "UDV.h"

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
    Eigen::HouseholderQR<Eigen::MatrixXd> qrHH(m);
    R = qrHH.matrixQR().triangularView<Eigen::Upper>();
    
    return qrHH.householderQ();
}

Eigen::MatrixXd UDV::getD()
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(m.rows(), m.rows());
    D.diagonal() = R.diagonal();
    return D;
}

Eigen::MatrixXd UDV::getV()
{
    int a;
    int b;
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


