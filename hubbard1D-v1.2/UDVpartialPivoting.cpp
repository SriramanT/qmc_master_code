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

UDVpartialPivoting::UDVpartialPivoting(Eigen::MatrixXd toDecompose)
{
    m = toDecompose;
}

void UDVpartialPivoting::printMatrixToDecompose()
{
    std::cout << "\nMatrixToDecompose\n" << m << std::endl;
}

Eigen::MatrixXd UDVpartialPivoting::QR_and_getU()
{
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrPartial(m);
    R = qrPartial.matrixQR().template triangularView<Eigen::Upper>();
    P = qrPartial.colsPermutation();
    return qrPartial.householderQ();
}

Eigen::MatrixXd UDVpartialPivoting::getD()
{
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(NSITES, NSITES);
    D.diagonal() = R.diagonal();
    return D;
}

Eigen::MatrixXd UDVpartialPivoting::getV()
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Identity(NSITES, NSITES);
    for (i = 0; i < NSITES; i++)
    {
        for (p = i + 1; p < NSITES; p++) // unit upper triangular -> loop starts in i + 1
        {
            V(i, p) = R(i, p) / R(i, i);
        }
    }
    return V;
}

Eigen::MatrixXd UDVpartialPivoting::getP()
{
    return P;
}
