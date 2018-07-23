#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <typeinfo>
#include <cmath>
#include <random>
#include <time.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <fstream>

//#include "QRP.h"
#include "SVD.h"

//  Try out the possibilities in Bai2010

int main()
{
    const int N = 4;
    Eigen::Matrix<double, N, N> A = Eigen::Matrix<double, N, N>::Random();
    Eigen::Matrix<double, N, N> B = Eigen::Matrix<double, N, N>::Random();
    std::cout << A << std::endl << std::endl;
    std::cout << B << std::endl << std::endl;
//    QRP< N > dec1;
//    Eigen::Matrix<double, N, N> q1 = dec1.QR_and_getQ(A);
//    Eigen::Matrix<double, N, N> d1 = dec1.getD();
//    Eigen::Matrix<double, N, N> v1 = dec1.getT();
//    QRP< N > dec2;
//    Eigen::Matrix<double, N, N> q2 = dec2.QR_and_getQ(B * q1  * d1);
//    Eigen::Matrix<double, N, N> d2 = dec2.getD();
//    Eigen::Matrix<double, N, N> v2 = dec2.getT();
    
    SVD< N > dec1;
    dec1.doSVD(A);
    Eigen::Matrix<double, N, N> u1 = dec1.getU();
    Eigen::Matrix<double, N, N> s1 = dec1.getS();
    Eigen::Matrix<double, N, N> v1 = dec1.getV();
    SVD< N > dec2;
    dec2.doSVD(B * u1 * s1);
    Eigen::Matrix<double, N, N> u2 = dec2.getU();
    Eigen::Matrix<double, N, N> s2 = dec2.getS();
    Eigen::Matrix<double, N, N> v2 = dec2.getV();
    
    
    std::cout << u2 * s2 * (v2 * v1) << std::endl << std::endl;
    std::cout << B * A << std::endl << std::endl;
    
    return 0;
}
