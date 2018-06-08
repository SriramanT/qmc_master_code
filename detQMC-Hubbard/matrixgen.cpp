//
//  createHoppingMatrix.cpp
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>

void genBmatrix(Eigen::MatrixXd* Bs, bool spin, double nu, int N, int L, Eigen::MatrixXd h, Eigen::MatrixXd BpreFactor, double dt, double mu)
{
    //  For a single B matrix do genBmatrix(&Bup[0], true, nu, N, L, h, BpreFactor, dt, mu); (to get the memory address of a given B-matrix)
    int i;
    int l;
    for (l = 0; l < L; l++)
    {
        Eigen::MatrixXd B = Eigen::MatrixXd::Identity(N, N);
        
        if (spin == true)
        {
            for (i = 0; i < N; i++)
            {
                B(i, i) = nu * h(l, i) + dt * mu;
            }
        }
        else
        {
            for (i = 0; i < N; i++)
            {
                B(i, i) = -1 * nu * h(l, i) + dt * mu;
            }
        }
        Bs[l] = BpreFactor * B.exp();
    }
}

Eigen::VectorXd uSigma(int N, Eigen::MatrixXd Green, int i)
{
    Eigen::MatrixXd e_i;
    e_i = Eigen::MatrixXd::Zero(N , 1) ;
    e_i(i, 0) = 1;
    return (Eigen::MatrixXd::Identity(N, N) - Green) * e_i;
}

Eigen::VectorXd wSigma(int N, Eigen::MatrixXd Green, int i)
{
    Eigen::MatrixXd e_i;
    e_i = Eigen::MatrixXd::Zero(N , 1) ;
    e_i(i, 0) = 1;
    return Green.transpose() * e_i;
}



