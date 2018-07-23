//
//  createHoppingMatrix.cpp
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

Eigen::MatrixXd genHoppingMatrix(int N)
{
    //  Generate the hopping matrix
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(N, N);
    
    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    K(0, 1) += 1;
    K(0, N - 1) += 1;
    K(N - 1, 0) += 1;
    K(N - 1, N - 2) += 1;
    
    //  Set the remaining ones
    int i;
    for (i = 1; i < N - 1; i++)
    {
        K(i, i - 1) += 1;
        K(i, i + 1) += 1;
    }
    
    return K;
}

Eigen::MatrixXd genHsMatrix(int L, int N)
{
    //  Generate the HS field matrix
    int l;
    int i;
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(L,N);
    for (l = 0; l < L; l++)
    {
        for (i = 0; i < N; i++)
        {
            if ( h(l, i) < 0 )
            {
                h(l, i) = -1;
            }
            else
            {
                h(l, i) = 1;
            }
        }
    }
    return h;
}

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

Eigen::MatrixXd regenB(bool spin, double nu, int N, Eigen::RowVectorXd h_l, Eigen::MatrixXd BpreFactor, double dt, double mu)
{
    int i;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(N, N);
    if (spin == true)
    {
        for (i = 0; i < N; i++)
        {
            B(i, i) = nu * h_l(i) + dt * mu;
        }
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            B(i, i) = - nu * h_l(i) + dt * mu;
        }
    }
    return BpreFactor * B.exp();
}




