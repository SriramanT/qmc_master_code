//
//  matrixgen.h
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#ifndef matrixgen_h
#define matrixgen_h

template<int N>
class Hoppings
{
    Eigen::Matrix<double, N, N> K;
public:
    Eigen::Matrix<double, N, N> OneDimensionalChainPBC();
};

template <int N>
Eigen::Matrix<double, N, N> Hoppings<N>::OneDimensionalChainPBC()
{
    int i, j;
    for (j = 0 ; j < N; j++)
    {
        K(0, j) = 0;
        K(N - 1, j) = 0;
    }
    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    K(0, 1) += 1;
    K(0, N - 1) += 1;
    K(N - 1, 0) += 1;
    K(N - 1, N - 2) += 1;
    
    //  Set the remaining ones
    for (i = 1; i < N - 1; i++)
    {
        K(i, i - 1) += 1; K(i, i + 1) += 1;
        for (j = 0; j < i - 1; j++)
        {
            K(i, j) = 0;
        }
        K(i, i) = 0;
        for (j = i + 2; j < N; j++)
        {
            K(i, j) = 0;
        }
    }
    return K;
}

template<int L, int N>
class HSfield
{
    Eigen::Matrix<double, L, N> h;
public:
    Eigen::Matrix<double, L, N> genHsMatrix();
};

template<int L, int N>
Eigen::Matrix<double, L, N> HSfield<L, N>::genHsMatrix()
{
    //  Generate the HS field matrix
    int l;
    int i;
    h = Eigen::Matrix<double, L, N>::Random(L,N);
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

template<int N, int L>
class OneParticlePropagators
{
    Eigen::Matrix<double, N, N> B[L];
public:
    void fillArr(bool spin, double nu, double dt, double mu, Eigen::Matrix<double, L, N> h, Eigen::Matrix<double, N, N> BpreFactor);
    Eigen::Matrix<double, N, N> getB(int slice);
};

template<int N, int L>
void OneParticlePropagators<N, L>::fillArr(bool spin, double nu, double dt, double mu, Eigen::Matrix<double, L, N> h, Eigen::Matrix<double, N, N> BpreFactor)
{
    int i;
    int l;
    for (l = 0; l < L; l++)
    {
        B[l] = Eigen::Matrix<double, N, N>::Zero();

        if (spin == true)
        {
            for (i = 0; i < N; i++)
            {
                B[l](i, i) = nu * h(l, i) + dt * mu;
            }
        }
        else
        {
            for (i = 0; i < N; i++)
            {
                B[l](i, i) = -1 * nu * h(l, i) + dt * mu;
            }
        }
        B[l] = BpreFactor * ( B[l] ).exp();
    }
}

template<int N, int L>
Eigen::Matrix<double, N, N> OneParticlePropagators<N, L>::getB(int slice)
{
    return B[slice];
}

Eigen::VectorXd uSigma(int N, Eigen::MatrixXd Green, int i);
Eigen::VectorXd wSigma(int N, Eigen::MatrixXd Green, int i);

#endif /* matrixgen_h */
