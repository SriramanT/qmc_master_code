//
//  matrixgen.h
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#ifndef matrixgen_h
#define matrixgen_h

template<int N>
class Geometry
{
    Eigen::Matrix<double, N, N> HoppingMatrix;
public:
    void oneDimensionalChainPBC();
    void nanoribbon(int Nx); //  Nx = width of the ribbon
    Eigen::Matrix<double, N, N> matrix();
};

template <int N>
void Geometry<N>::oneDimensionalChainPBC()
{
    HoppingMatrix = Eigen::Matrix<double, N, N>::Zero();
    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    HoppingMatrix(0, 1) += 1.;
    HoppingMatrix(0, N - 1) += 1.;
    HoppingMatrix(N - 1, 0) += 1.;
    HoppingMatrix(N - 1, N - 2) += 1.;
    //  Set the remaining ones
    for (int i = 1; i < N - 1; i++)
    {
        HoppingMatrix(i, i - 1) += 1; HoppingMatrix(i, i + 1) += 1;
    }
}

template <int N>
void Geometry<N>::nanoribbon(int Ny)
{   //  Nx = width of the ribbon
    HoppingMatrix = Eigen::Matrix<double, N, N>::Zero();
    int Nx = N / Ny / 2;
    for (int x = 0; x < Nx; x++)
    {
        for (int y = 0; y < Ny; y++)
        {
            //  SUBLATTICE A
            if (y == Ny - 1)
            {
                if (x == 0)
                {
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1;
                }
                else
                {
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1;
                }
            }
            else
            {
                if (x == 0)
                {
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + Nx - 1) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * (y + 1) + Nx - 1) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + Nx - 1, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * (y + 1) + Nx - 1, Nx * y + x) = 1;
                }
                else
                {
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + x - 1) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * ( y + 1 ) + x - 1) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x - 1, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * ( y + 1 ) + x - 1, Nx * y + x) = 1;
                }
            }
            //  SUBLATTICE B
            if (y == 0)
            {
                if (x == Nx - 1)
                {
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1;
                    HoppingMatrix(Nx * y, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1;
                }
                else
                {
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1;
                }
            }
            else
            {
                if (x == Nx - 1)
                {
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + Nx - 1) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) ) = 1;
                    HoppingMatrix(Nx * y, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * y + Nx - 1, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * ( y - 1 ) , Nx * Ny + Nx * y + x) = 1;
                }
                else
                {
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + x) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * y + x + 1) = 1;
                    HoppingMatrix(Nx * Ny + Nx * y + x, Nx * ( y - 1 ) + x + 1) = 1;
                    HoppingMatrix(Nx * y + x, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * y + x + 1, Nx * Ny + Nx * y + x) = 1;
                    HoppingMatrix(Nx * ( y - 1 ) + x + 1, Nx * Ny + Nx * y + x) = 1;
                }
            }
        }
    }
}

template <int N>
Eigen::Matrix<double, N, N> Geometry<N>::matrix()
{
    return HoppingMatrix;
}

template<int L, int N>
class Configuration
{
    Eigen::Matrix<double, L, N> HSfield;
public:
    void genHsMatrix();
    Eigen::Matrix<double, L, N> matrix();
    double get(int x, int y);
    void flip(int l, int i);
};

template<int L, int N>
void Configuration<L, N>::genHsMatrix()
{
    //  Generate the HS field matrix
    int l;
    int i;
    HSfield = Eigen::Matrix<double, L, N>::Random(L,N);
    for (l = 0; l < L; l++)
    {
        for (i = 0; i < N; i++)
        {
            if ( HSfield(l, i) < 0 )
            {
                HSfield(l, i) = -1;
            }
            else
            {
                HSfield(l, i) = 1;
            }
        }
    }
}

template<int L, int N>
Eigen::Matrix<double, L, N> Configuration<L, N>::matrix()
{
    return HSfield;
}

template<int L, int N>
double Configuration<L, N>::get(int x, int y)
{
    return HSfield(x, y);
}

template<int L, int N>
void Configuration<L, N>::flip(int l, int i)
{
    HSfield(l, i) *= -1;
}

template<int N, int L>
class OneParticlePropagators
{
    Eigen::Matrix<double, N, N> B[L];
public:
    void fillMatrices(bool spin, double nu, Eigen::Matrix<double, L, N> h, Eigen::Matrix<double, N, N> BpreFactor);
    Eigen::Matrix<double, N, N> matrix(int l);
    Eigen::Matrix<double, N, N> * list();
    void update(int l, int i, double alpha);
};

template<int N, int L>
void OneParticlePropagators<N, L>::fillMatrices(bool spin, double nu, Eigen::Matrix<double, L, N> h, Eigen::Matrix<double, N, N> BpreFactor)
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
                B[l](i, i) = exp ( nu * h(l, i) ) ;
            }
        }
        else
        {
            for (i = 0; i < N; i++)
            {
                B[l](i, i) = exp ( - 1. * nu * h(l, i) );
            }
        }
        B[l] = BpreFactor * B[l] ;
    }
}

template<int N, int L>
void OneParticlePropagators<N, L>::update(int l, int i, double alpha)
{
    B[l].col(i) *= ( alpha + 1 );
}

template<int N, int L>
Eigen::Matrix<double, N, N> OneParticlePropagators<N, L>::matrix(int l)
{
    return B[l];
}

template<int N, int L>
Eigen::Matrix<double, N, N> * OneParticlePropagators<N, L>::list()
{
    return B;
}

#endif /* matrixgen_h */
