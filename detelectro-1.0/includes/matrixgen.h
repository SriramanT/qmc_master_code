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
    Eigen::Matrix<double, N, N> B;
public:
    void oneDimensionalChainPBC(int t, double dt, double mu);
    void oneDimensionalChainOBC(int t, double dt, double mu);
    void twoDimensionalRectanglePBC(int Nx, int t, double dt, double mu);
    void twoDimensionalRectangleOBC(int Nx, int t, double dt, double mu);
    void nanoribbon(int Nx, int t, double dt, double mu); //  Nx = width of the ribbon
    Eigen::Matrix<double, N, N> BpreFactor();
};

template<int N>
void Geometry<N>::oneDimensionalChainPBC(int t, double dt, double mu)
{
//    Eigen::Matrix<double, N, N> HoppingMatrix = Eigen::Matrix<double, N, N>::Zero();
//    //  Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
//    HoppingMatrix(0, 1) += 1.;
//    HoppingMatrix(0, N - 1) += 1.;
//    HoppingMatrix(N - 1, 0) += 1.;
//    HoppingMatrix(N - 1, N - 2) += 1.;
//    //  Set the remaining ones
//    for (int i = 1; i < N - 1; i++)
//    {
//        HoppingMatrix(i, i - 1) += 1; HoppingMatrix(i, i + 1) += 1;
//    }
    //  Compute the exponential
    Eigen::Matrix<double, N, N> exp_k1 = Eigen::Matrix<double, N, N>::Zero();
    Eigen::Matrix<double, N, N> exp_k2 = Eigen::Matrix<double, N, N>::Zero();
    for (int i = 1; i < N / 2; i++)
    {
        exp_k1(2 * i, 2 * i) = cosh( t * dt );
        exp_k1(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1(0, 0) = cosh( t * dt );
    exp_k1(1, 1) = cosh( t * dt );
    exp_k1(0, 1) = sinh( t * dt );
    exp_k1(1, 0) = sinh( t * dt );
    exp_k2(0, 0) = cosh( t * dt / 2 );
    exp_k2(N - 1, N - 1) = cosh( t * dt / 2 );
    exp_k2(0, N - 1) = sinh( t * dt / 2 );
    exp_k2(N - 1, 0) = sinh( t * dt / 2 );
    B = exp(dt * mu) * Eigen::Matrix<double, N, N>::Identity() * exp_k2 * exp_k1 * exp_k2;
}

template<int N>
void Geometry<N>::oneDimensionalChainOBC(int t, double dt, double mu)
{
//    Eigen::Matrix<double, N, N> HoppingMatrix = Eigen::Matrix<double, N, N>::Zero();
//    //  Set the elements of the hopping matrix that define OBC corresponding to the ends of the 1D chain
//    HoppingMatrix(0, 1) += 1.;
//    //  HoppingMatrix(0, N - 1) += 1.; //   This hopping only occurs for PBC
//    //  HoppingMatrix(N - 1, 0) += 1.; //   So does this one
//    HoppingMatrix(N - 1, N - 2) += 1.;
//    //  Set the remaining ones
//    for (int i = 1; i < N - 1; i++)
//    {
//        HoppingMatrix(i, i - 1) += 1; HoppingMatrix(i, i + 1) += 1;
//    }
//      //Compute the exponential
//    B = exp(dt * mu) * Eigen::Matrix<double, N, N>::Identity() * (dt * t * HoppingMatrix).exp();
    Eigen::Matrix<double, N, N> exp_k1 = Eigen::Matrix<double, N, N>::Zero();
    Eigen::Matrix<double, N, N> exp_k2 = Eigen::Matrix<double, N, N>::Zero();
    for (int i = 1; i < N / 2; i++)
    {
        exp_k1(2 * i, 2 * i) = cosh( t * dt );
        exp_k1(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1(0, 0) = cosh( t * dt );
    exp_k1(1, 1) = cosh( t * dt );
    exp_k1(0, 1) = sinh( t * dt );
    exp_k1(1, 0) = sinh( t * dt );
    exp_k2(0, 0) = 1;
    exp_k2(N - 1, N - 1) = 1;
    B = exp(dt * mu) * Eigen::Matrix<double, N, N>::Identity() * exp_k2 * exp_k1 * exp_k2;
}

template<int N>
void Geometry<N>::twoDimensionalRectanglePBC(int Nx, int t, double dt, double mu)
{
    int Ny = N / Nx;
    //  Compute the exponential
    Eigen::MatrixXd exp_k1x = Eigen::MatrixXd::Zero(Nx, Nx);
    Eigen::MatrixXd exp_k2x = Eigen::MatrixXd::Zero(Nx, Nx);
    for (int i = 1; i < Nx / 2; i++)
    {
        exp_k1x(2 * i, 2 * i) = cosh( t * dt );
        exp_k1x(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1x(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1x(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2x(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2x(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1x(0, 0) = cosh( t * dt );
    exp_k1x(1, 1) = cosh( t * dt );
    exp_k1x(0, 1) = sinh( t * dt );
    exp_k1x(1, 0) = sinh( t * dt );
    exp_k2x(0, 0) = cosh( t * dt / 2 );
    exp_k2x(Nx - 1, Nx - 1) = cosh( t * dt / 2 );
    exp_k2x(0, Nx - 1) = sinh( t * dt / 2 );
    exp_k2x(Nx - 1, 0) = sinh( t * dt / 2 );
    
    Eigen::MatrixXd exp_k1y = Eigen::MatrixXd::Zero(Ny, Ny);
    Eigen::MatrixXd exp_k2y = Eigen::MatrixXd::Zero(Ny, Ny);
    for (int i = 1; i < Ny / 2; i++)
    {
        exp_k1y(2 * i, 2 * i) = cosh( t * dt );
        exp_k1y(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1y(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1y(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2y(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2y(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1y(0, 0) = cosh( t * dt );
    exp_k1y(1, 1) = cosh( t * dt );
    exp_k1y(0, 1) = sinh( t * dt );
    exp_k1y(1, 0) = sinh( t * dt );
    exp_k2y(0, 0) = cosh( t * dt / 2 );
    exp_k2y(Ny - 1, Ny - 1) = cosh( t * dt / 2 );
    exp_k2y(0, Ny - 1) = sinh( t * dt / 2 );
    exp_k2y(Ny - 1, 0) = sinh( t * dt / 2 );
    
    exp_k1x = exp_k2x * exp_k1x * exp_k2x;
    exp_k1y = exp_k2y * exp_k1y * exp_k2y;
    
    B = exp(dt * mu) * Eigen::Matrix<double, N, N>::Identity() * Eigen::kroneckerProduct(exp_k1y , exp_k1x);
}

template<int N>
void Geometry<N>::twoDimensionalRectangleOBC(int Nx, int t, double dt, double mu)
{
    int Ny = N / Nx;
    //  Compute the exponential
    Eigen::MatrixXd exp_k1x = Eigen::MatrixXd::Zero(Nx, Nx);
    Eigen::MatrixXd exp_k2x = Eigen::MatrixXd::Zero(Nx, Nx);
    for (int i = 1; i < Nx / 2; i++)
    {
        exp_k1x(2 * i, 2 * i) = cosh( t * dt );
        exp_k1x(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1x(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1x(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2x(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2x(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2x(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1x(0, 0) = cosh( t * dt );
    exp_k1x(1, 1) = cosh( t * dt );
    exp_k1x(0, 1) = sinh( t * dt );
    exp_k1x(1, 0) = sinh( t * dt );
    exp_k2x(0, 0) = 1;
    exp_k2x(Nx - 1, Nx - 1) = 1;
    
    Eigen::MatrixXd exp_k1y = Eigen::MatrixXd::Zero(Ny, Ny);
    Eigen::MatrixXd exp_k2y = Eigen::MatrixXd::Zero(Ny, Ny);
    for (int i = 1; i < Ny / 2; i++)
    {
        exp_k1y(2 * i, 2 * i) = cosh( t * dt );
        exp_k1y(2 * i + 1, 2 * i + 1) = cosh( t * dt );
        exp_k1y(2 * i, 2 * i + 1) = sinh( t * dt );
        exp_k1y(2 * i + 1, 2 * i) = sinh( t * dt );
        exp_k2y(2 * i, 2 * i) = cosh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i - 1) = cosh( t * dt / 2 );
        exp_k2y(2 * i, 2 * i - 1) = sinh( t * dt / 2 );
        exp_k2y(2 * i - 1, 2 * i) = sinh( t * dt / 2 );
    }
    exp_k1y(0, 0) = cosh( t * dt );
    exp_k1y(1, 1) = cosh( t * dt );
    exp_k1y(0, 1) = sinh( t * dt );
    exp_k1y(1, 0) = sinh( t * dt );
    exp_k2y(0, 0) = 1;
    exp_k2y(Ny - 1, Ny - 1) = 1;
    
    exp_k1x = exp_k2x * exp_k1x * exp_k2x;
    exp_k1y = exp_k2y * exp_k1y * exp_k2y;
    
    B = exp(dt * mu) * Eigen::Matrix<double, N, N>::Identity() * Eigen::kroneckerProduct(exp_k1y , exp_k1x);
}

template<int N>
void Geometry<N>::nanoribbon(int Ny, int t, double dt, double mu)
{   //  Nx = width of the ribbon
    Eigen::Matrix<double, N, N> HoppingMatrix = Eigen::Matrix<double, N, N>::Zero();
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
    //  Compute the exponential
    B = exp(dt * mu) * Eigen::Matrix<double, N, N>::Identity() * (t * dt * HoppingMatrix).exp();
}

template<int N>
Eigen::Matrix<double, N, N> Geometry<N>::BpreFactor()
{
    return B;
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
