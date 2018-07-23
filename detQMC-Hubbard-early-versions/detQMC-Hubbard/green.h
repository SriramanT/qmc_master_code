//
//  green.h
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#ifndef green_h
#define green_h

#include "UDV.h"
#include "VDU.h"

template<int N, int L, int Lbda>
class Green
{
    Eigen::Matrix<double, N, N> M;
    Eigen::Matrix<double, N, N> G;
    Eigen::Matrix<double, N, 1> u;
    Eigen::Matrix<double, 1, N> w;
    Eigen::Matrix<double, N, N> U;
    Eigen::Matrix<double, N, N> D;
    Eigen::Matrix<double, N, N> V;
    Eigen::Matrix<double, 2 * N, 2 * N> Udouble;
    Eigen::Matrix<double, 2 * N, 2 * N> Ddouble;
    Eigen::Matrix<double, 2 * N, 2 * N> Vdouble;
    Eigen::Matrix<double, N, N> Gforward;
    Eigen::Matrix<double, N, N> Gbackward;
    Eigen::Matrix<double, N, N> Gzero;
    //  ALLOCATE MEMORY TO STORE THE PARTIAL PRODUCTS INVOLVED IN
    //  SPEEDING UP THE LOW TEMPERATURE STABILIZATION.
    Eigen::Matrix<double, N, N> Us[Lbda];
    Eigen::Matrix<double, N, N> Ds[Lbda];
    Eigen::Matrix<double, N, N> Vs[Lbda];
public:
    void update(double alpha, double d, int i);
    void wrap(Eigen::Matrix<double, N, N> B);
    void computeGreenNaive(Eigen::Matrix<double, N, N>* Bs, int l);
    void computeStableGreenNaiveR(Eigen::Matrix<double, N, N>* Bs, int l);
    void computeStableGreenNaiveL(Eigen::Matrix<double, N, N>* Bs, int l);
    void computeGreenFromVDU();
    void initializeUneqs();
    void storeVDU(Eigen::Matrix<double, N, N>* Bs);
    void computeStableGreen(int l, int greenAfreshFreq);
    void computeBlockOfGreens(int l, int greenAfreshFreq);
    void storeUDV(Eigen::Matrix<double, N, N>* Bs, int l, int greenAfreshFreq);
    Eigen::Matrix<double, N, N> matrix();
    double uneqForward(int x, int y);
    double uneqBackward(int x, int y);
    double zero(int x, int y);
    double get(int x, int y);
    Eigen::Matrix<double, N, N> getM();
};

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::update(double alpha, double d, int i)
{
    u = ( Eigen::Matrix<double, N, N>::Identity() - G ).col(i);
    w = G.row(i);
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < N; y++)
        {
            G(x, y) -= alpha / d * u(x) * w(y);
        }
    }
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::wrap(Eigen::Matrix<double, N, N> B)
{
    G = B * G * B.inverse();
    Gforward = B * Gforward;
    Gbackward = Gbackward * B.inverse();
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeGreenNaive(Eigen::Matrix<double, N, N>* Bs, int l)
{
    M = Eigen::Matrix<double, N, N>::Identity();
    int slice;
    for (slice = l; slice >= 0; slice--)
    {
        M = M * Bs[slice];
    }
    for (slice = L - 1; slice > l; slice--)
    {
        M = M * Bs[slice];
    }
    M += Eigen::Matrix<double, N, N>::Identity();
    G = M.inverse();
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeStableGreenNaiveR(Eigen::Matrix<double, N, N>* Bs, int l)
{
    U = Eigen::Matrix<double, N, N>::Identity();
    D = Eigen::Matrix<double, N, N>::Identity();
    V = Eigen::Matrix<double, N, N>::Identity();
    
    int site;
    int slice;
    int sliceCounter = 0;
    
    //  INITIALIZE PARTIAL PRODUCTS.
    Eigen::Matrix<double, N, N> partialProdBs = Eigen::Matrix<double, N, N>::Identity();  //  an array of partial products
    
    for (slice = l + 1; slice < L; slice++)
    {
        partialProdBs = Bs[slice] * partialProdBs;
        sliceCounter += 1;
        if (sliceCounter == L / Lbda)
        {
            //  Householder QR applied to (previous partial product) * U * D
            //  where U, and D are from the current partial product's decomposition
            UDV udvRight( ( partialProdBs * U ) * D);
            U = udvRight.QR_and_getU();
            D = udvRight.getD();
            V = udvRight.getV() * V;
            sliceCounter = 0;
            partialProdBs = Eigen::Matrix<double, N, N>::Identity();
        }
    }
    
    for (slice = 0; slice <= l; slice++)
    {
        partialProdBs = Bs[slice] * partialProdBs;
        sliceCounter += 1;
        if (sliceCounter == L / Lbda)
        {
            //  Householder QR applied to (previous partial product) * U * D
            //  where U, and D are from the current partial product's decomposition
            UDV udvRight( ( partialProdBs * U ) * D);
            U = udvRight.QR_and_getU();
            D = udvRight.getD();
            V = udvRight.getV() * V;
            sliceCounter = 0;
            partialProdBs = Eigen::Matrix<double, N, N>::Identity();
        }
    }
    
    //  COMPUTE GREEN'S FUNCTION.
    //  See eq. (4.3) of "Stable simulations of models of interacting electrons"
    //  by Loh Jr and Gubernatis
    UDV middleMatrix(U.inverse() * V.inverse() + D);
    U = U * middleMatrix.QR_and_getU(); //  U U'
    D = middleMatrix.getD();
    //  Compute inverse of the diagonal matrix in the decomposition
    for (site = 0; site < N; site++)
    {
        D(site, site) = 1 / D(site, site);
    }
    V = middleMatrix.getV() * V;    //  V' V
    //  Final form of eq. (4.3)
    G = V.inverse() * D * U.inverse();
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeStableGreenNaiveL(Eigen::Matrix<double, N, N>* Bs, int l)
{
    U = Eigen::Matrix<double, N, N>::Identity();
    D = Eigen::Matrix<double, N, N>::Identity();
    V = Eigen::Matrix<double, N, N>::Identity();
    
    int site;
    int slice;
    int sliceCounter = 0;
    
    //  INITIALIZE PARTIAL PRODUCTS.
    Eigen::Matrix<double, N, N> partialProdBs = Eigen::Matrix<double, N, N>::Identity();  //  an array of partial products

    for (slice = l; slice >= 0; slice--)
    {
        partialProdBs *= Bs[slice];
        sliceCounter += 1;
        if (sliceCounter == L / Lbda)
        {
            //  Householder QR applied to (previous partial product) * U * D
            //  where U, and D are from the current partial product's decomposition
            VDU vduLeft(D * U * partialProdBs);
            U = vduLeft.QR_and_getU();
            D = vduLeft.getD();
            V = V * vduLeft.getV();
            sliceCounter = 0;
            partialProdBs = Eigen::Matrix<double, N, N>::Identity();
        }
    }
    
    for (slice = L - 1; slice > l; slice--)
    {
        partialProdBs *= Bs[slice];
        sliceCounter += 1;
        if (sliceCounter == L / Lbda)
        {
            //  Householder QR applied to (previous partial product) * U * D
            //  where U, and D are from the current partial product's decomposition
            VDU vduLeft(D * U * partialProdBs);
            U = vduLeft.QR_and_getU();
            D = vduLeft.getD();
            V = V * vduLeft.getV();
            sliceCounter = 0;
            partialProdBs = Eigen::Matrix<double, N, N>::Identity();
        }
    }

    //  COMPUTE GREEN'S FUNCTION.
    //  See eq. (4.3) of "Stable simulations of models of interacting electrons"
    //  by Loh Jr and Gubernatis
    VDU middleMatrix(V.inverse() * U.inverse() + D);
    U = middleMatrix.QR_and_getU() * U; //  U' U
    D = middleMatrix.getD();
    //  Compute inverse of the diagonal matrix in the decomposition
    for (site = 0; site < N; site++)
    {
        D(site, site) = 1 / D(site, site);
    }
    V = V * middleMatrix.getV();    //  V V'
    //  Final form of eq. (4.3)
    G = U.inverse() * D * V.inverse();
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeGreenFromVDU()
{
    U = Us[Lbda - 1];
    D = Ds[Lbda - 1];
    V = Vs[Lbda - 1];
    //  COMPUTE GREEN'S FUNCTION.
    //  See eq. (4.3) of "Stable simulations of models of interacting electrons"
    //  by Loh Jr and Gubernatis
    VDU middleMatrix(V.inverse() * U.inverse() + D);
    U = middleMatrix.QR_and_getU() * U; //  U' U
    D = middleMatrix.getD();
    //  Compute inverse of the diagonal matrix in the decomposition
    for (int site = 0; site < N; site++)
    {
        D(site, site) = 1 / D(site, site);
    }
    V = V * middleMatrix.getV();    //  V V'
    //  Final form of eq. (4.3)
    G = U.inverse() * D * V.inverse();
    Gzero = G;
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::initializeUneqs()
{
    Gforward = Gzero;
    Gbackward = Gzero - Eigen::Matrix<double, N, N>::Identity();
}


template<int N, int L, int Lbda>
void Green<N, L, Lbda>::storeVDU(Eigen::Matrix<double, N, N>* Bs)
{
    //  This is very similar to the function above. The only difference is that the
    //  Green's function is not computed, and this is to be used at the beginning of
    //  each imaginary time sweep, i.e. M = I + B_{L-1} B_{L-2}... B_{0}
    U = Eigen::Matrix<double, N, N>::Identity();
    D = Eigen::Matrix<double, N, N>::Identity();
    V = Eigen::Matrix<double, N, N>::Identity();
    
    int slice;
    int sliceCounter = 0;
    int lbda = 0;
    
    //  INITIALIZE PARTIAL PRODUCTS.
    Eigen::Matrix<double, N, N> partialProdBs = Eigen::Matrix<double, N, N>::Identity();  //  an array of partial products
    
    //  COMPUTE UDVs.
    for (slice = L - 1; slice >= 0; slice--)
    {
        partialProdBs *= Bs[slice];
        sliceCounter += 1;
        if (sliceCounter == L / Lbda)
        {
            //  Householder QR applied to (previous partial product) * U * D
            //  where U, and D are from the current partial product's decomposition
            VDU vduLeft(D * U * partialProdBs);
            U = vduLeft.QR_and_getU();
            Us[lbda] = U;
            D = vduLeft.getD();
            Ds[lbda] = D;
            V = V * vduLeft.getV();
            Vs[lbda] = V;
            sliceCounter = 0;
            lbda += 1;
            partialProdBs = Eigen::Matrix<double, N, N>::Identity();
        }
    }
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::storeUDV(Eigen::Matrix<double, N, N>* Bs, int l, int greenAfreshFreq)
{
    //  NOTE THAT THE ARGUMENT l will never be zero, since at l = 0, we compute the VDUs from the right and the Green's function
    int lbda = (l + 1) / greenAfreshFreq - 1;
    int slice;
    Eigen::Matrix<double, N, N> partialProdBs = Eigen::Matrix<double, N, N>::Identity();  //  an array of partial products
    
    //  INITIALIZE PARTIAL PRODUCTS.
    if ( lbda == 0 )
    {
        U = Eigen::Matrix<double, N, N>::Identity();
        D = Eigen::Matrix<double, N, N>::Identity();
        V = Eigen::Matrix<double, N, N>::Identity();
    }
    else
    {
        U = Us[Lbda - lbda] ;
        D = Ds[Lbda - lbda] ;
        V = Vs[Lbda - lbda] ;
    }
    
    //  COMPUTE UDVs.
    for (slice = ( l - (greenAfreshFreq - 1) ) ; slice <= l ; slice++)
    {
        partialProdBs = Bs[slice] * partialProdBs;
    }

    //  Householder QR applied to (previous partial product) * U * D
    //  where U, and D are from the current partial product's decomposition
    UDV udvRight(partialProdBs * U * D);
    U = udvRight.QR_and_getU();
    D = udvRight.getD();
    V = udvRight.getV() * V;
    //  STORE THEM. NOTE THAT WE SHOULD NEVER REACH lbda = Lbda. At that point, we must store the VDUs to restart.
    //  Hence, this is safe! Otherwise, it wouldn't be.
    Us[Lbda - lbda - 1] = U;
    Ds[Lbda - lbda - 1] = D;
    Vs[Lbda - lbda - 1] = V;
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeStableGreen(int l, int greenAfreshFreq)
{
    //  NOTE THAT THE ARGUMENT l will never be zero, since at l = 0, we compute the VDUs from the right and the Green's function
    int lbda = (l + 1) / greenAfreshFreq - 1;
    UDV udvIllConditionedSum( Us[Lbda - lbda - 1].inverse() * Us[Lbda - lbda - 2].inverse()
                             + Ds[Lbda - lbda - 1] * Vs[Lbda - lbda - 1] * Vs[Lbda - lbda - 2] * Ds[Lbda - lbda - 2] );
    //  U_R.inv * U_L.inv + D_R * V_R * V_L * D_L
    U = udvIllConditionedSum.QR_and_getU();
    D = udvIllConditionedSum.getD();
    V = udvIllConditionedSum.getV();
    G = Us[Lbda - lbda - 2].inverse() * V.inverse() * D.inverse() * U.inverse() * Us[Lbda - lbda - 1].inverse();
    //  U_L.inv V.inv D.inv U.inv U_R.inv
}

template<int N, int L, int Lbda>
void Green<N, L, Lbda>::computeBlockOfGreens(int l, int greenAfreshFreq)
{
    //  Compute both equal time and unequal time Green's functions
    //  (see Quantum Monte Carlo Methods on Lattices: The Determinantal Approach by Fakher F. Assaad)
    int lbda = (l + 1) / greenAfreshFreq - 1;
    Eigen::Matrix<double, 2 * N, 2 * N> tempMatrix;
    
    tempMatrix.block(0, 0, N, N) = Vs[Lbda - lbda - 2].inverse() * Vs[Lbda - lbda - 1].inverse();
    tempMatrix.block(0, N, N, N) = Ds[Lbda - lbda - 2];
    tempMatrix.block(N, 0, N, N) = - Ds[Lbda - lbda - 1];
    tempMatrix.block(N, N, N, N) =  Us[Lbda - lbda - 1].inverse() * Us[Lbda - lbda - 2].inverse();
    
    UDV decomposition(tempMatrix);
    Udouble = decomposition.QR_and_getU();
    Ddouble = decomposition.getD();
    Vdouble = decomposition.getV();
    
    Eigen::Matrix<double, 2 * N, 2 * N> RightMatrix;
    RightMatrix.block(0, 0, N, N) = Vs[Lbda - lbda - 2].inverse();
    RightMatrix.block(0, N, N, N) = Eigen::Matrix<double, N, N>::Zero();
    RightMatrix.block(N, 0, N, N) = Eigen::Matrix<double, N, N>::Zero();
    RightMatrix.block(N, N, N, N) = Us[Lbda - lbda - 1].inverse();
    RightMatrix = Udouble.inverse() * RightMatrix;
    
    Eigen::Matrix<double, 2 * N, 2 * N> LeftMatrix;
    LeftMatrix.block(0, 0, N, N) = Vs[Lbda - lbda - 1].inverse();
    LeftMatrix.block(0, N, N, N) = Eigen::Matrix<double, N, N>::Zero();
    LeftMatrix.block(N, 0, N, N) = Eigen::Matrix<double, N, N>::Zero();
    LeftMatrix.block(N, N, N, N) = Us[Lbda - lbda - 2].inverse();
    LeftMatrix *= Vdouble.inverse();
    
    tempMatrix = LeftMatrix * Ddouble.inverse() * RightMatrix;
    G = tempMatrix.block(N, N, N, N);
    Gforward = tempMatrix.block(N, 0, N, N);
    Gbackward = tempMatrix.block(0, N, N, N);
    Gzero = tempMatrix.block(0, 0, N, N);
}

template<int N, int L, int Lbda>
Eigen::Matrix<double, N, N> Green<N, L, Lbda>::matrix()
{
    return G;
}

template<int N, int L, int Lbda>
double Green<N, L, Lbda>::get(int x, int y)
{
    return G(x, y);
}

template<int N, int L, int Lbda>
double Green<N, L, Lbda>::uneqForward(int x, int y)
{
    return Gforward(x, y);
}

template<int N, int L, int Lbda>
double Green<N, L, Lbda>::uneqBackward(int x, int y)
{
    return Gbackward(x, y);
}

template<int N, int L, int Lbda>
double Green<N, L, Lbda>::zero(int x, int y)
{
    return Gzero(x, y);
}

template<int N, int L, int Lbda>
Eigen::Matrix<double, N, N> Green<N, L, Lbda>::getM()
{
    return M;
}

#endif /* green_h */
