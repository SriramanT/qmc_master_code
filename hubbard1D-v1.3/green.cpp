//
//  green.cpp
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "green.h"
#include "UDV.h"
#include "VDU.h"

Green::Green(int nSites, int nSlices)
{
    N = nSites;
    L = nSlices;
}

void Green::computeGreenNaive(Eigen::MatrixXd* Bs, int l)
{
    M = Eigen::MatrixXd::Identity(N, N);
    int slice;
    
    //  SET UP THE ITERATOR
    
    std::vector<int> sliceVector;
    std::vector<int>::iterator sliceIterator;
    for (slice = l; slice >= 0; slice--)
    {
        sliceVector.push_back(slice);
    }
    
    for (slice = L - 1; slice > l; slice--)
    {
        sliceVector.push_back(slice);
    }
    
    //  COMPUTE THE PRODUCT
    
    for (sliceIterator = sliceVector.begin(); sliceIterator != sliceVector.end(); sliceIterator++ )
    {
        M = M * Bs[*sliceIterator];
    }
    M += Eigen::MatrixXd::Identity(N, N);
    G = M.inverse();
}

void Green::computeStableGreenNaiveR(Eigen::MatrixXd* Bs, int l, int Lbda)
{
    U = Eigen::MatrixXd::Identity(N, N);
    D = Eigen::MatrixXd::Identity(N, N);
    V = Eigen::MatrixXd::Identity(N, N);
    
    int site;
    int slice;
    int sliceCounter = 0;

    //  INITIALIZE PARTIAL PRODUCTS.
    Eigen::MatrixXd partialProdBs = Eigen::MatrixXd::Identity(N, N);  //  an array of partial products

    //  SET UP THE ITERATOR.
    std::vector<int> sliceVector;
    std::vector<int>::iterator sliceIterator;
    for (slice = l + 1; slice < L; slice++)
    {
        sliceVector.push_back(slice);
    }

    for (slice = 0; slice <= l; slice++)
    {
        sliceVector.push_back(slice);
    }

    //  COMPUTE UDVs
    for (sliceIterator = sliceVector.begin(); sliceIterator != sliceVector.end(); sliceIterator++ )
    {
        partialProdBs = Bs[*sliceIterator] * partialProdBs;
        sliceCounter += 1;
        if (sliceCounter == L / Lbda)
        {
            //  Householder QR applied to (previous partial product) * U * D
            //  where U, and D are from the current partial product's decomposition
            UDV udvRight(partialProdBs * U * D);
            U = udvRight.QR_and_getU();
            D = udvRight.getD();
            V = udvRight.getV() * V;
            sliceCounter = 0;
            partialProdBs = Eigen::MatrixXd::Identity(N, N);
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

void Green::computeStableGreenNaiveL(Eigen::MatrixXd* Bs, int l, int Lbda)
{
    U = Eigen::MatrixXd::Identity(N, N);
    D = Eigen::MatrixXd::Identity(N, N);
    V = Eigen::MatrixXd::Identity(N, N);
    
    int site;
    int slice;
    int sliceCounter = 0;
    
    //  INITIALIZE PARTIAL PRODUCTS.
    Eigen::MatrixXd partialProdBs = Eigen::MatrixXd::Identity(N, N);  //  an array of partial products
    
    //  SET UP THE ITERATOR.
    std::vector<int> sliceVector;
    std::vector<int>::iterator sliceIterator;
    for (slice = l; slice >= 0; slice--)
    {
        sliceVector.push_back(slice);
    }
    
    for (slice = L - 1; slice > l; slice--)
    {
        sliceVector.push_back(slice);
    }
    
    //  COMPUTE UDVs
    for (sliceIterator = sliceVector.begin(); sliceIterator != sliceVector.end(); sliceIterator++ )
    {
        partialProdBs *= Bs[*sliceIterator];
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
            partialProdBs = Eigen::MatrixXd::Identity(N, N);
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

void Green::computeGreenFromVDU(Eigen::MatrixXd VsLast, Eigen::MatrixXd DsLast, Eigen::MatrixXd UsLast)
{
    U = UsLast;
    D = DsLast;
    V = VsLast;
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
}


void Green::storeVDU(Eigen::MatrixXd* Bs, int Lbda, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs)
{
    //  This is very similar to the function above. The only difference is that the
    //  Green's function is not computed, and this is to be used at the beginning of
    //  each imaginary time sweep, i.e. M = I + B_{L-1} B_{L-2}... B_{0}
    U = Eigen::MatrixXd::Identity(N, N);
    D = Eigen::MatrixXd::Identity(N, N);
    V = Eigen::MatrixXd::Identity(N, N);
    
    int slice;
    int sliceCounter = 0;
    int lbda = 0;
    
    //  INITIALIZE PARTIAL PRODUCTS.
    Eigen::MatrixXd partialProdBs = Eigen::MatrixXd::Identity(N, N);  //  an array of partial products
    
    //  SET UP THE ITERATOR.
    std::vector<int> sliceVector;
    std::vector<int>::iterator sliceIterator;
    for (slice = L - 1; slice >= 0; slice--)
    {
        sliceVector.push_back(slice);
    }
    
    //  COMPUTE UDVs
    for (sliceIterator = sliceVector.begin(); sliceIterator != sliceVector.end(); sliceIterator++ )
    {
        partialProdBs *= Bs[*sliceIterator];
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
            partialProdBs = Eigen::MatrixXd::Identity(N, N);
        }
    }
}

void Green::storeUDV(Eigen::MatrixXd* Bs, int l, int Lbda, int greenAfreshFreq, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs)
{
    //  NOTE THAT THE ARGUMENT l will never be zero, since at l = 0, we compute the VDUs from the right and the Green's function
    int lbda = (l + 1) / greenAfreshFreq - 1;
    int slice;
    Eigen::MatrixXd partialProdBs = Eigen::MatrixXd::Identity(N, N);  //  an array of partial products
    
    //  INITIALIZE PARTIAL PRODUCTS.
    if ( lbda == 0 )
    {
        U = Eigen::MatrixXd::Identity(N, N);
        D = Eigen::MatrixXd::Identity(N, N);
        V = Eigen::MatrixXd::Identity(N, N);
    }
    else
    {
        U = Us[Lbda - lbda] ;
        D = Ds[Lbda - lbda] ;
        V = Vs[Lbda - lbda] ;
    }
    
    //  SET UP THE ITERATOR.
    std::vector<int> sliceVector;
    std::vector<int>::iterator sliceIterator;
    for (slice = ( l - (greenAfreshFreq - 1) ) ; slice <= l ; slice++)
    {
        sliceVector.push_back(slice);
    }
    //  COMPUTE UDVs.
    for (sliceIterator = sliceVector.begin(); sliceIterator != sliceVector.end(); sliceIterator++ )
    {
        partialProdBs = Bs[*sliceIterator] * partialProdBs;
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

void Green::computeStableGreen(int l, int Lbda, int greenAfreshFreq, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs)
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

void Green::computeBlockOfGreens(int l, int Lbda, int greenAfreshFreq, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs)
{
    int lbda = (l + 1) / greenAfreshFreq - 1;
    Eigen::MatrixXd tempMatrix(2 * N, 2 * N);
    tempMatrix.block(0, 0, N, N) = ( Vs[Lbda - lbda - 1] * Vs[Lbda - lbda - 2] ).inverse();
    tempMatrix.block(0, N, N, N) = Ds[Lbda - lbda - 2];
    tempMatrix.block(N, 0, N, N) = - Ds[Lbda - lbda - 1];
    tempMatrix.block(N, N, N, N) = ( Us[Lbda - lbda - 2] * Us[Lbda - lbda - 1] ).inverse();
    UDV decomposition(tempMatrix);
    U = decomposition.QR_and_getU();
    D = decomposition.getD();
    V = decomposition.getV();
    Eigen::MatrixXd RightMatrix(2 * N, 2 * N);
    RightMatrix.block(0, 0, N, N) = Vs[Lbda - lbda - 2].inverse();
    RightMatrix.block(0, N, N, N) = Eigen::MatrixXd::Zero(N, N);
    RightMatrix.block(N, 0, N, N) = Eigen::MatrixXd::Zero(N, N);
    RightMatrix.block(N, N, N, N) = Us[Lbda - lbda - 1].inverse();
    RightMatrix = U.inverse() * RightMatrix;
    Eigen::MatrixXd LeftMatrix(2 * N, 2 * N);
    LeftMatrix.block(0, 0, N, N) = Vs[Lbda - lbda - 1].inverse();
    LeftMatrix.block(0, N, N, N) = Eigen::MatrixXd::Zero(N, N);
    LeftMatrix.block(N, 0, N, N) = Eigen::MatrixXd::Zero(N, N);
    LeftMatrix.block(N, N, N, N) = Us[Lbda - lbda - 2].inverse();
    LeftMatrix *= V.inverse();
    tempMatrix = LeftMatrix * D.inverse() * RightMatrix;
    G = tempMatrix.block(N, N, N, N);
    Gforward = tempMatrix.block(N, 0, N, N);
    Gbackward = tempMatrix.block(0, N, N, N);
}

Eigen::MatrixXd Green::getG()
{
    return G;
}

Eigen::MatrixXd Green::getGforward()
{
    return Gforward;
}

Eigen::MatrixXd Green::getGbackward()
{
    return Gbackward;
}

Eigen::MatrixXd Green::getM()
{
    return M;
}

void Green::printGreen(bool spin)
{
    if (spin == true)
    {
        std::cout << "\nGreenUp\n\n" << G << std::endl;
    }
    else
    {
        std::cout << "\nGreenDown\n\n" << G << std::endl;
    }
}

void Green::printM(bool spin)
{
    if (spin == true)
    {
        std::cout << "\nMUp\n\n" << M << std::endl;
    }
    else
    {
        std::cout << "\nMDown\n\n" << M << std::endl;
    }
}
