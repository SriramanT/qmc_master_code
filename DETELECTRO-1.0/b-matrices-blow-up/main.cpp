//
//  main.cpp
//  
//
//  Created by Francisco Brito on 08/06/2018.
//
//  This program simulates the Hubbard model on a 1D chain
//  using auxiliary field (or determinant) Quantum Monte Carlo: in particular, the BSS algorithm.
//  The used notation is based on the lecture notes "Numerical Methods for Quantum Monte Carlo
//  Simulations of the Hubbard Model by Zhaojun Bai, Wenbin Chen, Richard Scalettar, and
//  Ichitaro Yamazaki (2009)
//

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
#include <iomanip>

#include "matrixgen.h"
#include "green.h"


int main()
{
    //  SET SIMULATION PARAMETERS.
    const int N = 4;  //  # sites
    const int dtInv = 8;    //  Inverse Trotter error
    const double dt = 1. / dtInv;  //  Trotter error, or time subinterval width. error scales as dt^2
    const int beta = 2;  //  inverse temperature
    const int t = 1;  //  hopping parameter. set to 1 by default.
    const double U = 4.;  //  on-site interaction
    const double mu = 0;  //  chemical potential
    const int greenAfreshFreq = 4;  //   how often to calculate Green's functions afresh (in # im-time slices)
    const int L = beta * dtInv;  //  # slices
    const double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;  //  HS transformation parameter
    const int Lbda = L / greenAfreshFreq;  //  Lbda = # intervals in which the product of B's is divided to stabilize.
    
    //  DEFINE MATRICES TYPE.
    typedef Eigen::Matrix<double, N, N> MatrixNN;
    
    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 1;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double decisionMaker;
    const int totalMCSweeps = 1;
    const int totalMCSteps = totalMCSweeps * N * L;

    
    // -- INITIALIZATION ---
    
    
    //  HOPPING MATRIX FOR ARBITRARY GEOMETRY
    Geometry< N > K; K.oneDimensionalChainPBC();
    
    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Configuration< L , N > h; h.genHsMatrix();
    
    //  COMPUTE MATRIX 'PREFACTOR' OF THE B-MATRICES B = e^(t dt K).
    const MatrixNN BpreFactor = ( t * dt * K.matrix() ).exp();

    //  GENERATE THE B-MATRICES.
    OneParticlePropagators< N, L > Bup; OneParticlePropagators< N, L > Bdown;
    Bup.fillMatrices(true, nu, dt, mu, h.matrix(), BpreFactor); Bdown.fillMatrices(false, nu, dt, mu, h.matrix(), BpreFactor);

    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green< N, L, Lbda> Gup; Green< N, L, Lbda> Gdown;
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
    //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    Gup.computeGreenNaive(Bup.list(), L - 1); Gdown.computeGreenNaive(Bdown.list(), L - 1);
//    Gup.storeVDU( Bup.list() ); Gdown.storeVDU( Bdown.list() );
//    Gup.computeGreenFromVDU(); Gdown.computeGreenFromVDU();

    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double weight = 1.;
    double meanSign = std::copysign( 1. , Gup.matrix().determinant() )
                    * std::copysign( 1. , Gdown.matrix().determinant() );
    double avCondNumber = 0.;
    double cond[L];

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE, THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0; int i = 0; int latticeSweepUntilAfresh = 0; int sweep = 0; int step;


    // --- MC LOOP ---


    for (step = 0; step < totalMCSteps; step++)
    {


        //  COMPUTE THE ACCEPTANCE RATIO.
        alphaUp = ( exp( -2 * h.get(l, i) * nu ) - 1 ); alphaDown = ( exp( 2 * h.get(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup.get(i, i) ) ); dDown = ( 1 + alphaDown  * ( 1 - Gdown.get(i, i) ) );
        accRatio = dUp * dDown; //  Regular Metropolis

        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
        decisionMaker = dis(gen);

        if (decisionMaker <= abs(accRatio) )
        {
            //  KEEP TRACK OF WEIGHT
            weight = dUp * dDown * weight;

            //  FLIP A SPIN
            h.flip(l, i);

            //  UPDATE Bs
            Bup.update(l, i, alphaUp); Bdown.update(l, i, alphaDown);

            //  RANK-ONE UPDATE -> O(N^2)
            Gup.update(alphaUp, dUp, i); Gdown.update(alphaDown, dDown, i);
        }

        //  KEEP TRACK OF THE SIGN (RUNNING AVERAGE).
        meanSign += ( std::copysign(1., weight) - meanSign ) / ( step + 2 );


        // --- COMPUTE WRAPPED GREEN'S FUNCTIONS. ---


        if (i < N - 1)
        {   //  CONTINUE LOOPING THROUGH THE SPATIAL LATTICE
            i += 1;
        }
        else
        {
            //  EITHER WRAP OR COMPUTE GREEN'S FUNCTIONS FROM SCRATCH.
            latticeSweepUntilAfresh += 1;

            //  DEAL WITH THE GREEN'S FUNCTIONS.
                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == greenAfreshFreq)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                //  Uncomment to compute the product in the naive, unstable manner
//                Gup.computeGreenNaive(Bup.list(), l); Gdown.computeGreenNaive(Bdown.list(), l);
                //  Uncomment to compute the product in the stabilized, but slightly inefficient way
//                Gup.computeStableGreenNaiveR(Bup.list(), l); Gdown.computeStableGreenNaiveR(Bdown.list(), l);
                //  Most efficient solution (storing decompositions)
                if (l != ( L - 1 ) )
                {
                    Gup.storeUDV(Bup.list(), l, greenAfreshFreq); Gdown.storeUDV(Bdown.list(), l, greenAfreshFreq);
                    //  This is the standard way described in "Stable simulations of models of interacting electrons"
                    Gup.computeStableGreen(l, greenAfreshFreq); Gdown.computeStableGreen(l, greenAfreshFreq);
                }
                else
                {
                    Gup.storeVDU( Bup.list() ); Gdown.storeVDU( Bdown.list() );
                    Gup.computeGreenFromVDU(); Gdown.computeGreenFromVDU();
                }
                latticeSweepUntilAfresh = 0;
            }
            else
            {   //  WRAPPING.
                Gup.wrap( Bup.matrix(l) ); Gdown.wrap( Bdown.matrix(l) );
            }
            if (l < L - 1)
            {
                l += 1; i = 0;
            }
            else
            {
                l = 0; i = 0;
            }
        }
    }   //  END OF MC LOOP.

    Eigen::Matrix<double, N, N> partialProdBs = Eigen::Matrix<double, N, N>::Identity();
    
    for (int slice = L - 1; slice >= 0; slice--)
    {
        partialProdBs *= Bup.matrix(slice);
        Eigen::JacobiSVD< Eigen::Matrix<double, N, N> > svd( partialProdBs );
        cond[slice] = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size()-1);
    }
    
    Eigen::Matrix<double, N * L, N * L> HubbardMat = Eigen::Matrix<double, N * L, N * L>::Zero() ;
    
    HubbardMat.block(0, 0 ,N, N) =  Eigen::Matrix<double, N, N>::Identity();
    HubbardMat.block(0, N * ( L - 1 ), N, N) = Bup.matrix(0);
    for (int slice = 1; slice < L; slice++)
    {
        HubbardMat.block(N * slice, N * slice ,N, N) =  Eigen::Matrix<double, N, N>::Identity();
        HubbardMat.block(N * slice, N * (slice - 1), N, N) = - Bup.matrix(slice);
    }

//    return avCondNumber;
    //  SAVE OUTPUT.
    std::ofstream file1("conds.txt");
    if (file1.is_open())
    {
        for (int slice = L - 1; slice >= 0; slice--)
            file1 << cond[slice] << '\n';
    }
    std::ofstream file2("hubbMat.txt");
    if (file2.is_open())
    {
        file2 << HubbardMat << '\n';
    }
    return 0;
}

//int main()
//{
//    std::cout << simulation() << std::endl;
//    return 0;
//}

