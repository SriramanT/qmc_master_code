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
    const int N = 10;  //  # sites
    const int dtInv = 8;    //  Inverse Trotter error
    const double dt = 1. / dtInv;  //  Trotter error, or time subinterval width. error scales as dt^2
    const int beta = 30;  //  inverse temperature
    const int t = 1;  //  hopping parameter. set to 1 by default.
    const double U = 4.;  //  on-site interaction
    const double mu = 0.;  //  chemical potential
    const int greenAfreshFreq = 2;  //   how often to calculate Green's functions afresh (in # im-time slices)
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
    const int totalMCSweeps = 100;
    const int totalMCSteps = totalMCSweeps * N * L;

    
    // -- INITIALIZATION ---
    
    
    //  HOPPING MATRIX FOR ARBITRARY GEOMETRY
    Geometry< N > K; K.oneDimensionalChainPBC();
    
    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Configuration< L , N > h; h.genHsMatrix();

    //  COMPUTE MATRIX 'PREFACTOR' OF THE B-MATRICES B = e^(t dt K).
    const MatrixNN BpreFactor =  exp(dt * mu) * Eigen::Matrix<double, N, N>::Identity()  * ( t * dt * K.matrix() ).exp();

    //  GENERATE THE B-MATRICES.
    OneParticlePropagators< N, L > Bup; OneParticlePropagators< N, L > Bdown;
    Bup.fillMatrices(true, nu, h.matrix(), BpreFactor); Bdown.fillMatrices(false, nu, h.matrix(), BpreFactor);

    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green< N, L, Lbda> Gup; Green< N, L, Lbda> Gdown;
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
    //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
//    Gup.computeGreenNaive(Bup.list(), L - 1); Gdown.computeGreenNaive(Bdown.list(), L - 1);
    Gup.storeVDU( Bup.list() ); Gdown.storeVDU( Bdown.list() );
    Gup.computeGreenFromVDU(); Gdown.computeGreenFromVDU();
    Gup.initializeUneqs(); Gdown.initializeUneqs();

    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double * weights = new double[totalMCSweeps * L];
    double * electronDensities = new double[totalMCSweeps];
    double * doubleOcs = new double[totalMCSweeps];
    MatrixNN * magCorrs = new MatrixNN[totalMCSweeps];
//    MatrixNN * uneqMagCorrs = new MatrixNN[totalMCSweeps];
    magCorrs[0] = Eigen::Matrix<double, N, N>::Zero();
//    uneqMagCorrs[0] = Eigen::Matrix<double, N, N>::Zero();
    double weight = 0.;
    double electronDensity;
    double doubleOc;
    MatrixNN magCorr;
    electronDensities[0] = 0.; doubleOcs[0] = 0.;
    double meanSign = 1;

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE, THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0; int i = 0; int latticeSweepUntilAfresh = 0; int sweep = 0; int step;


    // --- MC LOOP ---


    std::cout << "\nMC loop started. Progress:\n";
    for (step = 0; step < totalMCSteps; step++)
    {
        //  DISPLAY PROGRESS OF THE RUN.
        if ( (step + 1)  % (totalMCSteps/8) == 0 )
        {
            std::cout << (step + 1) * 1. / totalMCSteps * 100 << " %" << std::endl << std::endl;
            std::cout << "Average Sign: " << meanSign << std::endl << std::endl;
            std::cout << "Log Weight: " << weight << std::endl << std::endl;
        }

        //  COMPUTE THE ACCEPTANCE RATIO.
        alphaUp = ( exp( -2 * h.get(l, i) * nu ) - 1 ); alphaDown = ( exp( 2 * h.get(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup.get(i, i) ) ); dDown = ( 1 + alphaDown  * ( 1 - Gdown.get(i, i) ) );
        accRatio = dUp * dDown; //  Regular Metropolis
//        accRatio = dUp * dDown / ( 1 + dUp * dDown );  // Heat Bath

        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
        decisionMaker = dis(gen);

        if (decisionMaker <= abs(accRatio) )
        {
            //  KEEP TRACK OF WEIGHT
            weight += log( abs( dUp ) ) + log( abs ( dDown ) );

            //  FLIP A SPIN
            h.flip(l, i);

            //  UPDATE Bs
            Bup.update(l, i, alphaUp); Bdown.update(l, i, alphaDown);

            //  RANK-ONE UPDATE -> O(N^2)
            Gup.update(alphaUp, dUp, i); Gdown.update(alphaDown, dDown, i);
        }

        //  KEEP TRACK OF THE SIGN (RUNNING AVERAGE).
        meanSign += ( std::copysign(1., dUp * dDown) - meanSign ) / ( step + 2 );


        // --- COMPUTE WRAPPED GREEN'S FUNCTIONS. ---


        if (i < N - 1)
        {   //  CONTINUE LOOPING THROUGH THE SPATIAL LATTICE
            i += 1;
        }
        else
        {
            //  EITHER WRAP OR COMPUTE GREEN'S FUNCTIONS FROM SCRATCH.
            latticeSweepUntilAfresh += 1;


            //  --- MEASUREMENTS ---


            //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
            weights[sweep * L + l] = weight;
            
            //  STORE ELECTRON DENSITY, DOUBLE OCCUPANCY, AND SPIN-SPIN CORRELATIONS.
            electronDensity = 0.; doubleOc = 0.;
            for (int x = 0; x < N; x++)
            {
                electronDensity -= ( Gup.get(x, x) + Gdown.get(x, x) );
                doubleOc = doubleOc - Gup.get(x, x) - Gdown.get(x, x) + Gup.get(x, x) * Gdown.get(x, x);
                magCorr(x, x) = 3 * ( Gup.get(x, x) + Gdown.get(x, x) ) - 6 * Gup.get(x, x) * Gdown.get(x, x);
//                if (l == 0)
//                {
//                    uneqMagCorrs[sweep](x, x) += ( 1 - Gup.zero(x, x) ) * ( 1 - Gup.zero(x, x) )
//                    + ( 1 - Gup.zero(x, x) ) * ( 1 - Gdown.zero(x, x) ) + ( 1 - Gdown.zero(x, x) ) * ( 1 - Gup.zero(x, x) )
//                    + ( 1 - Gdown.zero(x, x) ) * ( 1 - Gdown.zero(x, x) )
//                    - Gup.zero(x, x) * Gup.zero(x, x) - Gdown.zero(x, x) * Gdown.zero(x, x)
//                    - 4 * Gdown.zero(x, x) * Gup.zero(x, x);
//                }
//                else
//                {
//                    uneqMagCorrs[sweep](x, x) += ( 1 - Gup.get(x, x) ) * ( 1 - Gup.zero(x, x) )
//                    + ( 1 - Gup.get(x, x) ) * ( 1 - Gdown.zero(x, x) ) + ( 1 - Gdown.get(x, x) ) * ( 1 - Gup.zero(x, x) )
//                    + ( 1 - Gdown.get(x, x) ) * ( 1 - Gdown.zero(x, x) )
//                    - Gup.uneqBackward(x, x) * Gup.uneqForward(x, x) - Gdown.uneqBackward(x, x) * Gdown.uneqForward(x, x)
//                    - 2 * Gdown.uneqBackward(x, x) * Gup.uneqForward(x, x) - 2 * Gup.uneqBackward(x, x) * Gdown.uneqForward(x, x);
//                }
                for (int y = 0; y < x; y++)
                {
                    magCorr(x, y) = - Gup.get(x, y) * ( 2 * Gdown.get(y, x) + Gup.get(y, x) )
                    - Gdown.get(x, y) * ( 2 * Gup.get(y, x) + Gdown.get(y, x) )
                    + ( Gup.get(x, x) - Gdown.get(x, x) ) * ( Gup.get(y, y) - Gdown.get(y, y) );
                    magCorr(y, x) = magCorr(x, y);
//                    if (l == 0)
//                    {
//                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup.zero(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gup.zero(x, x) ) * ( 1 - Gdown.zero(y, y) ) + ( 1 - Gdown.zero(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gdown.zero(x, x) ) * ( 1 - Gdown.zero(y, y) )
//                        - Gup.zero(y, x) * Gup.zero(x, y) - Gdown.zero(y, x) * Gdown.zero(x, y)
//                        - 4 * Gdown.zero(y, x) * Gup.zero(x, y);
//                    }
//                    else
//                    {
//                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup.get(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gup.get(x, x) ) * ( 1 - Gdown.zero(y, y) ) + ( 1 - Gdown.get(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gdown.get(x, x) ) * ( 1 - Gdown.zero(y, y) )
//                        - Gup.uneqBackward(y, x) * Gup.uneqForward(x, y) - Gdown.uneqBackward(y, x) * Gdown.uneqForward(x, y)
//                        - 2 * Gdown.uneqBackward(y, x) * Gup.uneqForward(x, y) - 2 * Gup.uneqBackward(y, x) * Gdown.uneqForward(x, y);
//                    }
                }
//                for (int y = x + 1; y < N; y++)
//                {
//                    if (l == 0)
//                    {
//                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup.zero(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gup.zero(x, x) ) * ( 1 - Gdown.zero(y, y) ) + ( 1 - Gdown.zero(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gdown.zero(x, x) ) * ( 1 - Gdown.zero(y, y) )
//                        - Gup.zero(y, x) * Gup.zero(x, y) - Gdown.zero(y, x) * Gdown.zero(x, y)
//                        - 4 * Gdown.zero(y, x) * Gup.zero(x, y);
//                    }
//                    else
//                    {
//                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup.get(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gup.get(x, x) ) * ( 1 - Gdown.zero(y, y) ) + ( 1 - Gdown.get(x, x) ) * ( 1 - Gup.zero(y, y) )
//                        + ( 1 - Gdown.get(x, x) ) * ( 1 - Gdown.zero(y, y) )
//                        - Gup.uneqBackward(y, x) * Gup.uneqForward(x, y) - Gdown.uneqBackward(y, x) * Gdown.uneqForward(x, y)
//                        - 2 * Gdown.uneqBackward(y, x) * Gup.uneqForward(x, y) - 2 * Gup.uneqBackward(y, x) * Gdown.uneqForward(x, y);
//                    }
//                }
            }
            electronDensity /= N; electronDensity += 2;
            doubleOc /= N; doubleOc += 1;

            electronDensities[sweep] += ( electronDensity - electronDensities[sweep] ) / ( l + 1 ) ;
            doubleOcs[sweep] +=  ( doubleOc - doubleOcs[sweep] ) / ( l + 1 ) ;
            magCorrs[sweep] += (magCorr - magCorrs[sweep] ) / ( l + 1 );



            //  DEAL WITH THE GREEN'S FUNCTIONS.


                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == greenAfreshFreq)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                //  Uncomment to compute the product in the naive, unstable manner
//                Gup.computeGreenNaive(Bup.list(), l); Gdown.computeGreenNaive(Bdown.list(), l);
                //  Uncomment to compute the product in the stabilized, but slightly inefficient way
                Gup.computeStableGreenNaiveR(Bup.list(), l); Gdown.computeStableGreenNaiveR(Bdown.list(), l);
                //  Most efficient solution (storing decompositions)
//                if (l != ( L - 1 ) )
//                {
//                    Gup.storeUDV(Bup.list(), l, greenAfreshFreq); Gdown.storeUDV(Bdown.list(), l, greenAfreshFreq);
//                    //  This is the standard way described in "Stable simulations of models of interacting electrons"
////                    Gup.computeStableGreen(l, greenAfreshFreq); Gdown.computeStableGreen(l, greenAfreshFreq);
//                    //  Using the BlockOfGreens Method, we can obtain time-displaced Green's as well
//                    Gup.computeBlockOfGreens(l, greenAfreshFreq); Gdown.computeBlockOfGreens(l, greenAfreshFreq);
//                }
//                else
//                {
//                    Gup.storeVDU( Bup.list() ); Gdown.storeVDU( Bdown.list() );
//                    Gup.computeGreenFromVDU(); Gdown.computeGreenFromVDU();
//                    Gup.initializeUneqs(); Gdown.initializeUneqs();
//                }
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
                if (sweep != ( totalMCSweeps - 1 ) )
                {
                    //  MOVE SWEEP COUNTER
                    sweep += 1;
                    electronDensities[sweep] = 0.; doubleOcs[sweep] = 0.;
                    magCorrs[sweep] = Eigen::Matrix<double, N, N>::Zero();
//                    uneqMagCorrs[sweep] = Eigen::Matrix<double, N, N>::Zero();
                }
                l = 0; i = 0;
            }
        }
    }   //  END OF MC LOOP.

    std::cout << "Simulation ended" << std::endl << std::endl;
    std::cout << "Average Sign: " << meanSign << std::endl << std::endl;

    //  SAVE OUTPUT.
    std::ofstream file1("plots/simulationParameters.txt");
    if (file1.is_open())
    {
        file1 << "N\t" << N << '\n';
        file1 << "dt\t" << dt << '\n';
        file1 << "beta\t" << beta << '\n';
        file1 << "L\t" << L << '\n';
        file1 << "t\t" << t << '\n';
        file1 << "U\t" << U << '\n';
        file1 << "mu\t" << mu << '\n';
        file1 << "totalMCSweeps\t" << totalMCSweeps << '\n';
        file1 << "greenAfreshFreq\t" << greenAfreshFreq << '\n';
        file1 << "Lbda\t" << Lbda << '\n';
    }
    file1.close();
    //  STORE MEASUREMENTS
    std::ofstream file2("plots/controlMeasurements.txt");
    std::ofstream file3("plots/measurementsScalars.txt");
    std::ofstream file4("plots/SpinSpinCorrelations.txt");
//    std::ofstream file5("plots/UneqTimeSpinSpinCorrelations.txt");
    if ( file2.is_open() and file3.is_open() and file4.is_open() )
    {
        file2 << std::left << std::setw(25) << "Configuration weight";
        file2 << std::left << std::setw(25) << "Weight sign" << '\n';
        file3 << std::left << std::setw(25) << "Electron density <n>";
        file3 << std::left << std::setw(25) << "Double occupancy <n+ n->" << '\n';
        file4 << "Spin-spin correlation function <S_i S_j>" << '\n';
        file4 << std::left;
//        file5 << "Unequal time spin-spin correlation function <S_i S_j>" << '\n';
//        file5 << std::left;
        for (int s = 0; s < totalMCSweeps; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file2 << std::left << std::setw(25) << weights[s * L + slice];
                file2 << std::left << std::setw(25) << std::copysign( 1. , weights[s * L + slice] ) << '\n';
            }
            file3 << std::left << std::setw(25) << std::setprecision(10) << electronDensities[s];
            file3 << std::left << std::setw(25) << std::setprecision(10) << doubleOcs[s] << '\n';
            file4 << std::setprecision(10) << magCorrs[s] << std::endl << std::endl;
//            file5 << std::setprecision(10) << uneqMagCorrs[s] << std::endl << std::endl;
        }
        file2 << '\n';
        file3 << '\n';
        file4 << '\n';
//        file5 << '\n';
    }
    file2.close();
    file3.close();
    file4.close();
//    file5.close();

    delete[] weights; delete[] doubleOcs; delete[] electronDensities; delete[] magCorrs;
//    delete[] uneqMagCorrs;

    return 0;
}

