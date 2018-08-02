//
//  main.cpp
//
//
//  Created by Francisco Brito on 08/06/2018.
//
//  This program simulates the Hubbard model for an arbitrary geometry lattice
//  using auxiliary field (or determinant) Quantum Monte Carlo: in particular, the BSS algorithm.
//  The used notation is based on the lecture notes "Numerical Methods for Quantum Monte Carlo
//  Simulations of the Hubbard Model by Zhaojun Bai, Wenbin Chen, Richard Scalettar, and
//  Ichitaro Yamazaki (2009)
//

//  DEFAULT SIMULATION PARAMETERS FOR MINIMAL EXAMPLE.

#ifndef NSITES
#define NSITES 64 //  # sites
#endif

#ifndef DT_INV
#define DT_INV 10 //  Inverse Trotter error
#endif

#ifndef BETA
#define BETA 2  //  inverse temperature
#endif

#ifndef HOPPING
#define HOPPING 1 //  hopping parameter. set to 1 by default.
#endif

#ifndef GREEN_AFRESH_FREQ
#define GREEN_AFRESH_FREQ 4  //   how often to calculate Green's functions afresh (in # im-time slices)
#endif

#ifndef GEOM
#define GEOM 1  //   geometry (for example, 1 stands for a 1D chain with PBCs - see makefile)
#endif

#ifndef NY
#define NY 4    //   geometry parameter to define width of the simulated sample
#endif

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <random>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/KroneckerProduct>
#include <fstream>
#include <iomanip>
#include "matrixgen.h"
#include "green.h"


int main(int argc, char **argv)
{
    if ( argc != 4) //  U, mu, # sweeps
    {
        return -1;
    }

    double U = atof(argv[1]);  //  on-site interaction
    double mu = atof(argv[2]);  //  chemical potential
    int totalMCSweeps = atof(argv[3]);  //  number of sweeps

    double dt = 1. / DT_INV;  //  Trotter error, or time subinterval width. error scales as dt^2
    const int L = BETA * DT_INV;  //  # slices
    const int Lbda = L / GREEN_AFRESH_FREQ;  //  Lbda = # intervals in which the product of B's is divided to stabilize.
    double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;  //  HS transformation parameter

    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 1;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double decisionMaker; int totalMCSteps = totalMCSweeps * NSITES * L;


    // -- INITIALIZATION ---


    //  HOPPING MATRIX
    Geometry< NSITES > K;
    //  1D CHAIN PBC
    if (GEOM == 1)
    {
        K.oneDimensionalChainPBC(HOPPING, dt, mu);
    }
    //  1D CHAIN OBC
    if (GEOM == 2)
    {
        K.oneDimensionalChainOBC(HOPPING, dt, mu);
    }
    //  SQUARE LATTICE PBC
    if (GEOM == 3)
    {
        K.twoDimensionalRectanglePBC(sqrt(NSITES), HOPPING, dt, mu);
    }
    //  SQUARE LATTICE OBC
    if (GEOM == 4)
    {
        K.twoDimensionalRectangleOBC(sqrt(NSITES), HOPPING, dt, mu);
    }
    //  RECTANGULAR LATTICE PBC
    if (GEOM == 5)
    {
        K.twoDimensionalRectanglePBC(NY, HOPPING, dt, mu);
    }
    //  RECTANGULAR LATTICE OBC
    if (GEOM == 6)
    {
        K.twoDimensionalRectangleOBC(NY, HOPPING, dt, mu);
    }
    //  HONEYCOMB LATTICE PBC
    if (GEOM == 7)
    {
        std::cout << "Honeycomb" << std::endl;
    }
    //  NANORIBBON
    if (GEOM == 8)
    {
        K.nanoribbon(NY, HOPPING, dt, mu);
    }
    //  DOT
    if (GEOM == 9)
    {
        std::cout << "Dot" << std::endl;
    }

    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Configuration< L , NSITES > h; h.genHsMatrix();

    //  GENERATE THE B-MATRICES.
    OneParticlePropagators< NSITES, L > Bup; OneParticlePropagators< NSITES, L > Bdown;
    Bup.fillMatrices( true, nu, h.matrix(), K.BpreFactor() );
    Bdown.fillMatrices( false, nu, h.matrix(), K.BpreFactor() );

    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green< NSITES, L, Lbda> * Gup = new Green< NSITES, L, Lbda>;
    Green< NSITES, L, Lbda> * Gdown = new Green< NSITES, L, Lbda>;
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
    //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
//    Gup->computeGreenNaive(Bup.list(), L - 1); Gdown->computeGreenNaive(Bdown.list(), L - 1);
    Gup->storeVDU( Bup.list() ); Gdown->storeVDU( Bdown.list() );
    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();
    Gup->initializeUneqs(); Gdown->initializeUneqs();

    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double * weights = new double[totalMCSweeps * L];
    int * signs = new int[totalMCSweeps * L];
    double * electronDensities = new double[totalMCSweeps];
    double * doubleOcs = new double[totalMCSweeps];
    Eigen::Matrix<double, NSITES, NSITES> * magCorrs = new Eigen::Matrix<double, NSITES, NSITES>[totalMCSweeps];
    Eigen::Matrix<double, NSITES, NSITES> * uneqMagCorrs = new Eigen::Matrix<double, NSITES, NSITES>[totalMCSweeps];
    magCorrs[0] = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    uneqMagCorrs[0] = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    double LOGweight = 0.;
    double sign = std::copysign(1, Gup->matrix().determinant() * Gdown->matrix().determinant() );
    double electronDensity;
    double doubleOc;
    Eigen::Matrix<double, NSITES, NSITES> magCorr;
    electronDensities[0] = 0.; doubleOcs[0] = 0.;
    double meanSign = sign;

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
            std::cout << "Log Weight: " << LOGweight << std::endl << std::endl;
        }

        //  COMPUTE THE ACCEPTANCE RATIO.
        alphaUp = ( exp( -2 * h.get(l, i) * nu ) - 1 ); alphaDown = ( exp( 2 * h.get(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup->get(i, i) ) ); dDown = ( 1 + alphaDown  * ( 1 - Gdown->get(i, i) ) );
        //  SAMPLING: METROPOLIS OR HEAT BATH
        accRatio = abs( dUp * dDown ); // accRatio = dUp * dDown / ( 1 + dUp * dDown );
        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
        decisionMaker = dis(gen);

        if (decisionMaker <= accRatio )
        {
            //  KEEP TRACK OF WEIGHT
            LOGweight += log( fabs( dUp ) ) + log( fabs ( dDown ) );
            sign *= std::copysign(1, dUp * dDown );
            //  FLIP A SPIN
            h.flip(l, i);
            //  UPDATE Bs
            Bup.update(l, i, alphaUp); Bdown.update(l, i, alphaDown);
            //  RANK-ONE UPDATE -> O(N^2)
            Gup->update(alphaUp, dUp, i); Gdown->update(alphaDown, dDown, i);
        }
        //  KEEP TRACK OF THE SIGN.
        meanSign += ( sign - meanSign ) / ( step + 2 );


        // --- COMPUTE WRAPPED GREEN'S FUNCTIONS. ---


        if (i < NSITES - 1)
        {   //  CONTINUE LOOPING THROUGH THE SPATIAL LATTICE.
            i += 1;
        }
        else
        {
            //  EITHER WRAP OR COMPUTE GREEN'S FUNCTIONS FROM SCRATCH.
            latticeSweepUntilAfresh += 1;


            //  --- MEASUREMENTS ---


            //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
            weights[sweep * L + l] = LOGweight; signs[sweep * L + l] = std::copysign(1., dUp * dDown);

            //  STORE ELECTRON DENSITY, DOUBLE OCCUPANCY, AND SPIN-SPIN CORRELATIONS.
            electronDensity = 0.; doubleOc = 0.;
            for (int x = 0; x < NSITES; x++)
            {
                electronDensity -= ( Gup->get(x, x) + Gdown->get(x, x) );
                doubleOc += - Gup->get(x, x) - Gdown->get(x, x) + Gup->get(x, x)
                * Gdown->get(x, x);
                magCorr(x, x) = 3 * ( Gup->get(x, x) + Gdown->get(x, x) ) - 6 * Gup->get(x, x) * Gdown->get(x, x);
                if (l == 0)
                {
                    uneqMagCorrs[sweep](x, x) += ( 1 - Gup->zero(x, x) ) * ( 1 - Gup->zero(x, x) )
                    + ( 1 - Gup->zero(x, x) ) * ( 1 - Gdown->zero(x, x) ) + ( 1 - Gdown->zero(x, x) ) * ( 1 - Gup->zero(x, x) )
                    + ( 1 - Gdown->zero(x, x) ) * ( 1 - Gdown->zero(x, x) )
                    - Gup->zero(x, x) * Gup->zero(x, x) - Gdown->zero(x, x) * Gdown->zero(x, x)
                    - 4 * Gdown->zero(x, x) * Gup->zero(x, x);
                }
                else
                {
                    uneqMagCorrs[sweep](x, x) += ( 1 - Gup->get(x, x) ) * ( 1 - Gup->zero(x, x) )
                    + ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->zero(x, x) ) + ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->zero(x, x) )
                    + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->zero(x, x) )
                    - Gup->uneqBackward(x, x) * Gup->uneqForward(x, x) - Gdown->uneqBackward(x, x) * Gdown->uneqForward(x, x)
                    - 2 * Gdown->uneqBackward(x, x) * Gup->uneqForward(x, x) - 2 * Gup->uneqBackward(x, x) * Gdown->uneqForward(x, x);
                }
                for (int y = x + 1; y < NSITES; y++)
                {
                    // magCorr(x, y) = - Gup->get(x, y) * ( 2 * Gdown->get(y, x) + Gup->get(y, x) )
                    // - Gdown->get(x, y) * ( 2 * Gup->get(y, x) + Gdown->get(y, x) )
                    // + ( Gup->get(x, x) - Gdown->get(x, x) ) * ( Gup->get(y, y) - Gdown->get(y, y) );
                    // 2 * ( 1 - Gdown->get(y, x) ) * Gup->get(x, y)
                    // + 2 * ( 1 - Gup->get(y, x) ) * Gdown->get(x, y)
                    magCorr(x, y) =
                    - ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                    - ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->get(y, y) )
                    + ( 1 - Gup->get(x, x) ) * ( 1 - Gup->get(y, y) )
                    + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                    + ( 1 - Gup->get(y, x) ) * Gup->get(x, y)
                    + ( 1 - Gdown->get(y, x) ) * Gdown->get(x, y);
                    magCorr(y, x) = magCorr(x, y);
                    if (l == 0)
                    {
                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup->zero(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gup->zero(x, x) ) * ( 1 - Gdown->zero(y, y) ) + ( 1 - Gdown->zero(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gdown->zero(x, x) ) * ( 1 - Gdown->zero(y, y) )
                        - Gup->zero(y, x) * Gup->zero(x, y) - Gdown->zero(y, x) * Gdown->zero(x, y)
                        - 4 * Gdown->zero(y, x) * Gup->zero(x, y);
                    }
                    else
                    {
                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup->get(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->zero(y, y) ) + ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->zero(y, y) )
                        - Gup->uneqBackward(y, x) * Gup->uneqForward(x, y) - Gdown->uneqBackward(y, x) * Gdown->uneqForward(x, y)
                        - 2 * Gdown->uneqBackward(y, x) * Gup->uneqForward(x, y) - 2 * Gup->uneqBackward(y, x) * Gdown->uneqForward(x, y);
                    }
                }
                for (int y = 0; y < x; y++)
                {
                    if (l == 0)
                    {
                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup->zero(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gup->zero(x, x) ) * ( 1 - Gdown->zero(y, y) ) + ( 1 - Gdown->zero(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gdown->zero(x, x) ) * ( 1 - Gdown->zero(y, y) )
                        - Gup->zero(y, x) * Gup->zero(x, y) - Gdown->zero(y, x) * Gdown->zero(x, y)
                        - 4 * Gdown->zero(y, x) * Gup->zero(x, y);
                    }
                    else
                    {
                        uneqMagCorrs[sweep](x, y) += ( 1 - Gup->get(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->zero(y, y) ) + ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->zero(y, y) )
                        + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->zero(y, y) )
                        - Gup->uneqBackward(y, x) * Gup->uneqForward(x, y) - Gdown->uneqBackward(y, x) * Gdown->uneqForward(x, y)
                        - 2 * Gdown->uneqBackward(y, x) * Gup->uneqForward(x, y) - 2 * Gup->uneqBackward(y, x) * Gdown->uneqForward(x, y);
                    }
                }
            }
            electronDensity /= NSITES; electronDensity += 2;
            doubleOc /= NSITES; doubleOc += 1;

            electronDensities[sweep] += ( electronDensity - electronDensities[sweep] ) / ( l + 1 ) ;
            doubleOcs[sweep] +=  ( doubleOc - doubleOcs[sweep] ) / ( l + 1 ) ;
            magCorrs[sweep] += (magCorr - magCorrs[sweep] ) / ( l + 1 );



            //  DEAL WITH THE GREEN'S FUNCTIONS.


                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == GREEN_AFRESH_FREQ)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                //  Uncomment to compute the product in the naive, unstable manner
//                Gup->computeGreenNaive(Bup.list(), l); Gdown->computeGreenNaive(Bdown.list(), l);
                //  Uncomment to compute the product in the stabilized, but slightly inefficient way
//                Gup->computeStableGreenNaiveR(Bup.list(), l); Gdown->computeStableGreenNaiveR(Bdown.list(), l);
                //  Most efficient solution (storing decompositions)
                if (l != ( L - 1 ) )
                {
                    Gup->storeUDV(Bup.list(), l, GREEN_AFRESH_FREQ); Gdown->storeUDV(Bdown.list(), l, GREEN_AFRESH_FREQ);
                    //  This is the standard way described in "Stable simulations of models of interacting electrons"
//                    Gup->computeStableGreen(l, GREEN_AFRESH_FREQ); Gdown->computeStableGreen(l, GREEN_AFRESH_FREQ);
//                    //  Using the BlockOfGreens Method, we can obtain time-displaced Green's as well
                    Gup->computeBlockOfGreens(l, GREEN_AFRESH_FREQ); Gdown->computeBlockOfGreens(l, GREEN_AFRESH_FREQ);
                }
                else
                {
                    Gup->storeVDU( Bup.list() ); Gdown->storeVDU( Bdown.list() );
                    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();
                    Gup->initializeUneqs(); Gdown->initializeUneqs();
                }
                latticeSweepUntilAfresh = 0;
            }
            else
            {   //  WRAPPING.
                Gup->wrap( Bup.matrix(l) ); Gdown->wrap( Bdown.matrix(l) );
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
                    magCorrs[sweep] = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                    uneqMagCorrs[sweep] = Eigen::Matrix<double, NSITES, NSITES>::Zero();
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
        file1 << "N\t" << NSITES << '\n'; file1 << "dt\t" << dt << '\n';
        file1 << "BETA\t" << BETA << '\n'; file1 << "L\t" << L << '\n';
        file1 << "t\t" << HOPPING << '\n'; file1 << "U\t" << U << '\n';
        file1 << "mu\t" << mu << '\n'; file1 << "totalMCSweeps\t" << totalMCSweeps << '\n';
        file1 << "GREEN_AFRESH_FREQ\t" << GREEN_AFRESH_FREQ << '\n'; file1 << "Lbda\t" << Lbda << '\n';
    } file1.close();
    //  STORE MEASUREMENTS
    std::ofstream file2("plots/controlMeasurements.txt");
    std::ofstream file3("plots/measurementsScalars.txt");
    std::ofstream file4("plots/SpinSpinCorrelations.txt");
    std::ofstream file5("plots/UneqTimeSpinSpinCorrelations.txt");
    if ( file2.is_open() and file3.is_open() and file4.is_open() )
    {
        file2 << std::left << std::setw(25) << "Configuration weight";
        file2 << std::left << std::setw(25) << "Weight sign" << '\n';
        file3 << std::left << std::setw(25) << "Electron density <n>";
        file3 << std::left << std::setw(25) << "Double occupancy <n+ n->" << '\n';
        file4 << "Spin-spin correlation function <S_i S_j>" << '\n';
        file4 << std::left;
        file5 << "Unequal time spin-spin correlation function <S_i S_j>" << '\n';
        file5 << std::left;
        for (int s = 0; s < totalMCSweeps; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file2 << std::left << std::setw(25) << weights[s * L + slice];
                file2 << std::left << std::setw(25) << signs[s * L + slice] << '\n';
            }
            file3 << std::left << std::setw(25) << std::setprecision(10) << electronDensities[s];
            file3 << std::left << std::setw(25) << std::setprecision(10) << doubleOcs[s] << '\n';
            file4 << std::setprecision(10) << magCorrs[s] << std::endl << std::endl;
            file5 << std::setprecision(10) << uneqMagCorrs[s] << std::endl << std::endl;
        }
        file2 << '\n'; file3 << '\n'; file4 << '\n'; file5 << '\n';
    }
    file2.close(); file3.close(); file4.close(); file5.close();

    delete[] weights; delete[] doubleOcs; delete[] electronDensities; delete[] magCorrs; delete[] uneqMagCorrs;
    delete Gup; delete Gdown;

    return 0;
}
