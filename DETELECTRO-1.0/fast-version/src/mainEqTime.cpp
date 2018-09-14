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
//  Compare with J. E. Hirsch - Phys Rev B 28 7, 1983
//  For U = 4, we get <n_up n_dw > -> 0.1384 (exact)

#ifndef NSITES
#define NSITES 2 //  # sites
#endif

#ifndef DT_INV
#define DT_INV 16 //  Inverse Trotter error
#endif

#ifndef BETA
#define BETA 2  //  inverse temperature
#endif

#ifndef GREEN_AFRESH_FREQ
#define GREEN_AFRESH_FREQ 4  //   how often to calculate Green's functions afresh (in # im-time slices)
#endif

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include "unsupported/Eigen/KroneckerProduct"
#include "matrixgen.h"
#include "green.h"


int main(int argc, char **argv)
{
    if ( argc != 9) //  U, mu, # sweeps, # warm-up sweeps, geometry
    {
        return -1;
    }

    double t = atof(argv[1]);  //  tight binding parameter
    double U = atof(argv[2]);  //  on-site interaction
    double mu = atof(argv[3]);  //  chemical potential
    int geom = atof(argv[4]);  //   geometry (for example, 1 stands for a 1D chain with PBCs - see makefile)
    int Ny = atof(argv[5]);    //   geometry parameter to define width of the simulated sample
    int totalMCSweeps = atof(argv[6]);  //  number of sweeps
    int W = atof(argv[7]);  //  number of warm-up sweeps
    int A = atof(argv[8]);  //  number of auto-correlation sweeps

    double dt = 1. / DT_INV;  //  Trotter error, or time subinterval width. error scales as dt^2
    const int L = BETA * DT_INV;  //  # slices
    //  Lbda = # intervals in which the product of B's is divided to stabilize.
    const int Lbda = L / GREEN_AFRESH_FREQ;
    //  HS transformation parameter (to order dtau^2)
    double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;

    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 1;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double decisionMaker; int totalMCSteps = totalMCSweeps * NSITES * L;


    // -- INITIALIZATION ---


    //  HOPPING MATRIX
    Geometry< NSITES > K;
    //  1D CHAIN PBC
    if (geom == 1)
    {
        K.oneDimensionalChainPBC(t, dt, mu);
    }
    //  1D CHAIN OBC
    if (geom == 2)
    {
        K.oneDimensionalChainOBC(t, dt, mu);
    }
    //  SQUARE LATTICE PBC
    if (geom == 3)
    {
        K.twoDimensionalRectanglePBC(sqrt(NSITES), t, dt, mu);
    }
    //  SQUARE LATTICE OBC
    if (geom == 4)
    {
        K.twoDimensionalRectangleOBC(sqrt(NSITES), t, dt, mu);
    }
    //  RECTANGULAR LATTICE PBC
    if (geom == 5)
    {
        K.twoDimensionalRectanglePBC(Ny, t, dt, mu);
    }
    //  RECTANGULAR LATTICE OBC
    if (geom == 6)
    {
        K.twoDimensionalRectangleOBC(Ny, t, dt, mu);
    }
    //  TRIANGULAR LATTICE PBC
    if (geom == 6)
    {

    }
    //  HONEYCOMB LATTICE PBC
    if (geom == 9)
    {
        K.hcPBC();
        K.computeExponential(t, dt);
    }
    //  NANORIBBON
    if (geom == 10)
    {
        K.hcNanoribbon(Ny);
        K.computeExponential(t, dt);
    }
    //  3 ORBITAL TIGHT BIDING MODEL ON THE M-ATOM TRIANGULAR LATTICE
    //  MODEL OF A TMD NANORIBBON (SEE Liu2013)
    //  MoS2
    if (geom == 14)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {1.046, 2.104, 0.401, 0.507, 0.218, 0.338, 0.057};
        K.setParamsThreeOrbitalTB(params);
        K.tmdPBC();

        // FOR DEBUGGING: TEST MATRIX CREATED BY THE PROGRAM IN OTHER CODES
        std::ofstream TMDhopping("temp-data/tmd-hopping.csv");
        if (TMDhopping.is_open())
        {
            TMDhopping << K.getB() << '\n';
        }
        TMDhopping.close();
        K.computeExponential(t, dt);
    }
    if (geom == 15)
    {
        if (NSITES % 3 != 0)
        {
            std::cout << "Invalid number of sites (real + orbital spaces)." << std::endl;
            return -1;
        }
        double params[] = {1.046, 2.104, -0.184, 0.401, 0.507, 0.218, 0.338, 0.057};
        K.setParamsThreeOrbitalTB(params);
        K.tmdNanoribbon(Ny);

        // FOR DEBUGGING: TEST MATRIX CREATED BY THE PROGRAM IN OTHER CODES
        std::ofstream TMDhoppingNano("temp-data/tmd-hopping-nanoribbon.csv");
        if (TMDhoppingNano.is_open())
        {
            TMDhoppingNano << K.getB() << '\n';
        }
        TMDhoppingNano.close();
        K.computeExponential(t, dt);
    }

    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Configuration< L , NSITES > * h = new Configuration< L , NSITES >;
    h->genHsMatrix();

    //  GENERATE THE B-MATRICES.
    OneParticlePropagators< NSITES, L > * Bup =
      new OneParticlePropagators< NSITES, L >;
    OneParticlePropagators< NSITES, L > * Bdown=
      new OneParticlePropagators< NSITES, L >;
    Bup->fillMatrices( true, nu, h->matrix(), K.BpreFactor(dt, mu) );
    Bdown->fillMatrices( false, nu, h->matrix(), K.BpreFactor(dt, mu) );

    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green< NSITES, L, Lbda> * Gup = new Green< NSITES, L, Lbda>;
    Green< NSITES, L, Lbda> * Gdown = new Green< NSITES, L, Lbda>;
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
    //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    // Gup->computeGreenNaive(Bup->list(), L - 1);
    // Gdown->computeGreenNaive(Bdown->list(), L - 1);
    Gup->storeVDU( Bup->list() ); Gdown->storeVDU( Bdown->list() );
    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();

    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double * weights = new double[W * L];
    double * signs = new double[(totalMCSweeps - W) / A];
    double electronDensities = 0;
    double doubleOcs = 0;
    // double energies = 0;
    double zzMags = 0;
    Eigen::Matrix<double, NSITES, NSITES>  magCorrs =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    double LOGweight = 0.;
    double sign = std::copysign(1, Gup->matrix().determinant()
      * Gdown->matrix().determinant() );
    double electronDensity;
    double doubleOc;
    // double energy;
    double zzMag;
    Eigen::Matrix<double, NSITES, NSITES> magCorr;
    double nEl = 0;
    double nUp_nDw = 0;
    // double Hkin = 0;
    double zzAFstFactor = 0;
    Eigen::Matrix<double, NSITES, NSITES> SiSj =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    double meanSign = 0.;

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE,
    //  THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0; int i = 0; int latticeSweepUntilAfresh = 0; int sweep = 0; int step;


    // --- MC LOOP ---


    std::cout << "\nMC loop started. Progress:\n";
    for (step = 0; step < totalMCSteps; step++)
    {
        //  DISPLAY PROGRESS OF THE RUN.
        if ( (step + 1)  % (totalMCSteps/8) == 0 )
        {
            std::cout << (step + 1) * 1. / totalMCSteps * 100 << " %"
              << std::endl << std::endl;
            std::cout << "Average Sign: " << meanSign
              << std::endl << std::endl;
            std::cout << "Log Weight: " << LOGweight
              << std::endl << std::endl;
        }

        //  COMPUTE THE ACCEPTANCE RATIO.
        alphaUp = ( exp( -2 * h->get(l, i) * nu ) - 1 );
        alphaDown = ( exp( 2 * h->get(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup->get(i, i) ) );
        dDown = ( 1 + alphaDown  * ( 1 - Gdown->get(i, i) ) );
        //  SAMPLING: METROPOLIS OR HEAT BATH
        accRatio = fabs( dUp * dDown );
        // accRatio = fabs( dUp * dDown / ( 1 + dUp * dDown ) );

        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
        decisionMaker = dis(gen);

        if (decisionMaker <= accRatio )
        {
            //  KEEP TRACK OF WEIGHT
            LOGweight += log( fabs( dUp ) ) + log( fabs ( dDown ) );
            sign *= std::copysign(1, dUp * dDown );
            //  FLIP A SPIN
            h->flip(l, i);
            //  UPDATE Bs
            Bup->update(l, i, alphaUp); Bdown->update(l, i, alphaDown);
            //  RANK-ONE UPDATE -> O(N^2)
            Gup->update(alphaUp, dUp, i); Gdown->update(alphaDown, dDown, i);
        }


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
            if ( sweep < W )
            {
              //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
              weights[sweep * L + l] = LOGweight;
            }
            //  STORE ELECTRON DENSITY, DOUBLE OCCUPANCY, AND SPIN-SPIN CORRELATIONS.
            electronDensity = 0.; doubleOc = 0.; zzMag = 0.; energy = 0.;
            for (int x = 0; x < NSITES; x++)
            {
                electronDensity -= ( Gup->get(x, x) + Gdown->get(x, x) );
                doubleOc += - Gup->get(x, x) - Gdown->get(x, x) + Gup->get(x, x)
                * Gdown->get(x, x);
                magCorr(x, x) = ( Gup->get(x, x) + Gdown->get(x, x) )
                  - 2 * Gup->get(x, x) * Gdown->get(x, x);
                zzMag += magCorr(x, x);
                for (int y = 0; y < x; y++)
                {
                    magCorr(x, y) =
                      - ( 1 - Gup->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                      - ( 1 - Gdown->get(x, x) ) * ( 1 - Gup->get(y, y) )
                      + ( 1 - Gup->get(x, x) ) * ( 1 - Gup->get(y, y) )
                      + ( 1 - Gdown->get(x, x) ) * ( 1 - Gdown->get(y, y) )
                      - Gup->get(y, x) * Gup->get(x, y)
                      - Gdown->get(y, x) * Gdown->get(x, y);
		                magCorr(y, x) = magCorr(x, y);
                    if ( geom == 1 or geom == 2 )
                    {
                      zzMag += 2 * pow(-1, x - y ) * magCorr(x, y);
                    }

                    if ( geom == 3 or geom == 4 or geom == 5 or geom == 6 )
                    {
                		    if ( ( x + ( ( x - x % int (sqrt(NSITES)) )
                          / int (sqrt(NSITES)) ) % 2 ) % 2
                          == ( y + ( ( y - y % int (sqrt(NSITES)) )
                          / int (sqrt(NSITES)) ) % 2 ) % 2 )
                		    {
                            zzMag += 2 * magCorr(x, y);
                		    }
                		    else
                		    {
                		        zzMag -= 2 * magCorr(x, y);
                		    }
                    }
                  }
            }
            electronDensity /= NSITES; electronDensity += 2;
            doubleOc /= NSITES; doubleOc += 1;
            zzMag /= NSITES;

            //  KEEP TRACK OF THE SIGN OF THE SWEEP
            meanSign += ( sign - meanSign ) / ( l + 1 );
            electronDensities +=
              ( electronDensity * sign - electronDensities ) / ( l + 1 ) ;
            doubleOcs +=
              ( doubleOc * sign - doubleOcs ) / ( l + 1 ) ;
            magCorrs +=
              (magCorr * sign - magCorrs ) / ( l + 1 );
            zzMags +=
              (zzMag * sign - zzMags ) / ( l + 1 );

            //  DEAL WITH THE GREEN'S FUNCTIONS.


            //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == GREEN_AFRESH_FREQ)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                //  Uncomment to compute the product in the naive, unstable manner
                // Gup->computeGreenNaive(Bup->list(), l);
                // Gdown->computeGreenNaive(Bdown->list(), l);
                //  Uncomment to compute the product in the stabilized,
                //  but slightly inefficient way
                // Gup->computeStableGreenNaiveR(Bup->list(), l);
                // Gdown->computeStableGreenNaiveR(Bdown->list(), l);
                //  Most efficient solution (storing decompositions)
                if (l != ( L - 1 ) )
                {
                    Gup->storeUDV(Bup->list(), l, GREEN_AFRESH_FREQ);
                    Gdown->storeUDV(Bdown->list(), l, GREEN_AFRESH_FREQ);
                    //  This is the standard way described in
                    //  "Stable simulations of models of interacting electrons"
                    Gup->computeStableGreen(l, GREEN_AFRESH_FREQ);
                    Gdown->computeStableGreen(l, GREEN_AFRESH_FREQ);
                }
                else
                {
                    Gup->storeVDU( Bup->list() ); Gdown->storeVDU( Bdown->list() );
                    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();
                }
                latticeSweepUntilAfresh = 0;
            }
            else
            {   //  WRAPPING.
                Gup->wrap( Bup->matrix(l) ); Gdown->wrap( Bdown->matrix(l) );
            }
            if (l < L - 1)
            {
                l += 1; i = 0;
            }
            else
            {
                if ( (sweep >= W) )
                {
                    signs[sweep - W] = meanSign;
                    if ( sweep % A == 0 )
                    {
                      nEl += ( electronDensities / meanSign - nEl )
                       / ( (sweep - W)/A + 1 ) ;
                      nUp_nDw += ( doubleOcs / meanSign - nUp_nDw )
                       / ( (sweep - W)/A + 1 ) ;
                      SiSj += ( magCorrs / meanSign - SiSj )
                       / ( (sweep - W)/A + 1 ) ;
                      zzAFstFactor += ( zzMags / meanSign - zzAFstFactor )
                       / ( (sweep - W)/A + 1 ) ;
                    }
                    electronDensities = 0.; doubleOcs = 0.;
                    magCorrs = Eigen::Matrix<double, NSITES, NSITES>::Zero();
                    zzMags = 0.;

                }
                meanSign = 0.;
                //  MOVE SWEEP COUNTER
                sweep += 1;
                l = 0; i = 0;
            }
        }
    }   //  END OF MC LOOP.

    std::cout << "Simulation ended" << std::endl << std::endl;
    std::cout << "nEl: " << nEl << std::endl << std::endl;
    std::cout << "nUp_nDw: " << nUp_nDw << std::endl << std::endl;

    //  SAVE OUTPUT.
    std::ofstream file0("temp-data/simulationParameters.csv");
    if (file0.is_open())
    {
      file0 << std::left << std::setw(50) << "Number of sites," << NSITES << '\n';
      file0 << std::left << std::setw(50) << "dt," << dt << '\n';
      file0 << std::left << std::setw(50) << "beta," << BETA << '\n';
      file0 << std::left << std::setw(50) << "L," << L << '\n';
      file0 << std::left << std::setw(50) << "t," << t << '\n';
      file0 << std::left << std::setw(50) << "U," << U << '\n';
      file0 << std::left << std::setw(50) << "mu," << mu << '\n';
      file0 << std::left << std::setw(50) << "totalMCSweeps," << totalMCSweeps << '\n';
      file0 << std::left << std::setw(50) << "Frequency of recomputing G,"
        << GREEN_AFRESH_FREQ << '\n';
      file0 << std::left << std::setw(50)
        << "Number of multiplied Bs after stabilization," << Lbda << '\n';
      file0 << std::left << std::setw(50) << "Geometry," << geom << '\n';
      file0 << std::left << std::setw(50) << "Ny," << Ny << '\n';
    } file0.close();
    //  STORE MEASUREMENTS
    std::ofstream file1("temp-data/Log-weights.csv");
    std::ofstream file2("temp-data/Local-av-sign.csv");
    std::ofstream file3("temp-data/MeasurementsScalars.csv");
    std::ofstream file4("temp-data/EqTimeSzCorrelations.csv");
    if ( file1.is_open() and file2.is_open() and file3.is_open() and file4.is_open() )
    {
        file1 << std::left << std::setw(50) << "Configuration log weight" << '\n';
        file2 << std::left << std::setw(50) << "Local average sign" << '\n';
        for (int s = 0; s < W; s++)
        {
            for (int slice = 0; slice < L; slice++)
            {
                file1 << std::left << std::setw(50) << weights[s * L + slice] << '\n';
            }
        }
        for (int s = 0; s < (totalMCSweeps - W) / A ; s++)
        {
            file2 << std::left << std::setw(50) << signs[s] << '\n';
        }
        file3 << std::left << std::setw(50) << "Electron density <n>,";
        file3 << std::left << std::setw(50) << "Double occupancy <n+ n->,";
        file3 << std::left << std::setw(50) << "ZZ AF Structure Factor" << '\n';
        file3 << std::left << std::setw(50) << std::setprecision(10)
        << nEl << ",";
        file3 << std::left << std::setw(50) << std::setprecision(10)
        << nUp_nDw << ",";
        file3 << std::left << std::setw(50) << std::setprecision(10)
        << zzAFstFactor << '\n';
        file4 << std::left << std::setw(50) << "<Sz_i Sz_j >" << '\n';
        file4 << std::setprecision(10) << SiSj << '\n';
        file1 << '\n';
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();

    delete[] weights; delete[] signs;
    delete Gup; delete Gdown; delete h; delete Bup; delete Bdown;

    return 0;
}
