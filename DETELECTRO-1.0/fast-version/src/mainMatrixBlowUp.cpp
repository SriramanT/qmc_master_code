//
//  mainUneqTime.cpp
//
//
//  Created by Francisco Brito on 17/09/2018.
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

    Eigen::IOFormat CleanFmt(10, 0, ", ", "\n", "", "");

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
        //  See Python notebook to implement
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
            TMDhopping << K.matrix() << '\n';
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
            TMDhoppingNano << K.matrix() << '\n';
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
    Gup->initializeUneqs(); Gdown->initializeUneqs();

    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double * weights = new double[W * L];
    double LOGweight = 0.;

    double electronDensities = 0;
    double doubleOcs = 0;
    double energies = 0;
    double zzMags = 0;
    Eigen::MatrixXd magCorrZZs =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenFunctionUps =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenFunctionDowns =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd magCorrXXs =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd uneqMagCorrZZs =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd uneqMagCorrXXs =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();

    // double sign = std::copysign(1, Gup->matrix().determinant()
    //   * Gdown->matrix().determinant() );
    double sign = 1;
    double meanSign = 0;

    double electronDensity;
    double doubleOc;
    double energy;
    double zzMag;
    Eigen::MatrixXd magCorrZZ = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd magCorrXX = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd uneqMagCorrZZ = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd uneqMagCorrXX = Eigen::Matrix<double, NSITES, NSITES>::Zero();

    double nEl = 0;
    double nUp_nDw = 0;
    double Hkin = 0;
    double zzAFstFactor = 0;
    Eigen::MatrixXd SiSjZ =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd SiSjX =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd intSiTSjZ =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd intSiTSjX =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenUp = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenDown = Eigen::Matrix<double, NSITES, NSITES>::Zero();

    double nElSq = 0;
    double nUp_nDwSq = 0;
    double HkinSq = 0;
    double zzAFstFactorSq = 0;
    Eigen::MatrixXd SiSjZSq =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd SiSjXSq =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd intSiTSjZSq =
      Eigen::Matrix<double, NSITES, NSITES>::Zero();
    // Eigen::MatrixXd intSiTSjXSq =
    //   Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenUpSq = Eigen::Matrix<double, NSITES, NSITES>::Zero();
    Eigen::MatrixXd GreenDownSq = Eigen::Matrix<double, NSITES, NSITES>::Zero();

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
            //  DEAL WITH THE GREEN'S FUNCTIONS.


                //  DECIDE WHETHER TO COMPUTE GREEN'S FUNCTIONS AFRESH OR TO WRAP.
            if (latticeSweepUntilAfresh == GREEN_AFRESH_FREQ)
            {   //  COMPUTE SPIN-UP AND SPIN-DOWN GREEN'S FUNCTIONS AFRESH.
                //  Uncomment to compute the product in the naive, unstable manner
//                Gup->computeGreenNaive(Bup->list(), l);
//                Gdown->computeGreenNaive(Bdown->list(), l);
                //  Uncomment to compute the product in the stabilized, but slightly inefficient way
//                Gup->computeStableGreenNaiveR(Bup->list(), l);
//                Gdown->computeStableGreenNaiveR(Bdown->list(), l);
                //  Most efficient solution (storing decompositions)
                if (l != ( L - 1 ) )
                {
                    Gup->storeUDV(Bup->list(), l, GREEN_AFRESH_FREQ);
                    Gdown->storeUDV(Bdown->list(), l, GREEN_AFRESH_FREQ);
                    //  This is the standard way described in "Stable simulations
                    //  of models of interacting electrons"
//                    Gup->computeStableGreen(l, GREEN_AFRESH_FREQ);
//                    Gdown->computeStableGreen(l, GREEN_AFRESH_FREQ);
                    //  Using the BlockOfGreens Method, we can obtain
                    //  time-displaced Green's as well
                    Gup->computeBlockOfGreens(l, GREEN_AFRESH_FREQ);
                    Gdown->computeBlockOfGreens(l, GREEN_AFRESH_FREQ);
                }
                else
                {
                    Gup->storeVDU( Bup->list() ); Gdown->storeVDU( Bdown->list() );
                    Gup->computeGreenFromVDU(); Gdown->computeGreenFromVDU();
                    Gup->initializeUneqs(); Gdown->initializeUneqs();
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

                //  MOVE SWEEP COUNTER
                sweep += 1;
                l = 0; i = 0;
            }
        }
    }   //  END OF MC LOOP.

    //  Normalize to mean sign
    nEl /= meanSign; nUp_nDw /= meanSign; SiSjZ /= meanSign; zzAFstFactor /= meanSign;
    Hkin /= meanSign; GreenUp /= meanSign; GreenDown /= meanSign;
    nElSq /= meanSign; nUp_nDwSq /= meanSign; SiSjZSq /= meanSign; zzAFstFactorSq /= meanSign;
    HkinSq /= meanSign; GreenUpSq /= meanSign; GreenDownSq /= meanSign;

    std::cout << "Simulation ended" << std::endl << std::endl;
    file14.close();

    delete[] weights;
    delete Gup; delete Gdown; delete h; delete Bup; delete Bdown;

    return 0;
    }
