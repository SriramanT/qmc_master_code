//
//  hubbard-1d-finiteT.cpp
//  
//
//  Created by Francisco Brito on 09/05/2018.
//
//  This program simulates the Hubbard model on a 1D chain
//  using auxiliary field (or determinant) Quantum Monte Carlo.
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
#include "prints.h"
#include "green.h"


int main()
{
    //  TOOGLE DEBUGGING MODE.
    const bool debug = false;
    
    //  SET DEFAULT PARAMETERS.
    const int N = 4;  //  # sites
    const double dt = 0.125;  //  Trotter error, or time subinterval width. error scales as dt^2
    const double beta = 6.;  //  inverse temperature
    const double t = 1.;  //  hopping parameter. For 2 sites w/ PBCs use t = 2, corresponding to t = 1.
    const double U = 4.;  //  on-site interaction
    const double mu_U = 0.4;   //  chemical potential divided by interaction
    const int greenAfreshFreq = 4 ;  //   how often to calculate Green's functions afresh (in # im-time slices)
    const double mu = mu_U * U;  //  chemical potential
    const int L = beta / dt;  //  # slices
    const double nu = pow( (U * dt), 0.5) + pow( (U * dt), 1.5) / 12;  //  HS transformation parameter
    const int Lbda = L / greenAfreshFreq;  //  Lbda = # intervals in which the product of B's is divided to stabilize.

    
    //  RANDOM NUMBER GENERATION AND MONTE CARLO-RELATED VARIABLES.
    const int seed = 1;
    std::mt19937 gen(seed);  //  mt19937 algorithm to generate random numbers
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double decisionMaker;
        //  to accept or not to accept, hence the question.
    const int totalMCSweeps = 100;
        //  number of measurements will be totalMCSweeps * L ,
        //  i.e. we make a measurement every slice, then find correlations, etc. in post-processing
    const int totalMCSteps = totalMCSweeps * N * L;

    
    // -- INITIALIZATION ---
    
     
    //  HOPPING MATRIX FOR 1D CHAIN W/ PBCs.
    const Eigen::MatrixXd K = genHoppingMatrix(N);
    //  INITIALIZE THE HS MATRIX WITH +1 AND -1 RANDOMLY.
    Eigen::MatrixXd h = genHsMatrix(L, N);  //    HS field h_{l = 1,...,L ; i = 1,...,N}
    //  COMPUTE MATRIX 'PREFACTOR' OF THE B-MATRICES B = e^(t dt K).
    const Eigen::MatrixXd BpreFactor = (t * dt * K).exp();
    
    initialPrint(N, dt, beta, L, t, U, mu, totalMCSteps, debug, K, BpreFactor, h);
    
    //  GENERATE THE B-MATRICES.
    Eigen::MatrixXd Bup[L];
    Eigen::MatrixXd BupInv;
    genBmatrix(Bup, true, nu, N, L, h, BpreFactor, dt, mu);  //   this line inserts the B-matrices in the array above
    Eigen::MatrixXd Bdown[L];
    Eigen::MatrixXd BdownInv;
    genBmatrix(Bdown, false, nu, N, L, h, BpreFactor, dt, mu);  //   this line inserts the B-matrices in the array above
    
    //  ALLOCATE MEMORY TO STORE THE PARTIAL PRODUCTS INVOLVED IN
    //  SPEEDING UP THE LOW TEMPERATURE STABILIZATION.

    Eigen::MatrixXd UsUp[Lbda]; Eigen::MatrixXd DsUp[Lbda]; Eigen::MatrixXd VsUp[Lbda];
    Eigen::MatrixXd UsDown[Lbda]; Eigen::MatrixXd DsDown[Lbda]; Eigen::MatrixXd VsDown[Lbda];

    //  GENERATE THE SPIN-UP AND SPIN-DOWN GREEN FUNCTIONS.
    Green GreenUp(N, L);
    GreenUp.storeVDU(Bup, Lbda, UsUp, DsUp, VsUp);
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
//    GreenUp.computeGreenNaive(Bup, L - 1);  //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    GreenUp.computeGreenFromVDU(VsUp[Lbda - 1], DsUp[Lbda - 1], UsUp[Lbda - 1]);
    Eigen::MatrixXd Gup = GreenUp.getG();
    std::cout << Gup << std::endl;

    Green GreenDown(N, L);
    GreenDown.storeVDU(Bdown, Lbda, UsDown, DsDown, VsDown);
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
//    GreenDown.computeGreenNaive(Bdown, L - 1);  //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    GreenDown.computeGreenFromVDU(VsDown[Lbda - 1], DsDown[Lbda - 1], UsDown[Lbda - 1]);
    Eigen::MatrixXd Gdown = GreenDown.getG();
    
    //  ALLOCATE MEMORY FOR THE UNEQUAL TIME GREEN'S FUNCTIIONS
    Eigen::Matrix<double, N, N> * uneqGupForward = new Eigen::Matrix<double, N, N>[L];
    Eigen::Matrix<double, N, N> * uneqGdownForward = new Eigen::Matrix<double, N, N>[L];
    Eigen::Matrix<double, N, N> * uneqGupBackward = new Eigen::Matrix<double, N, N>[L];
    Eigen::Matrix<double, N, N> * uneqGdownBackward = new Eigen::Matrix<double, N, N>[L];
    uneqGupForward[0] = Gup;
    uneqGdownForward[0] = Gdown;
    uneqGupBackward[0] = Gup - Eigen::MatrixXd::Identity(N, N);
    uneqGdownBackward[0] = Gdown - Eigen::MatrixXd::Identity(N, N);
    
    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    Eigen::VectorXd uUp; Eigen::VectorXd uDown; Eigen::VectorXd wUp; Eigen::VectorXd wDown;
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double * weights = new double[totalMCSweeps];
    double * electronDensities = new double[totalMCSweeps];
    double * doubleOcs = new double[totalMCSweeps];
    Eigen::MatrixXd * magCorrs = new Eigen::MatrixXd[totalMCSweeps];
    magCorrs[0] = Eigen::MatrixXd::Zero(N, N);
    double weight = GreenUp.getM().determinant() * GreenDown.getM().determinant();
    double electronDensity;
    double doubleOc;
    Eigen::MatrixXd magCorr(N, N);
    electronDensities[0] = 0.; doubleOcs[0] = 0.;
    double meanSign = std::copysign( 1. , weight );

    //  INITIALIZE (l, i) <- (0, 0). INITIATIALIZE SPATIAL SWEEP COUNTER.
    //  FOR EACH IMAGINARY TIME SLICE l, LOOP OVER ALL SPATIAL LATTICE, THEN CHANGE SLICE, AND SO ON UNTIL l=L. REPEAT.
    int l = 0; int i = 0;
    int latticeSweepUntilAfresh = 0; int sweep = 0; int step;
    
    
    // --- MC LOOP ---


    std::cout << "\n\nMC loop started. Progress:\n\n";
    for (step = 0; step < totalMCSteps; step++)
    {
        //  DISPLAY PROGRESS OF THE RUN.
        if ( (step + 1)  % (totalMCSteps/8) == 0 )
        {
            std::cout << (step + 1)*1. / totalMCSteps * 100 << " %" << std::endl;
        }

        //  COMPUTE THE ACCEPTANCE RATIO.
        alphaUp = ( exp( -2 * h(l, i) * nu ) - 1 ); alphaDown = ( exp( 2 * h(l, i) * nu ) - 1 );
        dUp = ( 1 + alphaUp  * ( 1 - Gup(i, i) ) ); dDown = ( 1 + alphaDown  * ( 1 - Gdown(i, i) ) );
        accRatio = dUp * dDown; //  Regular Metropolis
//        accRatio = dUp * dDown / ( 1 + dUp * dDown );  // Heat Bath

        //  DECIDE WHETHER OR NOT TO ACCEPT THE STEP.
        decisionMaker = dis(gen);

        if (decisionMaker <= abs(accRatio) )
        {
            //  KEEP TRACK OF WEIGHT
            weight = dUp * dDown * weight;
            //  FLIP A SPIN
            h(l, i) *= -1;
            //  UPDATE Bs
            Bup[l].col(i) *= ( alphaUp + 1 );
            Bdown[l].col(i) *= ( alphaDown + 1 );

            //  RANK-ONE UPDATE -> O(N^2)
            uUp = uSigma(N, Gup, i); uDown = uSigma(N, Gdown, i); wUp = wSigma(N, Gup, i); wDown = wSigma(N, Gdown, i);
            for (int x = 0; x < N; x++)
            {
                for (int y = 0; y < N; y++)
                {
                    Gup(x, y) -= alphaUp / dUp * uUp(x) * wUp(y);
                    Gdown(x, y) -= alphaDown / dDown * uDown(x) * wDown(y);
                }
            }
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


            //  --- MEASUREMENTS ---


            //  STORE ELECTRON DENSITY, DOUBLE OCCUPANCY, AND SPIN-SPIN CORRELATIONS.
            electronDensity = 0.; doubleOc = 0.;
            for (int x = 0; x < N; x++)
            {
                electronDensity -= ( Gup(x, x) + Gdown(x, x) );
                doubleOc = doubleOc - Gup(x, x) - Gdown(x, x) + Gup(x, x) * Gdown(x, x);
                magCorr(x, x) = 3 * ( Gup(x, x) + Gdown(x, x) ) - 6 * Gup(x, x) * Gdown(x, x);
                for (int y = 0; y < x; y++)
                {
                    magCorr(x, y) = - Gup(x, y) * ( 2 * Gdown(y, x) + Gup(y, x) )
                    - Gdown(x, y) * ( 2 * Gup(y, x) + Gdown(y, x) )
                    + ( Gup(x, x) - Gdown(x, x) ) * ( Gup(y, y) - Gdown(y, y) );
                    magCorr(y, x) = magCorr(x, y);
                }
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
                Green GreenUp(N, L); Green GreenDown(N, L);
                //  Uncomment to compute the product in the naive, unstable manner
//                GreenUp.computeGreenNaive(Bup, l);
//                GreenDown.computeGreenNaive(Bdown, l);
                //  Uncomment to compute the product in the stabilized, but slightly inefficient way
//                GreenUp.computeStableGreenNaiveR(Bup, l, Lbda);
//                GreenDown.computeStableGreenNaiveR(Bdown, l, Lbda);
                //  Most efficient solution (storing decompositions)
                if (l != ( L - 1 ) )
                {
                    GreenUp.storeUDV(Bup, l, Lbda, greenAfreshFreq, UsUp, DsUp, VsUp);
                    GreenDown.storeUDV(Bdown, l, Lbda, greenAfreshFreq, UsDown, DsDown, VsDown);
                    //  Using the BlockOfGreens Method, we can obtain time-displaced Green's as well
//                    GreenUp.computeStableGreen(l, Lbda, greenAfreshFreq, UsUp, DsUp, VsUp);
//                    GreenDown.computeStableGreen(l, Lbda, greenAfreshFreq, UsDown, DsDown, VsDown);
                    GreenUp.computeBlockOfGreens(l, Lbda, greenAfreshFreq, UsUp, DsUp, VsUp);
                    GreenDown.computeBlockOfGreens(l, Lbda, greenAfreshFreq, UsDown, DsDown, VsDown);
                    Gup = GreenUp.getG(); Gdown = GreenDown.getG();
                    uneqGupForward[l] = GreenUp.getGforward(); uneqGdownForward[l] = GreenDown.getGforward();
                    uneqGupBackward[l] = GreenUp.getGbackward(); uneqGdownBackward[l] = GreenDown.getGbackward();
                }
                else
                {
                    GreenUp.storeVDU(Bup, Lbda, UsUp, DsUp, VsUp);
                    GreenDown.storeVDU(Bdown, Lbda, UsDown, DsDown, VsDown);
                    GreenUp.computeGreenFromVDU(VsUp[Lbda - 1], DsUp[Lbda - 1], UsUp[Lbda - 1]);
                    GreenDown.computeGreenFromVDU(VsDown[Lbda - 1], DsDown[Lbda - 1], UsDown[Lbda - 1]);
                    Gup = GreenUp.getG(); Gdown = GreenDown.getG();
                }
                latticeSweepUntilAfresh = 0;
            }
            else
            {   //  WRAPPING.
                BupInv = Bup[l].inverse();
                BdownInv = Bdown[l].inverse();
                Gup = Bup[l] * Gup * BupInv; Gdown = Bdown[l] * Gdown * BdownInv;
                uneqGupForward[l + 1] = Bup[l] * uneqGupForward[l]; uneqGdownForward[l + 1] = Bdown[l] * uneqGdownForward[l];
                uneqGupBackward[l + 1] = uneqGupBackward[l] * BupInv; uneqGdownBackward[l + 1] = uneqGdownBackward[l] * BdownInv;
            }
            if (l < L - 1)
            {
                l += 1; i = 0;
            }
            else
            {
                //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
                weights[sweep] = weight;
                if (sweep != ( totalMCSweeps - 1 ) )
                {
                    //  MOVE SWEEP COUNTER
                    sweep += 1;
                    electronDensities[sweep] = 0.; doubleOcs[sweep] = 0.;
                    magCorrs[sweep] = Eigen::MatrixXd::Zero(N, N);
                }
                l = 0; i = 0;
            }
        }
    }   //  END OF MC LOOP.

    std::cout << "Average Sign: " << meanSign << std::endl;

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
    std::ofstream file2("plots/measurementsScalars.txt");
    if ( file2.is_open() )
    {
        file2 << std::left << std::setw(25) << "Configuration weight";
        file2 << std::left << std::setw(25) << "Weight sign";
        file2 << std::left << std::setw(25) << "Electron density <n>";
        file2 << std::left << std::setw(25) << "Double occupancy <n+ n->" << '\n';
        for (int s = 0; s < totalMCSweeps; s++)
        {
            file2 << std::left << std::setw(25) << weights[s];
            file2 << std::left << std::setw(25) << std::copysign( 1. , weights[s] );
            file2 << std::left << std::setw(25) << std::setprecision(10) << electronDensities[s];
            file2 << std::left << std::setw(25) << std::setprecision(10) << doubleOcs[s] << '\n';
        }
    }
    file2.close();
    std::ofstream file3("plots/measurementsCorrelations.txt");
    if ( file3.is_open() )
    {
        file3 << "Spin-spin correlation function <S_i S_j>" << '\n';
        file3 << std::left;
        for (int s = 0; s < totalMCSweeps; s++)
        {
            file3 << std::setprecision(10) << magCorrs[s] << std::endl << std::endl;
        }
        file3 << '\n';
    }
    file3.close();
    
    delete[] weights; delete[] electronDensities; delete[] magCorrs; delete[] doubleOcs;
    delete[] uneqGupForward; delete[] uneqGdownForward; delete[] uneqGupBackward; delete[] uneqGdownBackward;
    
    return 0;
}

