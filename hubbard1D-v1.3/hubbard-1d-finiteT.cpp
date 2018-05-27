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
//  Ichitaro Yamazaki (2006)
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

#include "matrixgen.h"
#include "prints.h"
#include "green.h"

int main()
{
    //  TOOGLE DEBUGGING MODE AND MANUAL PARAMETER ENTERING MODE
    const bool debug = false;
    const bool manual = false;
    
    //  SET DEFAULT PARAMETERS.
    int N = 2;  //  # sites
    double dt = 0.125;  //  time subinterval width. error scales as dt^2
    double beta = 6.;  //  inverse temperature
    double t = 2.;  //  hopping matrix. For 2 sites w/ PBCs use t = 2.
    double U = 4.;  //  on-site interaction
    double mu_U = 0.6;   //  chemical potential divided by interaction
    int greenAfreshFreq = 4 ;  //   how often to calculate Green's functions afresh (in # im-time slices)
    
    if (manual == true)
    {
        std::cout << "Enter number of sites: " << std::endl;
        std::cin >> N;
        std::cout << "Enter Trotter error: " << std::endl;
        std::cin >> dt;
        std::cout << "Enter inverse temperature: " << std::endl;
        std::cin >> beta;
        std::cout << "Enter hopping parameter: " << std::endl;
        std::cin >> t;
        std::cout << "Enter interaction parameter: " << std::endl;
        std::cin >> U;
        std::cout << "Enter ratio of chemical potential to interaction: " << std::endl;
        std::cin >> mu_U;
        std::cout << "Enter frequency of computation of Green's function afresh: " << std::endl;
        std::cin >> greenAfreshFreq;
    }
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
    const int totalMCSweeps = 2048 ;
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
    genBmatrix(Bup, true, nu, N, L, h, BpreFactor, dt, mu);  //   this line inserts the B-matrices in the array above
    Eigen::MatrixXd Bdown[L];
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

    Green GreenDown(N, L);
    GreenDown.storeVDU(Bdown, Lbda, UsDown, DsDown, VsDown);
    //  Uncomment to compute the Green's function naively instead of using the stored VDU
//    GreenDown.computeGreenNaive(Bdown, L - 1);  //  start at l = L - 1, i.e. G = (1 + B_{L-1} B_{L-2} ... B_{0})^(-1)
    GreenDown.computeGreenFromVDU(VsDown[Lbda - 1], DsDown[Lbda - 1], UsDown[Lbda - 1]);
    Eigen::MatrixXd Gdown = GreenDown.getG();
    
    //  INITIALIZE RANK-ONE UPDATE-RELATED QUANTITIES AND ACCEPTANCE RATIO.
    Eigen::VectorXd uUp; Eigen::VectorXd uDown; Eigen::VectorXd wUp; Eigen::VectorXd wDown;
    double alphaUp; double alphaDown; double dUp; double dDown; double accRatio;

    //  INITIALIZE ARRAYS TO STORE MEASUREMENTS.
    double* weights = new double[totalMCSweeps * L];
    double* electronDensity = new double[totalMCSweeps * L];
    double* doubleOc = new double[totalMCSweeps * L];
    double weight = GreenUp.getM().determinant() * GreenDown.getM().determinant();
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
            
            
            //  STORE WEIGHT OF ACCEPTED CONFIGURATIONS
            weights[sweep] = weight;
            //  STORE ELECTRON DENSITY, AND DOUBLE OCCUPANCY
            electronDensity[sweep] = 0.; doubleOc[sweep] = 0.;
            for (int x = 0; x < N; x++)
            {
                electronDensity[sweep] -= ( Gup(x, x) + Gdown(x, x) );
                doubleOc[sweep] = doubleOc[sweep] - Gup(x, x) - Gdown(x, x) + Gup(x, x) * Gdown(x, x);
            }
            electronDensity[sweep] /= N; electronDensity[sweep] += 2;
            doubleOc[sweep] /= N; doubleOc[sweep] += 1;

            //  MOVE SWEEP COUNTER
            sweep += 1;
            
            //  DEAL WITH THE GREEN'S FUNCTIONS.
            
                //  REBUILD B-MATRICES.
            Bup[l] = regenB(true, nu, N, h.row(l), BpreFactor, dt, mu);
            Bdown[l] = regenB(false, nu, N, h.row(l), BpreFactor, dt, mu);

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
                    GreenUp.computeStableGreen(l, Lbda, greenAfreshFreq, UsUp, DsUp, VsUp);
                    GreenDown.computeStableGreen(l, Lbda, greenAfreshFreq, UsDown, DsDown, VsDown);
                }
                else
                {
                    GreenUp.storeVDU(Bup, Lbda, UsUp, DsUp, VsUp);
                    GreenDown.storeVDU(Bdown, Lbda, UsDown, DsDown, VsDown);
                    GreenUp.computeGreenFromVDU(VsUp[Lbda - 1], DsUp[Lbda - 1], UsUp[Lbda - 1]);
                    GreenDown.computeGreenFromVDU(VsDown[Lbda - 1], DsDown[Lbda - 1], UsDown[Lbda - 1]);
                }
                Gup = GreenUp.getG(); Gdown = GreenDown.getG();
                latticeSweepUntilAfresh = 0;
            }
            else
            {   //  WRAPPING.
                Gup = Bup[l] * Gup * Bup[l].inverse(); Gdown = Bdown[l] * Gdown * Bdown[l].inverse();
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
    //  THIS IS INEFFICIENT AND NEEDS IMPROVING
    std::ofstream file2("plots/weights.txt");
    if (file2.is_open())
    {
        file2 << "weights" << "\t\t\t" << "sign" << '\n';
        for (int s = 0; s < totalMCSweeps * L; s++)
        {
            file2 << weights[s] << "\t\t" << std::copysign( 1. , weights[s] ) << '\n';
        }
    }
    std::ofstream file3("plots/electronDensity.txt");
    if (file3.is_open())
    {
        file3 << "electronDensity" << '\n';
        for (int s = 0; s < totalMCSweeps * L; s++)
        {
            file3 << electronDensity[s] << '\n';
        }
    }
    std::ofstream file4("plots/doubleOc.txt");
    if (file4.is_open())
    {
        file4 << "doubleOccupancy" << '\n';
        for (int s = 0; s < totalMCSweeps * L; s++)
        {
            file4 << doubleOc[s] << '\n';
        }
    }
    
    
    return 0;
}

