//
//  afqmcHubbardOneDimension.cpp
//  
//
//  Created by Francisco Brito on 09/03/2018.
//
//  This code constitutes a first attempt at simulating the Hubbard model on a 1D chain
//  using auxiliary field quantum Monte Carlo.
//  The notation in use is based on the lecture notes "Numerical Methods for Quantum Monte Carlo
//  Simulations of the Hubbard Model by Zhaojun Bai, Wenbin Chen, Richard Scalettar, and
//  Ichitaro Yamazaki
//

#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <typeinfo>
#include <cmath>
#include <random>
#include <time.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <fstream>

#define NSITES 50

using namespace Eigen;
using namespace std;

// Auxiliary functions

MatrixXd createHoppingMatrix(int nSites)
{
    MatrixXd K = MatrixXd::Zero(NSITES, NSITES) ;
    // Set the elements of the hopping matrix that define PBC corresponding to the ends of the 1D chain
    K(0, 1) = 1;
    K(0, nSites - 1) = 1;
    K(nSites - 1, 0) = 1;
    K(nSites - 1, nSites - 2) = 1;
    
    // Set the remaining ones
    int i;
    for (i = 1; i < nSites - 1; i++)
    {
        K(i, i - 1) = 1;
        K(i, i + 1) = 1;
    }
    
    return K;
}

MatrixXd generateHSmatrix(int L, int nSites)
{
    // Generate the HS field matrix

    int l;
    int i;
    
    MatrixXd h = MatrixXd::Random(L,nSites);
    
    for (l = 0; l < L; l++)
    {
      
      for (i = 0; i < NSITES; i++)
      {
          if ( h(l, i) < 0 )
          {
              h(l, i) = -1;
          }
          else
          {
              h(l, i) = 1;
          }
      }
    }
    
    return h;
}
    
int main()
{
    // Toggle prints
    
    int printsOn = 0;                                               // 0 - NO PRINTS, 1 - PRINTS

    cout << "\nAuxiliary field QMC for 1D Hubbard chain\n" << endl;
    cout << "\nNumber of sites: " << NSITES << endl;
    
    //  Usual general use variables
    
    int i;
    int j;
    int l;
    int m = 0;
    int step;
    int i_chosen;
    int l_chosen;

    //  Things to do with random numbers throughout the code
    
    const static int seed = 12345;                              //  set a seed
    mt19937 gen(seed);                                          //  mt19937 algorithm to generate random numbers
    uniform_real_distribution<> dis(0.0, 1.0);                  //  set the distribution from which to draw numbers
    double decisionMaker;                                       //  to accept or not to accept, hence the question.
    
    //  Monte Carlo specific variables
    
    const static int totalMCSteps = 10;
    const static int W = 3000;                                  //  warm-up steps
    const static int autoCorrTime = 500;                        //  auto-correlation time
    const static int M = (totalMCSteps - W)/autoCorrTime;       //  number of measurements
    
    //  Physical parameters
    
    const static int L = 5;                                     //  number of imaginary time subintervals
    cout << "L: " << L << endl;
    const static double beta = 1.;                              //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
    cout << "beta: " << beta << endl;
    const static double dt = beta/L;                            //  time subinterval width
    cout << "dt: " << dt << endl;
    const static double t = 1.;                                 //  hopping parameters
    cout << "t: " << t << endl;
    const static double U = 10.;                                //  interaction energy
    cout << "U: " << U << endl;
    double nu = acosh( exp( U * dt / 2 ) );                     //  Hubbard Stratonovich transformation parameter
    
    // Initialize the HS field at all time slices in [0, beta] with +1 and -1 randomly
    
    MatrixXd h;
    h = generateHSmatrix(L, NSITES); // HS field h_{l = 1,... L, i = 1, ..., NSITES}
    
    // Hopping matrix 1D chain w/ PBC
    
    MatrixXd K;
    K = createHoppingMatrix(NSITES);
    
    // Compute matrix 'prefactor' of the B-matrices e^{t\delta\tau K}
    
    const static MatrixXd B_preFactor = (t * dt * K).exp();
  
    // Some prints for debugging. Toggle them if you want at the top.
    if (printsOn == 1)
    {
        cout <<"\n K = \n\n" << K << "\n\n";
        cout <<"exp (t * dt * K) = \n\n" << B_preFactor << "\n\n";
        cout << "Hubbard Stratonovich Field\n\n h = \n\n" << h << "\n\n";
    }
    
    // Build the B-matrices. We need a copy of each to perform the updates (naively)
    
    MatrixXd BpOld[L]; // B plus
    MatrixXd BmOld[L]; // B minus
    MatrixXd BpNew[L]; // B plus
    MatrixXd BmNew[L]; // B minus
    
    for (l = 0; l < L; l++)
    {
        BpOld[l] =  MatrixXd::Identity(NSITES,NSITES);
        BpOld[l].diagonal() = nu * h.row(l);
        BpOld[l] = B_preFactor * BpOld[l].exp();
        BmOld[l] =  MatrixXd::Identity(NSITES,NSITES);
        BmOld[l].diagonal() = (-nu) * h.row(l);
        BmOld[l] = B_preFactor * BmOld[l].exp();
        // The new matrices are initially set as equal to the old ones
        BpNew[l] = BpOld[l];
        BmNew[l] = BmOld[l];
        if (printsOn ==1)
        {
            cout << "\n\nHere is a B_plus matrix\n\n" << BpOld[l] << endl;
            cout << "\n\nHere is a B_minus matrix\n\n" << BmOld[l] << endl;
        }
    }
    
    // Build the M-matrices and Green's function (matrix)
    
    MatrixXd MPlusOld = MatrixXd::Identity(NSITES,NSITES);
    MatrixXd MMinusOld = MatrixXd::Identity(NSITES,NSITES);
    MatrixXd MPlusNew;
    MatrixXd MMinusNew;
    
    for (l = 0; l < L; l++)
    {
        MPlusOld *= BpOld[L - 1 - l];
        MMinusOld *= BmOld[L - 1 - l];
    }
    
    MPlusOld += MatrixXd::Identity(NSITES,NSITES);
    MMinusOld += MatrixXd::Identity(NSITES,NSITES);
    
    MatrixXd GreenPlus;
    MatrixXd GreenMinus;
    
    GreenPlus = MPlusOld.inverse();
    GreenMinus = MMinusOld.inverse();
    
    if (printsOn == 1)
    {
        cout << "\n\nHere is a M_plus matrix\n\n" << MPlusOld << endl;
        cout << "\n\nHere is a M_minus matrix\n\n" << MMinusOld << endl;
        cout << "\n\nHere is the M_plus determinant\n\n" << MPlusOld.determinant() << endl;
        cout << "\n\nHere is the M_minus determinant\n\n" << MMinusOld.determinant() << endl;
    }
    
    //  Initialize determinants and acceptance ratio
    
    double detOldPlus = MPlusOld.determinant();
    double detOldMinus = MMinusOld.determinant();
    double detsProdOld = abs( detOldPlus * detOldMinus );
    double detsProdNew = detsProdOld;
    double detNewPlus;
    double detNewMinus;
    double acceptanceRatio;
    double smartAccept;
    double alphaPlus;
    double alphaMinus;
    
    VectorXd weights(totalMCSteps);
//    VectorXd density(M);
    
    //  inititialize chosen entry of HS field matrix to (0, 0)
    
    l_chosen = 0;
    i_chosen = 0;
    
    for (step = 0; step < totalMCSteps; step++)
    {
        if (step % (totalMCSteps/10) == 0)
        {
            cout << "\nstep: " << step + 1 << "/" << totalMCSteps << endl; // print step
        }
        
        if (printsOn == 1)
        {
            cout << "Chosen l: " << l_chosen << endl;
            cout << "Chosen i: " << i_chosen << endl;
        }
        
//        //  Measurements
//
//        if ((step > W) and (step % autoCorrTime == 0))
//        {
//            density(m) = ( ((MPlusNew.inverse()).diagonal()).trace() + ((MMinusNew.inverse()).diagonal()).trace() ) / ( 2 * NSITES);
//        }
        
        // save weight of the configuration to see convergence
        
        weights(step) = detsProdNew;
        
        // flip
        h(l_chosen, i_chosen) *= -1;
        
        // rebuild B-matrices

        BpNew[l_chosen] =  MatrixXd::Identity(NSITES,NSITES);
    
        BpNew[l_chosen].diagonal() = nu * h.row(l_chosen);
        
        BpNew[l_chosen] = B_preFactor * BpNew[l_chosen].exp();
        
        if (printsOn == 1)
        {
            cout << "\n\nOld B_plus matrix\n\n" << BpOld[l_chosen] << endl;
            cout << "\n\nNew B_plus matrix\n\n" << BpNew[l_chosen] << endl;
        }

        BmNew[l_chosen] =  MatrixXd::Identity(NSITES,NSITES);

        BmNew[l_chosen].diagonal() = (-nu) * h.row(l_chosen);

        BmNew[l_chosen] = B_preFactor * BmNew[l_chosen].exp();
        
        if (printsOn == 1)
        {
            cout << "\n\nOld B_minus matrix\n\n" << BmOld[l_chosen] << endl;
            cout << "\n\nNew B_minus matrix\n\n" << BmNew[l_chosen] << endl;
        }
        
        // rebuild M-matrices
        
        MPlusNew = MatrixXd::Identity(NSITES,NSITES);
        MMinusNew = MatrixXd::Identity(NSITES,NSITES);
        
        for (l = 0; l < L; l++)
        {
            MPlusNew *= BpNew[L - 1 - l];
            MMinusNew *= BmNew[L - 1 - l];
        }

        MPlusNew += MatrixXd::Identity(NSITES,NSITES);
        MMinusNew += MatrixXd::Identity(NSITES,NSITES);
            
        if (printsOn == 1)
        {
            cout << "\n\nOld M_plus matrix\n\n" << MPlusOld << endl;
            cout << "\n\nOld M_minus matrix\n\n" << MMinusOld << endl;
            cout << "\n\nNew M_plus matrix\n\n" << MPlusNew << endl;
            cout << "\n\nNew M_minus matrix\n\n" << MMinusNew << endl;
            cout << "\n\nNew M_plus determinant\n\n" << MPlusNew.determinant() << endl;
            cout << "\n\nNew M_minus determinant\n\n" << MMinusNew.determinant() << endl;
        }

        detNewPlus = MPlusNew.determinant();
        detNewMinus = MMinusNew.determinant();
        detsProdNew = abs(detNewPlus * detNewMinus);
        acceptanceRatio = detsProdNew / detsProdOld;
        
        alphaPlus = ( exp( -2 * (-1) * h(l_chosen,i_chosen) * nu ) - 1 );
        alphaMinus = ( exp( 2 * (-1) * h(l_chosen,i_chosen) * nu ) - 1 );
        
        smartAccept = ( 1 + alphaPlus  * ( 1 - GreenPlus(i_chosen, i_chosen) ) ) * ( 1 + alphaMinus  * ( 1 - GreenMinus(i_chosen, i_chosen) ) );
        
        if (printsOn == 1)
        {
            cout << "\n\nacceptance ratio: " << acceptanceRatio << endl;
        }
        
        cout << "\n\nacceptance ratio: " << acceptanceRatio << endl;
        cout << "\n\nacceptance ratio: " << smartAccept << endl;
        
        // draw random number to decide
        
        decisionMaker = dis(gen);
        
        if (decisionMaker <= acceptanceRatio)
        {
            //update everything
            BpOld[l_chosen] = BpNew[l_chosen];
            BmOld[l_chosen] = BmNew[l_chosen];
            MPlusOld = MPlusNew;
            MMinusOld = MMinusNew;
            //h(l_chosen, i_chosen) = h(l_chosen, i_chosen);
            detsProdOld = detsProdNew;
            
            if (printsOn == 1)
            {
                cout << "Accepted\n" << endl;
            }
            
            GreenPlus = MPlusNew.inverse();
            GreenMinus = MMinusNew.inverse();
            
//            if ( l_chosen == i_chosen )
//            {
//                GreenPlus(i_chosen, i_chosen) -= alphaPlus/smartAccept * ( 1 - GreenPlus(i_chosen, i_chosen) ) * GreenPlus(i_chosen, i_chosen);
//                GreenMinus(i_chosen, i_chosen) -= alphaMinus/smartAccept * ( 1 - GreenMinus(i_chosen, i_chosen) ) * GreenMinus(i_chosen, i_chosen);
//            }
//            else
//            {
//                GreenPlus(i_chosen, i_chosen) -= alphaPlus/smartAccept * ( - GreenPlus(i_chosen, i_chosen) ) * GreenPlus(i_chosen, i_chosen);
//                GreenMinus(i_chosen, i_chosen) -= alphaMinus/smartAccept * ( - GreenMinus(i_chosen, i_chosen) ) * GreenMinus(i_chosen, i_chosen);
//            }
            
            
            
        }
        else
        {
            if (printsOn == 1)
            {
                cout << "Rejected\n" << endl;
            }
            // revert changes (switch new -> old wrt the above commands)
            BpNew[l_chosen] = BpOld[l_chosen];
            BmNew[l_chosen] = BmOld[l_chosen];
            MPlusNew = MPlusOld;
            MMinusNew = MMinusOld;
            h(l_chosen, i_chosen) *= -1;
            detsProdNew = detsProdOld;
        }
        
        // 'pick' a site (spatial + HS field)
        
        if (i_chosen < NSITES - 1)
        {
            i_chosen += 1;
            //l_chosen = l_chosen;
        }
        
        else
        {
            if (l_chosen < L - 1)
            {
                l_chosen += 1;
                i_chosen = 0;
            }
            else
            {
                l_chosen = 0;
                i_chosen = 0;
            }
        }
        
        
    }
    
    // Save weights of accepted configurations to file

    std::ofstream file("test.txt");
    if (file.is_open())
    {
        file << weights << '\n';
    }

//    // Save densities to file
//
//    std::ofstream file1("testDensities.txt");
//    if (file1.is_open())
//    {
//        file1 << density << '\n';
//    }
    
    
    
    return 0;
}
