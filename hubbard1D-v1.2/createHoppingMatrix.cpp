//
//  createHoppingMatrix.cpp
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#include <Eigen/Dense>


Eigen::MatrixXd createHoppingMatrix(int nSites)
{
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nSites, nSites) ;
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
