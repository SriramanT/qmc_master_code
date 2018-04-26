//
//  createHoppingMatrix.cpp
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#include <Eigen/Dense>

Eigen::MatrixXd generateHSmatrix(int L, int nSites)
{
    // Generate the HS field matrix
    
    int l;
    int i;
    
    Eigen::MatrixXd h = Eigen::MatrixXd::Random(L,nSites);
    
    for (l = 0; l < L; l++)
    {
        
        for (i = 0; i < nSites; i++)
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
