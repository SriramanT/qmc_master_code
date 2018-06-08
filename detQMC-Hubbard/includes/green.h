//
//  green.h
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#ifndef green_h
#define green_h

class Green
{
    Eigen::MatrixXd M;
    Eigen::MatrixXd G;
public:
    Green(int nSites);  // initialize Green's matrix as Id
    void computeGreenNaive(Eigen::MatrixXd* Bs, int L);
    void computeWrappedGreenNaive(Eigen::MatrixXd* Bs, int L, int l_);
    Eigen::MatrixXd getG();
    Eigen::MatrixXd getM();
    void printGreen(bool spin);
    void printM(bool spin);
};


#endif /* green_h */
