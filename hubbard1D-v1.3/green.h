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
    int N;
    int L;
    Eigen::MatrixXd M;
    Eigen::MatrixXd G;
public:
    Green(int nSites, int nSlices);  // initialize Green's matrix as Id
    void computeGreenNaive(Eigen::MatrixXd* Bs, int l);
    void computeStableGreenNaive(Eigen::MatrixXd* Bs, int l, int k);
    Eigen::MatrixXd getG();
    Eigen::MatrixXd getM();
    void printGreen(bool spin);
    void printM(bool spin);
};


#endif /* green_h */
