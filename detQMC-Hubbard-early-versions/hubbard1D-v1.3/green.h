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
    Eigen::MatrixXd U;
    Eigen::MatrixXd D;
    Eigen::MatrixXd V;
    Eigen::MatrixXd Gforward;
    Eigen::MatrixXd Gbackward;
public:
    Green(int nSites, int nSlices);  // initialize Green's matrix as Id
    void computeGreenNaive(Eigen::MatrixXd* Bs, int l);
    void computeStableGreenNaiveR(Eigen::MatrixXd* Bs, int l, int Lbda);
    void computeStableGreenNaiveL(Eigen::MatrixXd* Bs, int l, int Lbda);
    void computeGreenFromVDU(Eigen::MatrixXd VsLast, Eigen::MatrixXd DsLast, Eigen::MatrixXd UsLast);
    void storeVDU(Eigen::MatrixXd* Bs, int Lbda, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs);
    void computeStableGreen(int l, int Lbda, int greenAfreshFreq, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs);
    void computeBlockOfGreens(int l, int Lbda, int greenAfreshFreq, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs);
    void storeUDV(Eigen::MatrixXd* Bs, int l, int Lbda, int greenAfreshFreq, Eigen::MatrixXd* Us, Eigen::MatrixXd* Ds, Eigen::MatrixXd* Vs);
    Eigen::MatrixXd getG();
    Eigen::MatrixXd getGforward();
    Eigen::MatrixXd getGbackward();
    Eigen::MatrixXd getM();
    void printGreen(bool spin);
    void printM(bool spin);
};


#endif /* green_h */
