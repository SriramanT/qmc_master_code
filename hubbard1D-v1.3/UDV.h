//
//  UDV.h
//  
//
//  Created by Francisco Brito on 08/05/2018.
//

#ifndef UDV_h
#define UDV_h

class UDV
{
    Eigen::MatrixXd m;
    Eigen::MatrixXd R;
public:
    UDV(Eigen::MatrixXd toDecompose);
    void printMatrixToDecompose();
    Eigen::MatrixXd QR_and_getU();
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getV();
};

Eigen::MatrixXd computeGreen(Eigen::MatrixXd Bs[], int L, int k, int nSites);
Eigen::MatrixXd computeGreenNaive(Eigen::MatrixXd Bs[], int L, int nSites);

#endif /* UDV_h */
