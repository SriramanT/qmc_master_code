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
    Eigen::MatrixXd P;
public:
    UDV(Eigen::MatrixXd toDecompose);
    void printMatrixToDecompose();
    Eigen::MatrixXd QR_and_getU();
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getV();
};

#endif /* UDV_h */
