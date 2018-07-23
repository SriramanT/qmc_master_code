//
//  UDVpartialPivoting.h
//  
//
//  Created by Francisco Brito on 08/05/2018.
//

#ifndef UDVpartialPivoting_h
#define UDVpartialPivoting_h

class UDVpartialPivoting
{
    Eigen::MatrixXd m;
    Eigen::MatrixXd R;
    Eigen::MatrixXd P;
public:
    UDVpartialPivoting(Eigen::MatrixXd toDecompose);
    void printMatrixToDecompose();
    Eigen::MatrixXd QR_and_getU();
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getV();
    Eigen::MatrixXd getP();
};

#endif /* UDVpartialPivoting_h */
