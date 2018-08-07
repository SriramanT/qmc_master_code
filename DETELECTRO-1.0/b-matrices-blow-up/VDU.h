//
//  VDU.h
//  
//
//  Created by Francisco Brito on 21/05/2018.
//

#ifndef VDU_h
#define VDU_h

class VDU
{
    Eigen::MatrixXd m;
    Eigen::MatrixXd R;
public:
    VDU(Eigen::MatrixXd toDecompose);
    void printMatrixToDecompose();
    Eigen::MatrixXd QR_and_getU();
    Eigen::MatrixXd getD();
    Eigen::MatrixXd getV();
};

#endif /* VDU_h */
