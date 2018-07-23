//
//  matrixgen.h
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#ifndef matrixgen_h
#define matrixgen_h

Eigen::MatrixXd genHoppingMatrix(int N, double mu);
Eigen::MatrixXd genHsMatrix(int L, int N);
void genBmatrix(Eigen::MatrixXd* Bs, bool spin, double nu, int N, int L, Eigen::MatrixXd h, Eigen::MatrixXd BpreFactor); // pass an array of Bs
Eigen::MatrixXd regenB(bool spin, double nu, int N, Eigen::RowVectorXd h_l, Eigen::MatrixXd BpreFactor);
Eigen::VectorXd uSigma(int N, Eigen::MatrixXd Green, int i);
Eigen::VectorXd wSigma(int N, Eigen::MatrixXd Green, int i);

#endif /* matrixgen_h */
