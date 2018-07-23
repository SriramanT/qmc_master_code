//
//  matrixgen.h
//  
//
//  Created by Francisco Brito on 09/05/2018.
//

#ifndef matrixgen_h
#define matrixgen_h

Eigen::MatrixXd genHoppingMatrix(int N);
Eigen::MatrixXd genHsMatrix(int L, int N);
void genBmatrix(Eigen::MatrixXd* Bs, bool spin, double nu, int N, int L, Eigen::MatrixXd h, Eigen::MatrixXd BpreFactor, double dt, double mu); // pass an array of Bs
Eigen::MatrixXd regenB(bool spin, double nu, int N, Eigen::RowVectorXd h_l, Eigen::MatrixXd BpreFactor, double dt, double mu);
Eigen::VectorXd uSigma(int N, Eigen::MatrixXd Green, int i);
Eigen::VectorXd wSigma(int N, Eigen::MatrixXd Green, int i);

#endif /* matrixgen_h */
