//
//  matrices.h
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#ifndef matrices_h
#define matrices_h

Eigen::MatrixXd build_Bmatrix(bool sign, double nu, int nSites, Eigen::RowVectorXd h_l, Eigen::MatrixXd B_preFactor, double dt, double mu);

Eigen::VectorXd uSigma(int nSites, Eigen::MatrixXd Green, int i_chosen);

Eigen::VectorXd wSigma(int nSites, Eigen::MatrixXd Green, int i_chosen);

#endif /* matrices_h */
