//
//  prints.h
//  
//
//  Created by Francisco Brito on 27/04/2018.
//

#ifndef prints_h
#define prints_h

void printWelcome(int nSites);
void printMC(int totalMCSteps, int nSites, int L);
void printStartingMatrices(Eigen::MatrixXd K, Eigen::MatrixXd B_preFactor, Eigen::MatrixXd h, bool yesOrNo);

#endif /* prints_h */
