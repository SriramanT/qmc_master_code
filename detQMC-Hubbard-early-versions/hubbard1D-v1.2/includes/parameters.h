//
//  parameters.h
//  
//
//  Created by Francisco Brito on 23/04/2018.
//

#ifndef parameters_h
#define parameters_h

//  Physical parameters
extern double dt;                                   //  time subinterval width. error scales as dt^2
extern double beta;                                 //  imaginary time interval or equivalent maximum temperature of the (d+1) classical system
extern int L;                                       //  number of imaginary time subintervals
extern double t;                                    //  hopping parameters
extern double U;                                    //  interaction energy
extern double nu;                                   //  Hubbard Stratonovich transformation parameter
extern double mu;                                   //  chemical potential

void printParameters();

#endif /* parameters_h */
