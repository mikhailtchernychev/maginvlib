//
// @(#)l2estimator.h  1.1  misha-16may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  16may103  by  misha@misha.local
// Version:  16may103  by  misha@misha.local
//

#ifndef L2ESTIMATOR_H
#define L2ESTIMATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct _z29rls_params {

    double eps;                   // parameter to discard small singular numbers for z29rls
    double epsm;                  // parameter to check matrix V (if less them epsm, V[i,j]=0) (z29rls);
    double dm;                    // parameter to find regularization (z29rls)
    double am;                    // estimated reg. parameter after z29rls
    
    _z29rls_params();

} Z29RLS_PARAMS;


int Estimate_L2_rls(int n, int m, double *a, double *y, double *x, void * param,  double * synt = NULL, 
		    double * fit = NULL, double * max_fit = NULL, double *st_err = NULL, char * message = NULL);

#endif
