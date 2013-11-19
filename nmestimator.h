//
// @(#)nmestimator.h  1.1  misha-20aug2012
//
// Copyright (c) 2012 of Mikhail Tchernychev
// All rights reserved.
//

#ifndef NMESTIMATOR_H
#define NMESTIMATOR_H

#include "inversion.h"

typedef struct _nelder_mead {
  int n;                   // number of points
  int m;                   // number of variables;
  double eps;              // parameters for Nelder Mead
  double scale; 
  int n_iterations;
  CInversion * pInversion; 
  double * fvec;           // vector with function values;
  double * x_low;
  double * x_up;

  _nelder_mead () { n=m=n_iterations=0; eps=1.e-8; scale =1.; pInversion = NULL; fvec = x_low = x_up = NULL; }

} NELDER_MEAD;

double *  Estimator_Nelder_Mead(void * func_param, void * class_param,  int n, int m, double *x, 
			        double *x_up, double *x_low, double * fit, double *st_err, double *max_fit, char * message);

#endif


