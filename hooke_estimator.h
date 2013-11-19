//
// @(#)hooke_estimator.h  1.1  misha-20aug2012
//
// Copyright (c) 2012 of Mikhail Tchernychev
// All rights reserved.
//

#ifndef HOOKE_ESTIMATOR_H
#define HOOKE_ESTIMATOR_H

#include "inversion.h"

typedef struct _hooke_params {

  CInversion * pInversion; 
  double hooke_rho;
  double hooke_eps;
  int    hooke_niter;
  double * fvec;
  int m;
  int n;

 _hooke_params() {
	hooke_rho   = 0.2;
	hooke_eps   = 1.e-5;
	hooke_niter = 1000;
	pInversion  = NULL; 
	fvec        = NULL;
	n           = 0;
	m           = 0;
 }


} HOOKE_PARAMS;


typedef double (*HOOKE_FUNC)(double *x, int n, void * param);

int hooke(int nvars, double * startpt, double * endpt, double rho, double epsilon, int itermax, HOOKE_FUNC func, void *param,
	      PROGRESS_FUNC progress_func);

double *  Estimator_Hooke(void * func_param, void * class_param,  int n, int m, double *x, 
			  double *x_up, double *x_low, double * fit, double *st_err, double *max_fit, char * message);

#endif


