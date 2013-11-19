//
// @(#)nmestimator.cpp  1.1  misha-20aug2012
//
// Copyright (c) 2012 of Mikhail Tchernychev
// All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "nmestimator.h"

extern "C" {
double simplex(double (*objfunc)(double * x, void * pFuncParams), double start[],int n, double EPSILON, double scale, void (*constrain)(double *x,int n, void *pFuncParams),
				void *pFuncParams);

}

double objfunc(double *x, void *pParams) {

	if(!pParams) return 0.;
	NELDER_MEAD * pSimplexParams = (NELDER_MEAD *)pParams;
	
	pSimplexParams->pInversion->SetAllVariables(x, pSimplexParams->m);
	pSimplexParams->pInversion->ComputeField(pSimplexParams->fvec, pSimplexParams->n, 1);

   double discrepancy = 0;
   for(int i=0; i<pSimplexParams->n; i++) {
	   discrepancy += pSimplexParams->fvec[i]*pSimplexParams->fvec[i];
  }

   discrepancy = sqrt(discrepancy/pSimplexParams->n);

   PROGRESS_FUNC pFunc = pSimplexParams->pInversion->GetProgressFunc();
   if(pFunc) (*pFunc)((void *)pSimplexParams->pInversion, pSimplexParams->n_iterations++, discrepancy);

   return discrepancy;

}


void constrain(double *x, int n, void *pParams) {

	if(!pParams) return;
	NELDER_MEAD * pSimplexParams = (NELDER_MEAD *)pParams;

	if(!pSimplexParams->x_low && !pSimplexParams->x_up) return;

	for(int i=0; i<n; i++) {
		if(pSimplexParams->x_low && x[i]<pSimplexParams->x_low[i]) { x[i] = pSimplexParams->x_low[i]; continue; }
		if(pSimplexParams->x_up  && x[i]>pSimplexParams->x_up[i])  { x[i] = pSimplexParams->x_up[i]; continue;  }
	}
}

// func_param   -  parameters for Nelder Mead routine; if NULL assigns automatically
//                 if not NULL, has to be initialized by the caller
// class_param  -  CInversion class
// n            -  number of data points
// m            -  number of variables
// x            -  varaibles vector
// x_up         -  upper limits
// x_low        -  lower limits
// fit          -  rms for final model fit
// st_err       -  standard error for estimates
// max_fit      -  max. difference between observed and comoputed fields
// message      -  message with exit condition 


double *  Estimator_Nelder_Mead(void * func_param, void * class_param,  int n, int m, double *x, 
	double *x_up, double *x_low, double * fit, double *st_err, double *max_fit, char * message) {

		NELDER_MEAD * pSimplexParams = (NELDER_MEAD *)func_param;

		if(!func_param) {
			pSimplexParams = new NELDER_MEAD;
			if(!pSimplexParams) return NULL;
		}

		if(pSimplexParams->fvec) delete [] pSimplexParams->fvec;

		pSimplexParams->fvec       =  new double [n];
		pSimplexParams->x_low      = x_low;
		pSimplexParams->x_up       = x_up;
		pSimplexParams->pInversion = (CInversion *)class_param;
		pSimplexParams->n          = n;
		pSimplexParams->m          = m;

		*fit = simplex(objfunc, x, m, pSimplexParams->eps, pSimplexParams->scale, constrain, (void *)pSimplexParams);

		*max_fit = fabs(pSimplexParams->fvec[0]);
		for(int i=1; i<n; i++) {
			double mfit = fabs(pSimplexParams->fvec[i]);
			if(mfit > *max_fit) *max_fit = mfit;
		}

		if(!func_param)          delete pSimplexParams;
		return pSimplexParams->fvec;

err_return:
		if(pSimplexParams->fvec) delete [] pSimplexParams->fvec;
		if(!func_param)          delete pSimplexParams;
		return NULL;
}
