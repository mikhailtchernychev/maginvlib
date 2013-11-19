//
// @(#)hooke_estimator.cpp  1.1  misha-20aug2012
//
// Copyright (c) 2012 of Mikhail Tchernychev
// All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "hooke_estimator.h"

double hooke_func(double *x, int m, void * pParams) {

	if(!pParams) return 0.;
	HOOKE_PARAMS * pHookeParams = (HOOKE_PARAMS *)pParams;

	pHookeParams->pInversion->SetAllVariables(x, pHookeParams->m);
	pHookeParams->pInversion->ComputeField(pHookeParams->fvec, pHookeParams->n, 1);

	double discrepancy = 0;
   for(int i=0; i<pHookeParams->n; i++) {
	   discrepancy += pHookeParams->fvec[i]*pHookeParams->fvec[i];
  }

   discrepancy = sqrt(discrepancy/pHookeParams->n);

   return discrepancy;

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


double *  Estimator_Hooke(void * func_param, void * class_param,  int n, int m, double *x, 
	double *x_up, double *x_low, double * fit, double *st_err, double *max_fit, char * message) {

		HOOKE_PARAMS * pHookeParams = (HOOKE_PARAMS *)func_param;

		if(!func_param) {
			pHookeParams = new HOOKE_PARAMS;
			if(!pHookeParams) return NULL;
		}

		if(pHookeParams->fvec) delete [] pHookeParams->fvec;

		pHookeParams->fvec       =  new double [n];
		pHookeParams->pInversion = (CInversion *)class_param;
		pHookeParams->n          = n;
		pHookeParams->m          = m;

	    PROGRESS_FUNC pFunc = pHookeParams->pInversion->GetProgressFunc();
		double *x_out = new double[m];

		int ret = hooke(m, x, x_out, pHookeParams->hooke_rho, pHookeParams->hooke_eps, pHookeParams->hooke_niter, hooke_func, 
					    (void *)pHookeParams, pFunc);

		for(int i=0; i<m; i++) x[i] = x_out[i];
		delete [] x_out;

		*max_fit = fabs(pHookeParams->fvec[0]);
		*fit     = pow(pHookeParams->fvec[0],2);
		for(int i=1; i<n; i++) {
			double mfit = fabs(pHookeParams->fvec[i]);
			if(mfit > *max_fit) *max_fit = mfit;
			*fit += pow(pHookeParams->fvec[i],2);
		}

		*fit = sqrt(*fit/(n-m));

		if(!func_param)   delete pHookeParams;
		return pHookeParams->fvec;

err_return:
		if(pHookeParams->fvec) delete [] pHookeParams->fvec;
		if(!func_param)          delete pHookeParams;
		return NULL;
}
