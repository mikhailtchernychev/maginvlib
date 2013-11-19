//
// @(#)tn_estimator.cpp  1.1  misha-22sep103
//
// Copyright (c) 2003 Geometrics
//
// Created:  22sep103  by  misha@misha2.geometrics.com
// Version:  22sep103  by  misha@misha2.geometrics.com
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "tn_estimator.h"

int FCN(int *iflag, int *m, int *n1, double *params, 
	double *fvec, double *fjack, int *ldjack, void *p_in);


int tn_fcn(int *n, double *x, double *f, double *g, void *p_in) {


    FCN_PARAMS * pParam = (FCN_PARAMS *)p_in;
    TRUNCATED_NEWTON_STORAGE * p = (TRUNCATED_NEWTON_STORAGE *)pParam->pData;

    int i, iflag;
    int ldjack;

    iflag = 1;

    int m =  p->n_data;
    FCN(&iflag, &m, n, x, p->T_synt, NULL, &ldjack, pParam);

    *f = 0.;

    for(i=0; i<m; i++) {
      *f += p->T_synt[i]*p->T_synt[i];
    }

    iflag = 0;
    FCN(&iflag, &m, n, x, p->T_synt, NULL, &ldjack, pParam);


    int n1=1;

    for(i=0; i<*n; i++) g[i] = 0.;

    for(i=0; i<m; i++) {
	ldjack = i+1;

	if(p->is_fdiff) {
	  iflag = 1;
	  p->n_pos = i;
	  dfdjc3_(FCN, &n1, n, x, &p->T_synt[i], p->fjack, &n1, &iflag, &p->epsfcn, p->wa, pParam);  
	}
	else {
	  iflag = 3;
	  FCN(&iflag, &m, n, x, NULL, p->fjack, &ldjack, pParam);
	}
 
	for(int j=0; j<*n; j++) {
	  g[j] += 2.*(p->T_synt[i])*(p->fjack[j]);
	}

    }

    p->n_pos=-1;

    return 1;
}
   
  



// func_param   -  parameters for dnls1 routine; if NULL assigns automatically
//                 if not NULL, has to be initialized by the caller
// class_param  -  CInversion class
// n            -  number of data points
// m            -  number of variables
// x            -  varaibles vector
// x_up         -  upper limits (not used here)
// x_low        -  lower limits (not used here)
// fit          -  rms for final model fit
// st_err       -  standard error for estimates
// max_fit      -  max. difference between observed and comoputed fields
// message      -  message with exit condition from dnsl1


double *  Estimator_tn(void * func_param, void * class_param, int n, int m,  double *x, 
			  double * x_up, double *x_low, double * fit, 
			  double *st_err, double *max_fit, char * message) {

  if(!class_param) return NULL; 

  CInversion * pInversion = (CInversion *)class_param;

  TRUNCATED_NEWTON_STORAGE * p = NULL;

  double * z = new double [n];
  if(!z) return NULL;

  if(func_param)  
    p = (TRUNCATED_NEWTON_STORAGE *)func_param;
  else 
    p = new TRUNCATED_NEWTON_STORAGE;


   p->init(m, n);

    if(p->err) {
      delete p;
      delete [] z;
      return NULL;
      }
 
    p->sfun = tn_fcn;
    p->T_synt =  z;

    tn_general_driver(p, class_param, x, x_up, x_low, fit, st_err, max_fit, message);

    p->T_synt =  NULL;

    pInversion->ComputeField(z, n);

    if(!func_param) delete p;
    return z;

}














