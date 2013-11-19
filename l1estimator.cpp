//
// @(#)l1estimator.cpp  1.1  misha-19may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  19may103  by  misha@misha.local
// Version:  19may103  by  misha@misha.local
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#ifndef TRUE
   #define TRUE 1
#endif

#ifndef FALSE
   #define FALSE 0
#endif


#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) (((x) < (y)) ? (y) : (x))
#endif

//#include "magfunc.h"
#include "l1estimator.h"

extern "C" {
int l1_driver_d(int m, int n, double * ct,  double *f, double *a,
	      double prec, double eps, int *irank, int * iter,
	      double *r, double *z, int *ind, int (*cancel_func)());
}

_l1_params::_l1_params() {

    prec = 1.e-16;
    eps = 1.e-11;
    irank = 0;
    iter  = 0;
    ind   = 0;
    z     = 0.;

	m_pCancelFunc = NULL;
}



int Estimate_L1(int n, int m, double *a_matrix, double *y, double *x, void * param,
		double * synt, double * fit, double * max_fit, double *st_err, char * message) {

    double prec, eps, z;
    int irank, iter, ind;
    int ret = 0;
	int (*cancel_func)() = NULL;

    prec = 1.e-16;
    eps = 1.e-11;
    irank = 0;
    iter  = 0;
    ind   = 0;
    z     = 0.;

    if(param) {
	L1_PARAMS * p = (L1_PARAMS *)param;
	prec  = p->prec;
	eps   = p->eps;
	irank = p->irank;
	iter  = p->iter;
	ind   = p->ind;
	z     = p->z;
	cancel_func = p->m_pCancelFunc;
    }


    double  *r, * a_old_matrix;

    r = a_old_matrix = NULL;

    r = (double *)malloc(n*sizeof(double)); if(!r) goto err_return;

    a_old_matrix  = (double *)malloc(n*m*sizeof(double)); if(!a_old_matrix) goto err_return;

    memcpy(a_old_matrix, a_matrix, n*m*sizeof(double));


    l1_driver_d(m, n, a_matrix,  y,  x, prec, eps, &irank, &iter, r, &z, &ind, cancel_func);

    ret = 1;

    if(message) {
      switch(ind) {
      case 0:
	sprintf(message,"L1: THE SOLUTION VECTOR IS UNIQUE. RANK=%d ITERATIONS=%d", irank, iter); break;
      case 1:
	sprintf(message,"L1: THE SOLUTION VECTOR IS MOST PROBABLY NOT UNIQUE. RANK=%d ITERATIONS=%d", irank, iter);
	break;
      default:
	sprintf(message,"L1: PREMATURE TERMINATION OF THE CALCULATION DUE TO VERY ILL-CONDITIONING OF MATRIX");
      }
    }


    if(param) {
	L1_PARAMS * p = (L1_PARAMS *)param;
	p->irank = irank;
	p->iter  = iter;
	p->ind   = ind;
	p->z     = z;
    }

    memcpy(a_matrix, a_old_matrix,  n*m*sizeof(double));

    if(synt) {
      double d =0.;
      double d_max = 0.;

      for(int i=0; i<n; i++) {
	  d += (r[i]*r[i]);
	  d_max = MAX(d_max, fabs(r[i]));
	  synt[i] = r[i] + y[i];
	}

        *fit     = sqrt(d/n);
	*max_fit = d_max;

       if(st_err) for(int i=0; i<m; i++) st_err[i] = 0;

      } // end synt

 err_return:
    if(r) { free(r); r = NULL; }
    if(a_old_matrix) { free(a_old_matrix); a_old_matrix = NULL; }


    return ret;

}






