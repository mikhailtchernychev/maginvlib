//
// @(#)tn_dr.h  1.1  misha-043may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  04may103  by  misha@MISHA
// Version:  04mar103  by  misha@MISHA
//

#ifndef TN_DR_H
#define TN_DR_H

#include "fcnbase.h"

typedef /* Subroutine */ int (*tn_fp)(int *n, double *x, double *f, double *g, void *pParam);

typedef struct _truncated_newton_storage : public _fcnbase {

    int ifail;
    int n;
    int n_data;
    int is_fdiff;  // use finite difference to appoximate jacobian
    double * x, f, *g, *w, *T_synt, *fjack, *wa;
    int lw;
    tn_fp sfun;
    double * low, *up;
    int * ipivot;
    int msglvl;
    int maxit;
    double eta;
    double stepmx;
    double accrcy;
    double xtol;
    int err;

    _truncated_newton_storage();
    ~_truncated_newton_storage() { destroy();};
    void destroy();
    void init(int n_in, int ndata_in);

} TRUNCATED_NEWTON_STORAGE; 


// functions

extern "C" {

int lmqnbc_(int *ifail, int *n, double *x, 
	    double *f, double *g, double *w, int *lw, tn_fp sfun, 
	    double *low, double *up, int *ipivot, int *msglvl, 
	    int *maxit, int *maxfun, double *eta, double *stepmx, 
	    double *accrcy, double *xtol, void *pParam);

}

int tn_general_driver(TRUNCATED_NEWTON_STORAGE * p, void * pParam, 
		      double * x, double *x_up, double * x_low, double * fit, double *st_err,
		      double * max_diff, char * message) ;



#endif










