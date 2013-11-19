//
// @(#)snls1m_dr.h  1.1  misha-13mar103
//
// Copyright (c) 2003 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  13mar103  by  misha@MISHA
// Version:  13mar103  by  misha@MISHA
//

#ifndef DNLS1M_DR_H
#define DNLS1M_DR_H

#include "fcnbase.h"

typedef /* Subroutine */ int (*S_fp)(int *, int *, int *, double *, double *, double *, int *, void *);


typedef struct _dnls1_storage : public _fcnbase {

  // initial values for parameters:

  double *x;

  // error indicator;

 // storage arrays and parameters for dnls1_

  int * ipvt;
  double * fvec, * fjac, * diag, * qtf, * wa1, * wa2, * wa3, * wa4;

  int n;
  int m;
  int iopt,ldfjac,  mode, nprint, info, nfev, njev, nwrite;
  double ftol, xtol, gtol, factor, fnorm;
  int err;

  int compute_only_stdev;

  // FCN function for non-linear fit

  S_fp fcn;

  //int n_iterations; // number of FCN calls for printing

  _dnls1_storage() {
    ipvt = NULL;
    fvec = fjac = diag = qtf = wa1 = wa2 = wa3 = wa4 = x = NULL;
    n = m = 0;
    iopt = ldfjac =  mode   = nprint = info = nfev = njev = nwrite = -1;
    ftol =  xtol  = gtol   = factor = fnorm  =  epsfcn = -1.;
    err = 1;
    fcn = NULL;
	compute_only_stdev = 0;
 }

  void destroy();
  void init(int m, int n);
  ~_dnls1_storage() { destroy();}

}  DNLS1_STORAGE;


// functions


int dnls1_general_driver(DNLS1_STORAGE * pStorage, void * pGeop, double * x, double * fit,
			 double * st_err, double *max_diff, char * message);

#define real    double
#define integer int

int dnls1_(S_fp fcn, integer *iopt, integer *m, integer *n,
    real *x, real *fvec, real *fjac, integer *ldfjac, real *ftol, real *
    xtol, real *gtol, integer *maxfev, real *epsfcn, real *diag, integer *
    mode, real *factor, integer *nprint, integer *info, integer *nfev,
    integer *njev, integer *ipvt, real *qtf, real *wa1, real *wa2, real *
    wa3, real *wa4, void *pParam);

    double r1mach_(int *);
    double d1mach_(int *);

int dfdjc3_(S_fp fcn, int *m, int *n, double *x, double *
	    fvec, double *fjac, int *ldfjac, int *iflag, double *epsfcn, double
	    *wa, void *pParam);

int dcov_(S_fp fcn, int *iopt, int *m, integer *n,
	  double *x, double *fvec, double *r, int *ldr, int *info,
	  double *wa1, double *wa2, double *wa3, double * Wa4, void *pParam);

extern "C" {
int sgetri_(int *n, double *a, int *lda, int *ipiv,
	    double *work, int *lwork, int *info);
int sgetrf_(int *m, int *n, double *a, int *lda,
	    int *ipiv, int *info);
}

int dz29svd_(int *nm, int *m, int *n, double *a,
	double *w, long int *matu, double *u, long int *matv, double *v, int *
	ierr, double *rv1);

int compute_standard_errors(double * fjack, int m, int n, int * ipiv, double *work,
			    int pack, double fit, double * st_err);

int compute_standard_errors_svd(double * fjack, int m, int n,
			                     double fit, double * st_err);

#endif
