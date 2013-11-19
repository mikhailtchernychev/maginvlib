//
// @(#)fcnbase.h  1.1  misha-22sep103
//
// Copyright (c) 2003 Mikhail Tchernychev
//
// Created:  22sep103  by  misha@misha2
// Version:  22sep103  by  misha@misha2
//


#ifndef FCNBASE_H
#define FCNBASE_H

#pragma pack (1)

typedef struct _fcnbase {
  
  int n_iterations; // number of FCN calls for printing 
  int maxfev;       // max. iterations or fynction evaluations
  int n_pos;        // position where to compute field (if only one point is computed)
  double epsfcn;    

  _fcnbase () { n_iterations = 0; maxfev = -1; n_pos=-1; epsfcn=-1.;}

} FCNBASE;


typedef struct _fcn_pass_params {
  FCNBASE * pData;
  void * pFcnParams;
} FCN_PARAMS;


// part general to dnls1 - tn

int FCN(int *iflag, int *m, int *n1, double *params, 
	double *fvec, double *fjack, int *ldjack, void *p_in);

#endif
