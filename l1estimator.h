//
// @(#)l1estimator.h  1.1  misha-19may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  19may103  by  misha@misha.local
// Version:  19may103  by  misha@misha.local
//


#ifndef L1ESTIMATOR_H
#define L1ESTIMATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#pragma pack(1)

typedef struct _l1_params {

    double prec;  // THE ROUND-OFF LEVEL OF THE COMPUTER. FOR THE IBM
                  // 360/67 COMPUTER, PREC IS ABOUT 1.E-6 AND 1.E-16 FOR
                  // SINGLE AND DOUBLE PRECISION RESPECTIVELY.
    double eps;   // A SPECIFIED TOLERANCE SUCH THAT A CALCULATED NUMBER
                  // X IS CONSIDERED ZERO IF  ABS(X) <EPS. FOR THE IBM 360/67
                  // COMPUTER, EPS IS USUALLY TAKEN 1.E-4 AND 1.E-11 RESPECTIVELY.  

    int irank;    // RANK  THE CALCULATED RANK OF MATRIX C.
    int iter;     // THE NUMBER OF ITERATIONS WHICH THE SOLUTION NEEDED.
    int ind;      // A RETURN INDICATOR. IND=0 INDICATES THAT THE
                  // SOLUTION VECTOR A* IS UNIQUE. IND=1, INDICATES THAT
                  // A* IS MOST PROBABLY NOT UNIQUE. IND=-1 INDICATES
                  // PREMATURE TERMINATION OF THE CALCULATION DUE TO
                  // VERY ILL-CONDITIONING OF MATRIX C.

    double z;     // THE MINIMUM L1 NORM OF THE RESIDUAL VECTOR R.

	int (*m_pCancelFunc)();

    _l1_params();
} L1_PARAMS;



int Estimate_L1(int n, int m, double *a, double *y, double *x, void * param, double * synt = NULL, 
		double * fit = NULL, double * max_fit = NULL, double *st_err = NULL, char * message = NULL);


#endif
