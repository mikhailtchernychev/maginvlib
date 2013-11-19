//
// @(#)tn_estimator.h  1.1  misha-22sep103
//
// Copyright (c) 2003 Geometrics
//
// Created:  22sep103  by  misha@misha2.geometrics.com
// Version:  22sep103  by  misha@misha2.geometrics.com
//


#ifndef TN_ESTIMATOR_H
#define TN_ESTIMATOR_H

#include "inversion.h"
#include "tn_dr.h"
#include "dnls1m_dr.h"

double *  Estimator_tn(void * func_param, void * class_param,  int n, int m, double *x, 
			  double *x_up, double *x_low, double * fit, double *st_err, double *max_fit, char * message);

#endif
