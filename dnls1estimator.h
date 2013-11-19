//
// @(#)dnls1estimator.h  1.1  misha-08aug103
//
// Copyright (c) 2003 Mikhail Tchernychev
//
// Created:  08aug103  by  misha@misha2
// Version:  08aug103  by  misha@misha2
//


#ifndef DNLS1ESTIMATOR_H
#define DNLS1ESTIMATOR_H

#include "inversion.h"
#include "dnls1m_dr.h"

double *  Estimator_dnls1(void * func_param, void * class_param,  int n, int m, double *x, 
			  double *x_up, double *x_low, double * fit, double *st_err, double *max_fit, char * message);

#endif
