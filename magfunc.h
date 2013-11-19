//
// @(#)magfunc.h  1.1  misha-16may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  16may103  by  misha@misha.local
// Version:  16may103  by  misha@misha.local
//

#ifndef MAGFUNC_H
#define MAGFUNC_H


#include <vector>
#include <algorithm>

// linear estimator
typedef  int (*LINEAR_ESTIMATOR_FUNC)(int n, int m, double *a, double *y, double *x, void * param, double * synt,
				      double * fit,  double * max_fit, double *st_err, char * message);

// non-linear (or all parameters) estimator
typedef  double * (*NONLINEAR_ESTIMATOR_FUNC)(void * func_param, void *maginv_param, int n, int m, double *x,
					      double *x_up, double *x_low, double *fit, double *st_err,
					      double *max_fit, char * message);

// progress function
typedef int (*PROGRESS_FUNC)(void * data, int percent, double discrepancy);

// time functions

double date2sec(char * date_txt, char * time_txt, int format);
int sec2date (double dsec, char * date_txt, char * time_txt, int format);

void * LoadSurfer7Grid(char * file_name, long * nx, long *ny,
			 double * x1, double *y1,
			 double *dx, double *dy,
		         double *zmin, double *zmax, int is_float);


void MatMultVect(double *a, double *b, int n, int m, double *c);
void MatMultVectRow(double *a, double *b, int n, int m, double *c);
void compute_fit(int n, double * y, double * synt, double * fit, double * max_fit);

// mdedian

void CalculateMode(double *dataset, int numrow, double *mode);

// spline smooth

int splik(int n, double * x, double  * u,  double  dp,
          double * a, double * b, double  * c, double  * d, float *ex=NULL);

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


#ifndef NO_DATA
	#ifdef wx_x
	#define NO_DATA 1.e35
	#else
	#define NO_DATA 1.70141e38
	#endif
#endif


// attention - this sorts the vector

template <class T> T median_from_vector(std::vector<T> v) {

  int ss = v.size();

  if(!v.size()) return 0;

  std::sort(v.begin(), v.end());

  if(v.size()%2) return v[v.size()/2];


  return ( v[v.size()/2-1] + v[v.size()/2]) / 2.;

}



typedef double (*HOOKE_FUNC)(double *x, int n, void * param);

extern "C" {
int  hooke(int nvars, double * startpt, double * endpt, double rho, double epsilon, 
	   int itermax, HOOKE_FUNC func, void *param, PROGRESS_FUNC progress_func);
}



#endif
