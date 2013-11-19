//
// @(#)l2estimator.cpp  1.1  misha-16may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  16may103  by  misha@misha.local
// Version:  16may103  by  misha@misha.local
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "magfunc.h"
#include "l2estimator.h"

int dz29rls_(double *a, double *b, double *x, double *u, double *v,
   	    double *s, double *rv, double *vu, int *mn, int *m, int *n,
            double *eps, double *epsm, double *dm, double *am, int *k);

_z29rls_params::_z29rls_params() {

    eps=0.;
    epsm=0.;
    dm=1.e+20;
    am=1.e-20;
}


int Estimate_L2_rls(int n, int m, double *a, double *y, double *x, void * param,  double * synt, double * fit,
		    double * max_fit, double *st_err, char * message) {

    double eps, epsm, dm, am;
    int k;

    eps=0.;
    epsm=0.;
    dm=1.e+20;
    am=1.e-20;

    int ret = 0;

    if(param) {
	Z29RLS_PARAMS * p = (Z29RLS_PARAMS *)param;
	eps  = p->eps;
	epsm = p->epsm;
	dm   = p->dm;
	am   = p->am;
    }

    double *u, *v, *vu,  *s, *rv;

    int mn_rls = n;
    int m_rls  = n;
    int n_rls  = m;

    u = v = vu =  s = rv =  NULL;

    u      = (double *)malloc(n*m*sizeof(double));  if(!u)      goto err_return;
    v      = (double *)malloc(n*m*sizeof(double));  if(!v)      goto err_return;
    vu     = (double *)malloc(n*m*sizeof(double));  if(!vu)     goto err_return;
    s      = (double *)malloc(m*sizeof(double));    if(!s)      goto err_return;
    rv     = (double *)malloc(m*sizeof(double));    if(!rv)     goto err_return;

    // for(int i=0; i<n; i++) {
    // 	for(int j=0; j<m; j++) fprintf(stdout, "%le ", a[n*j+i]);
    //	fprintf(stdout,"\n%lf\n", y[i]);
    //	}

    dz29rls_(a, y, x, u,v,s,rv, vu, &mn_rls, &m_rls,
            &n_rls, &eps, &epsm, &dm, &am, &k);

   //for(int i=0; i<m; i++) fprintf(stdout, "s[%d]=%lf rv[%d]=%lf\n", i, s[i],i,rv[i]);

    // compute synthetic field and related values

    if(synt) {
      MatMultVect(a, x, n, m, synt);

      if(fit && max_fit) {
	compute_fit(n, y, synt, fit, max_fit);

	if(st_err) {
	  int i, j, ll;
	  double b1 = 0.;
	  for(i=0; i<m; i++)
	    for(j=0; j<m; j++) {
              b1 = 0.;
	      for(ll=0; ll<k; ll++) {
		b1 =  b1 + v[m*ll+i]*pow(1./s[ll], 2)*v[m*ll+j];
	      }
	      u[m*j+i] = b1;
	    }

	  for(i=0; i<m; i++) st_err[i] = (*fit)*sqrt(u[i*m+i]);

	} // end st_err
      } // end fit & max_fit
    } // end compute synt. field


    ret = 1;


 err_return:

 if(u)       { free(u);      u      = NULL; }
 if(v)       { free(v);      v      = NULL; }
 if(vu)      { free(vu);     vu     = NULL; }
 if(s)       { free(s);      s      = NULL; }
 if(rv)      { free(rv);     rv     = NULL; }

 return ret;

}
