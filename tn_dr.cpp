//
// @(#)tn_dr.cpp  1.1  misha-043may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  04may103  by  misha@MISHA
// Version:  04mar103  by  misha@MISHA
//

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include <string.h>

#include <math.h>

#include "tn_dr.h"


_truncated_newton_storage::_truncated_newton_storage() {
        x = g = w = low = up = T_synt = fjack = wa = NULL;
	ipivot = NULL;
	f = 1.;
	lw = 0;
	eta = 0.25;
	stepmx = 10.;
	accrcy = 1.e-15;
	xtol = sqrt(accrcy);
	msglvl = 1;
	err = 1;
	n = 0;
	n_data = 0;
	is_fdiff = 0;

	maxit =  10;
	maxfev = 100000;
}


void _truncated_newton_storage::destroy() {

    if(x)       { free(x);             x = NULL; }
    if(g)       { free(g);             g = NULL; }
    if(w)       { free(w);             w = NULL; }
    if(low)     { free(low);         low = NULL; }
    if(up)      { free(up);           up = NULL; }
    if(fjack)   { free(fjack);     fjack = NULL; }
    if(wa)      { free(wa);           wa = NULL; }
    if(ipivot)  { free(ipivot);   ipivot = NULL; }

    n  =  0;
    lw =  0;
    err = 1;
}

void _truncated_newton_storage::init(int n_in, int n_data_in) {

    int ret =0;
    int lw_in = 14*n_in;
    int i;


    if(n_in > n || err) {
	destroy();
         x   = (double *) malloc(n_in*sizeof(double));      if(!x)     goto error;
         g   = (double *) malloc(n_in*sizeof(double));      if(!g)     goto error;
         low = (double *) malloc(n_in*sizeof(double));      if(!low)   goto error;
         up  = (double *) malloc(n_in*sizeof(double));      if(!up)    goto error;
           w = (double *) malloc(lw_in*sizeof(double));     if(!w)     goto error;
       fjack = (double *) malloc(n_in*sizeof(double));      if(!fjack) goto error;
          wa = (double *) malloc(n_in*sizeof(double));      if(!wa)    goto error;

      ipivot = (int *) malloc(n_in*sizeof(int));            if(!ipivot) goto error;
    }

    n  =  n_in;
    lw = lw_in;
    
    err = 0;
    
    for(i=0; i<n; i++) {
	low[i] = -1.E38;
	up[i]  =  1.E38;
    }
    
    n_data = n_data_in;
    return;

 error:
    destroy();
}


int tn_general_driver(TRUNCATED_NEWTON_STORAGE * p, void * pParam, double * x, double *x_up, double *x_low, 
		                                 double * fit, double *st_err,
						 double * max_diff, char * message) {

    if(!p)      return 0;
    if(p->err)  return 0;
    if(!pParam) return 0;
    
     FCN_PARAMS params;
     params.pData = p;
     params.pFcnParams = pParam;
    
    *fit = -1;
    *max_diff = -1;
    

    for(int i=0; i<p->n; i++) st_err[i] = -1;
    
    memcpy(p->x,   x,      sizeof(double)*p->n);
    memcpy(p->up,  x_up,   sizeof(double)*p->n);
    memcpy(p->low, x_low,  sizeof(double)*p->n);


    /*p->x[0]=8.450000;
    p->x[1]=10.230000;
    p->x[2]=3.000000;
    p->x[3]=1441.249987;
    p->x[4]=-37.460175;
    p->x[5]=1228.875079;
    p->x[6]=-0.502384;
    p->x[7]=0.413681;
    p->x[8]=-4.928084;



    for(int i=0; i<p->n; i++) fprintf(stderr,"limits %e %e\n", p->low[i], p->up[i]);

    fprintf(stderr,"p->n=%d\n", p->n);    
    for(int i=0; i<p->n; i++) {
      fprintf(stderr,"p->x[%d]=%lf\n", i, p->x[i]);
      fprintf(stderr,"p->low[%d]=%le\n", i, p->low[i]);
      fprintf(stderr,"p->up[%d]=%le\n", i, p->up[i]);
    }

    fprintf(stderr,"p->j=%lf\n", p->f);
    fprintf(stderr,"p->lw=%d\n", p->lw);     
    fprintf(stderr,"p->msglvl=%d\n", p->msglvl);     
    fprintf(stderr,"p->maxit=%d\n", p->maxit);     
    fprintf(stderr,"p->maxfun=%d\n", p->maxfev);     
    fprintf(stderr,"p->eta=%lf\n", p->eta);     
    fprintf(stderr,"p->stepmx=%lf\n", p->stepmx);  
    fprintf(stderr,"p->accrcy=%le\n", p->accrcy);  
    fprintf(stderr,"p->xtol=%le\n", p->xtol);  

    */

lmqnbc_(&p->ifail, &p->n, p->x,  &p->f, p->g, p->w, &p->lw, p->sfun, 
	p->low, p->up, p->ipivot, &p->msglvl, &p->maxit, &p->maxfev, &p->eta, &p->stepmx, 
	&p->accrcy, &p->xtol, &params);

  *fit = sqrt(p->f/(p->n_data - p->n));

  *max_diff = fabs(p->T_synt[0]);
  for(int i=0; i<p->n_data; i++) {
	  double t_abs = fabs(p->T_synt[i]);
	  if(*max_diff < t_abs) * max_diff = t_abs;
  }


 switch(p->ifail) {
   case 0:
     strcpy(message, "TN:  NORMAL RETURN"); break;
   case 2:
     strcpy(message, "TN: MORE THAN MAXFUN EVALUATIONS OR CANCELED"); break;
   case 3:
     strcpy(message, "TN:  LINE SEARCH FAILED TO FIND LOWER POINT (MAY NOT BE SERIOUS)"); break;
 default:
     strcpy(message,"TN: ERROR IN INPUT PARAMETERS");
 } 

    memcpy(x, p->x, sizeof(double)*p->n);


 return 1;

}





