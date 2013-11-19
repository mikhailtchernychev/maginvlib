//
// @(#)dnls1m_dr.cpp  1.1  misha-13mar103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  13mar103  by  misha@MISHA
// Version:  13mar103  by  misha@MISHA
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>

#include "dnls1m_dr.h"

//#define USE_SIMPLE_C_CODE
//#define OLD_COVARIANCE

double discrep(double *a, int n, int n_vars) {
  int i;
  double d;

  for(i=0, d=0.; i<n; i++) d += (a[i]*a[i]);
  return sqrt(d/(double)(n-n_vars));
}


void MatPrint ( double *m1 , int m , int n , char *t )
{
    int i , j ;

    printf ( "%s\n\n" , t ) ;
    for ( i = 0 ; i < m ; ++i )
    {
        for ( j = 0 ; j < n ; ++j ) printf( "%10.3f  " , m1[i*n+j] ) ;
        printf ( "\n" ) ;
    }
    printf ( "\n" ) ;
} /* MathPrint */

#ifdef USE_SIMPLE_C_CODE
void MatPrintD ( double *m1 , int m , int n , char *t )
{
    int i , j ;

    printf ( "%s\n\n" , t ) ;
    for ( i = 0 ; i < m ; ++i )
    {
        for ( j = 0 ; j < n ; ++j ) printf( "%le  " , m1[i*n+j] ) ;
        printf ( "\n" ) ;
    }
    printf ( "\n" ) ;
} /* MathPrint */


void matswap ( double *s1 , double *s2 )
{
    double t ;

    t = *s1 ;
    *s1 = *s2 ;
    *s2 = t ;
} /* matswap */


void MatInvert ( double  *m , int c , double *d , double *a )
{
    int  *pl ;
    char *pc ;
    int i , j , k , l , l1 , lr , lc ;
    double  p , t , v;

    pl = (int *)calloc (1000 , 2 ) ;
    pc = (char *)calloc ( 500 , 2 ) ;
    *d = 1 ;
    for ( i = 0 ; i < c ; ++i )
    {
        pc[i] = 0 ;
        for ( j = 0 ; j < c ; ++j ) a[i*c+j] = m[i*c+j] ;
    }
    for ( i = 0 ; i < c ; ++i )
    {
        v = 0.0 ;
        for ( j = 0 ; j < c ; ++j )
        {
            if ( ! ( pc[j] ) )
            {
                for ( k = 0 ; k < c ; ++k )
                {
                    if ( ! ( pc[k] ) )
                    {
                        if ( fabs ( a[j*c+k] ) > v )
                        {
                            lr = j ;
                            lc = k ;
                            v = fabs ( a[j*c+k] ) ;
                        }
                    }
                }
            }
        }
        pc[lc] = 1 ;
        pl[i<<1] = lr ;
        pl[i*2+1] = lc ;
        if ( lr != lc )
        {
            *d = - *d ;
            for ( l = 0 ; l < c ; ++l ) matswap ( &a[lr*c+l] , &a[lc*c+l] ) ;
        }
        p = a[lc*c+lc] ;
        *d = *d * p ;
        if ( *d > 1.0e+30 ) *d = 1 ;
        a[lc*c+lc] = 1.0 ;
        for ( l = 0 ; l < c ; ++l ) a[lc*c+l] = a[lc*c+l] / p ;
        for ( l1 = 0 ; l1 < c ; ++l1 )
        {
            if ( l1 != lc )
            {
                t = a[l1*c+lc] ;
                a[l1*c+lc] = 0 ;
                for ( l = 0 ; l < c ; ++l )
                 a[l1*c+l] = a[l1*c+l] - a[lc*c+l] * t ;
            }
        }
    }
    for ( i = 0 ; i < c ; ++i )
    {
        l = c - i - 1 ;
        if ( pl[l*2] != pl[l*2+1] )
        {
            lr = pl[l*2] ;
            lc = pl[l*2+1] ;
            for ( k = 0 ; k < c ; ++k ) matswap ( &a[k*c+lr] , &a[k*c+lc] ) ;
        }
    }
    free ( pl ) ;
    free ( pc ) ;
} /* MatInvert */

#endif




void _dnls1_storage::destroy() {

  if (ipvt) { free(ipvt); ipvt = NULL; }
  if (fvec) { free(fvec); fvec = NULL; }
  if (fjac) { free(fjac); fjac = NULL; }
  if (diag) { free(diag); diag = NULL; }
  if (qtf)  { free(qtf);  qtf  = NULL; }
  if (wa1)  { free(wa1);  wa1  = NULL; }
  if (wa2)  { free(wa2);  wa2  = NULL; }
  if (wa3)  { free(wa3);  wa3  = NULL; }
  if (wa4)  { free(wa4);  wa4  = NULL; }
  if (x)    { free(x);    x    = NULL; }

  n = 0;
  m = 0;
  err = 1;
  n_iterations = 0;

}

void _dnls1_storage::init(int n_in, int m_in) {

	if(n_in > n || m_in > m)  {
		destroy();
		ipvt = (int *)malloc(n_in*sizeof(int));            if (!ipvt) goto error;
		fvec = (double *)malloc(m_in*sizeof(double));      if (!fvec) goto error;
		fjac = (double *)malloc(m_in*n_in*sizeof(double)); if (!fjac) goto error;
		diag = (double *)malloc(n_in*sizeof(double));      if (!diag) goto error;
		qtf  = (double *)malloc(n_in*sizeof(double));      if (!qtf)  goto error;
		wa1  = (double *)malloc(n_in*sizeof(double));      if (!wa1)  goto error;
		wa2  = (double *)malloc(n_in*sizeof(double));      if (!wa2)  goto error;
		wa3  = (double *)malloc(n_in*sizeof(double));      if (!wa3)  goto error;
		wa4  = (double *)malloc(m_in*sizeof(double));      if (!wa4)  goto error;
		x    = (double *)malloc(m_in*sizeof(double));      if (!x)    goto error;
	}

  n = n_in;
  m = m_in;
  err = 0;
  n_iterations = 0;
  ldfjac = m;

  compute_only_stdev = 0;

  return;

 error:
  destroy();

}





int dnls1_general_driver(DNLS1_STORAGE * pStorage, void * pGeop, double * x, double * fit, double *st_err,
						 double * max_diff, char * message) {

  int iw = 4, i, err = 0;
  double dummy_fit = 0.;

  if(!pGeop)    return 0;
  if(!pStorage) return 0;

  if(pStorage->err) return -1;


  if(pStorage->iopt==-1)   pStorage->iopt=1;
  if(pStorage->ftol==-1.)  pStorage->ftol= sqrt(d1mach_(&iw));
  if(pStorage->xtol==-1.)  pStorage->xtol= sqrt(d1mach_(&iw));
  if(pStorage->gtol==-1.)  pStorage->gtol= 0.;
  if(pStorage->maxfev==-1) pStorage->maxfev= 40000;
  if(pStorage->epsfcn==-1) pStorage->epsfcn= 0.;
  if(pStorage->mode  ==-1) pStorage->mode= 1;
  if(pStorage->factor==-1) pStorage->factor = 100;
  if(pStorage->nprint==-1) pStorage->nprint = 1;


  pStorage->ldfjac       = pStorage->m;
  pStorage->n_iterations = 0;

  memcpy(pStorage->x, x, sizeof(double)*pStorage->n);

  FCN_PARAMS params;
  params.pData = pStorage;
  params.pFcnParams = pGeop;


/*  if(pStorage->iopt==1) { // forward differencing for jackobian
    int iflag=1;
    memset(pStorage->fjac, sizeof(double)*pStorage->n*pStorage->m, 0);

   (*pStorage->fcn)(&iflag, &pStorage->m, &pStorage->n, pStorage->x,  // need to compute fvec
			  pStorage->fvec, pStorage->fjac, &pStorage->ldfjac, (void *)&params);


    dfdjc3_(pStorage->fcn, &pStorage->m, &pStorage->n, pStorage->x,
    		pStorage->fvec, pStorage->fjac, &pStorage->ldfjac, &iflag, &pStorage->epsfcn,
    		pStorage->wa4, (void *)&params);
    		fprintf(stderr,"%d %d %d\n", pStorage->m, pStorage->n, pStorage->ldfjac);
        }
     else {
       int iflag=2;
       int ldfjac = 1;
 	   (*pStorage->fcn)(&iflag, &pStorage->m, &pStorage->n, pStorage->x,
			  pStorage->fvec, pStorage->fjac, &pStorage->ldfjac, (void *)&params);

       }



       fprintf(stderr,"pStorage->n=%d pStorage->m=%d\n", pStorage->n, pStorage->m);

        for(i=0; i<pStorage->m; i++) {
                 fprintf(stdout,"%d " ,i);
                 for(int j=0; j<pStorage->n; j++)
                      fprintf(stdout,"%lg ", pStorage->fjac[pStorage->m*j+i]);
                  fprintf(stdout,"\n");
              }

  exit(0);  */

  /*fprintf(stderr,"Before---------------------\n");
  fprintf(stderr,"pStorage->ldfjac=%d\n",  pStorage->ldfjac);
  fprintf(stderr,"pStorage->ftol=%lg\n",   pStorage->ftol);
  fprintf(stderr,"pStorage->xtol=%lg\n",   pStorage->xtol);
  fprintf(stderr,"pStorage->gtol=%lg\n",   pStorage->gtol);
  fprintf(stderr,"pStorage->maxfev=%d\n",  pStorage->maxfev);
  fprintf(stderr,"pStorage->epsfcn=%lg\n", pStorage->epsfcn);
  fprintf(stderr,"pStorage->mode=%d\n",    pStorage->mode);
  fprintf(stderr,"pStorage->factor=%lg\n", pStorage->factor);
  fprintf(stderr,"pStorage->nprint=%d\n",  pStorage->nprint);*/

  if(!pStorage->compute_only_stdev) {

	  dnls1_(pStorage->fcn, &pStorage->iopt, &pStorage->m, &pStorage->n, pStorage->x,
		  pStorage->fvec, pStorage->fjac, &pStorage->ldfjac, &pStorage->ftol, &pStorage->xtol,
		  &pStorage->gtol, &pStorage->maxfev, &pStorage->epsfcn, pStorage->diag,
		  &pStorage->mode, &pStorage->factor, &pStorage->nprint, &pStorage->info,
		  &pStorage->nfev, &pStorage->njev, pStorage->ipvt, pStorage->qtf,
		  pStorage->wa1, pStorage->wa2, pStorage->wa3, pStorage->wa4, (void *)&params);

	  if(message) {
		  switch(pStorage->info) {
		  case 0: strcpy(message, "DNLS1: improper input parameters."); break;

		  case 1: sprintf(message, "DNLS1: both actual and predicted relative reductions in the sum of squares are at most FTOL=%lg",pStorage->ftol);
			  break;

		  case 2: sprintf(message, "DNLS1: relative error between two consecutive iterates is at most XTOL=%lg",pStorage->xtol);
			  break;

		  case 3: sprintf(message, "DNLS1:  both actual and predicted relative reductions in the sum of squares are at most FTOL=%lg and relative error between two consecutive iterates is at most XTOL=%lg", pStorage->ftol, pStorage->xtol);
			  break;

		  case 4: sprintf(message, "DNLS1: the cosine of the angle between FVEC and any column of the Jacobian is at most GTOL=%lg in absolute value.",pStorage->gtol);
			  break;

		  case 5: strcpy(message, "DNLS1: number of calls to FCN for function evaluation has reached its limit or executaion has been cancelled");
			  break;

		  case 6: sprintf(message, "DNLS1: FTOL=%lg is too small.  No further reduction in the sum of squares is possible.", pStorage->ftol);
			  break;

		  case 7: sprintf(message, "DNLS1: XTOL=%lg is too small.  No further improvement in the approximate solution X is possible",pStorage->xtol);
			  break;

		  case 8: sprintf(message, "DNLS1: GTOL=%lg is too small.  FVEC is orthogonal to the columns of the Jacobian to machine precision.", pStorage->xtol);
			  break;

		  default: strcpy(message, "DNLS1: No information is available.");

		  }
	  }

	  *fit = discrep(pStorage->fvec,pStorage->m, pStorage->n);

	  // find max. difference

	  *max_diff = fabs(pStorage->fvec[0]);
	  double am;
	  for(i=0; i<pStorage->m; i++) {
		  am = fabs(pStorage->fvec[i]);
		  if(am > *max_diff) *max_diff = am;
	  }
  }
  else { // only estimate std. dev
	//*fit      = 0.;
	//*max_diff = 0.;
  }

   memcpy(x, pStorage->x, sizeof(double)*pStorage->n);


   //fprintf(stdout,"-----------------scov matrix ---------------\n");
   //for(i=0; i<pStorage->n; i++)  {
   //	   for(int j=0; j<pStorage->n; j++) {
   //			   fprintf(stdout,"%10e ", r[ldr*j+i]/((*fit)*(*fit)));
   //	   }
   //	   fprintf(stdout,"\n");
   //}
   //fprintf(stdout,"-----------------end of scov matrix ---------------\n");


// old way how to compute covariance

#ifndef OLD_COVARIANCE

  // compute standard covariance matrix

   int ldr =  (pStorage->iopt == 3) ? pStorage->n :  pStorage->m;
   int info;

   int iflag=1;


   dcov_(pStorage->fcn, &pStorage->iopt, &pStorage->m, &pStorage->n,
		 pStorage->x, pStorage->fvec, pStorage->fjac, &ldr, &info,
		 pStorage->wa1, pStorage->wa2, pStorage->wa3, pStorage->wa4, (void *)&params);


   if(pStorage->compute_only_stdev) { // only could find discr AFTER dcov_
	 dummy_fit = discrep(pStorage->fvec,pStorage->m, pStorage->n);
	 for(i=0; i<pStorage->n; i++) st_err[i] =  sqrt(fabs(pStorage->fjac[ldr*i+i]))/(dummy_fit);
   }
   else
	 for(i=0; i<pStorage->n; i++) {
		 //fprintf(stderr,"Jc %lf\n", pStorage->fjac[ldr*i+i]);
		 st_err[i] =  sqrt(fabs(pStorage->fjac[ldr*i+i]));
	 }


   for(i=0; i<pStorage->n; i++) {
     for(int j=0; j<pStorage->n; j++) {
         //fprintf(stderr, "Jc %d %d %6.2lf ", i, j, pStorage->fjac[ldr*i+j]);
         pStorage->fjac[ldr*i+j] /= (st_err[i]*st_err[j]);
     }
     //fprintf(stderr, "\n");
   }

#else

  // compute jackobian


  int pack = 1;

  if(pStorage->iopt==1) { // forward differencing for jackobian
    int iflag=1;
    pack =1;

  (*pStorage->fcn)(&iflag, &pStorage->m, &pStorage->n, pStorage->x,  // need to compute fvec
			  pStorage->fvec, pStorage->fjac, &pStorage->ldfjac, (void *)&params);

    dfdjc3_(pStorage->fcn, &pStorage->m, &pStorage->n, pStorage->x,
    		pStorage->fvec, pStorage->fjac, &pStorage->ldfjac, &iflag, &pStorage->epsfcn,
    		pStorage->wa4, (void *)&params);
        }
     else {

       int iflag=2;
       int ldfjac = 1;
 	   (*pStorage->fcn)(&iflag, &pStorage->m, &pStorage->n, pStorage->x,
			  pStorage->fvec, pStorage->fjac, &ldfjac, (void *)&params);

       }


  compute_standard_errors(pStorage->fjac, pStorage->m, pStorage->n, pStorage->ipvt,
   			  pStorage->wa1, pack, *fit, st_err);

  //compute_standard_errors_svd(pStorage->fjac, pStorage->m, pStorage->n, *fit,
  //                            st_err);

#endif


  return 1;

}



// input: fjack - jacobian
// m - number of data
// n - number of variables ( m > n);
// ipiv - service array dim. n
// work - service array dim. n
// st_err = standard errors

int compute_standard_errors(double * fjack, int m, int n, int * ipiv, double *work,
			    int pack, double fit, double * st_err) {


#ifdef USE_SIMPLE_C_CODE
  double * xTx_inv = new double[n*n];
  if(!xTx_inv) return -1;
  double * xTx     = new double[n*n];
#else
  double * xTx = new double[n*n];
#endif

  int i,j,k;

  if(!xTx) return -1;

  for(i=0; i<n*n; i++) xTx[i] = 0.;

  if(pack==0) {
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
	for(k=0; k<m; k++) {
	  xTx[n*j+i] += fjack[n*k+i]*fjack[n*k+j];
	}
  }
  else {
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
	for(k=0; k<m; k++) {
	  xTx[n*j+i] += fjack[m*i+k]*fjack[m*j+k];
	}
  }

  /*fprintf(stderr,"-----------------xTx ----------------\n");

   for(i=0; i<n; i++) {
      for(j=0; j<n; j++)
        fprintf(stderr,"%lg ",xTx[n*j+i]);
      fprintf(stderr,"\n");
      }*/

  // invert matrix

#ifdef USE_SIMPLE_C_CODE
   double d;
   MatInvert(xTx , n , &d , xTx_inv );
   if(d==0.) return -1;
   fprintf(stderr,"det=%lf\n", d);
   memcpy(xTx, xTx_inv, n*n*sizeof(double));
   MatPrintD(xTx_inv, n, n, "covar. matrix\n");
   delete [] xTx_inv;
#else
   int lwork = n;
   int info;

   fprintf(stderr,"compute_standard_errors\n");

   sgetrf_(&n, &n, xTx, &n, ipiv,  &info);
   if( info != 0) {
     delete [] xTx;
     fprintf(stderr,"Point 1 %d\n", info);
     return -1;
   }

   sgetri_(&n,xTx,&n,ipiv,work,&n,&info);
   if( info != 0) {
     delete [] xTx;
     fprintf(stderr,"Point 2\n");
     return -1;
   }
#endif

   for(i=0; i<n; i++) {
     st_err[i] = sqrt(xTx[n*i+i])*fit;
     fprintf(stderr,"st_err[%d]=%lg\n", i, st_err[i]);
   }

  delete [] xTx;
  return 1;

}


// input: fjack - jacobian
// m - number of data
// n - number of variables ( m > n);
// st_err = standard errors

// column packed (FORTRAN default)

int compute_standard_errors_svd(double * fjack, int m, int n,
			                     double fit, double * st_err) {


   int ret = -1;
   double * u, *v, *w, *rv1;
   int i;

   u = v = w = rv1 = NULL;

   u   = new double[m*n];
   v   = new double[m*n];
   w   = new double[n];
   rv1 = new double[n];

   if(!u || !v || !w ||  !rv1) goto mret;

   long int matu, matv;
   int ierr;

    matu = matv = 1;

    dz29svd_(&m, &m, &n, fjack, w, &matu, u, &matv, v, &ierr, rv1);

    if(ierr) goto mret;

    fprintf(stderr,"singular numbers:\n");
    for(i=0; i<n; i++)
        fprintf(stderr,"Sing[%d]=%lg\n", i, w[i]);

mret:
    if(u)   delete [] u;
    if(v)   delete [] v;
    if(w)   delete [] w;
    if(rv1) delete [] rv1;

    return ret;


}

