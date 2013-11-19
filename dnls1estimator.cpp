//
// @(#)dnls1estimator.cpp  1.1  misha-08aug103
//
// Copyright (c) 2003 Mikhail Tchernychev
//
// Created:  08aug103  by  misha@misha2
// Version:  08aug103  by  misha@misha2
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "dnls1estimator.h"
#include "inversion.h"


// part general to dnls1

int FCN(int *iflag, int *m, int *n1, double *params, 
	            double *fvec, double *fjack, int *ldjack, void *p_in) {

  // this is general for each driver function

  assert(iflag);
  assert(m);
  assert(n1);
  assert(params);
  if(*iflag <= 1) assert(fvec);
  if(*iflag >  1) assert(fjack);
  assert(ldjack);
  assert(p_in);

  FCN_PARAMS * pParams = (FCN_PARAMS *)p_in;
  FCNBASE * pData = (FCNBASE *)pParams->pData;

  CInversion * pInversion = (CInversion *)pParams->pFcnParams;
  assert(pInversion);

  pInversion->SetAllVariables(params, *n1);

  switch(*iflag) { 

  case 0: // ----------------------- show progress only -------------------------------------
         { 
         PROGRESS_FUNC pFunc = pInversion->GetProgressFunc();
        if(pFunc) {
	  pData->n_iterations++;
	  // compute RMS of residual
	  double discrepancy = 0;
	  int n = *m;
	  if(n>0) {
		  for(int i=0; i<n; i++) discrepancy += fvec[i]*fvec[i];
		  discrepancy = sqrt(discrepancy/n);
		  //for(int i=0; i<*n1; i++) fprintf(stdout,"%lf ", params[i]); fprintf(stdout,"\n");
		  }
          if((*pFunc)((void *)pInversion, pData->n_iterations, discrepancy)) pData->maxfev = 0; // to cancel execution
	  //if(pData->n_iterations == 5) pData->maxfev = 0;
        }
    break;
    }
    
   case 1: { // --------------------- compute function values -----------------------------------

     int compute_difference = 1;
	 if(pInversion->GetTheoryStdDevFlag()) compute_difference = 0;

     if(pData->n_pos !=-1 && *m == 1) 
	   	fvec[0] = pInversion->ComputeFieldAtPosition( pData->n_pos, compute_difference);
     else 
   		pInversion->ComputeField(fvec, *m, compute_difference);
      } 

      //for(int i=0; i<*m; i++) fprintf(stderr,"fvec[%d]=%lf\n", i, fvec[i]);

	  break;

  case 2: { // --------------------- compute full Jacobian
        pInversion->GetWholeInversionMatrix(fjack,1);
        int n_total = (*m)*(*n1);
        for(int i=0; i<n_total; i++) fjack[i] = -fjack[i];
       /*for(int i=0; i< *m; i++) {
	   for(int j=0; j< *n1; j++) fprintf(stdout,"%le ", fjack[j*(*m)+i]);
	   fprintf(stdout,"\n");
	   }
	   fprintf(stdout,"----------------------------\n");
	 //exit(0);*/
        break;
      } 

  case 3: { // --------------------- compute one string of Jacobian
       pInversion->GetWholeInversionMatrixRow(fjack, (*ldjack)-1, *m);
       for(int i=0; i<*n1; i++) fjack[i] = -fjack[i]; 
       /*for(int i=0; i< *m; i++) {
	  pInversion->GetWholeInversionMatrixRow(fjack, i, *m);
	  for(int i=0; i<*n1; i++) { fjack[i] = -fjack[i];  fprintf(stdout,"%lf ", fjack[i]);}
	  fprintf(stdout, "\n");
	}
	exit(0);*/
        break;
      } 


    
  }

 
  return 1;
}



void check_size() {

 int dnls_size =  sizeof(DNLS1_STORAGE);
 int col_size   = sizeof(CMagDataCollection);

 return; 
}


// func_param   -  parameters for dnls1 routine; if NULL assigns automatically
//                 if not NULL, has to be initialized by the caller
// class_param  -  CInversion class
// n            -  number of data points
// m            -  number of variables
// x            -  varaibles vector
// x_up         -  upper limits (not used here)
// x_low        -  lower limits (not used here)
// fit          -  rms for final model fit
// st_err       -  standard error for estimates
// max_fit      -  max. difference between observed and comoputed fields
// message      -  message with exit condition from dnsl1


double *  Estimator_dnls1(void * func_param, void * class_param, int n, int m,  double *x, 
			  double * x_up, double *x_low, double * fit, 
			  double *st_err, double *max_fit, char * message) {

  if(!class_param) return NULL; 

  CInversion * pInversion = (CInversion *)class_param;

  DNLS1_STORAGE * p = NULL;

  int dnls_size =  sizeof(DNLS1_STORAGE);
  int fp_size =  sizeof(S_fp);

  double * z = new double[n];
  if(!z) return NULL;

  if(func_param)  
    p = (DNLS1_STORAGE *)func_param;
  else 
    p = new DNLS1_STORAGE;

   p->init(m, n);
   p->compute_only_stdev =  pInversion->GetTheoryStdDevFlag();

    if(p->err) {
      delete p;
      delete [] z;
      return NULL;
      }
 
    p->fcn = FCN;

    int iab = pInversion->GetTotalAbilities();
   
    if(iab==1 && iab < p->iopt)   p->iopt = iab;

    dnls1_general_driver(p, class_param, x, fit, st_err, max_fit, message);

    if(p->info == 0) {
         free(z);
         if(!func_param) delete p;
         return NULL;
    }

    pInversion->GetDataVector(z);

	if(!pInversion->GetTheoryStdDevFlag()) {
		for(int i=0; i<n; i++) z[i] -= p->fvec[i];
	}
	else {
	   for(int i=0; i<n; i++) z[i] = p->fvec[i];
	}

    if(!func_param) delete p;
    return z;

}


void print_jack(double * fjack, int m, int n) {
   for(int i=0; i< m; i++) {
     for(int j=0; j< n; j++) fprintf(stdout,"%lf ", fjack[j*m+i]);
     fprintf(stdout,"\n");
    }
       exit(0);
}









