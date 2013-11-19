//
// @(#)magtest.cpp  1.1  misha-10may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  10may103  by  misha@misha.local
// Version:  10may103  by  misha@misha.local
//

#include <stdio.h>
#include <stdlib.h>

#include "magobject.h"
#include "magdipole.h"
#include "magcollection.h"
#include "maggriddata.h"
#include "maginv.h"
#include "magfunc.h"
#include "l2estimator.h"
#include "l1estimator.h"
#include "dnls1estimator.h"
#include "tn_estimator.h"
#include "magbackground.h"

#include <signal.h>

static int ibreak = 0;

void ctrl_c(int sig) {
  ibreak = 1;
  fprintf(stderr,"signal\n");
}

int progress_func(void * data, int percent, double discr) { // progress / cancel function

  fprintf(stderr,"iter=%d discrepancy=%lf\n", percent, discr); fflush(stderr);
  return ibreak;

}

main(int argc, char ** argv) {

  long nx, ny;
  double x1, y1, dx, dy, zmin, zmax;
  double params[100];
  double std_err[100];
  double vars[100];

  signal(SIGINT,ctrl_c);

  double * z = (double *) LoadSurfer7Grid(argv[1], &nx, &ny, &x1, &y1, &dx, &dy, &zmin, &zmax, 0);
  fprintf(stderr,"x1=%lf y1=%lf\n", x1, y1);
  CMagGridData * data = new CMagGridData(nx, ny, x1, y1, 0., dx, dy, z);

  CMagDataCollection data_collection;
  data_collection.AddData(data);

  int is_induced = 0;

  CMagCollection col;

  while(!feof(stdin)) {
    fscanf(stdin,"%lf%lf%lf\n", &params[0], &params[1], &params[2]);
    CMagDipole *pDipole = new CMagDipole(is_induced);
    //pDipole->FixNonLinearParam(0, TRUE, params[0]); 
    //pDipole->FixNonLinearParam(1, TRUE, params[1]); 
    //pDipole->FixNonLinearParam(2, TRUE, params[2]); 

    pDipole->SetNonLinearParams(params, 3);
    col.AddMagObject(pDipole);
    }


  double I = 67.8879;
  double D = -0.4387;
  double Az = 90.;

  // background

  CMagBackGround * back = new CMagBackGround;
  back->FixLinearParam(2, TRUE, 0.);

  col.AddMagObject(back);



  CMagInv inversion;
  inversion.SetProgressFunc(progress_func);


  inversion.SetMagData(&data_collection);
  inversion.SetMagCollection(&col);
  inversion.SetEarthFieldDirections(I,D,Az);
  int n_total_data, n_linear_vars, n_nonline_vars;
  inversion.GetTotalParameters(&n_total_data, &n_linear_vars, &n_nonline_vars);


  fprintf(stderr,"data points %d linear params %d non linear params %d\n",
	  n_total_data, n_linear_vars, n_nonline_vars);


  double fit;
  Z29RLS_PARAMS rls_params;
  rls_params.dm = 1.e+15;

  inversion.EstimateLinearParameters(Estimate_L2_rls, &fit, 1, &rls_params);
  fprintf(stderr,"fit=%lf\n", fit);

  col.ResetPosition();
  int nonlin_start, lin_start;
  CMagObject * object = col.GetCurrentMagObject(&nonlin_start, &lin_start);
 
  inversion.DumpObjects(stderr);
  //exit(0);
  //double pars[4]; pars[0] =pars[1] = 0; pars[2] = 0.;
  //back->SetLinearParams(pars, 4, 0);

  
  // my CVS test
  
   while(object) {
  	int m1 = object->GetTotalLinearParams();
  	object->GetLinearParams(vars, 100);
	for(int j=0; j<m1; j++) fprintf(stderr,"%le ", vars[j]);
	fprintf(stderr,"\n");
	object = col.GetCurrentMagObject(&nonlin_start, &lin_start);
	}


   // now non-linear part

   DNLS1_STORAGE dstorage;
   dstorage.iopt =  3;
   TRUNCATED_NEWTON_STORAGE tn;
   tn.is_fdiff = 0;
   char message[1024];
   //inversion.EstimateAllParameters(Estimator_dnls1, &dstorage, message);
   inversion.EstimateAllParameters(Estimator_tn, &tn, message);

   col.ResetPosition();
   object = col.GetCurrentMagObject(&nonlin_start, &lin_start);

   fprintf(stderr,"------------ non-linear vars -----------------\n");
   while(object) {
  	int m1 = object->GetTotalNonLinearParams();
  	object->GetNonLinearParams(vars, 100);
	for(int j=0; j<m1; j++) fprintf(stderr,"%le ", vars[j]);
	fprintf(stderr,"\n");
	object = col.GetCurrentMagObject(&nonlin_start, &lin_start);
	}

   inversion.DumpObjects(stderr);

   fprintf(stderr,"%s\n", message);

}









