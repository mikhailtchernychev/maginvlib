//
// @(#)dipole1dinv.h  1.1  misha-12oct108
//
// Copyright (c) 2008 of Geometrics.
// All rights reserved.
//
// Created:  12oct108  by  misha
// Version:  12oct108  by  misha
//

#ifndef DIPOLE1DINV_H
#define DIPOLE1DINV_H


#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "magobject.h"
#include "magcollection.h"
#include "maggriddata.h"
#include "maginv.h"
#include "magfunc.h"
#include "magbackground.h"
#include "magbackground2.h"
#include "magdipole.h"


#include "l2estimator.h"
#include "l1estimator.h"
#include "dnls1estimator.h"

#include <vector>

#ifdef NO_DATA
#undef NO_DATA
#endif

#ifdef wx_x
#define NO_DATA 1.e35
#else
#define NO_DATA 1.70141e38
#endif




template <class T> int EstimateDipole1D(T *x, T *y, T *z, T * alt,
				       T *Td, T *Tsynt, int n,
				       double x0, double z0,
				       double * rms_fit,  double * max_fit,
				       double * synt_max, double * synt_min,
				       double * position_estimate = NULL,
				       double * linear_estimate   = NULL,
				       double * trend_estimate    = NULL,
				       int      n_extend          = 0,
				       int      max_iterations    = 2000,
				       PROGRESS_FUNC pFunc = NULL,
				       char * inv_message  = NULL) {


    char message[1024];
    double fit;

    CMagDataCollection DataCollection;
    CMagCollection     ObjectCollection;
    CMagInv inversion;

    Z29RLS_PARAMS rls_params;
    DNLS1_STORAGE dstorage;
    rls_params.dm = 1.e+20;

    dstorage.iopt   =  1;
    dstorage.maxfev =  max_iterations;

   // data
   CMagData1Dt<T> *pData =  new CMagData1Dt<T>;
   if(!pData) {
	return 0;
   }

   pData->SetArrays(x, y, z, alt, Td, n);
   pData->SetWindow(0, n-1);

   // dipole
   CMagDipole1D *pDipole = new CMagDipole1D;
   if(!pDipole) {
       delete pData;
       return 0;
   }

   pDipole->FixNonLinearParam(X_OFFSET, FALSE, x0);
   pDipole->FixNonLinearParam(Y_OFFSET, FALSE, z0);

   // background
   CMagBackGround * pBackground = new CMagBackGround;
   if(!pBackground)  {
       delete pData;
       delete pDipole;
       return 0;
   }

   pBackground->FixLinearParam(Y_OFFSET, TRUE, 0.);
   pBackground->FixLinearParam(Z_OFFSET, TRUE, 0.);

   DataCollection.AddData(pData);
   ObjectCollection.AddMagObject(pDipole);
   ObjectCollection.AddMagObject(pBackground);

   if(pFunc) inversion.SetProgressFunc(pFunc);
   inversion.SetMagData(&DataCollection);
   inversion.SetMagCollection(&ObjectCollection);


   // linear estination
   inversion.EstimateLinearParameters(Estimate_L2_rls, &fit, 1, &rls_params);

   // non-linear estimation
   inversion.EstimateAllParameters(Estimator_dnls1, &dstorage, message);

   *max_fit = inversion.GetNonLinearMaxDiff();
   *rms_fit = inversion.GetNonLinearDiscrepancy();

   if(inv_message) strcpy(inv_message, message);

   DataCollection.Dump(stdout);
   inversion.DumpObjects(stderr);

   // synt min / max

   ObjectCollection.RemoveObject(pBackground);
   inversion.ComputeField();
   pData->GetMinMax(1, synt_min, synt_max);
   ObjectCollection.AddMagObject(pBackground);


   // position estimate x0, z0 x_dev, z_dev
   if(position_estimate) {
       position_estimate[0] = pDipole->GetNonLinearParam(0, 0);
       position_estimate[2] = pDipole->GetNonLinearParam(1, 0);
       position_estimate[3] = position_estimate[4] = pDipole->GetNonLinearParam(0, 1);
       position_estimate[5] = pDipole->GetNonLinearParam(1, 1);

       pData->From1Dto2D(position_estimate[0], &position_estimate[0], &position_estimate[1]);
   }


   if(linear_estimate) {
       pDipole->GetLinearParams(linear_estimate,     3, 0);
       pDipole->GetLinearParams(&linear_estimate[3], 3, 1);
   }

   if(trend_estimate) {
       pBackground->GetLinearParams(trend_estimate,     2, 0);
       pBackground->GetLinearParams(&trend_estimate[2], 2, 1);
   }


   //DataCollection.ReplaceData(pData, pData2);
   //delete pData;
   //inversion.ComputeField();

   //pData2->GetSyntheticField(*Tsynt, n+2*n_extend);

   return 1;
}



#endif
