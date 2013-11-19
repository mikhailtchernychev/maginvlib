//
// @(#)arrayinv.h  1.1  misha1-07oct108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  07oct108  by  misha
// Version:  07oct108  by  misha
//


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
#include "monopole.h"

#include "magarrayt.h"

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

#include <iostream>
#include <fstream>

template <class T> int ArrayInversion(
				      CMagCollection & ObjectCollection,
					  CMagDataCollection & DataCollection,
				      std::vector<T**> vX,    // one vector element per data array.
				      std::vector<T**> vY,
				      std::vector<T**> vZ,
				      std::vector<T**> vA,
				      std::vector<T**> vT,
				      std::vector<T**> vTsynt,
				      std::vector<double *> vTime,
				      std::vector<int> vDataTypes,
				      std::vector<double> vSeparation,
				      std::vector<std::string> vLineNames,
				      std::vector<int> vCh,
				      std::vector<int> vN,
				      std::vector<int> vN1,
				      std::vector<int> vN2,
				      std::vector<double> & vX0, // initial dipole positions.
				      std::vector<double> & vY0, // vector size == number of dipoles
				      std::vector<double> & vZ0, // on exit - new dipoles positions + std. dev (2*n_dipoles)
   				      std::vector<double> & vJ,  // on exit - new dipoles moments   + std. dev (2*n_dipoles)
				      double * rms_fit,  double * max_fit,  // model fit
				      double * obs_max, double * obs_min,
				      double I, double D, double Az,   // Inclination, declination, azimuth
				      int object_type,
					  double * median_water_depth,
					  double * median_sensor_depth,
					  std::vector<std::string> vIDs,
				      int is_induced = 0,
				      double depth_dm =  -1,
				      int      max_iterations    = 2000,
				      PROGRESS_FUNC pFunc = NULL,
				      char * inv_message  = NULL,
				      int background = 1,
				      FILE * objects_out = stderr,
				      FILE * data_out    = stdout,
				      std::ostream & stream_out = std::cout,
					  const char * protocol_folder =  NULL,
					  double xshift = 0., double yshift = 0.,
					  int is_fixed_depth = 0,
					  int reported_depth = 2) {


   char message[1024];
   double fit;
   int i, j, k;
   double middle_x, middle_y, middle_z;
   int  n_dipoles;
   double mag_min, mag_max, mag_ampl;

   //CMagDataCollection DataCollection;
   //CMagCollection     ObjectCollection;
   CMagInv inversion;
   int ret = 0, nb, jj;
   int  n, m, nlin;
   int n_data_points = 0;

   std::vector<double> vlfMedianWDepth;
   vlfMedianWDepth.clear();

   m = n = nlin = 0;

   inversion.SetEarthFieldDirections(I,D,Az);
   inversion.SetInvID(vIDs[0]);

   Z29RLS_PARAMS rls_params;
   DNLS1_STORAGE dstorage;
   rls_params.dm = 1.e+20;

	int dnls_size =  sizeof(DNLS1_STORAGE);
	int fp_size =  sizeof(S_fp);

   dstorage.iopt   =  2;
   dstorage.maxfev =  max_iterations;

   // check if z is const
   int z_first    =  1;
   int is_z_const =  1;
   double z_const =  0.;

   for(i=0; i<vZ.size() && is_z_const; i++) {
	   for(j=0; j<vCh[i] && is_z_const; j++) {
		   for(k=vN1[i]; k<=vN2[i] && is_z_const; k++) {
			   if(z_first) {
				   z_const = vZ[i][j][k];
				   z_first = 0;
			   }
			   else {
				   if(z_const!=vZ[i][j][k]) is_z_const = 0;
			   }
		   }
	   }
   }


    // data collection
    middle_x = middle_y = middle_z = 0.;

   for(i=0; i<vX.size(); i++) {
     CMagArray<T> * array = new CMagArray<T>;
     if(!array) return 0;
     array->SetObjectName((char*)vLineNames[i].c_str());
     array->SetArrays(vCh[i], vN[i], vX[i], vY[i], vT[i], vTsynt[i], vZ[i], vA[i], vTime[i]);
     if(vN1[i]>= vN2[i]) {
       delete array;
       continue;
     }

     array->SetWindow(vN1[i], vN2[i]);
     n_data_points += (vN2[i] - vN1[i] + 1);

     switch(vDataTypes[i]) {
       case 2: // vert. gradient
	 array->SetFixedDirection(0.,0.,1);
	 array->SetFieldType(totalmaggradient);
	 break;
        case 3: // finite vert. grad
	 array->SetFixedDirection(0.,0.,1);
	 array->SetFiniteSeparation(vSeparation[i]);
         array->SetFieldType(totalfinitemaggrad);
	 break;
     default:
         array->SetFieldType(totalmagfield);
	 }

     DataCollection.AddData(array);
     middle_x += array->GetMiddleX();
     middle_y += array->GetMiddleY();
     middle_z += array->GetMiddleZ();
     if(depth_dm > 0) array->DepthSplineSmooth(depth_dm);
   }

   middle_x /= double(vX.size());
   middle_y /= double(vX.size());
   middle_z /= double(vX.size());


   // compute median sensor depth

   vlfMedianWDepth.clear();
   vlfMedianWDepth.reserve(n_data_points);
   for(i=0; i<vX.size(); i++) {
	   int n1 = vN1[i];
	   int n2 = vN2[i];

	   for(j=0; j<vCh[i]; j++) {
		   for(k=n1; k<=n2; k++) vlfMedianWDepth.push_back(vZ[i][j][k]);
	   }

   }

   if(vlfMedianWDepth.size() > 0) {
	   *median_sensor_depth = median_from_vector<double>(vlfMedianWDepth);
  	   vlfMedianWDepth.clear();
	   inversion.AddDoubleProperty("MEDIAN_SENSOR_DEPTH", *median_sensor_depth);
   }
   else *median_sensor_depth = NO_DATA;

   // compute median altitude and water depth if any

   vlfMedianWDepth.reserve(n_data_points);

   int is_alt   = 0;


   for(i=0; i<vX.size(); i++) {
	   int n1 = vN1[i];
	   int n2 = vN2[i];

   	   for(j=0; j<vCh[i]; j++) {
		   for(k=n1; k<=n2; k++) {
			   if(vA[i][j][n1]!=vA[i][j][k]) {
				   is_alt = 1;
				   break;
			   }
		   }
	   }

	   for(j=0; j<vCh[i]; j++) {
		   for(k=n1; k<=n2 && is_alt; k++) vlfMedianWDepth.push_back(vA[i][j][k]);
	   }

   }

   if(vlfMedianWDepth.size() > 0 && *median_sensor_depth != NO_DATA) {
	   *median_water_depth = median_from_vector<double>(vlfMedianWDepth) + *median_sensor_depth;
  	   vlfMedianWDepth.clear();
	   inversion.AddDoubleProperty("MEDIAN_WATER_DEPTH", *median_water_depth);
   }
   else {
	   if(is_alt==0) *median_water_depth = vA[0][0][vN1[0]] + *median_sensor_depth;
	            else *median_water_depth = NO_DATA;
   }




    // dipoles

   n_dipoles = vX0.size();

    for(i=0; i<vX0.size(); i++) {
		double pos[3];

		if(vX0[i] !=NO_DATA && vY0[i]!=NO_DATA) {
			pos[0] = vX0[i]; pos[1] = vY0[i]; //pos[2] = middle_z + vZ0[i];
		}
		else {
			pos[0] = middle_x + xshift;
			pos[1] = middle_y + yshift;
			//pos[2] = middle_z + vZ0[i];
		}

		switch(reported_depth) {
		case 0:
			if(*median_water_depth == NO_DATA) pos[2]  = vZ0[i];
			else pos[2] = (*median_water_depth) + vZ0[i];
			break;
		case 1:
			pos[2] = vZ0[i];
			break;
		default:
			pos[2]  = *median_sensor_depth + vZ0[i];
			break;
		}

			switch(object_type) {
		case 0: {
		  CMagDipole * pDipole = new CMagDipole(is_induced);
		  if(!pDipole) return 0;
		  pDipole->SetNonLinearParams(pos, 3);
		  switch(is_fixed_depth) {
		  case 1:
			  pos[2] = *median_water_depth;
			  pDipole->FixNonLinearParam(2, 1, pos[2]);
			  break;
		  case 2:
			  pDipole->FixNonLinearParam(2, 1, pos[2]);
			  break;
		  }
		  pDipole->SetID(vIDs[i]);
		  ObjectCollection.AddMagObject(pDipole);
		  break;
		}
		case 1:
		  CMagMonopole * pMonopole =  new CMagMonopole();
		  if(!pMonopole) return 0;
		  pMonopole->SetID(vIDs[i]);
		  pMonopole->SetNonLinearParams(pos, 3);
		  switch(is_fixed_depth) {
		  case 1:
			  pos[2] = *median_water_depth;
			  pMonopole->FixNonLinearParam(2, 1, pos[2]);
			  break;
		  case 2:
			  pMonopole->FixNonLinearParam(2, 1, pos[2]);
			  break;
		  }
		  ObjectCollection.AddMagObject(pMonopole);
		  break;
		}
    }


   CMagObject *  pBackground;

   switch(background) {

     case 0:
         pBackground = new CMagBackGround;
 	 if(pBackground)  {
	     pBackground->FixLinearParam(X_OFFSET, TRUE, 0.);
	     pBackground->FixLinearParam(Y_OFFSET, TRUE, 0.);
   	     pBackground->FixLinearParam(Z_OFFSET, TRUE, 0.);
	 }
	 break;

     case 2:
     pBackground = new CMagBackGround2;
     break;

     default:
       pBackground = new CMagBackGround;
       break;
   }

    if(!pBackground)  return 0;
	if(is_z_const)	pBackground->FixLinearParam(Z_OFFSET, TRUE, 0.);

    ObjectCollection.AddMagObject(pBackground);

    // now check if there are enough data points

    if(ObjectCollection.GetTotalLinearParams() + ObjectCollection.GetTotalNonLinearParams() >= n_data_points) {
      return 0;
    }



     if(pFunc) inversion.SetProgressFunc(pFunc);
     inversion.SetMagData(&DataCollection);
     inversion.SetMagCollection(&ObjectCollection);
     inversion.SetEarthFieldDirections(I,D,Az);

     // linear estination
     inversion.EstimateLinearParameters(Estimate_L2_rls, &fit, 1, &rls_params);

     // non-linear estimation
     inversion.EstimateAllParameters(Estimator_dnls1, &dstorage, message);

     *max_fit = inversion.GetNonLinearMaxDiff();
     *rms_fit = inversion.GetNonLinearDiscrepancy();

     if(inv_message) strcpy(inv_message, message);

     if(!stream_out) ; else inversion.DumpToXMLStream(stream_out);

     // compute model field to estimate anomaly amplitude

     inversion.GetTotalParameters(&n, &m, &nlin);
     if(n>0) {
       double *T_synt = new double[n];
       if(T_synt) {
	 nb = pBackground->GetAllLinearParams();
	 for(jj=0; jj<nb; jj++) {
           pBackground->FixLinearParam(jj, TRUE, 0.);
	 }

	 mag_min =  1.e+38;
	 mag_max = -1.e+38;
	 inversion.ComputeField(T_synt, n);
	 for(jj=0; jj<n; jj++) {
           mag_min = min(mag_min, T_synt[jj]);
           mag_max = max(mag_max, T_synt[jj]);
	 }
	 mag_ampl = fabs(mag_min-mag_max);
	 if(objects_out) fprintf(objects_out,"mag_min = %lf mag_max = %lf\n", mag_min, mag_max);
         delete [] T_synt;
       }
     }


	 // damping by sections
	 DataCollection.ResetDataPos();

	 CMagData * pData = DataCollection.GetNextData();
	 while(pData && data_out) {
		 	int n_channels = pData->GetNumberOfChannels();
			int n_data     = pData->GetNDataPerChannel();
			pData->Dump(data_out);
			pData = DataCollection.GetNextData();
	 }

	 //DataCollection.Dump(data_out);
	 if(objects_out) inversion.DumpObjects(objects_out);

     vX0.clear();
     vY0.clear();
     vZ0.clear();
     vJ.clear();

     // values

     for(i=0; i<n_dipoles; i++) {
       CMagObject * pObject = ObjectCollection.GetObjectAt(i);
       if(pObject) {
	 vX0.push_back(pObject->GetNonLinearParam(0, 0));
	 vY0.push_back(pObject->GetNonLinearParam(1, 0));
	 vZ0.push_back(pObject->GetNonLinearParam(2, 0));
	 for(i=0; i<3; i++) vJ.push_back(pObject->GetLinearParam(i,0));
       }
     }

     // std. deviations
     for(i=0; i<n_dipoles; i++) {
       CMagObject * pObject = ObjectCollection.GetObjectAt(i);
       if(pObject) {
	 vX0.push_back(pObject->GetNonLinearParam(0, 1));
	 vY0.push_back(pObject->GetNonLinearParam(1, 1));
	 vZ0.push_back(pObject->GetNonLinearParam(2, 1));
	 for(i=0; i<3; i++) vJ.push_back(pObject->GetLinearParam(i,1));
       }
     }

	 // observed min / max

	 DataCollection.GetMinMax(0, obs_min, obs_max);

	 // synt min / max

	 //ObjectCollection.RemoveObject(pBackground);
	 //inversion.ComputeField();
	 //DataCollection.Dump(data_out);
	 //DataCollection.GetMinMax(1, synt_min, synt_max);
     //ObjectCollection.AddMagObject(pBackground);


  return 1;

}
