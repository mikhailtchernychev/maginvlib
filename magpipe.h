//
// @(#)magpipe.h  1.1  misha-04nov103
//
// Copyright (c) 2003 Geometrics, Inc.
//
// Created:  04nov103  by  misha@misha2.geometrics.com
// Version:  04nov103  by  misha@misha2.geometrics.com
//


#ifndef MAGPIPE_H
#define MAGPIPE_H

#include "magobject.h"

// for CMagPipe:
// Non linear parameters: x0, y0, z0, a (half a size) 
// declination (angle in XOY plane) inclination (vertical angle)
// angles in radians


#define X0_POS    0
#define Y0_POS    1
#define Z0_POS    2
#define AHALF_POS 3
#define DECL_POS  4
#define INCL_POS  5

class CMagPipe : public CMagObject {

public:
  CMagPipe();
  CMagPipe(int is_induced);


  // x_obs, y_obs, z_obs - position of the observation point
  // A, B, C             - Earth's magnetic field direction coefficients

  virtual double ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
			              int n_pos);


  virtual int GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
				     double * pos_deriv,
				     double x_obs, double y_obs, double z_obs,
				     double x0,    double y0,    double z0,
				     double A,     double B,     double C,
			  	     int n_pos, int *n_params, int * start);

  virtual int GetTotalLinearParams();
  virtual void DumpParameters(FILE * dat);
  virtual int CanComputeField() { return 1; }



 // computational functions

 void GetCoefficients(double * add_params, double &x_obs, double &y_obs, double &z_obs,
		      double x0,    double y0,    double z0,
		      double A,     double B,     double C,
		      double A_tr[9], double n_vect[3], double r[3],
		      double r0[3], double n_earth[3], double J[3],  double * alen);


 double ComputeTotalField(double * add_params, double x_obs, double y_obs, double z_obs,
			  double x0,    double y0,    double z0,
			  double A,     double B,     double C);

 double ComputeGradient(double * add_params, double x_obs, double y_obs, double z_obs,
		        double x0,    double y0,    double z0,
			double A,     double B,     double C);

 double ComputeAbsGradient(double * add_params, double x_obs, double y_obs, double z_obs,
		        double x0,    double y0,    double z0,
			double A,     double B,     double C);

 int GetTotalFieldLinearDerivatives(double * add_params, double * deriv, int size,
				    double * pos_deriv,
				    double x_obs, double y_obs, double z_obs,
	                            double x0,    double y0,    double z0,
				    double A,     double B,     double C,
			            int * n_params, int * start = NULL);


 int GetGradientLinearDerivatives(double * add_params, double * deriv, int size,
				  double * pos_deriv,
				  double x_obs, double y_obs, double z_obs,
	                          double x0,    double y0,    double z0,
				  double A,     double B,     double C,
			          int * n_params, int * start = NULL);


};


#endif
