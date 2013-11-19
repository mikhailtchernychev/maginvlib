//
// @(#)dipole.h  1.1  misha-10may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  10may103  by  misha@misha.local
// Version:  10may103  by  misha@misha.local
//

#ifndef MAGDIPOLE_H
#define MAGDIPOLE_H

#include "magobject.h"

class CMagDipole : public CMagObject {


public:
  CMagDipole();
  virtual  ~CMagDipole();
  CMagDipole(int is_induced);

  void GetDipoleMoment(double A, double B, double C, double Az,
                       double * Jtotal, double * incl, double  *decl,
                       double * angle);

  // how to compute magnetic field and its derivatives

  // x_obs, y_obs, z_obs - position of the observation point
  // A, B, C             - Earth's magnetic field direction coefficients

  virtual double ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
			              int data_pos);

  virtual int GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				   double * pos_deriv,
				   double x_obs, double y_obs, double z_obs,
	                           double x0,    double y0,    double z0,
				   double A,     double B,     double C,
			           int data_pos, int * n_params, int * start = NULL);

  virtual int GetNonLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				      double * pos_deriv, 
				      double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
				      int data_pos, int * n_params, int * start = NULL);



  virtual int GetTotalLinearParams();


  virtual void DumpParameters(FILE * dat);


  virtual int CanComputeField() { return 1; }


  // computational functions

 double ComputeTotalField(double * add_params, double x_obs, double y_obs, double z_obs,
			  double x0,    double y0,    double z0,
			  double A,     double B,     double C);

 double ComputeGradient(double * add_params, double x_obs, double y_obs, double z_obs,
		        double x0,    double y0,    double z0,
			double A,     double B,     double C);

 double ComputeAbsGradient(double * add_params, double x_obs, double y_obs, double z_obs,
		        double x0,    double y0,    double z0,
			double A,     double B,     double C);



 double ComputeFiniteGradient(double * add_params, double x_obs, double y_obs, double z_obs,
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

int GetTotalFieldNonLinearDerivatives(double * add_params, double * deriv, int size,
				    double * pos_deriv,
				    double x_obs, double y_obs, double z_obs,
	                            double x0,    double y0,    double z0,
				    double A,     double B,     double C,
			            int * n_params, int * start = NULL);

int GetFiniteGradientLinearDerivatives(double * add_params, double * deriv, int size,
				       double * pos_deriv,
				       double x_obs, double y_obs, double z_obs,
				       double x0,    double y0,    double z0,
				       double A,     double B,     double C,
				       int * n_params, int * start = NULL);


int GetAbsGradientLinearDerivatives(double * add_params, double * deriv, int size, 
				    double * pos_deriv,
				    double x_obs, double y_obs, double z_obs,
				    double x0,    double y0,    double z0,
				    double A,     double B,     double C,
				    int *n_params, int * start);

int GetFiniteGradientNonLinearDerivatives(double * add_params, double * deriv, int size,
				       double * pos_deriv,
				       double x_obs, double y_obs, double z_obs,
				       double x0,    double y0,    double z0,
				       double A,     double B,     double C,
				       int * n_params, int * start = NULL);


 int GetGradientNonLinearDerivatives(double * add_params, double * deriv, int size,
				  double * pos_deriv,
				  double x_obs, double y_obs, double z_obs,
	                          double x0,    double y0,    double z0,
				  double A,     double B,     double C,
			          int * n_params, int * start = NULL);


int GetAbsGradientNonLinearDerivatives(double * add_params, double * deriv, int size, 
				       double * pos_deriv,
				       double x_obs, double y_obs, double z_obs,
				       double x0,    double y0,    double z0,
				       double A,     double B,     double C,
				       int *n_params, int * start = NULL);

void GetAbsFinalGradSecondPoint(int axis, double s, double x1, double y1, double z1,
                                double * x2, double *y2, double *z2); 

double ComputeAbsFiniteGradient(double * add_params, double x_obs, double y_obs, double z_obs,
			      double x0,    double y0,    double z0,
			      double A,     double B,     double C);
			      
int GetAbsFiniteGradientLinearDerivatives(double * add_params, double * deriv, int size, 
			       double * pos_deriv,
			       double x_obs, double y_obs, double z_obs,
			       double x0,    double y0,    double z0,
			       double A,     double B,     double C,
			       int *n_params, int * start = NULL);		      

int GetAbsFiniteGradientNonLinearDerivatives(double * add_params, double * deriv, int size, 
			       double * pos_deriv,
			       double x_obs, double y_obs, double z_obs,
			       double x0,    double y0,    double z0,
			       double A,     double B,     double C,
			       int *n_params, int * start = NULL);	

 int GetInfoHeader(char * str, int size);
 int GetDataHeader(char * str, int size);


 void DumpToStream(std::ostream & str, std::string comment="#", 
						   char * close_str="/", char * delimit=",",
						   std::string offset = "  ");

 virtual void DumpToXMLStream(std::ostream & str);

};






#endif

