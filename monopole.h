//
// @(#)monopole.h  1.1  misha1-09oct108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  09oct108  by  misha
// Version:  09oct108  by  misha
//


#ifndef MONOPOLE_H
#define MONOPILE_H


#include "magobject.h"

class CMagMonopole : public CMagObject {
	
	
public:
	CMagMonopole();
	virtual  ~CMagMonopole();
	
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
	
		
	
	virtual void DumpParameters(FILE * dat);
	
	
	virtual int CanComputeField() { return 1; }
	
	
	
	// computational functions
	
	
	double ComputeTotalField(double * add_params, double x_obs, double y_obs, double z_obs,
		double x0,    double y0,    double z0,
		double A,     double B,     double C);
	
	
	int GetTotalFieldLinearDerivatives(double * add_params, double * deriv, int size,
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
		int * n_params, int * start);

	void DumpToXMLStream(std::ostream & str);
	
};









#endif


