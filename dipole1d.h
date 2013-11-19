//
// @(#)dipole1d.h  1.1  misha1-29sep108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  29sep108  by  misha
// Version:  29sep108  by  misha
//


#include "magobject.h"

#ifndef MAGDIPOLE1D_H
#define MAGDIPOLE1D_H


class CMagDipole1D : public CMagObject {
	
	
public:
	CMagDipole1D();
	virtual  ~CMagDipole1D();
	
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
	
	
	
	virtual int GetTotalLinearParams() { return 3; }
	
	
	virtual void DumpParameters(FILE * dat);
	
	
	virtual int CanComputeField() { return 1; }
	
	
};


#endif


