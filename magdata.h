//
// @(#)magdata.h  1.1  misha-13may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  13may103  by  misha@misha.local
// Version:  13may103  by  misha@misha.local
//

// base class for different types of magnetic data

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifndef MAGDATA_H
#define MAGDATA_H

#include "magobject.h"
#include "fieldtype.h"

class CMagData: public CMagObject  {

 protected:
 
  double m_Nx[3]; // local X axis direction
  double m_Ny[3]; // local Y axis direction
  double m_Nz[3]; // local Z axis direction

  int m_iIsFixedGradDirection; // use fixed direction (empty m_vNx ...)

  double m_lfNx;   // fixed direction
  double m_lfNy;   
  double m_lfNz;

  double m_lfSeparation; // separation for finite-base gradient

public:
    
    CMagData();
    CMagData(FIELDTYPE type);
    virtual ~CMagData();

    virtual FIELDTYPE GetFieldType() { return m_iFieldType;}
    virtual void SetFieldType(FIELDTYPE type) { m_iFieldType = type;} 

    virtual int GetNonLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				      double * pos_deriv, 
				      double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
				      int data_pos, int * n_params, int * start = NULL);

    virtual double ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
			              int data_pos) { return 0; }

    virtual int get_data_origin(double *x, double *y, double *z) { return 0;}
    virtual int reset(){return 0;};                         // resets data pointer to start 
    virtual int get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				  double *add_parameters = NULL, int size_in=0, int *size_out=NULL) 
	{ return 0;}

    virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
				  double *add_parameters = NULL, int size_in=0, int *size_out=NULL) 
	{ return 0;}

    virtual int get_total_data_points() { return 0;}
    
    virtual int set_synthetic_data(double *data, int replace=0) { return 0;}
    virtual int is_valid() {return 0;}
    virtual int get_position(int n, double *xd, double *yd, double *zd) { return 0;} 
    virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2) { return 0; }  

    virtual void Dump(FILE * dat) {};
    virtual int GetMinMax(int type, double * data_min, double * data_max) { return 0; }
    virtual int GetAllPositionParams(double * params, int size);

	virtual int GetSyntheticField(double *Tsynt, int n) { return 0; }

    // local coordinate system

    void SetLocalSystem(double * nx, double * ny, double * nz); 
    void RecomputeIntoGlobal(double *x);
    virtual int ComputeLocalCoordinates() {return 0; }

	// for multi-channel data

	virtual int GetNumberOfChannels() { return 0;}
	virtual int GetNDataPerChannel()  { return 0;}

	// nearest point in the data

	virtual int GetNearestDataPoint(double x0, double y0, double * d, double *x, double *y,
		double * dtime, std::string & line_name) { return -1; }

	// fixed and finite gradient
	virtual  int SetVariableDirection() {m_iIsFixedGradDirection = 0; return 1; }
    virtual  int SetFixedDirection(double nx, double ny, double nz);
    int IsFixedDirection() { return m_iIsFixedGradDirection; }
    double GetFiniteSeparation() { return m_lfSeparation; };
    void   SetFiniteSeparation(double s) { m_lfSeparation = s; } 

};



#endif











