//
// @(#)magdipolethread.h  1.1  misha1-26jun107
//
// Copyright (c) 2007 of Geometrics.
// All rights reserved.
//
// Created:  26jun107  by  misha1@misha.geometrics.local
// Version:  26jun107  by  misha1@misha.geometrics.local
//


#ifndef MAGDIPOLETHREAD_H
#define MAGDIPOLETHREAD_H

#include "magobject.h"
#include "magdipole.h"
#include "magpipe.h"

// for CMagDiploeThreade:
// Non linear parameters: x0, y0, z0, a (half a size) 
// declination (angle in XOY plane) inclination (vertical angle)
// angles in radians
// Number of dipoles per unit length


#define X0_POS_TH    0
#define Y0_POS_TH    1
#define Z0_POS_TH    2
#define AHALF_POS_TH 3
#define DECL_POS_TH  4
#define INCL_POS_TH  5

class CMagDipoleThread : public CMagObject {

protected:
    CMagDipole * m_pDipole;    // use to compute field
    double m_lfDipolesPerLength;
    void LoadDipole(double dist);  

    void get_dipoles(double * dlen, int * n);

public:

   CMagDipoleThread(int is_induced=0, double diploles_per_unit=1);
   CMagDipoleThread(CMagPipe *pipe, double diploles_per_unit=1);  
   ~CMagDipoleThread();


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



};


#endif
