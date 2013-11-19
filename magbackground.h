
//
// @(#)magbackground.h  1.1  misha-14may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  14may103  by  misha@misha.local
// Version:  14may103  by  misha@misha.local
//

#ifndef MAGBACKGROUND_H
#define MAGBACKGROUND_H

#include "magobject.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class CMagBackGround : public CMagObject   {

 protected:

  double m_lfXFalse;
  double m_lfYFalse;
  double m_lfZFalse;


public:

    //  1  A*x + B*y + C*z + D;

    CMagBackGround(); // 
virtual    ~CMagBackGround(); 


  virtual double ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C, int data_pos);

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
				      int data_pos, int * n_params, int * start = NULL) { *n_params = 0; return 0;}

  virtual void DumpParameters(FILE * dat);

  virtual int CanComputeField() { return 1; }

  virtual void DumpToXMLStream(std::ostream & str);
};



#endif
