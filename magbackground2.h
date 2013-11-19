//
// @(#)magbackground2.h  1.1  misha-15nov104
//
// Copyright (c) 2004 of Geometrics.
// All rights reserved.
//
// Created:  15nov104  by  misha@MISHA
// Version:  15nov104  by  misha@MISHA
//


#ifndef MAGBACKGROUND2_H
#define MAGBACKGROUND2_H

#include "magobject.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class CMagBackGround2 : public CMagObject   {

 protected:

  double m_lfXFalse;
  double m_lfYFalse;
  double m_lfZFalse;


public:

    //  1  A*x + B*y + C*z + D +
    //     A2*x**2 + B2*y**2 + AB*x*y

    CMagBackGround2(); // 
    ~CMagBackGround2();


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
