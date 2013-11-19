//
// @(#)dipole1d.cpp  1.1  misha1-29sep108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  29sep108  by  misha
// Version:  29sep108  by  misha
//


#include "dipole1d.h"
#include "math.h"

CMagDipole1D::CMagDipole1D(): CMagObject(0, 3,2, "CMagDipole1D", "dipole1d") { 

  m_iObjectLinearAbilities    = 2; 
  m_iObjectNonLinearAbilities = 2;

}


CMagDipole1D::~CMagDipole1D() { 



}

double CMagDipole1D::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				  double x0,    double y0,    double z0,
				  double A,     double B,     double C,
				  int data_pos) {

  double x = (x_obs - (m_lfPos[0]-x0));
  double z = (m_lfPos[1]-z0);  // here is sign to match gnuplot
  //double z = -(m_lfPos[1]-z0);

  double R = sqrt(x*x+z*z);
  double R5 = R*R*R*R*R;

  return (m_lfJ[0]*x*z + m_lfJ[1]*x*x + m_lfJ[2]*z*z) / R5;
  

}


int CMagDipole1D::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				       double * pos_deriv,
				       double x_obs, double y_obs, double z_obs,
				       double x0,    double y0,    double z0,
				       double A,     double B,     double C,
				       int data_pos, int * n_params, int * start) {

  if(start) *start =  m_iLinearStart;


  double x = (x_obs - (m_lfPos[0]-x0));
  double z = (m_lfPos[1]-z0);
  //double z = -(m_lfPos[1]-z0);

  double R = sqrt(x*x+z*z);
  double R5 = R*R*R*R*R;

  deriv[0] = x*z/R5;
  deriv[1] = x*x/R5;
  deriv[2] = z*z/R5;

  *n_params = 3;

  return 1;

}


int  CMagDipole1D::GetNonLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
										   double * pos_deriv, 
										   double x_obs, double y_obs, double z_obs,
										   double x0,    double y0,    double z0,
										   double A,     double B,     double C,
										   int data_pos, int * n_params, int * start) {

  if(start) *start =  m_iNonLinearStart;	

  double x = (x_obs - (m_lfPos[0]-x0));
  double z = (m_lfPos[1]-z0);
  //double z = -(m_lfPos[1]-z0);

  double R = sqrt(x*x+z*z);
  double R5 = R*R*R*R*R;	
  double R7 = R5*R*R;

  *n_params = 2;

  /* for z = -(m_lfPos[1]-z0);

  deriv[0] = (5.*x*(z*z*m_lfJ[2]+x*x*m_lfJ[1]+x*z*m_lfJ[0]) / R7 ) +
	         (-2.*x*m_lfJ[1]-z*m_lfJ[0])/R5;
	
  deriv[1] = -(2.*z*m_lfJ[2]+x*m_lfJ[0])/R5 +
	         (5.*z*(z*z*m_lfJ[2] + x*x*m_lfJ[1] + x*z*m_lfJ[0])/R7);*/


  deriv[0] = (5.*x*(z*z*m_lfJ[2]+x*x*m_lfJ[1]+x*z*m_lfJ[0]) / R7 ) +
	         (-2.*x*m_lfJ[1]-z*m_lfJ[0])/R5;
	
  deriv[1] = (2.*z*m_lfJ[2]+x*m_lfJ[0])/R5 -
	         (5.*z*(z*z*m_lfJ[2] + x*x*m_lfJ[1] + x*z*m_lfJ[0])/R7);
  return 1;

}


void CMagDipole1D::DumpParameters(FILE * dat) {

  CMagObject::DumpParameters(dat);
  fprintf(dat,"1D dipole position:\n");
  
  fprintf(dat,"T = (A*(x-x0)*z + B*(x-x0)^2 + C*z^2)/R^5)\n");

  if(m_pPosfixed[0]) 
    fprintf(dat," X: %.3lf FIXED\n", m_lfPos[0]);
  else 
    fprintf(dat," X: %.3lf +/- %lg\n", m_lfPos[0], m_lfStdDevPos[0]);

 if(m_pPosfixed[1]) 
    fprintf(dat," Z: %.3lf FIXED\n", m_lfPos[1]);
  else 
    fprintf(dat," Z: %.3lf +/- %lg\n", m_lfPos[1], m_lfStdDevPos[1]);

 
    if(m_pJfixed[0]) 
      fprintf(dat," A: %lg FIXED\n",  m_lfJ[0]);
    else
      fprintf(dat," A: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);

   if(m_pJfixed[0]) 
      fprintf(dat," B: %lg FIXED\n",  m_lfJ[1]);
    else
      fprintf(dat," B: %lg +/- %lg\n",  m_lfJ[1], m_lfStdDevJ[1]);

   if(m_pJfixed[0]) 
      fprintf(dat," C: %lg FIXED\n",  m_lfJ[2]);
    else
      fprintf(dat," C: %lg +/- %lg\n",  m_lfJ[2], m_lfStdDevJ[2]);

}
