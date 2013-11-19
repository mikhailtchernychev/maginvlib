//
// @(#)magdata.cpp  1.1  misha-13may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  13may103  by  misha@misha.local
// Version:  13may103  by  misha@misha.local
//


#include "magdata.h"
#include "magfunc.h"

CMagData::CMagData():CMagObject(0, MAX_FIELD_PARAMS, MAX_FIELD_PARAMS, "data") { 

  m_iFieldType = totalmagfield;

  for(int i=0; i< MAX_FIELD_PARAMS; i++) {
    CMagObject::FixNonLinearParam(i, TRUE, 0.);
    CMagObject::FixLinearParam(i, TRUE, 0.);
  }


  m_Nx[0] = 1.;   m_Nx[1] = 0.;  m_Nx[2] = 0.;  // local X axis direction
  m_Ny[0] = 0.;   m_Ny[1] = 1.;  m_Ny[2] = 0.;  // local Y axis direction
  m_Nz[0] = 0.;   m_Nz[1] = 0.;  m_Nz[2] = 1.;  // local Z axis direction

  m_iIsFixedGradDirection = 0;
  m_lfNx = 0; 
  m_lfNy = 0;   
  m_lfNz = 1;
  m_lfSeparation = 1.;

}

CMagData::CMagData(FIELDTYPE type):CMagObject(0, MAX_FIELD_PARAMS, MAX_FIELD_PARAMS, "data") { 

  m_iFieldType = type;


  for(int i=0; i< MAX_FIELD_PARAMS; i++)
    CMagObject::FixNonLinearParam(i, TRUE, 0.);

  m_Nx[0] = 1.;   m_Nx[1] = 0.;  m_Nx[2] = 0.;  // local X axis direction
  m_Ny[0] = 0.;   m_Ny[1] = 1.;  m_Ny[2] = 0.;  // local Y axis direction
  m_Nz[0] = 0.;   m_Nz[1] = 0.;  m_Nz[2] = 1.;  // local Z axis direction

  m_iIsFixedGradDirection = 0;
  m_lfNx = 0; 
  m_lfNy = 0;   
  m_lfNz = 1;
  m_lfSeparation = 1.;

}

CMagData::~CMagData() { 

}


void CMagData::SetLocalSystem(double * nx, double * ny, double * nz) {

    int i;
    for(i=0; i<3; i++) {
	m_Nx[i] = nx[i];
	m_Ny[i] = ny[i];
	m_Nz[i] = nz[i];
    }

}



int CMagData::GetAllPositionParams(double * params,    int size) {

  if(!m_lfPos && !size) return 0;

  CMagObject::GetAllPositionParams(params, size);

  if(size >= 3) {
      params[0] = m_lfPos[0]*m_Nx[0] +  m_lfPos[1]*m_Ny[0] + m_lfPos[2]*m_Nz[0];
      params[1] = m_lfPos[0]*m_Nx[1] +  m_lfPos[1]*m_Ny[1] + m_lfPos[2]*m_Nz[1];
      params[2] = m_lfPos[0]*m_Nx[2] +  m_lfPos[1]*m_Ny[2] + m_lfPos[2]*m_Nz[2];
  }
  
  if(size>= FINITE_SEPARATION) {
	  params[DIRECTION_X]       = m_lfNx;
	  params[DIRECTION_Y]       = m_lfNy;
	  params[DIRECTION_Z]       = m_lfNz;
	  params[FINITE_SEPARATION] = m_lfSeparation;
  }
 

  return 1;
}

void CMagData::RecomputeIntoGlobal(double *x) {

    double params[3];

   params[0] = x[0]*m_Nx[0] +  x[1]*m_Ny[0] + x[2]*m_Nz[0];
   params[1] = x[0]*m_Nx[1] +  x[1]*m_Ny[1] + x[2]*m_Nz[1];
   params[2] = x[0]*m_Nx[2] +  x[1]*m_Ny[2] + x[2]*m_Nz[2];

   for(int i=0; i<3; i++) x[i] = params[i];

}


int CMagData::GetNonLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				      double * pos_deriv, 
				      double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
				      int data_pos, int * n_params, int * start) {

  int i, j;

  double deriv_local[3];

  deriv_local[0] = pos_deriv[0]*m_Nx[0] + pos_deriv[1]*m_Nx[1] + pos_deriv[2]*m_Nx[2];
  deriv_local[1] = pos_deriv[0]*m_Ny[0] + pos_deriv[1]*m_Ny[1] + pos_deriv[2]*m_Ny[2];
  deriv_local[2] = pos_deriv[0]*m_Nz[0] + pos_deriv[1]*m_Nz[1] + pos_deriv[2]*m_Nz[2];

 if(start) *start =  m_iNonLinearStart;
 
  j=0;
  *n_params = 0;
  for(i=0; i<3; i++) {
	  if(!m_pPosfixed[i]) { 
	    deriv[j] = -deriv_local[i];
	    j++;
	    (*n_params)++;
	  }
	}

  return 1;
}


int CMagData::SetFixedDirection(double nx, double ny, double nz) {

  double s = nx*nx + ny*ny + nz*nz;
  if(s == 0.) return 0;

  s = sqrt(s);

  nx /= s; ny /=s; nz /=s;

  m_lfNx = nx;
  m_lfNy = ny;   
  m_lfNz = nz;  

  m_iIsFixedGradDirection = 1;
  return 1;

}