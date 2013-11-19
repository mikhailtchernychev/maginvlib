//
// @(#)magdipolethread.cpp  1.1  misha1-26jun107
//
// Copyright (c) 2007 of Geometrics.
// All rights reserved.
//
// Created:  26jun107  by  misha1@misha.geometrics.local
// Version:  26jun107  by  misha1@misha.geometrics.local
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "magdipolethread.h"


CMagDipoleThread::CMagDipoleThread(int is_induced, double dipoles_per_length) : CMagObject(is_induced, 3, 6,
							     "dipole_thread", "CMagDipoleThread") {

    m_lfDipolesPerLength = dipoles_per_length;
    m_pDipole =  new CMagDipole(is_induced);
   
}


CMagDipoleThread::CMagDipoleThread(CMagPipe * pipe, double dipoles_per_length) : CMagObject(0, 3, 6,
							     "dipole_thread", "CMagDipoleThread") {

    m_lfDipolesPerLength = dipoles_per_length;
    m_pDipole =  new CMagDipole(0);
    
    *((CMagObject *)this) =  *((CMagObject *)pipe);
    
    double alen, scale;
    int n;
    get_dipoles(&alen, &n); 
    
    scale =  (2.*m_lfPos[AHALF_POS_TH])/n;  
    
    for(int i=0; i<3; i++) {
           m_lfJ[i] *= scale; 
    }
    
}

CMagDipoleThread::~CMagDipoleThread() {
    if(m_pDipole)      delete m_pDipole;
}

void CMagDipoleThread::get_dipoles(double *alen, int *n) {

    double len2  =  2*m_lfPos[AHALF_POS_TH];
    
    *n = int((len2/m_lfDipolesPerLength)+0.5);

    *alen = len2 / double((*n)-1);
    
}


void CMagDipoleThread::LoadDipole(double dist) {

    double params[3];

    double ny = cos(m_lfPos[INCL_POS_TH])*cos(m_lfPos[DECL_POS_TH]);
    double nx = cos(m_lfPos[INCL_POS_TH])*sin(m_lfPos[DECL_POS_TH]);
    double nz = sin(m_lfPos[INCL_POS_TH]);

    double x0 = m_lfPos[X0_POS_TH];
    double y0 = m_lfPos[Y0_POS_TH];
    double z0 = m_lfPos[Z0_POS_TH]; 


    double len  =  m_lfPos[AHALF_POS_TH];

    params[0] = x0  + nx*(dist - len);
    params[1] = y0  + ny*(dist - len);
    params[2] = z0  + nz*(dist - len);

    m_pDipole->SetNonLinearParams(params, 3);
    
    if(m_bInduced) {
    	params[0] = m_lfJ[0];
     	m_pDipole->SetLinearParams(params, 1);
    }
    else {
    	params[0] = m_lfJ[1];
     	params[1] = m_lfJ[0];
      	params[2] = m_lfJ[2];
       	m_pDipole->SetLinearParams(params, 3);
    }
}


double CMagDipoleThread::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
			       double x0,    double y0,    double z0,
			       double A,     double B,     double C,
			       int n_pos) {

    if(!m_pDipole) return 0;

    double T = 0;
    int n;
    double dlen;

    get_dipoles(&dlen, &n);

    for(int i=0; i<n; i++) {
	LoadDipole(dlen*i);
	T += m_pDipole->ComputeField(type, add_params, x_obs, y_obs, z_obs,  x0,  y0, z0, A, B, C, n_pos);
    }
	
    return T;

}



int CMagDipoleThread::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
				    double * pos_deriv,
				    double x_obs, double y_obs, double z_obs,
				    double x0,    double y0,    double z0,
				    double A,     double B,     double C,
				    int n_pos, int *n_params, int * start) {

  if(!m_pDipole) return 0;

  int n_params_tmp = 0;
  int size_tmp, i, j;

  int n;
  double dlen;

  get_dipoles(&dlen, &n);

  if(m_bInduced) size_tmp = 1;
            else size_tmp = 3;

  if(start) *start =  m_iLinearStart;
  *n_params = 0;

  double lfDeriv[3], lfDerivAll[3];

  for(i=0; i<3; i++) lfDerivAll[i] = 0.;

  for(i=0; i<n; i++) {
	LoadDipole(dlen*i);
	m_pDipole->GetLinearDerivatives(type, add_params, lfDeriv, size_tmp, 
					 pos_deriv, x_obs, y_obs, z_obs,
					 x0, y0, z0,  A, B, C, n_pos, &n_params_tmp, NULL);

	for(j=0; j<3; j++) lfDerivAll[j] += lfDeriv[j];
   }

   *n_params  = 0;
    j=0;
    for(i=0; i<3; i++) {
    if(!m_pJfixed[i]) { 
	    deriv[j] = lfDerivAll[i];
	    j++;
	    (*n_params)++;
	    }
    }

    return 1;

}


int CMagDipoleThread::GetTotalLinearParams() {
    if(m_bInduced) return 1;
    else return m_nLinearParams;

}




void CMagDipoleThread:: DumpParameters(FILE * dat) {


  CMagObject::DumpParameters(dat);

  fprintf(dat,"Position:\n"); 
  if(m_pPosfixed[X0_POS_TH]) 
    fprintf(dat," X: %lg FIXED\n", m_lfPos[X0_POS_TH]);
  else 
    fprintf(dat," X: %lg +/- %lg\n", m_lfPos[X0_POS_TH], m_lfStdDevPos[X0_POS_TH]);

 if(m_pPosfixed[Y0_POS_TH]) 
    fprintf(dat," Y: %lg FIXED\n", m_lfPos[Y0_POS_TH]);
  else 
    fprintf(dat," Y: %lg +/- %lg\n", m_lfPos[Y0_POS_TH], m_lfStdDevPos[Y0_POS_TH]);

 if(m_pPosfixed[Z0_POS_TH]) 
    fprintf(dat," Z: %lg FIXED\n", m_lfPos[Z0_POS_TH]);
  else 
    fprintf(dat," Z: %lg +/- %lg\n", m_lfPos[Z0_POS_TH], m_lfStdDevPos[Z0_POS_TH]);


 if(m_pPosfixed[AHALF_POS_TH]) 
    fprintf(dat," Half length: %lg FIXED\n", m_lfPos[AHALF_POS_TH]);
  else 
    fprintf(dat," Half length: %lg +/- %lg\n", m_lfPos[AHALF_POS_TH], m_lfStdDevPos[AHALF_POS_TH]);


  double rad=M_PI/180.;

 if(m_pPosfixed[DECL_POS_TH]) 
    fprintf(dat," Direction declination: %lg FIXED\n", m_lfPos[DECL_POS_TH]/rad);
  else 
    fprintf(dat," Direction declination: %lg +/- %lg\n", m_lfPos[DECL_POS_TH]/rad, m_lfStdDevPos[DECL_POS_TH]/rad);


 if(m_pPosfixed[INCL_POS_TH]) 
    fprintf(dat," Direction inclination: %lg FIXED\n", m_lfPos[INCL_POS_TH]/rad);
  else 
    fprintf(dat," Direction inclination: %lg +/- %lg\n", m_lfPos[INCL_POS_TH]/rad, m_lfStdDevPos[INCL_POS_TH]/rad);



  if(m_bInduced) {
    fprintf(dat,"Mag. Moment - induced only:\n");
    if(m_pJfixed[0]) 
      fprintf(dat," Jtotal: %lg FIXED\n",  m_lfJ[0]);
    else
      fprintf(dat," Jtotal: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);
   }
  else {
    fprintf(dat,"Mag. Moment:\n");

    if(m_pJfixed[0]) 
      fprintf(dat," Jx: %lg FIXED\n",  m_lfJ[0]);
    else
      fprintf(dat," Jx: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);

    if(m_pJfixed[1]) 
      fprintf(dat," Jy: %lg FIXED\n",  m_lfJ[1]);
    else
      fprintf(dat," Jy: %lg +/- %lg\n",  m_lfJ[1], m_lfStdDevJ[1]);

    if(m_pJfixed[2]) 
      fprintf(dat," Jz: %lg FIXED\n",  m_lfJ[2]);
    else
      fprintf(dat," Jz: %lg +/- %lg\n",  m_lfJ[2], m_lfStdDevJ[2]);
  }

}






