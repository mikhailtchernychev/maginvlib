//
// @(#)magbackground.cpp  1.1  misha-14may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  14may103  by  misha@misha.local
// Version:  14may103  by  misha@misha.local
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "magbackground.h"


CMagBackGround::CMagBackGround():CMagObject() {

  m_nNonLinearParams = 0;

  strcpy(m_sObjectType, "CMagBackGround");
  m_nLinearParams = 4; 
  strcpy(m_sObjectName, "A*x+B*y+C*z+D");

  m_lfJ = new double[m_nLinearParams];
  if(!m_lfJ) return;

  m_pJfixed = new int[m_nLinearParams];
  if(!m_pJfixed) return;

  m_lfStdDevJ = new double[m_nLinearParams];
  if(!m_lfStdDevJ) return;


  m_lfJUp = new double[m_nLinearParams];
  if(!m_lfJUp) return;

  m_lfJLow = new double[m_nLinearParams];
  if(!m_lfJLow) return;
 
  m_lfJMax   = new double[m_nLinearParams];
  m_lfJMin   = new double[m_nLinearParams];   
  
  for(int i=0; i<m_nLinearParams; i++) {
    m_lfJ[i] = m_lfStdDevJ[i] = 0.;
    m_pJfixed[i] = 0;
    m_lfJUp[i]   =  1.E38;
    m_lfJLow[i]  = -1.E38;
    m_lfJMin[i]  =  1.E38;
    m_lfJMax[i]  = -1.E38;         
  }

  m_lfXFalse = 0.;
  m_lfYFalse = 0.;
  m_lfZFalse = 0.;

  m_bIsValid         = 1;

  m_iObjectLinearAbilities    = 2;
  m_iObjectNonLinearAbilities = 2; 
}



CMagBackGround::~CMagBackGround() {
}


double CMagBackGround::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				            double x0,    double y0,    double z0,
					    double A,     double B,     double C, int data_pos) {
  if(!m_bIsValid) return 0.;
 
 if(m_iFieldType != anyfield && m_iFieldType != type) return 0.;

  m_lfXFalse = -x0;
  m_lfYFalse = -y0;
  m_lfZFalse = -z0;

  if(m_iFieldType != anyfield && m_iFieldType != type) { // if it is nt requested for this field type, return 0.
    return 0.;
  }

  return m_lfJ[0]*x_obs + m_lfJ[1]*y_obs + m_lfJ[2]*z_obs +  m_lfJ[3];

}


void CMagBackGround::DumpParameters(FILE * dat) {

  CMagObject::DumpParameters(dat);

  fprintf(stderr,"Type: A*(x+xfalse) + B*(y+yfalse) + C*(z+zfalse) + D\n"); 

  switch(m_iFieldType) {
      case totalmagfield:   fprintf(dat,"Compute for total magnetic field only\n"); break; 
      case totalmaggradient:   fprintf(dat,"Compute for magnetic field gradient only\n"); break; 
  default:
    fprintf(dat,"Compute for all fields\n");
  }


  fprintf(dat,"xfalse = %lg\nyfalse = %lg\nzfalse = %lg\n", m_lfXFalse, m_lfYFalse, m_lfZFalse);

  if(m_pJfixed[0])
	 fprintf(dat," A=%lg FIXED\n", m_lfJ[0]);
  else 
	 fprintf(dat," A=%lg +/- %lg\n", m_lfJ[0], m_lfStdDevJ[0]);

  if(m_pJfixed[1])
	 fprintf(dat," B=%lg FIXED\n", m_lfJ[1]);
  else 
	 fprintf(dat," B=%lg +/- %lg\n", m_lfJ[1], m_lfStdDevJ[1]);

  if(m_pJfixed[2])
	 fprintf(dat," C=%lg FIXED\n", m_lfJ[2]);
  else 
	 fprintf(dat," C=%lg +/- %lg\n", m_lfJ[2], m_lfStdDevJ[2]);

  if(m_pJfixed[3])
	 fprintf(dat," D=%lg FIXED\n", m_lfJ[3]);
  else 
	 fprintf(dat," D=%lg +/- %lg\n", m_lfJ[3], m_lfStdDevJ[3]);

}



int CMagBackGround::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
					 double * pos_deriv, 
					 double x_obs, double y_obs, double z_obs,
					 double x0,    double y0,    double z0,
					 double A,     double B,     double C, 
					 int data_pos, int * n_params, int * start) {

  if(m_iFieldType != anyfield && m_iFieldType != type) { // if it is nt requested for this field type, return 0.
    deriv[0] =  deriv[1] = deriv[2] = deriv[3] = 0.;
    return 1;
  }

  pos_deriv[0] = pos_deriv[1] = pos_deriv[2] = 0.; 

  if(!m_bIsValid) return 0;

  if(start) *start =  m_iLinearStart;

  m_lfXFalse = -x0;
  m_lfYFalse = -y0;
  m_lfZFalse = -z0;

  int j = 0;
  if(!m_pJfixed[0]) { deriv[j] = x_obs; j++; }
  if(!m_pJfixed[1]) { deriv[j] = y_obs; j++; }
  if(!m_pJfixed[2]) { deriv[j] = z_obs; j++; }
  if(!m_pJfixed[3]) { deriv[j] = 1; j++; }

  *n_params = j;

  return 1;
}





void CMagBackGround::DumpToXMLStream(std::ostream & str) {

   char delimit = ' ';
   char quot = '\"';
   str << "<" << m_sObjectType << ">" << std::endl;

   str << "<formula formula=" << quot << "A*(x+xfalse) + B*(y+yfalse) + C*(z+zfalse) + D" << quot << "/>" << std::endl;

   str << "<params ";

   str << "xfalse=" << quot << m_lfXFalse << quot << delimit;
   str << "yfalse=" << quot << m_lfYFalse << quot << delimit;
   str << "zfalse=" << quot << m_lfZFalse << quot << delimit;
   
   str << "A=" << quot << m_lfJ[0] << quot << delimit;
   str << "B=" << quot << m_lfJ[1] << quot << delimit;
   str << "C=" << quot << m_lfJ[2] << quot << delimit; 
   str << "D=" << quot << m_lfJ[3] << quot << delimit;

   str << "A_std=" << quot << m_lfStdDevJ[0] << quot << delimit;
   str << "B_std=" << quot << m_lfStdDevJ[1] << quot << delimit;
   str << "C_std=" << quot << m_lfStdDevJ[2] << quot << delimit; 
   str << "D_std=" << quot << m_lfStdDevJ[3] << quot << delimit;

   str << "/>" << std::endl;
   str << "</" << m_sObjectType << ">" << std::endl;

}






