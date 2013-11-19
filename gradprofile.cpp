//
// @(#)gradprofile.cpp  1.1  misha-06oct103
//
// Copyright (c) 2003 Geometrics
//
// Created:  06oct103  by  misha@misha2.geometrics.com
// Version:  06oct103  by  misha@misha2.geometrics.com
//


#include "gradprofile.h"

void CGradProfile::clear_all() {

   CTestProfile::clear_all();

   m_vNx.clear();
   m_vNy.clear();
   m_vNz.clear();
  
   m_iIsFixedGradDirection = 1;

   m_lfNx = 0.;
   m_lfNy = 0.;   
   m_lfNz = 1.;  

   m_lfSeparation = 1.;
}





void CGradProfile::Add(double x, double y, double time, double field, double z, 
	               double nx, double ny, double nz) {

  CTestProfile::Add(x, y, time, field, z);

  if(!m_iIsFixedGradDirection) {
      m_vNx.push_back(nx);
      m_vNy.push_back(ny);
      m_vNz.push_back(nz);
  }

}


int CGradProfile::GetGradDirection(int n, double * dir) {
  
  if(m_iIsFixedGradDirection) {
    dir[0] = m_lfNx;
    dir[1] = m_lfNy;
    dir[2] = m_lfNz;
    return 1;
  }

  if(n >= int(m_vNz.size())) return 0;

  dir[0] = m_vNx[n];
  dir[1] = m_vNy[n];
  dir[2] = m_vNz[n];

  return 1;
}


int CGradProfile::get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				    double *add_parameters, int size_in, int *size_out) {


  int n = m_nPos;
  if(!CTestProfile::get_next_position(x, y, z, T, n_pos, add_parameters, size_in, size_out)) {
    return 0;
  }

   double dir[3];

  if(add_parameters && GetGradDirection(n, dir) ) {
    add_parameters[DIRECTION_X] = dir[0];
    add_parameters[DIRECTION_Y] = dir[1];
    add_parameters[DIRECTION_Z] = dir[2];
    add_parameters[FINITE_SEPARATION] = m_lfSeparation;
  }


   return 1;

}

int CGradProfile::get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
				  double *add_parameters, int size_in, int *size_out) {

  if(!CTestProfile::get_at_position(n_pos, x, y, z, T, add_parameters, size_in, size_out)) return 0;


  double dir[3];

   if(add_parameters &&  GetGradDirection(n_pos, dir)) {
     add_parameters[DIRECTION_X] = dir[0];
     add_parameters[DIRECTION_Y] = dir[1];
     add_parameters[DIRECTION_Z] = dir[2];
     add_parameters[FINITE_SEPARATION] = m_lfSeparation;
   }


   return 1;

}



void CGradProfile::Dump(FILE * dat) {


    switch(m_iFieldType) {
	case  totalmagfield:
	    fprintf(dat,"X   Y    Z   T  Line Tsynt \n"); break;
	default:   
	    fprintf(dat,"X   Y    Z   T  Line Tsynt DIR_X DIR_Y DIR_Z \n"); 
    }

  int n_synt = GetSyntheticPoints();

   double x_off[3];

   x_off[0] = m_lfPos[X_OFFSET];
   x_off[1] = m_lfPos[Y_OFFSET];
   x_off[2] = m_lfPos[Z_OFFSET];

   RecomputeIntoGlobal(x_off);

  for(int i=0; i<GetTotalPoints(); i++) { 
      fprintf(dat,"%lf  %lf  %lf  %le %s",  GetX(i) + x_off[0], 
	      GetY(i) + x_off[1],
	      GetZ(i) + x_off[2],
	      GetField(i), GetName());  
      if(n_synt) fprintf(dat," %le", GetSyntField(i));
            else fprintf(dat," 0");

    switch(m_iFieldType) {
      case  totalmagfield:
	  fprintf(dat,"\n");
	  break;
      default: {	  
         if(m_iIsFixedGradDirection) 
	   fprintf(dat," %lf %lf %lf\n", m_lfNx, m_lfNy, m_lfNz); 
         else
	   fprintf(dat," %lf %lf %lf\n", m_vNx[i], m_vNy[i], m_vNz[i]); 
	  }
    }

    }

}



void CGradProfile::DumpParameters(FILE * dat){

    CTestProfile::DumpParameters(dat);

    if(m_lfSeparation !=0) 
	  fprintf(stderr,"  Finite gradient separation %lf\n", m_lfSeparation);    


}







