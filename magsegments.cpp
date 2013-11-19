//
// @(#)magsegments.cpp  1.1  misha-15feb10
//
// Copyright (c) 2006 Geometrics, Inc.
//
// Created:  15feb106  by  misha@misha2.geometrics.com
// Version:  15feb106  by  misha@misha2.geometrics.com
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "magsegments.h"


CMagSegments::CMagSegments(int is_induced, int n_seg) : CMagObject(is_induced, 3*n_seg, 3*n_seg+3,
							     "seg_pipe", "CMagSegments") {

    m_nSegments = n_seg;
    m_pSegment =  new CMagPipe(is_induced);

    m_pLinearDeriv = new double[3*n_seg];

}

CMagSegments::~CMagSegments() {
    if(m_pSegment)     delete m_pSegment;
    if(m_pLinearDeriv) delete [] m_pLinearDeriv;
}


void CMagSegments::LoadSegment(int n) {

    double params[6];

    if(n>=m_nSegments) return;

    params[0] = (m_lfPos[X_SEG_POS(n)] + m_lfPos[X_SEG_POS(n+1)])/2.;
    params[1] = (m_lfPos[Y_SEG_POS(n)] + m_lfPos[Y_SEG_POS(n+1)])/2.;
    params[2] = (m_lfPos[Z_SEG_POS(n)] + m_lfPos[Z_SEG_POS(n+1)])/2.;

    params[5] = 0;            // inclination

    double len = sqrt(pow(m_lfPos[X_SEG_POS(n)] - m_lfPos[X_SEG_POS(n+1)],2) + 
                      pow(m_lfPos[Y_SEG_POS(n)] - m_lfPos[Y_SEG_POS(n+1)],2) + 
                      pow(m_lfPos[Z_SEG_POS(n)] - m_lfPos[Z_SEG_POS(n+1)],2));

    double lenxy = sqrt(pow(m_lfPos[X_SEG_POS(n)] - m_lfPos[X_SEG_POS(n+1)],2) + 
		       pow(m_lfPos[Y_SEG_POS(n)] - m_lfPos[Y_SEG_POS(n+1)],2));

    params[3] = len / 2.;     // half length

    // azimuth (declination);
    params[4] = atan2(m_lfPos[X_SEG_POS(n+1)] - m_lfPos[X_SEG_POS(n)], 
		      m_lfPos[Y_SEG_POS(n+1)] - m_lfPos[Y_SEG_POS(n)]);

    // inclination
    params[5] = atan2(m_lfPos[Z_SEG_POS(n+1)] - m_lfPos[Z_SEG_POS(n)], lenxy); 


    m_pSegment->SetNonLinearParams(params, 6);

    if(m_bInduced) {
	params[0] = m_lfJ[n];
	m_pSegment->SetLinearParams(params, 1);
    }
    else {
	params[0] = m_lfJ[n*3];
	params[1] = m_lfJ[n*3+1];
	params[2] = m_lfJ[n*3+2];
	m_pSegment->SetLinearParams(params, 3);
    }

}


double CMagSegments::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
			       double x0,    double y0,    double z0,
			       double A,     double B,     double C,
			       int n_pos) {

    if(!m_pSegment) return 0;

    double T = 0;

    for(int i=0; i<m_nSegments; i++) {

	LoadSegment(i);
	T += m_pSegment->ComputeField(type, add_params, x_obs, y_obs, z_obs,  x0,  y0, z0, A, B, C, n_pos);
    }
	
    return T;

}



int CMagSegments::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
				    double * pos_deriv,
				    double x_obs, double y_obs, double z_obs,
				    double x0,    double y0,    double z0,
				    double A,     double B,     double C,
				    int n_pos, int *n_params, int * start) {

  if(!m_pSegment || !m_pLinearDeriv) return 0;

    int n_params_tmp = 0;
    int size_tmp, i, j;

    if(m_bInduced) size_tmp = 1;
              else size_tmp = 3;

    if(start) *start =  m_iLinearStart;
    *n_params = 0;

    for(i=0; i<size_tmp*m_nSegments; i++) m_pLinearDeriv[i] = 0.;

    for(i=0; i<m_nSegments; i++) {

	LoadSegment(i);
	m_pSegment->GetLinearDerivatives(type, add_params, &m_pLinearDeriv[size_tmp*i], size_tmp, 
					 pos_deriv, x_obs, y_obs, z_obs,
					 x0, y0, z0,  A, B, C, n_pos, &n_params_tmp, NULL); 
    }

     *n_params  = 0;
     j=0;
    for(i=0; i<size_tmp*m_nSegments; i++) {
    if(!m_pJfixed[i]) { 
	    deriv[j] = m_pLinearDeriv[i];
	    j++;
	    (*n_params)++;
	    }
    }

    return 1;

}


int CMagSegments::GetTotalLinearParams() {
    if(m_bInduced) return 1;
    else return m_nLinearParams;

}



int CMagSegments::GetNodes(double *pos, int size) {

    if(size < (m_nSegments+1)*3) return -1;
    for(int i=0; i<(m_nSegments+1)*3; i++) pos[i] = m_lfPos[i];

    return 1;
}



void CMagSegments:: DumpParameters(FILE * dat) {


  CMagObject::DumpParameters(dat);

  fprintf(dat,"Segmented pipe object:\n"); 

  if(!m_bInduced)  fprintf(dat,"#   X  Y  Z   Jx  Jy  Jz\n");
  else fprintf(dat,"#   X  Y  Z   Jtotal\n");

  for(int i=0; i<m_nSegments+1; i++) {

      fprintf(dat, "%lg %lg %lg ", m_lfPos[X_SEG_POS(i)],m_lfPos[Y_SEG_POS(i)], m_lfPos[Z_SEG_POS(i)]);

      if(m_bInduced) fprintf(dat,"%lg\n",m_lfJ[i]);
                else fprintf(dat,"%lg %lg %lg\n",m_lfJ[3*i], m_lfJ[3*i+1], m_lfJ[3*i+2]);

      if(i<m_nSegments) {
	  LoadSegment(i);
	  m_pSegment->DumpParameters(dat);
	  }

  }

}





