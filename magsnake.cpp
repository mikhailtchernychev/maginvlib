//
// @(#)magsnake.cpp  1.1  misha-13feb10
//
// Copyright (c) 2006 Geometrics, Inc.
//
// Created:  13feb106  by  misha@misha2.geometrics.com
// Version:  13feb106  by  misha@misha2.geometrics.com
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "magsnake.h"


CMagSnake::CMagSnake(int is_induced, int n_seg) : CMagObject(is_induced, 3*n_seg, 3*n_seg+3,
							     "segmented_pipe", "CMagSnake") {

    m_nSegments = n_seg;
    m_pSegment =  new CMagPipe(is_induced);

    m_x0 = m_y0 = m_z0 = 0.;
    m_nCurrentSegment = 0;

    m_pLinearDeriv = new double[3*n_seg];

}

CMagSnake::~CMagSnake() {
    if(m_pSegment)     delete m_pSegment;
    if(m_pLinearDeriv) delete [] m_pLinearDeriv;
}


void CMagSnake::LoadSegment(int n) {

    if(n>=m_nSegments) return;

    double n_vect[3], x1, y1, z1, params[6];

    n_vect[1] = cos(m_lfPos[INCL_N_POS(n)])*cos(m_lfPos[DECL_N_POS(n)]);
    n_vect[0] = cos(m_lfPos[INCL_N_POS(n)])*sin(m_lfPos[DECL_N_POS(n)]);
    n_vect[2] = sin(m_lfPos[INCL_N_POS(n)]);

    x1 = m_x0 + n_vect[0]*m_lfPos[LEN_N_POS(n)];
    y1 = m_y0 + n_vect[1]*m_lfPos[LEN_N_POS(n)];
    z1 = m_z0 + n_vect[2]*m_lfPos[LEN_N_POS(n)];

    params[0] = (x1+m_x0) / 2.;
    params[1] = (y1+m_y0) / 2.;
    params[2] = (z1+m_z0) / 2.;

    params[3] = m_lfPos[LEN_N_POS(n)] / 2.;
    params[4] = m_lfPos[DECL_N_POS(n)];
    params[5] = m_lfPos[INCL_N_POS(n)];

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

    m_x0 = x1;
    m_y0 = y1;
    m_z0 = z1;

}


double CMagSnake::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
			       double x0,    double y0,    double z0,
			       double A,     double B,     double C,
			       int n_pos) {

    if(!m_pSegment) return 0;

    m_x0 = m_lfPos[0];
    m_y0 = m_lfPos[1];
    m_z0 = m_lfPos[2];

    double T = 0;

    for(int i=0; i<m_nSegments; i++) {

	LoadSegment(i);
	T += m_pSegment->ComputeField(type, add_params, x_obs, y_obs, z_obs,  x0,  y0, z0, A, B, C, n_pos);
    }
	
    return T;

}



int CMagSnake::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
				    double * pos_deriv,
				    double x_obs, double y_obs, double z_obs,
				    double x0,    double y0,    double z0,
				    double A,     double B,     double C,
				    int n_pos, int *n_params, int * start) {

  if(!m_pSegment || !m_pLinearDeriv) return 0;

    m_x0 = m_lfPos[0];
    m_y0 = m_lfPos[1];
    m_z0 = m_lfPos[2];

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


int CMagSnake::GetTotalLinearParams() {
    if(m_bInduced) return 1;
    else return m_nLinearParams;

}



int CMagSnake::GetNodes(double *pos, int size) {

    int j=0;
    if(size < 3*(m_nSegments+1)) return -1;

    pos[j++] = m_x0 = m_lfPos[0];
    pos[j++] = m_y0 = m_lfPos[1];
    pos[j++] = m_z0 = m_lfPos[2];

    for(int i=0; i<m_nSegments; i++) { 
	LoadSegment(i);
	pos[j++] = m_x0;
	pos[j++] = m_y0;
	pos[j++] = m_z0;
    }
    return 1;
}



void CMagSnake:: DumpParameters(FILE * dat) {


   m_x0 = m_lfPos[0];
   m_y0 = m_lfPos[1];
   m_z0 = m_lfPos[2];

  CMagObject::DumpParameters(dat);

  fprintf(dat,"Segmented pipe object:\n"); 
  //for(int i=0; i<m_nNonLinearParams; i++) fprintf(dat,"#%d %lg\n", i, m_lfPos[i]);


  if(!m_bInduced)  fprintf(dat,"#   X  Y  Z   Jx  Jy  Jz\n");
  else fprintf(dat,"#   X  Y  Z   Jtotal\n");

   fprintf(dat, "%lg %lg %lg\n", m_x0, m_y0, m_z0); 

  for(int i=0; i<m_nSegments; i++) {

      fprintf(dat, "%lg %lg %lg ", m_x0, m_y0, m_z0);
      LoadSegment(i);

      if(m_bInduced) fprintf(dat,"%lg\n",m_lfJ[i]);
                else fprintf(dat,"%lg %lg %lg\n",m_lfJ[3*i], m_lfJ[3*i+1], m_lfJ[3*i+2]);

      m_pSegment->DumpParameters(dat);

      }

}





