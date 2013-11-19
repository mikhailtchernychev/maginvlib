//
// @(#)magchannel.cpp  1.1  misha-15nov104
//
// Copyright (c) 2004 of Geometrics.
// All rights reserved.
//
// Created:  15nov104  by  misha@MISHA
// Version:  15nov104  by  misha@MISHA
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "magchannel.h"



CMagChannel::CMagChannel() {
    m_nChannel       = 0;
    m_pMagArray      = NULL;
    m_nPos           = 0;
    m_nTotalChannels = 0;
  }


void CMagChannel::SetChannel(int channel, CMagProfileArray * array) {

  if(channel >= array->GetChannels()) return;

  m_pMagArray      = array;
  m_nChannel       = channel;
  m_nTotalChannels = array->GetChannels();

  char msg[100];
  sprintf(msg,"channel#%d", channel);
  SetObjectName(msg);

}


int CMagChannel::GetTotalPoints() {

  if(!m_pMagArray) return 0;
  return m_pMagArray->GetTotalPoints()/m_nTotalChannels;

}


double CMagChannel::GetX(int n) {
  if(!m_pMagArray) return 0;
  return m_pMagArray->GetX(m_nTotalChannels*n + m_nChannel);
}

double CMagChannel::GetY(int n) {
  if(!m_pMagArray) return 0;
  return m_pMagArray->GetY(m_nTotalChannels*n + m_nChannel);
}

double CMagChannel::GetZ(int n) {
  if(!m_pMagArray) return 0;
  return m_pMagArray->GetZ(m_nTotalChannels*n + m_nChannel);
}

double CMagChannel::GetField(int n) {
  if(!m_pMagArray) return 0;
  return m_pMagArray->GetField(m_nTotalChannels*n + m_nChannel);
}

double CMagChannel::GetSyntField(int n) {
  if(!m_pMagArray) return 0;
  return m_pMagArray->GetSyntField(m_nTotalChannels*n + m_nChannel);
}


int CMagChannel::GetSyntheticPoints() {
  if(!m_pMagArray) return 0;
  return m_pMagArray->GetSyntheticPoints()/m_nTotalChannels;
}


// implementation for magdata collection

int CMagChannel::get_data_origin(double *x, double *y, double *z){
  if(!m_pMagArray) return 0;
  return m_pMagArray->get_data_origin(x, y, z);
}


int CMagChannel::reset() {
  m_nPos = 0;
  return 1;
}

int CMagChannel::get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				   double *add_parameters, int size_in, int *size_out) {

  if(!m_pMagArray) return 0;
  if(m_nPos >=  GetTotalPoints()) return 0;

  *x = GetX( m_nPos);
  *y = GetY( m_nPos);
  *z = GetZ( m_nPos);
  *T = GetField(m_nPos);
  *n_pos = m_nPos;
  GetAllPositionParams(add_parameters, size_in);

  m_nPos++;

  return 1;
}


int CMagChannel::get_at_position(int n_pos, double *x, double *y, double *z, double * T,
				  double *add_parameters, int size_in, int *size_out) {

  if(!m_pMagArray) return 0;
  if(n_pos <0 || n_pos >= GetTotalPoints()) return 0;

  *x = GetX(n_pos);
  *y = GetY(n_pos);
  *z = GetZ(n_pos);
  *T = GetField(n_pos);

  GetAllPositionParams(add_parameters, size_in);

  return 1;

}



int CMagChannel::get_total_data_points() {

  if(!m_pMagArray) return 0;
  return m_pMagArray->get_total_data_points()/m_nTotalChannels;
}


int CMagChannel::set_synthetic_data(double *data, int replace) {

  if(!m_pMagArray) return 0;
  m_pMagArray->set_synthetic_data_channel(data, m_nChannel,replace);

  return 0;
}


int CMagChannel::get_position(int n, double *xd, double *yd, double *zd) {

  double T;
  if(!m_pMagArray) return 0;
  return m_pMagArray->get_at_position(n*m_nTotalChannels+m_nChannel, xd, yd, zd, &T);

}



int CMagChannel::GetBoundingBox(double *x1, double *x2, double *y1, double *y2,
				double *z1, double *z2) {

  if(!m_pMagArray) return 0;
  return m_pMagArray->GetBoundingBox(x1,x2, y1, y2, z1, z2);

}



int CMagChannel::GetMinMax(int type, double * data_min, double * data_max) {

	int i = 0;
	int n = GetTotalPoints();
	double f;

	for(i=0; i<n; i++) {
		switch(type) {
		case 1: 
			f = GetSyntField(i); break;
		default:
			f = GetField(i); break;
		}
	

	 if(!i) {
			 *data_min = f;
			 *data_max = f;
			 continue;
		 }

		 *data_max = MAX(f, *data_max);
		 *data_min = MIN(f, *data_min);

	 }

	return 1;
}

/*void CMagChannel::DumpParameters(FILE * dat){

  assert(dat);
  CTestProfile::DumpParameters(dat);

  if(m_pPosfixed[X_OFFSET])
    fprintf(dat,"X-shift = %lg FIXED\n", m_lfPos[X_OFFSET]);
  else
    fprintf(dat,"X-shift = %lg +/- %lg\n", m_lfPos[X_OFFSET], m_lfStdDevPos[X_OFFSET]);

  if(m_pPosfixed[Y_OFFSET])
    fprintf(dat,"Y-shift = %lg FIXED\n", m_lfPos[Y_OFFSET]);
  else
    fprintf(dat,"Y-shift = %lg +/- %lg\n", m_lfPos[Y_OFFSET], m_lfStdDevPos[Y_OFFSET]);


  if(m_pPosfixed[Z_OFFSET])
    fprintf(dat,"Z-shift = %lg FIXED\n", m_lfPos[Z_OFFSET]);
  else
    fprintf(dat,"Z-shift = %lg +/- %lg\n", m_lfPos[Z_OFFSET], m_lfStdDevPos[Z_OFFSET]);

}*/
