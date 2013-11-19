//
// @(#)magdata1d.cpp  1.1  misha1-29sep108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  29sep108  by  misha
// Version:  29sep108  by  misha
//


#include "magdata1d.h"



CMagData1D::CMagData1D () {

  m_pT       = NULL;
  m_pTsynt   = NULL;
  m_lfDx     = 0;
  m_nPoints  = 0;
  m_nCurrent = 0;
  m_nStart   = 0;
}

void CMagData1D::clean () {

  if(m_pT)      { delete [] m_pT;      m_pT     = NULL; }
  if(m_pTsynt)  { delete [] m_pTsynt;  m_pTsynt = NULL; }

  m_lfDx     = 0;
  m_nPoints  = 0;
  m_nCurrent = 0;
  m_nStart   = 0;

}

CMagData1D::~CMagData1D () {

  clean();
}


CMagData1D::CMagData1D (int n, double *x, double *y, double *z, double *T) {

  // allocate arrays;

  m_pT     = new double [n];
  m_pTsynt = new double [n];

  m_lfDx     = 0;
  m_nPoints  = 0;
  m_nCurrent = 0;
  m_nStart   = 0;

  if(!m_pT || !m_pTsynt) {
    clean();
    return;
  }


  // compute average step;

  double dx = 0., a_tmp;

  m_pT[0]     = T[0];
  m_pTsynt[0] = 0.;

  for(int i=1; i<n; i++) {

    m_pT[i]     =  T[i];
    m_pTsynt[i] =  0.;

    a_tmp  = (x[i]-x[i-1])*(x[i]-x[i-1]) + (y[i]-y[i-1])*(y[i]-y[i-1]);
    if(z) a_tmp += (z[i]-z[i-1])*(z[i]-z[i-1]);
    dx += sqrt(a_tmp);
  }

  m_lfDx = dx / float(n);

  if(m_lfDx == 0.) {
    clean();
    return;
  }

  m_nPoints  = n;
  m_nCurrent = 0;
  m_nStart   = n/2;


}


CMagData1D::CMagData1D(int n, double dx) {

  m_pT     = new double [n];
  m_pTsynt = new double [n];

  m_lfDx     = dx;
  m_nPoints  = 0;
  m_nCurrent = 0;
  m_nStart   = 0;

  if(!m_pT || !m_pTsynt) {
    clean();
    return;
  }

  m_nPoints  = n;
  m_nCurrent = 0;
  m_nStart   = n/2;

  for(int i=0; i<n; i++) {
	  m_pT[i] = m_pTsynt[i] = 0;
  }

}

int CMagData1D::get_data_origin(double *x, double *y, double *z) {

  *x = *y = *z = 0.;
  return 1;

}


int CMagData1D::get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				  double *add_parameters, int size_in, int *size_out) {


  if(!m_pT) return 0;
  if(m_nCurrent>= m_nPoints) return 0;	

  *n_pos =  m_nCurrent;
  *xd    = (m_nCurrent-m_nStart)*m_lfDx;
  *yd    = 0;
  *zd    = 0;
  *Td    = m_pT[m_nCurrent];
   
  GetAllPositionParams(add_parameters, size_in);

  m_nCurrent++;
  

  return 1;

}


int CMagData1D::get_at_position(int n_pos, double *x, double *y, double *z, double * T,
				double *add_parameters, int size_in, int *size_out) {

  if(!m_pT) return 0;
  assert(n_pos < m_nPoints && n_pos >= 0);

  *x    = (n_pos-m_nStart)*m_lfDx;
  *y    = 0;
  *z    = 0;
  *T    = m_pT[n_pos];

  GetAllPositionParams(add_parameters, size_in);

  return 1;

}


int CMagData1D::set_synthetic_data(double *data, int replace) {

    if(!m_pT) return 0;

    for(int i=0; i<m_nPoints; i++) {
      m_pTsynt[i] = data[i];
      if(replace) m_pT[i] = data[i];
    }

    return 1;

}


int CMagData1D::get_position(int n, double * xd, double * yd, double * zd) {

  if(!m_pT) return 0;
  assert(n < m_nPoints && n >= 0);

  *xd = (n-m_nStart)*m_lfDx;

  *yd = *zd = 0.;

  return 1;

}


int CMagData1D::GetBoundingBox(double *x1, double *x2, double *y1, double *y2,
			       double *z1, double *z2) {

  *x1 = -m_nStart*m_lfDx;
  *x2 = (m_nPoints-m_nStart-1)*m_lfDx;

  *y1 = *y2 = *z1 = *z2 = 0;

  return 1;
}


void CMagData1D::Dump(FILE * dat) {

  if(!m_pT) return;

  for(int i=0; i<m_nPoints; i++) {
    fprintf(dat,"%lf %lf %lf\n", (i-m_nStart)*m_lfDx, m_pT[i], m_pTsynt[i]);

  }

}


int CMagData1D::GetSyntheticField(double *Tsynt, int n) {

	if(!Tsynt || n<m_nPoints) return 0;

	  for(int i=0; i<m_nPoints; i++) {
		Tsynt[i] =  m_pTsynt[i];
	  }

	return 1;
}


int CMagData1D::GetMinMax(int type, double * data_min, double * data_max) {

  if(!m_pT) return 0;

  double dmax, dmin;
  double *pData = m_pT;

  if(type==1) pData = m_pTsynt;

  dmax = dmin = pData[0];
  for(int i=0; i<m_nPoints; i++) {
        dmax = MAX(dmax, pData[i]);
        dmin = MIN(dmin, pData[i]);
    }

    *data_min = dmin;
    *data_max = dmax;

    return 1;
}


// from external to internal coordinates
double  CMagData1D::RecalcCoord(double x) {

	return x - m_nStart*m_lfDx;


}

// here x and y should be the same as when created

int CMagData1D::From1Dto2D(double *x, double *y, int n, double x0,
						   double *x_out, double *y_out) {


int k = m_nStart + int(x0/m_lfDx);
if(n<=k+1) return 0;

double t = (x0 - (k-m_nStart)*m_lfDx)/m_lfDx;

*x_out = x[k]+t*(x[k+1]-x[k]);
*y_out = y[k]+t*(y[k+1]-y[k]);

return 1;

}
