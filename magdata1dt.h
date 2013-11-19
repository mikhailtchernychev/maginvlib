//
// @(#)magdata1dt.h  1.1  misha1-01oct108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  01oct108  by  misha
// Version:  01oct108  by  misha
//


#include "magdata.h"
#include <assert.h>

#ifndef MAGDATA1DT_H
#define MAGDATA1DT_H


template <class T = const float>
class CMagData1Dt : public CMagData {

protected:
  T * m_pX;              // positions of external data
  T * m_pY;              // do not remove
  T * m_pZ;
  T * m_pAlt;
  T * m_pT;

  int m_nPoints;         //  total data points in external arrays
  int m_nWindowStart;    //  point where data window starts
  int m_nWindowStop;     //  point where data window stops

  double * m_pTsynt;     //  synthetic magnetic field (also storage array)
  double * m_pDistance;  //  distance along the profile
  int m_nCurrent;        //  current position in array
  int m_nStart;          //  point with X=0 in current window

  double m_lfWindowWidth; // selected width of the window

  void clean() {
	  if(m_pTsynt)    { delete [] m_pTsynt;    m_pTsynt    = NULL; }
	  if(m_pDistance) { delete [] m_pDistance; m_pDistance = NULL; }
  }

public:
    CMagData1Dt(){
		m_pX = m_pY = m_pZ = m_pAlt = m_pT = NULL;
		m_pTsynt        = NULL;
		m_pDistance     = NULL;
		m_nPoints       = 0;
		m_nWindowStart  = 0;
		m_nWindowStop   = 0;
		m_nCurrent      = 0;
		m_nStart        = 0;
		m_lfWindowWidth = 0.;
	};

    virtual ~CMagData1Dt(){
		clean();
	};

	int SetArrays(T * x, T * y, T * z, T *alt, T * field, int n) {

		clean();

		m_pDistance = new double[n];
		m_pTsynt    = new double[n];

		if(!m_pDistance || !m_pTsynt) {
			clean();
			return 0;
		}

		m_pX      = x;
		m_pY      = y;
		m_pZ      = z;
		m_pAlt    = alt;
		m_pT      = field;
		m_nPoints = n;

		double d = 0;
		m_pDistance[0] = m_pTsynt[0] = 0.;

		for(int i=1; i<n; i++) {
			d += sqrt((m_pX[i-1]-m_pX[i])*(m_pX[i-1]-m_pX[i]) + (m_pY[i-1]-m_pY[i])*(m_pY[i-1]-m_pY[i]));
			m_pDistance[i] = d;
			m_pTsynt[i]    = 0.;
		}

		return 1;
	}

	void SetWindow(int n1, int n2) {
		assert(n1<m_nPoints && n1>=0);
		assert(n2<m_nPoints && n2>=0);
		assert(n1< n2);

		m_nWindowStart = n1;
		m_nWindowStop  = n2;

		m_nCurrent = 0;
		m_nStart   =  (n1+n2)/2;

	}


	int SetWindow(double w) {

		double x, y;
		int n1;

		m_nWindowStart = 0;
		m_nWindowStop  = m_nPoints-1;

		m_nCurrent =  0;
		m_nStart   =  0;

		if(!From1Dto2D(w, &x, &y, NULL, NULL, NULL, &n1)) return 0;

		SetWindow(0, n1+1);

		m_lfWindowWidth = w;
		return 1;

	}

	int Advance() {

		if(m_nWindowStop >= m_nPoints-1) return 0;

		int n1 = m_nWindowStart+1;
		int n2 = m_nWindowStop;

		for(;n2<m_nPoints;n2++) {
			if((m_pDistance[n2] - m_pDistance[n1])> m_lfWindowWidth) break;
		}

		if(n2>=m_nPoints) n2 = m_nPoints-1;

		SetWindow(n1, n2);
		return 1;
	}


	void DumpWindow(FILE * dat) {
		fprintf(dat,"%d %d %d %lf %lf\n", m_nWindowStop-m_nWindowStart, m_nWindowStart, m_nWindowStop,
			m_pDistance[m_nWindowStop]-m_pDistance[m_nWindowStart],
			m_pDistance[m_nWindowStop-2]-m_pDistance[m_nWindowStart]);
	}

    virtual int get_data_origin(double *x, double *y, double *z) {
		*x = *y = *z = 0;
	    return 1;
	}

    virtual int reset() {  m_nCurrent = 0; return 1;}


	virtual int get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
		double *add_parameters = NULL, int size_in=0, int *size_out=NULL) {


		if(!m_pT || !m_pDistance) return 0;

		*n_pos = m_nCurrent;
		*xd    = m_pDistance[m_nWindowStart + m_nCurrent]-m_pDistance[m_nStart];
		*yd    = 0;
		*zd    = 0;
		*Td    = m_pT[m_nWindowStart + m_nCurrent];

		m_nCurrent++;

		return 1;

	}


    virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * Td,
								double *add_parameters = NULL, int size_in=0, int *size_out=NULL) {

		if(!m_pT || &m_pDistance) return 0;
		assert(n_pos < (m_nWindowStop-m_nWindowStart+1) && n_pos >= 0);

		*x    = m_pDistance[m_nWindowStart + n_pos]-m_pDistance[m_nStart];
		*y    = 0;
		*z    = 0;
		*Td   = m_pT[m_nWindowStart + n_pos];

		return 1;

	}



    virtual int get_total_data_points() { return m_nWindowStop-m_nWindowStart+1; }

    virtual int set_synthetic_data(double *data, int replace = 0) {

		if(!m_pT) return 0;
		int n = m_nWindowStop-m_nWindowStart+1;

		for(int i=0; i<n; i++) {
			m_pTsynt[i+m_nWindowStart] = data[i];
		}

		return 1;
	}


    virtual int get_position(int n, double * xd, double * yd, double * zd) {

		if(!m_pDistance) return 0;
		assert(n < (m_nWindowStop-m_nWindowStart+1) && n >= 0);

		*xd = m_pDistance[m_nWindowStart + n]-m_pDistance[m_nStart];
		*yd = *zd = 0.;

		return 1;

	}

    virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2,
		double *z1, double *z2) {

		if(!m_pDistance) return 0;

		*x1 = m_pDistance[m_nWindowStart] - m_pDistance[m_nStart];
		*x2 = m_pDistance[m_nWindowStop]  - m_pDistance[m_nStart];

		*y1 = *y2 = *z1 = *z2 = 0;

		return 1;

	}

    virtual void Dump(FILE * dat) {

		if(!m_pT || !m_pDistance) return;

		fprintf(dat,"X Y DISTANCE T Tsynt\n");

		int n = m_nWindowStop - m_nWindowStart + 1;

		for(int i=0; i<n; i++) {
			fprintf(dat,"%lf %lf %lf %lf %lf\n",
				double(m_pX[m_nWindowStart+i]),
				double(m_pY[m_nWindowStart+i]),
				double(m_pDistance[m_nWindowStart+i]),
				double(m_pT[m_nWindowStart+i]),
				double(m_pTsynt[m_nWindowStart+i]));
		}

	}


	virtual int GetMinMax(int type, double * data_min, double * data_max) {

		if(!m_pT) return 0;

		double dmax, dmin;
		int n = m_nWindowStop - m_nWindowStart + 1;

		if(type==1) {
			dmax = dmin = m_pTsynt[m_nWindowStart];
			for(int i=0; i<n; i++) {
				dmax = MAX(dmax, m_pTsynt[m_nWindowStart+i]);
				dmin = MIN(dmin, m_pTsynt[m_nWindowStart+i]);
			}
		}
		else {
			dmax = dmin = (double)m_pT[m_nWindowStart];
			for(int i=0; i<n; i++) {
				dmax = MAX(dmax, (double)m_pT[m_nWindowStart]);
				dmin = MIN(dmin, (double)m_pT[m_nWindowStart]);
			}
		}

		*data_min = dmin;
		*data_max = dmax;

		return 1;

	}



	virtual int GetSyntheticField(double *Tsynt, int n) {

		int k = m_nWindowStop - m_nWindowStart + 1;

		if(!Tsynt || n<k) return 0;

		for(int i=0; i<k; i++) {
			Tsynt[i] =  m_pTsynt[m_nWindowStart+i];
		}

		return 1;

	}


	void GetWindow(int *n1, int *n2) {
		*n1 = m_nWindowStart;
		*n2 = m_nWindowStop;
	}


	int From1Dto2D(double x0, double *x_out, double *y_out,
				   double *dist_out=NULL, double *z_out=NULL, double *alt_out=NULL,
				   int * n_out = NULL) {

		if(!m_pDistance) return 0;

		x0 += m_pDistance[m_nStart];

		if(x0 < m_pDistance[m_nWindowStart] || x0 > m_pDistance[m_nWindowStop]) return 0;


		int n1 = m_nWindowStart;
		int n2 = m_nWindowStop;

		int k = (n1+n2)/2;

		while(1) {

			if(x0>=m_pDistance[n1]  && x0<= m_pDistance[n2] && n2<=n1+1) break;

			if(x0 < m_pDistance[k]) {
				n2 = k;
			}
			else {
				n1 = k;
			}

			k = (n2+n1)/2;
		}

		if(n_out) *n_out = n1;

		if(n1==n2) {
			*x_out = m_pX[n1];
			*y_out = m_pY[n1];
			return 1;
		}

		double t = (x0-m_pDistance[n1]) / (m_pDistance[n2]  - m_pDistance[n1]);

		*x_out = m_pX[n1]+t*(m_pX[n2]-m_pX[n1]);
		*y_out = m_pY[n1]+t*(m_pY[n2]-m_pY[n1]);


		if(dist_out)          *dist_out = x0; 
		if(m_pZ && z_out)     *z_out    = m_pZ[n1]   + t*(m_pZ[n2]-m_pZ[n1]);
		if(m_pAlt && alt_out) *alt_out  = m_pAlt[n1] + t*(m_pAlt[n2]-m_pAlt[n1]);

		return 1;
	}

};

#endif



