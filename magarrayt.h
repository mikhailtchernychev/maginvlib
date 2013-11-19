//
// @(#)magarrayt.h  1.1  misha1-06oct108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  06oct108  by  misha
// Version:  06oct108  by  misha
//


#include "magdata.h"
#include "magfunc.h"

#include <assert.h>
#include <math.h>

#ifndef MAGARRAY_H
#define MAGARRAY_H

#include <vector>
#include <iomanip>

template <class T = const float>
class CMagArray : public CMagData {

protected:

  std::vector<const T *> m_vX;      // positions of external data
  std::vector<const T *> m_vY;      // do not delete
  std::vector<const T *> m_vZ; 
  std::vector<const T *> m_vAlt;
  std::vector<const T *> m_vT;
  std::vector<T *> m_vTsynt;
  const double * m_pTime;        // do not delete

  int m_nPoints;          //  total data points in external arrays per channel
  int m_nWindowStart;     //  point where data window starts
  int m_nWindowStop;      //  point where data window stops

  double * m_pDistance;   //  distance along the profile (average of all channels, delete in dtor)
  int m_nCurrentPos;      //  current position in array
  int m_nCurrentChannel;  //  current channel in the array

  double m_lfWindowWidth; // selected width of the window

  double m_lfXmiddle;     // middle point
  double m_lfYmiddle;
  double m_lfZmiddle;
  double m_lfX1, m_lfX2, m_lfY1, m_lfY2, m_lfZ1, m_lfZ2; // data limits

  int m_iIsSmoothedDepth;

  void clean() {

    if(m_pDistance) { delete [] m_pDistance; m_pDistance = NULL; }

    m_pTime = NULL;

    m_vTsynt.clear();
    m_vX.clear();
    m_vY.clear();
    m_vAlt.clear();
    m_vT.clear();

    m_nPoints         = 0;
    m_nWindowStart     = 0;
    m_nWindowStop     = 0;
    m_nCurrentPos     = 0;
    m_nCurrentChannel = 0;
    m_lfWindowWidth   = 0.;
    m_lfXmiddle       = 0.;
    m_lfYmiddle       = 0.;
    m_lfZmiddle       = 0.;
    m_lfX1 =  m_lfX2 =  m_lfY1 = m_lfY2 = m_lfZ1 =  m_lfZ2 = 0.;

	if(m_iIsSmoothedDepth) {
		for(int i=0; i<m_vZ.size(); i++) {
			T * z = (T *)m_vZ[i];
			if(z) delete [] z;
		}
	}

	m_vZ.clear();

  }

public:
    CMagArray(){
      m_vX.clear();
      m_vY.clear();
      m_vZ.clear();
      m_vAlt.clear();
      m_vT.clear();
      m_vTsynt.clear();
      m_pDistance       = NULL;
      m_nPoints         = 0;
      m_nWindowStart    = 0;
      m_nWindowStop     = 0;
      m_nCurrentPos     = 0;
      m_nCurrentChannel = 0;
      m_lfWindowWidth   = 0.;
      m_lfXmiddle       = 0.;
      m_lfYmiddle       = 0.;
      m_lfZmiddle       = 0.;
      m_pTime           = NULL;
      m_lfX1 =  m_lfX2 =  m_lfY1 = m_lfY2 = m_lfZ1 =  m_lfZ2 = 0.;
	  m_iIsSmoothedDepth = 0;

    };

  virtual ~CMagArray(){
    clean();
  };

  int SetArrays(int n_channels, int n, T ** x, T ** y, T ** field, T ** synt,
	  T ** z = NULL, T **alt=NULL, double * pTime = NULL) {

	  int i,j;

	  clean();

	  m_pDistance = new double[n];
	  if(!m_pDistance) {
		  clean();
		  return 0;
	  }


	  double x_av1, y_av1, x_av2, y_av2;

	  x_av1 = y_av1 = x_av2 = y_av2 = 0.;
	  m_pDistance[0] = 0.;

	  for(i=0; i<n; i++) {

		  x_av2 = y_av2 = 0.;
		  for(j=0; j<n_channels; j++) {
			  x_av2 += x[j][i];
			  y_av2 += x[j][i];
		  }

		  x_av2 /= n_channels;
		  y_av2 /= n_channels;

		  if(i) {
			  m_pDistance[i] = m_pDistance[i-1] +
				  sqrt((x_av2-x_av1)*(x_av2-x_av1) + (y_av2-y_av1)*(y_av2-y_av1));
		  }

		  x_av1 = x_av2;
		  y_av1 = y_av2;

	  }

	  for(i=0; i<n_channels; i++) {
		  m_vX.push_back(x[i]);
		  m_vY.push_back(y[i]);
		  m_vT.push_back(field[i]);
		  m_vTsynt.push_back(synt[i]);

		  if(z)   m_vZ.push_back(z[i]);
		  if(alt) m_vAlt.push_back(alt[i]);
	  }

	  m_nPoints = n;
	  m_pTime   = pTime;

	  return 1;
  }



  void SetWindow(int n1, int n2) {
    assert(n1<m_nPoints && n1>=0);
    assert(n2<m_nPoints && n2>=0);
    assert(n1< n2);

    m_nWindowStart = n1;
    m_nWindowStop  = n2;

    m_nCurrentPos      = n1;
    m_nCurrentChannel  = 0;

    m_lfX1 =  m_lfX2 =  m_lfY1 = m_lfY2 = m_lfZ1 =  m_lfZ2 = 0.;

    int n_channels = m_vX.size();
    int n = n2-n1+1;
    int is_start = 1;

    m_lfXmiddle     = 0.;
    m_lfYmiddle     = 0.;
    m_lfZmiddle     = 0.;

    for(int i=0; i<n_channels; i++) {

		for(int j=n1; j<=n2; j++) {

		m_lfXmiddle += m_vX[i][j];
		m_lfYmiddle += m_vY[i][j];

		if(m_vZ.size()) m_lfZmiddle += m_vZ[i][j];

			if(is_start) {
				m_lfX1 = m_lfX2 = m_vX[i][j];
				m_lfY1 = m_lfY2 = m_vY[i][j];
				if(m_vZ.size()) m_lfZ1 = m_lfZ2 = m_vZ[i][j];
				is_start = 0;
			}
			else {
				m_lfX1 = MIN(m_lfX1, m_vX[i][j]);   m_lfX2 = MAX(m_lfX2,  m_vX[i][j]);
				m_lfY1 = MIN(m_lfY1, m_vY[i][j]);   m_lfY2 = MAX(m_lfY2,  m_vY[i][j]);
				if(m_vZ.size()) {
					m_lfZ1 = MIN(m_lfZ1, m_vZ[i][j]); m_lfZ2 = MAX(m_lfZ2,  m_vZ[i][j]);
				}
			}
		}
    }

    int n_total = (n2-n1+1)*n_channels;
    m_lfXmiddle /= double(n_total);
    m_lfYmiddle /= double(n_total);
    m_lfZmiddle /= double(n_total);


}


 double GetMiddleX() { return m_lfXmiddle; }
 double GetMiddleY() { return m_lfYmiddle; }
 double GetMiddleZ() { return m_lfZmiddle; }



int SetWindow(double w) {

	int n1;

	m_nWindowStart = 0;
	m_nWindowStop  = m_nPoints-1;

	m_nCurrentPos     = 0;
	m_nCurrentChannel = 0;


	if(!FromDtoN(w, &n1)) return 0;

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
     *x = m_lfXmiddle;
     *y = m_lfYmiddle;
     *z = m_lfZmiddle;
    return 1;
  }

 virtual int reset() {  m_nCurrentPos = m_nWindowStart; m_nCurrentChannel = 0; return 1;}

 virtual int GetAllPositionParams(double * params, int size) {

   CMagObject::GetAllPositionParams(params, size);

   if(params && size > FINITE_SEPARATION) {
     params[DIRECTION_X]       = m_lfNx;
     params[DIRECTION_Y]       = m_lfNy;
     params[DIRECTION_Z]       = m_lfNz;
     params[FINITE_SEPARATION] = m_lfSeparation;
     return 1;
   }
   return 0;
 }


 virtual int get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
			       double *add_parameters = NULL, int size_in=0, int *size_out=NULL) {


	if(m_vX.size()<=0) return 0;
	if(m_nCurrentChannel >= m_vX.size()) return 0;


	int n = m_nWindowStop - m_nWindowStart +1;

	*n_pos = m_nCurrentChannel*n+(m_nCurrentPos-m_nWindowStart);
	*xd    = m_vX[m_nCurrentChannel][m_nCurrentPos];
	*yd    = m_vY[m_nCurrentChannel][m_nCurrentPos];
	*Td    = m_vT[m_nCurrentChannel][m_nCurrentPos];

	if(m_vZ.size())
	  *zd    = m_vZ[m_nCurrentChannel][m_nCurrentPos];
	else
	  *zd    = 0;

	GetAllPositionParams(add_parameters, size_in);

	m_nCurrentPos++;
	if(m_nCurrentPos > m_nWindowStop ) {
	  m_nCurrentChannel++;
	  m_nCurrentPos=m_nWindowStart;
	}


	return 1;
 }


virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * Td,
								double *add_parameters = NULL, int size_in=0, int *size_out=NULL) {

	if(m_vX.size()<=0) return 0;

	int n_total   =  m_nWindowStop - m_nWindowStart +1;
	int n_channel =  n_pos/n_total;
	int n         =  n_pos - n_channel*n_total + m_nWindowStart;

	if(n_channel <0 || n_channel >= m_vX.size()) return 0;
	if(n<m_nWindowStart || n > m_nWindowStop) return 0;

	*x    = m_vX[n_channel][n];
	*y    = m_vY[n_channel][n];
	*Td   = m_vT[n_channel][n];

	if(m_vZ.size())
	  *z =  m_vZ[n_channel][n];
	else
	  *z = 0.;

	GetAllPositionParams(add_parameters, size_in);

	return 1;

}



    virtual int get_total_data_points() { return m_vX.size()*(m_nWindowStop-m_nWindowStart+1); }


    virtual int set_synthetic_data(double *data, int replace = 0) {

		int n_channels = m_vX.size();
		if(n_channels<=0) return 0;

		int n = m_nWindowStop-m_nWindowStart+1;
		int k = 0;

		for(int j=0; j<n_channels; j++) {
			for(int i=0; i<n; i++) {
				m_vTsynt[j][i+m_nWindowStart] = data[k];
				k++;
			}
		}
		return 1;
    }


    virtual int get_position(int n_pos, double * xd, double * yd, double * zd) {

	if(m_vX.size()<=0) return 0;

	int n_total   =  m_nWindowStop - m_nWindowStart +1;
	int n_channel =  n_pos/n_total;
	int n         =  n_pos - n_channel*n_total + m_nWindowStart;

	assert(n_channel >=0 && n_channel <= m_vX.size());
	assert(n >=m_nWindowStart && n <= m_nWindowStop);

	*xd    = m_vX[n_channel][n];
	*yd    = m_vY[n_channel][n];

	if(m_vZ.size())
	  *zd =  m_vZ[n_channel][n];
	else
	  *zd = 0.;

	return 1;
    }


    virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2,
		double *z1, double *z2) {

      *x1 = m_lfX1; *x2 = m_lfX2;
      *y1 = m_lfY1; *y2 = m_lfY2;
      *z1 = m_lfZ1; *z2 = m_lfZ2;

      return 1;
    }

    virtual void Dump(FILE * dat) {

		int i, j;
		char date_txt[20], time_txt[20];

		int n_channels = m_vX.size();
		if(n_channels<=0) return;

		for(i=0; i<n_channels; i++) {
			j = i+1;
			fprintf(dat,"X%d   Y%d  Z%d  ALT%d  T%d  Tsynt%d     ", j,j,j,j,j,j);
		}

		if(m_pTime) fprintf(dat,"DATE TIME LINE\n");
		else fprintf(dat,"LINE\n");

		for(i=m_nWindowStart; i<=m_nWindowStop; i++) {

			if(m_pTime) sec2date(m_pTime[i], date_txt, time_txt, 1);

			for(j=0; j<n_channels; j++) {

				fprintf(dat,"%lf  %lf  %lf  %lf  %lf %lf ",
					m_vX[j][i] + m_lfPos[X_OFFSET],
					m_vY[j][i] + m_lfPos[Y_OFFSET],
					(m_vZ.size()>0)   ? m_vZ[j][i]   +  m_lfPos[Z_OFFSET] : 0.,
					(m_vAlt.size()>0) ? m_vAlt[j][i] +  m_lfPos[Z_OFFSET] : 0.,
					m_vT[j][i],
					m_vTsynt[j][i]);
			}

				if(m_pTime) fprintf(dat,"%s %s %s\n", date_txt, time_txt, GetName());
				else fprintf(dat,"%s\n", date_txt, time_txt, GetName());
		}

    }


   virtual void DumpToXMLStream(std::ostream & str) {

    int i,j;
    char delimit = ' ';
    char date_txt[20], time_txt[20];
    std::string data_type;
    char quot = '\"';

    GetFieldTypeText(data_type);


    str << "<" << m_sObjectType << ">" << std::endl;

	str << "<params ";
    str << "data_type=" << quot << data_type << quot << delimit;
	str << "/>" << std::endl;

	str << "<header header=\"";

	for(i=0; i<m_vX.size(); i++) {
		str << "X"     << i+1 << delimit;
		str << "Y"     << i+1 << delimit;
		str << "Z"     << i+1 << delimit;
		str << "T"     << i+1 << delimit;
		str << "ALT"   << i+1 << delimit;
		str << "Tsynt" << i+1 << delimit;
	}

	if(m_pTime) str << "DATE TIME" << delimit;   
	str << "LINE\"/>" << std::endl;

	str << "<![CDATA[" << std::endl;

	int prec = str.precision();
	str << std::fixed << std::setprecision(3);

	for(i=m_nWindowStart; i<=m_nWindowStop; i++) {

		if(m_pTime) sec2date(m_pTime[i], date_txt, time_txt, 1);

		for(j=0; j<m_vX.size(); j++) {

			    str << m_vX[j][i] + m_lfPos[X_OFFSET] << delimit;
				str << m_vY[j][i] + m_lfPos[Y_OFFSET] << delimit;
				str << ((m_vZ.size()>0)   ? m_vZ[j][i]   +  m_lfPos[Z_OFFSET] : 0.) << delimit;
				str << ((m_vAlt.size()>0) ? m_vAlt[j][i] +  m_lfPos[Z_OFFSET] : 0.) << delimit;
				str << m_vT[j][i] << delimit;
				str << m_vTsynt[j][i] << delimit;
				if(m_pTime) str << date_txt << " " << time_txt << delimit;
				str << GetName() << std::endl;
		}
	}


	str << std::resetiosflags(ios_base::fixed);
    str << std::setprecision(prec);

	str << "]]>" << std::endl;
	str << "</" << m_sObjectType << ">" << std::endl;

	}




  virtual int GetMinMax(int type, double * data_min, double * data_max) {

		if(m_vT.size()<=0 || m_vTsynt.size()<=0) return 0;

		int n_channels = m_vX.size();
		double dmax, dmin;
		int n = m_nWindowStop - m_nWindowStart + 1;

		if(type==1) {
			dmax = dmin = m_vTsynt[0][m_nWindowStart];
			for(int j=0; j<n_channels; j++) {
			  for(int i=0; i<n; i++) {
			    dmax = MAX(dmax, m_vTsynt[j][m_nWindowStart+i]);
			    dmin = MIN(dmin, m_vTsynt[j][m_nWindowStart+i]);
			  }
			}
		}
		else {
			dmax = dmin = (double)m_vT[0][m_nWindowStart];
			for(int j=0; j<n_channels; j++) {
			  for(int i=0; i<n; i++) {
			    dmax = MAX(dmax, (double)m_vT[j][m_nWindowStart+i]);
			    dmin = MIN(dmin, (double)m_vT[j][m_nWindowStart+i]);
			  }
			}
		}

		*data_min = dmin;
		*data_max = dmax;

		return 1;

	}


	void GetWindow(int *n1, int *n2) {
		*n1 = m_nWindowStart;
		*n2 = m_nWindowStop;
	}


  int FromDtoN(double x0, int * n_out) {

    if(!m_pDistance) return 0;

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

     *n_out = n1;
     return 1;
  }


   const char * GetName() { return m_sObjectName; }


	virtual int GetNumberOfChannels() { return m_vX.size(); }
	virtual int GetNDataPerChannel()  { return m_nWindowStop-m_nWindowStart+1; }


    void DepthSplineSmooth(double dm) {

    int n, i, j, ret2;

	for(i=0; i<m_vZ.size(); i++) {
		T * z = (T *)m_vZ[i];	 
		T * znew = new T[m_nPoints];
		for(j=0; j<m_nPoints; j++) znew[j] = z[j];
		m_vZ[i] = znew;
	}


	m_iIsSmoothedDepth =  1;


	for(i=0;  i<m_vZ.size(); i++) {

		T * z = (T *)m_vZ[i];
		n = m_nWindowStop - m_nWindowStart + 1;
		
        double * x = new double[n];
        double * u = new double[n];
        double * a = new double[n];
        double * b = new double[n];
        double * c = new double[n];
        double * d = new double[n];
		
        if(!x || !u || !a || !b || !c || !d) goto m_ret;
		
        for(j=0; j<n; j++) {
			x[j] = double(j);
			u[j] = z[m_nWindowStart+j];
		}
		
		for(j=0; j<n; j++) x[i] = double(i);
		
		ret2 = splik(n, x, u, dm, a, b,c, d); 
		
		for(j=0; j<n; j++)  z[m_nWindowStart+j] = d[j];    
		
m_ret:
		if(x) delete [] x;
		if(u) delete [] u;      
		if(a) delete [] a;
		if(b) delete [] b;
		if(c) delete [] c;
		if(d) delete [] d;
		
		}


	}

	 int GetNearestDataPoint(double x0, double y0, double * d, double *x, double *y, double *dtime,
		 	              std::string & line_name) {


		int i, j, k;
		double x1, y1, dist1, dist2;

		int n_channels = m_vX.size();
		if(n_channels<=0) return -1;

		k = 0;

		for(i=m_nWindowStart; i<=m_nWindowStop; i++) {

			for(j=0; j<n_channels; j++) {

				x1 = m_vX[j][i] + m_lfPos[X_OFFSET];
				y1 = m_vY[j][i] + m_lfPos[Y_OFFSET];

				dist1 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));

				if(!k || dist1 < dist2) {
					*d = dist1; *x = x1; *y=y1; 
					if(m_pTime) *dtime =  m_pTime[i]; else *dtime =  0.;
					line_name = GetName();
				}

				dist2 = dist1;
				k++;
			}

		}

		return 1;
	 }
};

#endif



