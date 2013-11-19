//
// @(#)magarray.h  1.1  misha-11nov104
//
// Copyright (c) 2004 of Geometrics.
// All rights reserved.
//
// Created:  11nov104  by  misha@MISHA
// Version:  11nov104  by  misha@MISHA
//

#ifndef MAGARRAY_H
#define MAGARRAY_H

#include "profwrapper.h"
#include "magdata.h"

#include <vector>
#include <string>
#include <map>

using namespace std;

class CMagProfileArray: public CProfWrapper, public CMagData {

protected:

  vector<double> m_vX;
  vector<double> m_vY;
  vector<double> m_vZ;
  vector<double> m_vAlt;
  vector<double> m_vTime;
  vector<double> m_vField;
  vector<double> m_vComputedField;

  vector<double> m_vXmiddle;
  vector<double> m_vYmiddle;

  vector<std::string> m_vPositionNames;

  int m_nPos;
  double m_lfX1, m_lfX2, m_lfY1, m_lfY2, m_lfZ1, m_lfZ2;

  double m_lfXaverage;
  double m_lfYaverage;
  double m_lfZaverage;

  int m_nChannels;    // number of magnetometers in array;

  int m_iStart;       // limited area - start index;
  int m_iStop;        // limited area - stop index;
  int m_iTotalStart;
  int m_iTotalStop;

  // median depth / alt parameters for the area

  double m_lfMedianDepth;
  double m_llfMedianAlt;

  double m_lfNX;  // unit vector along line section
  double m_lfNY;

  double  m_lfMiddleX;
  double  m_lfMiddleY;

  double m_lfLength;

  void clear_all();

public:

  CMagProfileArray() : CProfWrapper(){
    m_iObjectLinearAbilities    =  1;
    m_iObjectNonLinearAbilities =  1;
    clear_all();
  }

  ~CMagProfileArray() { clear_all(); }


  // sample implementation

  void InitProfile(char * name);
  void SetChannels(int n_channels) { m_nChannels = n_channels;}
  int  GetChannels() { return m_nChannels; }

  void Add(double * x, double * y, double time, double * field, double * z = NULL,
	    double * alt = NULL, double x_middle=0., double y_middle=0., char * name = NULL);

  virtual void Correct();

  void RecomputePositions90(double * shifts90);

  const char * GetName() { return m_sObjectName; }

  // mag array specific

  int  SetStartStop(int start, int stop);
  void GetStartStop(int * start, int * stop) { *start = m_iStart; *stop = m_iStop; }
  int GetMiddlePointCoords(double *x, double *y, double *z, double *alt);
  int GetTotalTimePoints() { return m_vTime.size(); }
  double GetDistanceFromMiddle(double x, double y);

  double GetMedianDepth()  { return m_lfMedianDepth; }
  double GetMedianAlt()    { return m_llfMedianAlt;  }
  double GetAreaLength()   { return m_lfLength;      }
  void GetNxNy(double *nx, double *ny) { *nx=m_lfNX; *ny=m_lfNY;}
  double GetNearestTime(double x, double y);

  // real implementation for cross-analysis

  virtual int GetTotalPoints() { return (m_iStop - m_iStart)*m_nChannels; }
  virtual double GetX(int n);
  virtual double GetY(int n);
  virtual double GetZ(int n);
  virtual double GetAlt(int n);
  virtual double GetTime(int n);
  virtual double GetAbsoluteTime(int n) { return 0; };
  virtual double GetField(int n);
  virtual double GetSyntField(int n);

  int GetSyntheticPoints() { return GetTotalPoints(); }

 // implementation for magdata collection

  virtual int get_data_origin(double *x, double *y, double *z);
  virtual int reset() { m_nPos = 0; return 0; };  // resets data pointer to start

  virtual int get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL);


  virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T,
			      double *add_parameters = NULL, int size_in=0, int *size_out=NULL);

  virtual int get_total_data_points() { return GetTotalPoints(); }

  virtual int set_synthetic_data(double *data, int replace=0);
  int set_synthetic_data_channel(double *data, int channel, int replace=0);

  virtual int is_valid() {return (GetTotalPoints() > 0); }
  virtual int get_position(int n, double *xd, double *yd, double *zd);
  virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2,
			       double *z1, double *z2);



  void Dump(FILE * dat);
  void Dump(FILE * dat, char * head1, char * head2);
  void DumpBinary(FILE * dat);
  void DumpParameters(FILE * dat);
  int DepthSplineSmooth(double dm);
  int CheckBackgroundL1(double * res);
};

extern "C" {
int l1_driver_d(int m, int n, double * ct,  double *f, double *a,
	      double prec, double eps, int *irank, int * iter,
	      double *r, double *z, int *ind, int (*cancel_func)());
}

#endif

