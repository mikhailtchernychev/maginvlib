//
// @(#)testprofile.h  1.1  Owner-18jul103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  18jul103  by  misha@MISHA
// Version:  18jul103  by  misha@MISHA
//

#ifndef TESTPROFILE_H
#define TESTPROFILE_H

#include "profwrapper.h"
#include "magdata.h"

#include <vector>
#include <string>
#include <map>

using namespace std;


class CTestProfile: public CProfWrapper, public CMagData {

protected:
  
  vector<double> m_vX;
  vector<double> m_vY;
  vector<double> m_vZ;
  vector<double> m_vTime;
  vector<double> m_vField;
  vector<double> m_vComputedField;
 
  int m_nPos;
  double m_lfX1, m_lfX2, m_lfY1, m_lfY2, m_lfZ1, m_lfZ2;

  double m_lfXaverage;
  double m_lfYaverage;
  double m_lfZaverage;

  
  void clear_all();

public:
  
  CTestProfile() : CProfWrapper(){
    m_iObjectLinearAbilities    =  1;
    m_iObjectNonLinearAbilities =  1;
    clear_all();
  }

  ~CTestProfile() { clear_all(); }


  // sample implementation

  void InitProfile(char * name);

  void Add(double x, double y, double time, double field, double z = 0.);
  virtual void Correct();


  const char * GetName() { return m_sObjectName; }

  // real implementation for cross-analysis

  virtual int GetTotalPoints()   { return m_vTime.size(); }
  virtual double GetX(int n)     { return m_vX[n]; }
  virtual double GetY(int n)     { return m_vY[n];}
  virtual double GetZ(int n)     { return m_vZ[n];}
  virtual double GetTime(int n)  { return m_vTime[n];}
  virtual double GetAbsoluteTime(int n)  { return m_vTime[n]+GetTime0();}
  virtual double GetField(int n) { return m_vField[n];}
  virtual double GetSyntField(int n) { return m_vComputedField[n];}

  void SetZ(double z);
  
  int GetSyntheticPoints() { return m_vComputedField.size(); }
 
 // implementation for magdata collection

  virtual int get_data_origin(double *x, double *y, double *z);
  virtual int reset(){ m_nPos = 0; return 1;};  

                       // resets data pointer to start 
  virtual int get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 


  virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
			      double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual int get_total_data_points() { return GetTotalPoints(); }
    
  virtual int set_synthetic_data(double *data, int replace=0); 
  virtual int is_valid() {return (GetTotalPoints() > 0); }
  virtual int get_position(int n, double *xd, double *yd, double *zd);
  virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2);



  void Dump(FILE * dat);
  void DumpParameters(FILE * dat);
  
  // local coordinate system
  
  virtual int ComputeLocalCoordinates();
  double GetTimeOrFiducial(int n);


  // computing linear background

 // profile has ability to compute its linear background as b + a*t,
 // where b at pos 0 and a at pos 2
 // t is time of fiducial 

  virtual double ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
			              int data_pos);

  virtual int GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				   double * pos_deriv,
				   double x_obs, double y_obs, double z_obs,
	               double x0,    double y0,    double z0,
				   double A,     double B,     double C,
			       int data_pos, int * n_params, int * start = NULL);

  virtual int GetMinMax(int type, double * data_min, double * data_max);

  virtual int CanComputeField();

};

#endif







