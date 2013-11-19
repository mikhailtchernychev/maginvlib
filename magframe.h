//
// @(#)magframe.h  1.1  misha-12jul2012
//
// Copyright (c) 2012 of Mikhail Tchernychev
// All rights reserved.
//


#ifndef MAGFRAME_H
#define MAGFRAME_H

#include <map>
#include <vector>
#include <string>
#include "magdata.h"

typedef struct _rigid_sensor {
  double x;
  double y;
  double z;
  double F;
  double Fsynt;
  _rigid_sensor() {
    x = y = z = F =  NO_DATA;
	Fsynt = 0.;
  }
} RIGID_SENSOR;


class CMagFrame : public CMagData {

protected:
  std::map<int, RIGID_SENSOR> m_mSensorMap;
  std::vector<std::map<int, RIGID_SENSOR>::iterator> m_vPositions;

  double m_lfnX, m_lfnY, m_lfnZ; // directional cosines relatively to Earth's mag field
  int m_iIsValid;
  int m_iCurrentPos;

  double m_lfXShift;
  double m_lfYShift;
  double m_lfZShift;

public:
  CMagFrame();
  virtual ~CMagFrame();

  int ClearFrame();
  int SetPosition(int id, double x, double y, double z);
  int SetField(int id, double field);
  int SetFluxgates(double vx, double vy, double vz);
  int GetCosines(double * nx, double *ny, double *nz);
  int SetPositions(std::string s_file);

  void SetShifts(double x, double y, double z) {
	m_lfXShift = x;
    m_lfYShift = y;
    m_lfZShift = z;
  }


  virtual int get_data_origin(double *x, double *y, double *z);
  virtual int reset();
  virtual int get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				  double *add_parameters = NULL, int size_in=0, int *size_out=NULL);

  virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual int get_total_data_points();
  virtual int set_synthetic_data(double *data, int replace = 0);
  virtual int get_position(int n, double * xd, double * yd, double * zd); 
  int get_total_positions();

  virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2); 
           
  virtual void Dump(FILE * dat);
  virtual int GetMinMax(int type, double * data_min, double * data_max);		


};





#endif
