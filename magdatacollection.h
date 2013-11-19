//
// @(#)magdatacollection.h  1.1  misha-10sep103
//
// Copyright (c) 2003 Geometrics, Inc
//
// Created:  10sep103  by  misha@misha2.geometrics.com
// Version:  10sep103  by  misha@misha2.geometrics.com
//


#ifndef MAGDATACOLECTION_H
#define MAGDATACOLECTION_H

#include "magdata.h"

#include <vector>
#include <algorithm>

class CMagDataCollection {

protected:

  // bounding 3-D box

  double m_lfX1, m_lfX2, m_lfY1, m_lfY2, m_lfZ1, m_lfZ2;  

  // average point

  double m_lfXaverage;
  double m_lfYaverage;
  double m_lfZaverage;

  std::vector<CMagData *> m_vData;   

  int m_nTotalPoints;  
  int m_iCurrentPos;
  int m_iCurrentTotalData;

  int m_iCurrentData;

  void destroy();
  void init();

public:

  CMagDataCollection();
  virtual ~CMagDataCollection();

  int AddData(CMagData * data); 
  void ReScan();

  void ResetDataPos() { m_iCurrentData = 0; };
  CMagData * GetNextData(); 

  // now interface to general magdata class

  virtual int get_data_origin(double *x, double *y, double *z);
  virtual int reset();                         // resets data pointer to start 
  virtual CMagData * get_next_position(double *x, double *y, double *z, double * T, int *n_pos, int * data_pos,
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual CMagData * get_at_position(int n_pos, double *x, double *y, double *z, double * T, int * data_pos, 
			      double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual int get_total_data_points() { return  m_nTotalPoints; }
    
  virtual int set_synthetic_data(double *data, int replace=0);
  virtual int is_valid();
  virtual CMagData *  get_position(int n, double *xd, double *yd, double *zd);

  virtual int GetMinMax(int type, double * data_min, double * data_max);


  // debugging functions

  void Dump(FILE * dat);

  // removing elements 

  int ReplaceData(CMagData *p_old, CMagData *p_new);
  int RemoveData(CMagData * pData);

  virtual void DumpToXMLStream(std::ostream & str);
  int GetDataPoints(std::vector<int> & points);

  int GetNearestDataPoint(double x0, double y0, double * d, double *x, double *y, double *time,
		 	              std::string & line_name);
};

#endif




