//
// @(#)magchannel.h  1.1  misha-15nov104
//
// Copyright (c) 2004 of Geometrics.
// All rights reserved.
//
// Created:  15nov104  by  misha@MISHA
// Version:  15nov104  by  misha@MISHA
//


#ifndef MAGCHANNEL_H
#define MAGCHANNEL_H

#include "testprofile.h"
#include "magdata.h"
#include "magarray.h"

#include <vector>
#include <string>
#include <map>

using namespace std;

class CMagChannel: public CTestProfile {
        
private:
  CMagProfileArray * m_pMagArray;
  int m_nChannel;
  int m_nPos;
  int m_nTotalChannels;

public:

  CMagChannel(); 
  ~CMagChannel(){};

 // sample implementation

  void SetChannel(int channel, CMagProfileArray * array);
  
  // real implementation for cross-analysis

  virtual int GetTotalPoints();
  virtual double GetX(int n);
  virtual double GetY(int n);
  virtual double GetZ(int n);
  virtual double GetField(int n);
  virtual double GetSyntField(int n);
  
  int GetSyntheticPoints();
 
 // implementation for magdata collection

  virtual int get_data_origin(double *x, double *y, double *z);
  virtual int reset();

                       // resets data pointer to start 
  virtual int get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 


  virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
			      double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual int get_total_data_points();
    
  virtual int set_synthetic_data(double *data, int replace=0); 
  virtual int is_valid() {return (GetTotalPoints() > 0); }
  virtual int get_position(int n, double *xd, double *yd, double *zd);
  virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2);

  virtual int GetMinMax(int type, double * data_min, double * data_max);

  //virtual void DumpParameters(FILE * dat);

};




#endif
