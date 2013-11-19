//
// @(#)gradprofile.h  1.1  misha-06oct103
//
// Copyright (c) 2003 Geometrics
//
// Created:  06oct103  by  misha@misha2.geometrics.com
// Version:  06oct103  by  misha@misha2.geometrics.com
//

#ifndef GRADPROFILE_H
#define GRADPROFILE_H

#include "testprofile.h"

#include <vector>
#include <string>
#include <map>

using namespace std;

class CGradProfile: public CTestProfile {

protected:
  
  vector<double> m_vNx;  // gradient direction
  vector<double> m_vNy;  // in every data point
  vector<double> m_vNz;

  void clear_all();

public:
  CGradProfile() : CTestProfile() {
    clear_all();
  }

  virtual ~CGradProfile() { clear_all(); }

  void Add(double x, double y, double time, double field, double z = 0., 
	   double nx = 0., double ny =0., double nz = 1.);

  int GetGradDirection(int n, double * dir); 

  virtual int get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 
  
  virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
			      double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual void Dump(FILE * dat);
  virtual void DumpParameters(FILE * dat);


};

#endif

