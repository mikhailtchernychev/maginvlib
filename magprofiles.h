//
// @(#)magprofiles.h  1.1  misha-08sep103
//
// Copyright (c) 2003 Geometrics
//
//
// Created:  08sep103  by  misha@misha2.geometrics.com
// Version:  08sep103  by  misha@misha2.geometrics.com
//

// sample profile implementation

#ifndef MAGPROFILES_H
#define MAGPROFILES_H

#include "magdata.h"
#include "testprofile.h"

#include <map>

class CMagProfiles: public CMagData {

protected:

  // average point for all profiles

  double m_lfXaverage;
  double m_lfYaverage;
  double m_lfZaverage;

  int m_nCurrentPoints;

  // map with profile names and profiles themselves

  map<string, CProfWrapper *> m_mapProfiles;      // holds all profile data
  map<string, CProfWrapper *>::iterator prof_pos; // to iterate through profiles.
  
  int m_iCurrentPos;
  int m_iCurrentProfilePos;  

  void destroy();


public:

  CMagProfiles();
  virtual ~CMagProfiles();

  // Function specific to profiles

  void Add(char * name, double x, double y, double time, double field, double z); 


  // now interface to general magdata class

  virtual int get_data_origin(double *x, double *y, double *z);
  virtual int reset();                         // resets data pointer to start 
  virtual int get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
			      double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

  virtual int get_total_data_points();
    
  virtual int set_synthetic_data(double *data);
  virtual int is_valid();
  virtual int get_position(int n, double *xd, double *yd, double *zd);


  // debugging functions

  void Dump(FILE * dat);


};

#endif




