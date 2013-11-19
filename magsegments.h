//
// @(#)magsegments.h  1.1  misha-15feb10
//
// Copyright (c) 2006 Geometrics, Inc.
//
// Created:  15feb106  by  misha@misha2.geometrics.com
// Version:  15feb106  by  misha@misha2.geometrics.com
//

#ifndef MAGSEGMENTS_H
#define MAGSEGMENTS_H

#include "magpipe.h"

// for CMagSnake: - sequence of pipe segments. X, Y, Z at each segment

#define X_SEG_POS(N)  (3*(N))
#define Y_SEG_POS(N)  (3*(N) + 1) 
#define Z_SEG_POS(N)  (3*(N) + 2) 

class CMagSegments : public CMagObject {

protected:
    CMagPipe * m_pSegment; // use to compute field
    int m_nSegments;       // number of segments
   
    int m_nCurrentSegment;

    double * m_pLinearDeriv;

    void LoadSegment(int n);   

public:
    CMagSegments(int is_induced=0, int n_seg=2);
    ~CMagSegments();

    int GetNSegments() { return m_nSegments; }
    int GetNodes(double *pos, int size);


  // x_obs, y_obs, z_obs - position of the observation point
  // A, B, C             - Earth's magnetic field direction coefficients

  virtual double ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
			              int n_pos);


  virtual int GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
				     double * pos_deriv,
				     double x_obs, double y_obs, double z_obs,
				     double x0,    double y0,    double z0,
				     double A,     double B,     double C,
			  	     int n_pos, int *n_params, int * start);

  virtual int GetTotalLinearParams();
  virtual void DumpParameters(FILE * dat);
  virtual int CanComputeField() { return 1; }


};



#endif


