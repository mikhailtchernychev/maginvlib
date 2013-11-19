//
// @(#)magsnake.h  1.1  misha-13feb10
//
// Copyright (c) 2006 Geometrics, Inc.
//
// Created:  13feb106  by  misha@misha2.geometrics.com
// Version:  13feb106  by  misha@misha2.geometrics.com
//

#ifndef MAGSNAKE_H
#define MAGSHANE_H

#include "magpipe.h"

// for CMagSnake: - sequence of pipe segments
// Non linear parameters: x0, y0, z0 - start point
// for each segment:
//                declination (angle in XOY plane) inclination (vertical angle) lenght
//                angles in radians

#define X0_POS    0
#define Y0_POS    1
#define Z0_POS    2
#define SEGMENT_START (Z0_POS+1)
#define DECL_N_POS(NSEG)  (SEGMENT_START+(NSEG)*3)
#define INCL_N_POS(NSEG)  (SEGMENT_START+(NSEG)*3 + 1)
#define LEN_N_POS(NSEG)   (SEGMENT_START+(NSEG)*3 + 2)


class CMagSnake : public CMagObject {

protected:
    CMagPipe * m_pSegment; // use to compute field
    int m_nSegments;       // number of segments
   
    int m_nCurrentSegment;
    double m_x0, m_y0, m_z0;    // start of current segment

    double * m_pLinearDeriv;

    void LoadSegment(int n);   

public:
    CMagSnake(int is_induced=0, int n_seg=2);
    ~CMagSnake();

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


