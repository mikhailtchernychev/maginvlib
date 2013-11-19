//
// @(#)magdata1d.h  1.1  misha1-29sep108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  29sep108  by  misha
// Version:  29sep108  by  misha
//


#include "magdata.h"

#ifndef MAGDATA1D_H
#define MAGDATA1D_H


class CMagData1D : public CMagData {

protected:
  double * m_pT;       //  observed magnetic field 
  double * m_pTsynt;   //  synthetic magnetic field (also storage array) 
  double   m_lfDx;     //  step along the profile
  int m_nPoints;       //  total data points
  int m_nCurrent;      //  current position in array 
  int m_nStart;        //  point with X=0

  void clean();

public:
    CMagData1D();
    virtual ~CMagData1D();

    CMagData1D(int n, double *x, double *y, double *z, double *T);
	CMagData1D(int n, double dx);

    virtual int get_data_origin(double *x, double *y, double *z);
    virtual int reset() {  m_nCurrent = 0; return 1;} 
    virtual int get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				  double *add_parameters = NULL, int size_in=0, int *size_out=NULL);

    virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

    virtual int get_total_data_points() { return m_nPoints; }
    virtual int set_synthetic_data(double *data, int replace = 0);
    virtual int get_position(int n, double * xd, double * yd, double * zd); 

    virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2); 
           
    virtual void Dump(FILE * dat);
    virtual int GetMinMax(int type, double * data_min, double * data_max);		         
	virtual int GetSyntheticField(double *Tsynt, int n);


	double RecalcCoord(double x);
	double GetDX() { return m_lfDx; }
	int From1Dto2D(double *x, double *y, int n, double x0, 
			       double *x_out, double *y_out);

};

#endif



