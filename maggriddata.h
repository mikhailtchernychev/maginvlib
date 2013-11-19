//
// @(#)maggriddata.h  1.1  misha-13may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  13may103  by  misha@misha.local
// Version:  13may103  by  misha@misha.local
//


#include "magdata.h"

#ifndef MAGGRIDDATA_H
#define MAGGRIDDATA_H

class CMagGridData : public CMagData {

protected:
    int is_valid;                 //  if this data is valid
    int nx, ny;                   //  nx nodes along x nx along y for grid    
    double x1, y1, z1;            //  origin of grid and elevation
    double dx, dy;                //  steps along x and y
    double *T;                    //  observed magnetic field 
    double *T_synt;               //  synthetic magnetic field (also storage array) 
    int is_clean;                 //  clean data on destructor 
    int i_current,j_current;      //  current row and column
    int n_current;                //  current position in array 

public:
    CMagGridData();
    virtual ~CMagGridData();
    CMagGridData(int in_nx, int in_ny, double x_in, double y_in, double z_in,
	     double dx_in, double dy_in, double *T_in, int in_clean=1);
	CMagGridData(std::string s_param_file);

    virtual int get_data_origin(double *x, double *y, double *z);
    virtual int reset() { n_current = i_current =j_current = 0; return is_valid;}
    virtual int get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				  double *add_parameters = NULL, int size_in=0, int *size_out=NULL);

    virtual int get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
				double *add_parameters = NULL, int size_in=0, int *size_out=NULL); 

    virtual int get_total_data_points() { return nx*ny; }
    virtual int set_synthetic_data(double *data, int replace = 0);
    virtual int get_position(int n, double * xd, double * yd, double * zd); 

    virtual int GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2); 
           
    virtual void Dump(FILE * dat);
    virtual int GetMinMax(int type, double * data_min, double * data_max);		         
};

#endif






