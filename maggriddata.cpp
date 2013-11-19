//
// @(#)maggriddata.cpp  1.1  misha-13may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  13may103  by  misha@misha.local
// Version:  13may103  by  misha@misha.local
//


#include "maggriddata.h"

CMagGridData::CMagGridData() {

    is_valid = 0;     
    nx = ny = 0;      
    x1 = y1 = z1 = 0.;
    dx = dy = 0.;

    T        = NULL;       
    T_synt   = NULL;  
    is_clean = 0;
    n_current = i_current = j_current = 0;
}

CMagGridData::~CMagGridData() {
    if(is_clean) {
	if(T)      { free(T);           T = NULL; }
	if(T_synt) { free(T_synt); T_synt = NULL; }
    }
}


CMagGridData::CMagGridData(int in_nx, int in_ny, double x_in, double y_in, 
			   double z_in, double dx_in, double dy_in, double *T_in, 
			   int in_clean) {

    // default constructor part

    is_valid = 0;     
    nx = ny = 0;      
    x1 = y1 = z1 = 0.;
    dx = dy = 0.;
    n_current = i_current = j_current = 0;

    T        = NULL;       
    T_synt   = NULL;  
    is_clean = 0;


    if(!T_in) return; 

    nx = in_nx;
    ny = in_ny;
    x1 = x_in;
    y1 = y_in;
    z1 = z_in;
    dx = dx_in;
    dy = dy_in;
    T = T_in;
    is_clean = in_clean;
    is_valid = 1;
}


// create dummy or load grid - use magpick dummy file format

CMagGridData::CMagGridData(std::string s_param_file) {

	int is_dummy = 0;
	double x_min, x_max;
	int i = 0;

	is_valid = 0;     
    nx = ny = 0;      
    x1 = y1 = z1 = 0.;
    dx = dy = 0.;

    T        = NULL;       
    T_synt   = NULL;  
    is_clean = 0;
    n_current = i_current = j_current = 0;


	char buffer[1024];
	FILE * dat = fopen(s_param_file.c_str(),"rt");

	if(!dat) return;

	buffer[0] = 0;
	fgets(buffer, 1022, dat);

	if(strncmp(buffer,"DUMMY", 5) ==0) is_dummy = 1;
	else if(strncmp(buffer,"DSAA", 4) ==0) is_dummy = 0;
	else goto m_ret;

	if(!fgets(buffer, 1022, dat)) goto m_ret;
	sscanf(buffer,"%d%d", &nx, &ny);

	if(!fgets(buffer, 1022, dat)) goto m_ret;
	sscanf(buffer,"%lf%lf", &x_min, &x_max);

	x1 = x_min;	dx = (x_max-x_min)/(nx-1);

	if(!fgets(buffer, 1022, dat)) goto m_ret;
	sscanf(buffer,"%lf%lf", &x_min, &x_max);

	y1 = x_min;	dy = (x_max-x_min)/(ny-1);

	if(!fgets(buffer, 1022, dat)) goto m_ret;
	sscanf(buffer,"%lf%lf", &x_min, &x_max);

	T = (double *)malloc(sizeof(double)*nx*ny);

	if(!T) goto m_ret;

	T_synt = (double *)malloc(sizeof(double)*nx*ny);

	if(!T_synt) goto m_ret;

	for(i =0; i<nx*ny; i++) T[i] = T_synt[i] = 0;
	
	if(!is_dummy) {
		for(i =0; i<nx*ny; i++) fscanf(dat,"%lf", &T[i]);
	}

	is_valid = 1;
	is_clean = 1;

	reset();

m_ret:

	fclose(dat);

}



int CMagGridData::get_data_origin(double *x, double *y, double *z) {

    *x = x1;
    *y = y1;
    *z = z1;

    return 1;
}

int CMagGridData::get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				   double *add_parameters, int size_in, int *size_out) {


    if(!is_valid) return 0;
    if(n_current >= nx*ny) return 0;

    *xd = x1 + i_current*dx;
    *yd = y1 + j_current*dy;
    *zd = z1;
    *Td = T[n_current];
    *n_pos = n_current;	

    GetAllPositionParams(add_parameters, size_in);

    //if(add_parameters && size_in >=1 && size_out) {
    //	add_parameters[0] = T_synt[n_current];
    //	*size_out = 1;
    //}

    n_current++;
    
    i_current++;
    if(i_current >= nx) {
	i_current = 0;
	j_current++;
    }


    return 1;
}

int CMagGridData::get_at_position(int n_pos, double *xd, double *yd, double *zd, double * Td, 
				 double *add_parameters, int size_in, int *size_out) {


  assert(n_pos < nx*ny && n_pos >= 0);

  // figure out what is i, j based on *ldfjack

    int j = n_pos/nx;
    int k = n_pos - j*nx;

    *xd = x1 + k*dx;
    *yd = y1 + j*dy;
    *zd = z1;
    *Td = T[n_pos];

    GetAllPositionParams(add_parameters, size_in);

    //if(add_parameters && size_in >=1 && size_out) {
    //	add_parameters[0] = T_synt[n_pos];
    //	*size_out = 1;
    //}

   return 1;

}



int CMagGridData::set_synthetic_data(double *data, int replace) {

    if(T_synt) free(T_synt);
    int n = nx*ny;

    T_synt = (double *) malloc(sizeof(double)*n);
    if(!T_synt) return 0;
    
    if(replace) {
        if(T) free(T);
        T = (double *) malloc(sizeof(double)*n);
        if(!T) {
            free(T_synt);
            T_synt = NULL;
            return 0;
        }
    }        

    for(int i=0; i<n; i++) {
        T_synt[i] = data[i];
        if(replace) T[i] = data[i];
        }

    return 1;

}


int CMagGridData::get_position(int n, double *xd, double *yd, double *zd) {

    if( n < 0 || n>= nx*ny) return 0;

   // figure out what is i, j based on *ldfjack

    int j = n/nx;
    int k = n  - j*nx;

    *xd = x1 + k*dx;
    *yd = y1 + j*dy;
    *zd = z1;

    return 1;

}

int CMagGridData::GetBoundingBox(double *x1b, double *x2b, double *y1b, double *y2b, 
				 double *z1b, double *z2b) {


  *z1b = *z2b = z1;
  *x1b = x1;  *x2b = x1 + (nx-1)*dx;
  *y1b = y1;  *y2b = y1 + (ny-1)*dy;

  return 1;

}


int CMagGridData::GetMinMax(int type, double * data_min, double * data_max) {

    double dmax, dmin;
    
    long total = nx*ny;
        
    double *pData = T;    
    
    if(type==1) pData = T_synt;

    dmax = dmin = pData[0];
    for(long i=0; i<total; i++) {
        dmax = MAX(dmax, pData[i]);
        dmin = MIN(dmin, pData[i]);
    }    

    *data_min = dmin;
    *data_max = dmax;

    return 1;
}


    
void CMagGridData::Dump(FILE * dat) {

    double dmin, dmax;
    
    GetMinMax(1, &dmin, &dmax);
    fprintf(dat,"DSAA\n%d %d\n%lf %lf\n%lf %lf\n%lg %lg\n",
        nx, ny, x1, x1 + (nx-1)*dx, y1, y1 + (ny-1)*dy, dmin, dmax);
 
    int k = 0;
    for(int i=0; i<ny; i++) {           
        for(int j=0; j<nx; j++) fprintf(dat,"%lg ", T_synt[i*nx+j]);        
        fprintf(dat,"\n");
        }


}
