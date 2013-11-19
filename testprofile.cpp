//
// @(#)testprofile.cpp  1.1  misha-10sep103
//
// Copyright (c) 2003 Geometrics
//
// Created:  10sep103  by  misha@misha2.geometrics.com
// Version:  10sep103  by  misha@misha2.geometrics.com
//


#include "testprofile.h"
#include "l1estimator.h"

void CTestProfile::clear_all() {
    m_vX.clear();
    m_vY.clear();
    m_vZ.clear();
    m_vTime.clear();
    m_vField.clear();
    m_vComputedField.clear();
    m_nPos = 0;

    m_lfX1 = m_lfX2 = m_lfY1 = m_lfY2 = m_lfZ1 = m_lfZ2 = 0.;

    m_lfXaverage = 0.;
    m_lfYaverage = 0.;
    m_lfZaverage = 0.;

  }

void CTestProfile::InitProfile(char * name) {
   clear_all();
   SetObjectName(name);
}


void CTestProfile::Add(double x, double y, double time, double field, double z) {

  if(m_vField.size()) {
    m_lfX1 = MIN(m_lfX1, x);
    m_lfX2 = MAX(m_lfX2, x);
    m_lfY1 = MIN(m_lfY1, y);
    m_lfY2 = MAX(m_lfY2, y);
    m_lfZ1 = MIN(m_lfZ1, z);
    m_lfZ2 = MAX(m_lfZ2, z);
  }
  else {
    m_lfX1 = x;
    m_lfX2 = x;
    m_lfY1 = y;
    m_lfY2 = y;
    m_lfZ1 = z;
    m_lfZ2 = z;
  }


    m_vX.push_back(x);
    m_vY.push_back(y);
    m_vZ.push_back(z);
    m_vTime.push_back(time-GetTime0());
    m_vField.push_back(field);

    m_lfXaverage = (m_lfX1 + m_lfX2) / 2.;
    m_lfYaverage = (m_lfY1 + m_lfY2) / 2.;
    m_lfZaverage = (m_lfZ1 + m_lfZ2) / 2.;

 }


void CTestProfile::Correct() {
    int n = m_vTime.size();
    for(int i=0; i<n; i++) {
      m_vField[i] = m_vField[i] + m_vTime[i]* GetDrift() + GetDCShift(); 
    }
 }


void CTestProfile::SetZ(double z) {
    int n = m_vZ.size();
    for(int i=0; i<n; i++) {
      m_vZ[i] = z; 
    }
 }


int CTestProfile::get_data_origin(double *x, double *y, double *z) { 

    *x = m_lfXaverage;
    *y = m_lfYaverage;
    *z = m_lfZaverage;

    return 1;

}

int CTestProfile::get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				     double *add_parameters, int size_in, int *size_out){ 

  if(m_nPos >=  GetTotalPoints()) return 0;

  *x = GetX( m_nPos);
  *y = GetY( m_nPos);
  *z = GetZ( m_nPos);
  *T = GetField(m_nPos);
  *n_pos = m_nPos;
  GetAllPositionParams(add_parameters, size_in);
  
  m_nPos++;
 
  return 1;

}

int CTestProfile::get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
				  double *add_parameters, int size_in, int *size_out) {

  if(n_pos <0 || n_pos >= GetTotalPoints()) return 0;

  *x = GetX(n_pos);
  *y = GetY(n_pos);
  *z = GetZ(n_pos);
  *T = GetField(n_pos);

  GetAllPositionParams(add_parameters, size_in);

  return 1;

}

int CTestProfile::set_synthetic_data(double *data, int replace) {

  int n = m_vField.size();
  m_vComputedField.clear();
  for(int i=0; i<n; i++) m_vComputedField.push_back(data[i]);

  if(replace) {
     m_vField.clear();
     for(int i=0; i<n; i++) m_vField.push_back(data[i]);    
  }

  return 1;
} 


int CTestProfile::get_position(int n, double *xd, double *yd, double *zd) { 

  double T;
  return get_at_position(n, xd, yd, zd, &T); 

}


int CTestProfile::GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2) { 

  *x1 = m_lfX1;
  *x2 = m_lfX2;
  *y1 = m_lfY1;
  *y2 = m_lfY2;
  *z1 = m_lfZ1;
  *z2 = m_lfZ2;
 
  return 1;

} 

double CTestProfile::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
				  int data_pos){

    if(data_pos < 0) return 0.;

    return m_lfJ[0] + m_lfJ[1]*GetTimeOrFiducial(data_pos);


}

 int CTestProfile::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				   double * pos_deriv,
				   double x_obs, double y_obs, double z_obs,
	                           double x0,    double y0,    double z0,
				   double A,     double B,     double C,
				   int data_pos, int * n_params, int * start) {

     int j;
     *n_params = 0;
     j=0;
     if(data_pos < 0) return 0;

     if(start) *start =  m_iLinearStart;

     if(!m_pJfixed[0]) {
	 deriv[j] = 1.;
	 j++;
	 (*n_params)++;
     }

    if(!m_pJfixed[1]) {
	 deriv[j] = GetTimeOrFiducial(data_pos);
	 j++;
	 (*n_params)++;
     }

    return 1;
 }



void CTestProfile::Dump(FILE * dat) {

  fprintf(dat,"X   Y    Z   T  Line Tsynt\n"); 

  int n_synt = GetSyntheticPoints();

  double x_off[3];

  x_off[0] = m_lfPos[X_OFFSET];
  x_off[1] = m_lfPos[Y_OFFSET];
  x_off[2] = m_lfPos[Z_OFFSET];

  RecomputeIntoGlobal(x_off);

  for(int i=0; i<GetTotalPoints(); i++) { 
      fprintf(dat,"%lf  %lf  %lf  %lf %s",  GetX(i) + x_off[0], 
	      GetY(i) + x_off[1],
	      GetZ(i) + x_off[2],
	      GetField(i), GetName());  
      if(n_synt) fprintf(dat," %lf\n", GetSyntField(i));
      else fprintf(dat,"\n");
    }

}


void CTestProfile::DumpParameters(FILE * dat){

  assert(dat);  
  CMagObject::DumpParameters(dat);
   
  double rad=M_PI/180.; 
  double az = 90. - atan2(m_Nx[1], m_Nx[0])/rad;
    
  fprintf(dat,"Azimuth=%lg\n", az);
  fprintf(dat,"N vector %lg %lg %lg\n", m_Nx[0], m_Nx[1], m_Nx[2]);

  fprintf(dat,"Linear profile background: B+A*t (or fiducial)\n");

 if(m_pJfixed[0])
	 fprintf(dat," B=%lg FIXED\n", m_lfJ[0]);
  else 
	 fprintf(dat," B=%lg +/- %lg\n", m_lfJ[0], m_lfStdDevJ[0]);

  if(m_pJfixed[1])
	 fprintf(dat," A=%lg FIXED\n", m_lfJ[1]);
  else 
	 fprintf(dat," A=%lg +/- %lg\n", m_lfJ[1], m_lfStdDevJ[1]);

    
  if(m_pPosfixed[X_OFFSET])
    fprintf(dat,"X-shift = %lg FIXED\n", m_lfPos[X_OFFSET]); 
  else
    fprintf(dat,"X-shift = %lg +/- %lg\n", m_lfPos[X_OFFSET], m_lfStdDevPos[X_OFFSET]); 

  if(m_pPosfixed[Y_OFFSET])
    fprintf(dat,"Y-shift = %lg FIXED\n", m_lfPos[Y_OFFSET]); 
  else
    fprintf(dat,"Y-shift = %lg +/- %lg\n", m_lfPos[Y_OFFSET], m_lfStdDevPos[Y_OFFSET]); 


  if(m_pPosfixed[Z_OFFSET])
    fprintf(dat,"Z-shift = %lg FIXED\n", m_lfPos[Z_OFFSET]);
  else
    fprintf(dat,"Z-shift = %lg +/- %lg\n", m_lfPos[Z_OFFSET], m_lfStdDevPos[Z_OFFSET]); 

}


/*void CTestProfile::DumpParameters(FILE * dat){

  assert(dat);  
  CMagObject::DumpParameters(dat);

  double x_off[3], std_off[3];

  x_off[0] = m_lfPos[X_OFFSET];
  x_off[1] = m_lfPos[Y_OFFSET];
  x_off[2] = m_lfPos[Z_OFFSET];

  std_off[0] = m_lfStdDevPos[X_OFFSET];
  std_off[1] = m_lfStdDevPos[Y_OFFSET];
  std_off[2] = m_lfStdDevPos[Z_OFFSET];

  RecomputeIntoGlobal(x_off);
  RecomputeIntoGlobal(std_off);


  if(m_pPosfixed[X_OFFSET])
    fprintf(dat,"X-shift = %lg FIXED\n", x_off[0]); 
  else
    fprintf(dat,"X-shift = %lg +/- %lg\n", x_off[0],std_off[0]); 

  if(m_pPosfixed[Y_OFFSET])
    fprintf(dat,"Y-shift = %lg FIXED\n",  x_off[1]); 
  else
    fprintf(dat,"Y-shift = %lg +/- %lg\n",  x_off[1], std_off[1]); 


  if(m_pPosfixed[Z_OFFSET])
    fprintf(dat,"Z-shift = %lg FIXED\n", x_off[2]);
  else
    fprintf(dat,"Z-shift = %lg +/- %lg\n", x_off[2], std_off[2]); 

}*/

double CTestProfile::GetTimeOrFiducial(int k) {

   int n = GetTotalPoints();    
   if(m_vTime.size()==0 || m_vTime.at(0) == m_vTime.at(n-1)) return double(k);
   
   return m_vTime.at(k); 
}

int CTestProfile::ComputeLocalCoordinates() {
    
    int ret = 0;
    int i;
    char msg[120];
    double coeff_x[2], coeff_y[2];
    
    int n = GetTotalPoints();
       
       
    if(n<2) return ret;

    double * a_matrix = new double[2*n];
    double * y        = new double[n];

    double xx1, xx2, yy1, yy2, nx, ny, s;
 
    for(i=0; i<n; i++) {
             a_matrix[2*i]   = GetTimeOrFiducial(i);
             a_matrix[2*i+1] = 1;
             y[i]            = GetX(i);
    }
    
    ret = Estimate_L1(n, 2, a_matrix, y, coeff_x, NULL, NULL, NULL, NULL, NULL, msg);
    if(!ret) goto m_ret;
    
    
    for(i=0; i<n; i++) {
             a_matrix[2*i]   = GetTimeOrFiducial(i);
             a_matrix[2*i+1] = 1;
             y[i]            = GetY(i);
    }
    
    ret = Estimate_L1(n, 2, a_matrix, y, coeff_y, NULL, NULL, NULL, NULL, NULL, msg);
    if(!ret) goto m_ret;  
    
    xx1 = GetTimeOrFiducial(0)   * coeff_x[0] + coeff_x[1];
    xx2 = GetTimeOrFiducial(n-1) * coeff_x[0] + coeff_x[1];
    
    yy1 = GetTimeOrFiducial(0)   * coeff_y[0]   + coeff_y[1];
    yy2 = GetTimeOrFiducial(n-1) * coeff_y[0]   + coeff_y[1];
      
    nx = (xx2 - xx1);
    ny = (yy2 - yy1);
    
    s = sqrt(nx*nx + ny*ny);
    
    if(s==0.) { ret = 0; goto m_ret; }
    
    nx /= s; ny /= s;
      
     m_Nx[0] = nx;   m_Nx[1] =  ny;  m_Nx[2] = 0.;  // local X axis direction
     m_Ny[0] = ny;   m_Ny[1] = -nx;  m_Ny[2] = 0.;  // local Y axis direction
     m_Nz[0] = 0.;   m_Nz[1] =  0.;  m_Nz[2] = 1.;  // local Z axis direction
    
    // degug output
      
     /*for(i=0; i<n; i++) {
            fprintf(stdout,"%lf %lf %lf %lf\n",  GetX(i), GetY(i),
            GetTimeOrFiducial(i)*coeff_x[0]   + coeff_x[1], 
            GetTimeOrFiducial(i)*coeff_y[0]   + coeff_y[1]);
     }
    

     xx1 = GetTimeOrFiducial(n/2)*coeff_x[0]   + coeff_x[1];
     yy1 = GetTimeOrFiducial(n/2)*coeff_y[0]   + coeff_y[1];
     s = 0.1;
    
     fprintf(stdout,"\n%lf %lf\n", xx1, yy1);
     fprintf(stdout,"%lf %lf\n", xx1+s*m_Nx[0], yy1+s*m_Nx[1]);   
  
     fprintf(stdout,"\n%lf %lf\n", xx1, yy1);
     fprintf(stdout,"%lf %lf\n", xx1+s*m_Ny[0], yy1+s*m_Ny[1]);   */
     
 m_ret:
       if(a_matrix) { delete [] a_matrix; }
       if(y)        { delete [] y; }  
       return ret;   
}



int CTestProfile::GetMinMax(int type, double * data_min, double * data_max) {

	vector<double> * f;

	if(!type) f = &m_vField;
	     else f = &m_vComputedField;

	 for(int i=0; i<int(f->size()); i++) {
		 if(!i) {
			 *data_min = f->at(i);
			 *data_max = f->at(i);
			 continue;
		 }

		 *data_max = MAX(f->at(i), *data_max);
		 *data_min = MIN(f->at(i), *data_min);

	 }

	return 1;
}


int CTestProfile::CanComputeField() {

   // scale and bias are fixed and 0 - field is 0
   return !(m_pJfixed[0] && m_lfJ[0]==0. && m_pJfixed[1] && m_lfJ[1]==0.);

}
