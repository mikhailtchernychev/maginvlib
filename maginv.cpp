//
// @(#)maginv.cpp  1.1  misha-14may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  14may103  by  misha@misha.local
// Version:  14may103  by  misha@misha.local
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "maginv.h"
#include "dnls1estimator.h"


CMagInv::CMagInv () {

    // direction part
    I=90.; D=0.;
    Az = 90;
    SetEarthCoeff();

    // linear part
    //eps=0.; 
    //epsm=0.;
    //dm=1.e+10;
    //am=1.e-10;
    //linear_fit = -1.;

    // non-linear part
    //n_iter = 0;
    //iopt = 1;
    //maxfev = 10000;
    //fit = -1.;
    //max_diff = -1.;
    //info = -1;
    //std_err = NULL;

    m_lfLinearDiscrepancy     = -1.;
    m_lfLinearMaxDiff        = -1.;
    m_lfNonLinearDiscrepancy = -1.;
    m_lfNonLinearMaxDiff     = -1.;

    m_pData       = NULL;
    m_pCollection = NULL;

    m_ProgressFunc = NULL;

	m_sInvID="";

	m_vPropNames.clear();
	m_sPropVals.clear();

	m_iTheoryStdDevFlag = 0;

}


CMagInv::~CMagInv() {
}

void CMagInv::SetEarthCoeff() {

    double rad=M_PI/180.;
    D_rad= (90.-Az)+D;
    I_rad = I*rad; D_rad*=rad;

    A=cos(I_rad)*cos(D_rad);
    B=cos(I_rad)*sin(D_rad);
    C=sin(I_rad);	

}

void CMagInv::SetEarthCosines(double nx, double ny, double nz) {

	// warning: nx and ny changed!
    A = ny;
    B = nx;
	C = nz;

}


void CMagInv::GetEarthCoeff(double *Aout, double *Bout, double *Cout, 
                            double *Azout) {

    *Aout  = A;
    *Bout  = B;    	
    *Cout  = C;
    *Azout = Az;
}
  
void CMagInv::SetEarthFieldDirections(double Incl, double Decl, double Azimuth) {

     I  = Incl;
     D  = Decl;
     Az = Azimuth;

     SetEarthCoeff();

}

int CMagInv::GetTotalParameters(int *n_data, int *n_linear_vars, int *n_nonline_vars) {

    if(!m_pData || !m_pCollection) return 0;

    // count all data

    *n_data         = m_pData->get_total_data_points();
    *n_linear_vars  = m_pCollection->GetTotalLinearParams();
    *n_nonline_vars = m_pCollection->GetTotalNonLinearParams();

    return 1;

}

// get non-linear and linear  inversion matrix 

int CMagInv::GetWholeInversionMatrixRow(double * a, int n_pos, int size) {

  int n, m, nlin, m_total;

  if(!GetTotalParameters(&n, &m, &nlin)) return 0;
  m_total = nlin + m;
  if(m_total > size) return 0;
  for(int i=0; i<m_total; i++) a[i] = 0.;

  if(!GetGeneralInversionMatrixRow(n_pos, a, &CMagObject::GetLinearDerivatives,  0)) return 0;
  if(!GetGeneralInversionMatrixRow(n_pos, a, &CMagObject::GetNonLinearDerivatives, 0)) return 0;

  //FILE * out = fopen("A_row", "a+");
  //for(int j=0; j<m_total; j++) fprintf(out,"%10lf ", a[j]);
  //fprintf(out,"\n"); fclose(out);

  return 1;

 }


// get non-linear and linear  inversion matrix 

int CMagInv::GetWholeInversionMatrix(double * a, int packing) {

  int n, m, nlin;

  if(!GetTotalParameters(&n, &m, &nlin)) return 0;

  int n_all = n*(m+nlin);
  for(int i=0; i<n_all; i++) a[i] = 0.;
 

  if(!GetGeneralInversionMatrix(a, &CMagObject::GetLinearDerivatives, m+nlin, 0,  packing)) return 0;
  if(!GetGeneralInversionMatrix(a, &CMagObject::GetNonLinearDerivatives, m+nlin, 0,  packing)) return 0;

  //FILE * out = fopen("A_out", "wt+");
  //
  //for(int i = 0; i<n; i++) {
  //	  for(int j=0; j<m+nlin; j++) {
  //		  fprintf(out,"%10lf ", a[j*n+i]);
  //	  }
  //	  fprintf(out,"\n");
  //  }
  //
  //fclose(out);
  //
  return 1;

 }



// get nonlinear inversion matrix 

int CMagInv::GetNonLinearInversionMatrix(double * a, double *y_vect, int packing) {

  int n, m, nlin;

  if(!GetTotalParameters(&n, &m, &nlin)) return 0;

  int n_all = n*(m+nlin);
  for(int i=0; i<n_all; i++) a[i] = 0.;
  fprintf(stderr,"n=%d params=%d\n", n, m+nlin);

  if(!GetDataVector(y_vect)) return 0;

  if(!GetGeneralInversionMatrix(a, &CMagObject::GetNonLinearDerivatives, nlin, 0,  packing)) return 0;

  return 1;

 }



// get linear inversion matrix - note - it does not care about non-linear params,
// takes first block of the matrix

int CMagInv::GetLinearInversionMatrix(double * a, double *y_vect, int packing) {

  int n, m, nlin;

  if(!GetTotalParameters(&n, &m, &nlin)) return 0;

  if(!GetDataVector(y_vect)) return 0;

  int n_all = n*nlin;
  for(int i=0; i<n_all; i++) a[i] = 0.;

  if(!GetGeneralInversionMatrix(a, &CMagObject::GetLinearDerivatives, m, nlin,  packing)) return 0;

  return 1;

 }




// get general matrix for inversion

// a          - output matrix
// df         - function to compute derivatives for CMagObject
// m_vars     - total number of columns in the calculated matrix (including columns are not calculated here);
// lin_offset - offset to move columns from their original locations in the matrix


int CMagInv::GetGeneralInversionMatrix(double * a, GetDerivativesFunc df, int m_vars, 
                                       int lin_offset, int packing) {


    int n_linear, n_non_linear, n_total, n, data_pos;
    double * deriv = NULL;
  
    double add_params[MAX_FIELD_PARAMS];
    int size_in = MAX_FIELD_PARAMS;
    int size_out, i,j,j1;

    double pos_deriv[3], pos_deriv_tmp[3];
    for(i=0; i<3; i++) pos_deriv[i] = pos_deriv_tmp[i] = 0.;

    assert(m_pData);
    assert(m_pCollection);

    if(!GetTotalParameters(&n, &n_linear, &n_non_linear)) return 0;

    n_total = n_linear + n_non_linear; 

    deriv =  new double[n_total];
    if(!deriv) return 0;
    
    int nonlin_start, lin_start, n_pos, n_start, n_params;;
    double x, y, z, T;
    double x1, y1, z1;

    int m_total_objects = m_pCollection->GetTotalObjects();

    m_pData->reset();
    m_pData->get_data_origin(&x1, &y1, &z1);


    for(j=0; j<n; j++) {    // j loop over data points
    
     CMagData * pData = m_pData->get_next_position(&x, &y, &z, &T, &n_pos, &data_pos, add_params, size_in, &size_out);
     x -= x1; y -= y1; z -= z1;


     CMagObject * pDataObject = NULL;
     for(j1=0; j1<3; j1++) pos_deriv[j1] = 0.;

      m_pCollection->ResetPosition();

     for(i=0; i<m_total_objects; i++) {   // i loop over objects

       CMagObject * p = m_pCollection->GetCurrentMagObject(&nonlin_start, &lin_start);
       if(p == pData && p) pDataObject = p;

       if(!p || !p->CanComputeField()) continue;

       int ret = (p->*df)(pData->GetFieldType(), add_params, deriv, n_total, pos_deriv_tmp,
		x, y, z, x1, y1, z1, A, B, C, data_pos, &n_params, &n_start);

       if(!ret) continue;

           
         n_start -= lin_offset;

	     if(packing) 
		 for(int k=0; k<n_params; k++) a[n*n_start+k*n+n_pos] = deriv[k];  // columns packed
	     else
		 for(int k=0; k<n_params; k++) a[n_pos*m_vars+n_start+k] = deriv[k]; // rows packed
    
	     for(j1=0; j1<3; j1++) {
	       pos_deriv[j1]    += pos_deriv_tmp[j1];
	       pos_deriv_tmp[j1] = 0.;
	     }

     }  // end loop over objects

     // now treat data as object (for possible position corretion)


    if(pDataObject) {
               int ret = (pDataObject->*df)(pData->GetFieldType(), add_params, deriv, n_total, pos_deriv,  
		                  x, y, z, x1, y1, z1, A, B, C, data_pos, &n_params, &n_start);

	       if(!ret) continue;

	       n_start -= lin_offset;

	       if(packing) 
		   for(int k=0; k<n_params; k++) a[n*n_start+k*n+n_pos]    = deriv[k];  // columns packed
	       else
		   for(int k=0; k<n_params; k++) a[n_pos*m_vars+n_start+k] = deriv[k]; // rows packed
      }

    }  // end loop over data ponts 

  delete [] deriv;

  return 1;
}


int CMagInv::GetGeneralInversionMatrixRow(int n_pos, double * a, GetDerivativesFunc df, int lin_offset) {


    int n_linear, n_non_linear, n_total, n, data_pos;

    double add_params[MAX_FIELD_PARAMS];
    double pos_deriv[3], pos_deriv_tmp[3];

    int size_in = MAX_FIELD_PARAMS;
    int size_out, i, j, k;

    for(i=0; i<3; i++) pos_deriv[i] = pos_deriv_tmp[i] = 0.;

    assert(m_pData);
    assert(m_pCollection);

    if(!GetTotalParameters(&n, &n_linear, &n_non_linear)) return 0;

    n_total = n_linear + n_non_linear;
 
    double * deriv =  new double[n_total];
    if(!deriv) return 0;

    int nonlin_start, lin_start, n_start, n_params;;
    double x, y, z, T;
    double x1, y1, z1;

    m_pData->get_data_origin(&x1, &y1, &z1);
    m_pCollection->ResetPosition();

    int m_total_objects = m_pCollection->GetTotalObjects();

    CMagData *pData = m_pData->get_at_position(n_pos, &x, &y, &z, &T, &data_pos, add_params, size_in, &size_out);
    
    x -= x1; y -= y1; z -= z1;

    CMagObject * pDataObject = NULL;

    for(i=0; i<m_total_objects; i++) {   // i loop over objects

	CMagObject * p = m_pCollection->GetCurrentMagObject(&nonlin_start, &lin_start);
	
	if(p == pData) pDataObject = p;

        if(!p || !p->CanComputeField()) continue;

	    int ret = (p->*df)(pData->GetFieldType(), add_params, deriv, n_total, pos_deriv_tmp,  
		              x, y, z, x1, y1, z1, A, B, C, data_pos, &n_params, &n_start);

	    if(!ret) continue;

	    for(k=0; k<n_params; k++) {
                       a[n_start+k] = deriv[k];
                    } // columns packed

	    for(j=0; j<3; j++) {
	      pos_deriv[j]    += pos_deriv_tmp[j];
	      pos_deriv_tmp[j] = 0.;
	    }


   }  // end loop over objects

    // now treat data as object (for possible position corretion)

    if(pDataObject) {
      int ret = (pDataObject->*df)(pData->GetFieldType(), add_params, deriv, n_total, pos_deriv,  
		 x, y, z, x1, y1, z1, A, B, C, data_pos, &n_params, &n_start);

      for( k=0; k<n_params && ret; k++)  a[n_start+k] = deriv[k]; 
    }

   delete [] deriv;

  return 1;
}



int CMagInv::GetDataVector(double *y_vect) {

    int n_linear, n_non_linear, n, data_pos;

    assert(m_pData);
    assert(m_pCollection);

    if(!GetTotalParameters(&n, &n_linear, &n_non_linear)) return 0;

    m_pData->reset();
    double x, y, z, T;
    int n_pos;

    for(int j=0; j<n; j++) { 
	    m_pData->get_next_position(&x, &y, &z, &T, &n_pos, &data_pos);
	    y_vect[n_pos] = T;
	}  

    return 1;
}



int CMagInv::GetLinearVariables(double * vars, int size, int linear_only, int is_errors) {

   if(!m_pCollection) return 0;
   return  m_pCollection->GetSetVariables(0, vars, size, 0, linear_only, is_errors); 
	
}



int CMagInv::SetLinearVariables(double * vars, int size, int linear_only, int is_errors) {

   if(!m_pCollection) return 0;
   return m_pCollection->GetSetVariables(1, vars, size, 0, linear_only, is_errors); 

}



int CMagInv::GetNonLinearVariables(double * vars, int size, int is_errors) {

         if(!m_pCollection) return 0;
         return  m_pCollection->GetSetVariables(0, vars, size, 1, 0, is_errors); 
}



int CMagInv::SetNonLinearVariables(double * vars, int size, int is_errors) {

   if(!m_pCollection) return 0;
   return  m_pCollection->GetSetVariables(1, vars, size, 1, 0, is_errors); 

}


int CMagInv::GetAllVariables(double * vars, int size, int is_errors) {

   int n, m, nlin, m_total;
   if(!GetTotalParameters(&n, &m, &nlin)) return 0;

   m_total = m + nlin;

   if(m_total > size) return 0; 

    GetNonLinearVariables(vars, size, is_errors);
    GetLinearVariables(vars,size, 0, is_errors);

    return 1;

}

int CMagInv::SetAllVariables(double * vars, int size, int is_errors) {

   int n, m, nlin, m_total;
   if(!GetTotalParameters(&n, &m, &nlin)) return 0;

   m_total = m + nlin;

   if(m_total > size) return 0; 

    SetNonLinearVariables(vars, size, is_errors);
    SetLinearVariables(vars, size, 0, is_errors);
   
    return 1;

}


int CMagInv::EstimateLinearParameters(LINEAR_ESTIMATOR_FUNC func, double * fit, int packing, void * params,
		                     char * message, int store_results) {

	
    int n, m, nlin;
    double d = 0.;
    
    int ret = 0;

    double * a_matrix, * y_vector, * x_init, * T_synt, * st_err;	
    a_matrix = y_vector = x_init = st_err =  T_synt = st_err = NULL;	
    
    CMagDataCollection * pData = GetMagData();
    if(!pData) return 0; 
    if(!m_pCollection) return 0;
   
     m_pData->ReScan();

     m_pCollection->ReScanAllObjects(&m, &nlin);
    if(m<=0) return 0;


    if(!GetTotalParameters(&n, &m, &nlin)) return 0;

    a_matrix = new double[(m+nlin)*n]; // using whole matrix even for linear 
    if(!a_matrix) goto m_end;
    
    y_vector = new double[n];
    if(!y_vector) goto m_end;

    T_synt = new double[n];
    if(!T_synt) goto m_end;
    				
    x_init = new double [m];
    if(!x_init) goto m_end; 

    st_err = new double [m];
    if(!st_err) goto m_end;

    if(!GetLinearInversionMatrix(a_matrix, y_vector, packing)) goto m_end;
    if(!GetLinearVariables(x_init, m)) goto m_end;

    ret = (*func)(n, m, a_matrix, y_vector, x_init, params, T_synt, fit, &m_lfLinearMaxDiff, st_err, message);
    if(!ret) goto m_end;
    
    if(store_results) {
        m_lfLinearDiscrepancy = *fit;

        SetLinearVariables(x_init, m); 		
        SetLinearVariables(st_err, m, TRUE, 1); 		

        m_pData->set_synthetic_data(T_synt);
    }

m_end:
	if(a_matrix)       delete [] a_matrix;
	if(y_vector)       delete [] y_vector;
	if(T_synt)         delete [] T_synt;
	if(x_init)         delete [] x_init;
	if(st_err)         delete [] st_err;

	return ret;
}


// computes magnetic field for all structure.
// T_synt    - array for new magnetic field
// size      - size of array
// type == 0 - magnetic field
//         1 - difference observed - syntetic

int CMagInv::ComputeField(double *T_synt, int size, int type) {

   int i, j, n,m, nlin,  m_total_objects, lin_start, nonlin_start, n_pos, data_pos;
   double x, y, z, x1, y1, z1, T;

   double add_params[MAX_FIELD_PARAMS];
   int size_in = MAX_FIELD_PARAMS;
   int size_out;

   if(!GetTotalParameters(&n, &m, &nlin)) return 0;

   if(n > size) return 0;

    for(i=0; i<n; i++) T_synt[i] = 0.;
 
   // now calculate synthetic field

    m_pCollection->ResetPosition();
    m_pData->get_data_origin(&x1, &y1, &z1);

    m_total_objects = m_pCollection->GetTotalObjects();

    for(i=0; i<m_total_objects; i++) {   // i loop over objects


	CMagObject * p = m_pCollection->GetCurrentMagObject(&nonlin_start, &lin_start);
        if(!p || !p->CanComputeField()) continue;

	m_pData->reset();

	for(j=0; j<n; j++) {    // j loop over data points
	    CMagData * pData = m_pData->get_next_position(&x, &y, &z, &T, &n_pos, &data_pos, 
							  add_params, size_in, &size_out);
            x -= x1; y -= y1; z -= z1;
	    if(pData != p) data_pos = -1;
            T_synt[n_pos] += p->ComputeField(pData->GetFieldType(), add_params, x, y, z, x1, y1, z1, A, B, C, 
					     data_pos);

 	}  // end loop over data ponts 
   }  // end loop over objects

    if(type) { // compute difference
            m_pData->reset();
            for(int j=0; j<n; j++) {
              m_pData->get_next_position(&x, &y, &z, &T, &n_pos, &data_pos);
              T_synt[n_pos] = T-T_synt[n_pos];
            }
	}


return 1;
}



// computes magnetic field for all structure at position n_pos
// T_synt    - computed magnetic field
// type == 0 - magnetic field
//         1 - difference observed - syntetic

double CMagInv::ComputeFieldAtPosition(int n_pos, int type) {

   int n,m, nlin,  m_total_objects, lin_start, nonlin_start, data_pos;
   double x, y, z, x1, y1, z1, T;

   double add_params[MAX_FIELD_PARAMS];
   int size_in = MAX_FIELD_PARAMS;
   int size_out;

   double T_synt = 0.;

   if(!GetTotalParameters(&n, &m, &nlin)) return 0;

   // now calculate synthetic field

    m_pCollection->ResetPosition();
    m_pData->get_data_origin(&x1, &y1, &z1);
    m_total_objects = m_pCollection->GetTotalObjects();

    CMagData * pData = m_pData->get_at_position(n_pos, &x, &y, &z, &T, &data_pos, add_params, size_in, &size_out);
    if(!pData) return 0.;
    
    x -= x1; y -= y1; z -= z1;

    for(int i=0; i<m_total_objects; i++) {   // i loop over objects


	CMagObject * p = m_pCollection->GetCurrentMagObject(&nonlin_start, &lin_start);
        if(!p || !p->CanComputeField()) continue;

	if(p != pData) data_pos = -1;
         T_synt += p->ComputeField(pData->GetFieldType(), add_params, x, y, z, x1, y1, z1, A, B, C, data_pos);

   }  // end loop over objects


    if(type) T_synt = T-T_synt;

return T_synt;
}




int CMagInv::EstimateAllParameters(NONLINEAR_ESTIMATOR_FUNC func, 
    void * func_param, char * message, int store_data) {


  int n, m, nlin;
  double * x_vec        = NULL;
  double * x_up         = NULL; // upper limits
  double * x_low        = NULL; // lower limits
  double * standard_err = NULL;
  double * T_synt = NULL;

  int ret = 0;

  int inversion_size = sizeof(CInversion);

  CMagDataCollection * pData = GetMagData();
  if(!pData) return 0; 
  if(!m_pCollection) return 0;
  
  m_pData->ReScan();
  
  m_pCollection->ReScanAllObjects(&m, &nlin);
  if(m+nlin <=0) return 0;

  if(!GetTotalParameters(&n, &m, &nlin)) return 0;

   int m_total = m + nlin;

   x_vec = new double[m_total];
   if(!x_vec) goto m_end;

   x_up = new double[m_total];
   if(!x_up) goto m_end;

   x_low = new double[m_total];
   if(!x_low) goto m_end;


   standard_err = new double[m_total];
   if(!standard_err) goto m_end;

   if(!GetAllVariables(x_vec, m_total)) goto m_end;

   // get limits

   if(!GetAllVariables(x_up,  m_total, 2)) goto m_end;
   if(!GetAllVariables(x_low, m_total, 3)) goto m_end;

   
   T_synt = (*func)(func_param, (void *)this, n, m_total, x_vec, x_up, x_low,  &m_lfNonLinearDiscrepancy, standard_err, 
                    &m_lfNonLinearMaxDiff, message);
                                  
   if(!T_synt) goto m_end;

   ret = 1;
   
   SetAllVariables(x_vec, m_total);
   
   if(store_data) {
      SetAllVariables(standard_err, m_total, 1);
      pData->set_synthetic_data(T_synt);
      }
 

m_end:
   if(x_vec)         { delete [] x_vec;           x_vec = NULL; }
   if(x_up)          { delete [] x_up;            x_up  = NULL; }
   if(x_low)         { delete [] x_low;           x_low = NULL; }
   if(standard_err)  { delete [] standard_err;    standard_err = NULL; }
   if(T_synt)        { delete [] T_synt;          T_synt = NULL; }

  return ret;
}


int CMagInv::EstimateNonLinearParametersRange(LINEAR_ESTIMATOR_FUNC func, double fit, double step, double max_step,
                                         int packing, void * params){
                                             
                                             
    char message[1024];                                        
    int n, m, nlin;
    double d = 0., fit_out, max_fit_out;
    
    int ret = 0, i, j;

    double * x_params, *x_params_save,  *dx_min, *dx_max, *T_save, *T_synt;	
    x_params  = x_params_save = dx_min = dx_max =  T_save = T_synt = NULL;	
    
    CMagDataCollection * pData = GetMagData();
    if(!pData) return 0; 
    if(!m_pCollection) return 0;
 
    m_pData->ReScan();
  
    if(!GetTotalParameters(&n, &m, &nlin)) return 0;

    T_save = new double[n];
    if(!T_save) goto m_end;
 
    T_synt = new double[n];
    if(!T_synt) goto m_end;
    
    x_params = new double[nlin];
    if(!x_params) goto m_end;
  
    x_params_save    = new double[nlin];
    if(!x_params_save) goto m_end;
       				
    dx_min = new double [nlin];
    if(!dx_min) goto m_end; 

    dx_max = new double [nlin];
    if(!dx_max) goto m_end; 
 
    for(i=0; i<nlin; i++) dx_min[i] = dx_max[i] = 0.;                                                 

    if(!GetNonLinearVariables(x_params, nlin)) goto m_end;
                        
    for(i=0; i<nlin; i++) x_params_save[i] = x_params[i];                                                                                           
    
    ComputeField(T_save, n);
    
    // main loop                                                                                                                      
    for(i=0; i<nlin; i++) {
        
      SetNonLinearVariables(x_params_save, nlin);
      GetNonLinearVariables(x_params, nlin);
      double dx = 0.;
      double x = x_params[i];
      j=0;
      while(dx < max_step) { // increase parameter loop
         dx = j*step;   
         x_params[i] = x+dx;
         SetNonLinearVariables(x_params, nlin);
         EstimateLinearParameters(func, &fit_out, packing, params, message, 0);
         ComputeField(T_synt, n);
         compute_fit(n, T_save, T_synt, &fit_out, &max_fit_out);
         fprintf(stderr,"fit_out = %lf %lf\n", fit_out,x_params[i] );
         if(fit_out > fit) break;
         dx_max[i] = MAX(dx_max[i], dx); 
         j++;
      }     
      
      SetNonLinearVariables(x_params_save, nlin);
      GetNonLinearVariables(x_params, nlin);
      dx = 0.;
      x = x_params[i];
      j=0;
      while(dx < max_step) { // decrease parameter loop
         dx = j*step;   
         x_params[i] = x-dx;
         SetNonLinearVariables(x_params, nlin);
         EstimateLinearParameters(func, &fit_out, packing, params, message, 0);
         ComputeField(T_synt, n);
         compute_fit(n, T_save, T_synt, &fit_out, &max_fit_out);
         fprintf(stderr,"fit_out = %lf\n", fit_out);       
         if(fit_out > fit) break;
         dx_min[i] = MAX(dx_min[i], dx); 
         j++;
      }               
                                                              
    } // end main loop                                                                                                                    

    // print the results
    
    for(i = 0; i<nlin; i++)
        fprintf(stderr,"Value: %lf min %lf max %lf\n", x_params_save[i], dx_min[i], dx_max[i]);                                            
                                                                                                                                                                                                                                                                                                                                                                                                                                          
    ret = 1;
 
    SetNonLinearVariables(x_params_save, nlin);
    
m_end:
 
    if(x_params)      delete [] x_params;
    if(x_params_save) delete [] x_params_save;
    if(dx_min)        delete [] dx_min;
    if(dx_max)        delete [] dx_max;
    if(T_save)        delete [] T_save;
    if(T_synt)        delete [] T_synt;
        
return ret;                                             
}    


int CMagInv::EstimateNonLinearParametersRange2(NONLINEAR_ESTIMATOR_FUNC func, 
          double fit, double step, double max_step, void * params){
                                             
                                             
    char message[1024];                                        
    int n, m, nlin, m_total;
    double d = 0., fit_out, max_fit_out;
    
    int ret = 0, i, n_var, n_objects;

    double * x_params,  *dx_min, *dx_max, *T_save, *T_synt;	
    x_params  = dx_min = dx_max =  T_save = T_synt = NULL;	
    
    CMagDataCollection * pData = GetMagData();
    if(!pData) return 0; 
    if(!m_pCollection) return 0;
  
    m_pData->ReScan();
 
    if(!GetTotalParameters(&n, &m, &nlin)) return 0;
    m_total = m + nlin;
     
    m_pCollection->ReScanAllObjects(&m, &nlin);
    if(m+nlin <=0) return 0;
   
    T_save = new double[n];
    if(!T_save) goto m_end;
 
    T_synt = new double[n];
    if(!T_synt) goto m_end;
    
    x_params = new double[m_total];
    if(!x_params) goto m_end;
       				
    dx_min = new double [nlin];
    if(!dx_min) goto m_end; 

    dx_max = new double [nlin];
    if(!dx_max) goto m_end; 
 
    for(i=0; i<nlin; i++) dx_min[i] = dx_max[i] = 0.;                                                 
                                                                                   
    if(!GetAllVariables(x_params, m_total)) goto m_end;
    
    ComputeField(T_save, n);
    
    // object loop loop                                                                                                                      
    
    n_objects = m_pCollection->GetTotalObjects();
    n_var = 0;

    double x_vars[10];

    m_pCollection->InitMinMaxForObjects();
      
    for(i=0; i<n_objects; i++) {

        CMagObject * obj = m_pCollection->GetObjectAt(i);
        if(!obj) continue;
        
        int n_params = obj->GetAllNonLinearParams();

        for(int j=0; j<n_params; j++) {
            
               if(obj->IsNonLinearParamFixed(j)) continue;
                 
                double x_up, x_low;
                x_up =  obj->GetNonLinearParam(j, 2);
                x_low = obj->GetNonLinearParam(j, 3);
                               
                SetAllVariables(x_params, m_total);
                      
                double dx = 0.;
                double x = obj->GetNonLinearParam(j);
                           
                int j1 = 0;
                while(dx < max_step) { // increase parameter loop
                 dx = j1*step;   
                 SetAllVariables(x_params, m_total);
                 if(x+dx < x_low || x+dx > x_up) break;
                 obj->FixNonLinearParam(j, TRUE, x+dx);
                 if(!EstimateAllParameters(func, params, message, 0)) break;
                 ComputeField(T_synt, n);
                 compute_fit(n, T_save, T_synt, &fit_out, &max_fit_out);
                 obj->FixNonLinearParam(j, FALSE, x+dx);
                 m_pCollection->ReScanAllObjects(&m, &nlin);
                 GetAllVariables(x_vars, m_total);                              
                 if(max_fit_out > fit) break;
                 //fprintf(stdout,"%lf %lf %lf %lf\n", x_vars[0], x_vars[1], x_vars[2], fit_out); 
                 m_pCollection->UpdateMinMaxForObjects();
                 dx_max[n_var] = MAX(dx_max[n_var], dx); 
                 j1++; 
                 }
               
               j1 = 0;
               dx = 0.;
               while(dx < max_step) { // decrease parameter loop
                 dx = j1*step;   
                 SetAllVariables(x_params, m_total);
                 if(x-dx < x_low || x-dx > x_up) break;
                 obj->FixNonLinearParam(j, TRUE, x-dx);
                 if(!EstimateAllParameters(func, params, message, 0)) break;
                 ComputeField(T_synt, n);
                 compute_fit(n, T_save, T_synt, &fit_out, &max_fit_out);          
                 obj->FixNonLinearParam(j, FALSE, x-dx);
                 m_pCollection->ReScanAllObjects(&m, &nlin);
                 GetAllVariables(x_vars, m_total);
                 if(max_fit_out > fit) break;
                 m_pCollection->UpdateMinMaxForObjects();
                 //fprintf(stdout,"%lf %lf %lf %lf\n", x_vars[0], x_vars[1], x_vars[2], fit_out);
                 dx_min[n_var] = MAX(dx_min[n_var], dx); 
                 j1++;
               }
               
              n_var++;              
           }                
                                                                      
    } // end main loop                                                                                                                    

    // print the results
    
    for(i = 0; i<nlin; i++)
        fprintf(stderr,"Value: %lf min %lf max %lf\n", x_params[i], dx_min[i], dx_max[i]);                                            
                                                         
    m_pCollection->DumpMinMaxForObjects(stderr);
                                                                                                                                                                                                                                                                                                                                                                                     
    ret = 1;
    
    SetAllVariables(x_params, m_total);
    m_pCollection->ReScanAllObjects(&m, &nlin);
     
m_end:
 
    if(x_params)      delete [] x_params;
    if(dx_min)        delete [] dx_min;
    if(dx_max)        delete [] dx_max;
    if(T_save)        delete [] T_save;
    if(T_synt)        delete [] T_synt;
        
return ret;                                             
}    


void CMagInv::DumpObjects(FILE * dat) {

  assert(dat);
  if(m_pCollection) m_pCollection->DumpObjects(dat); 
}



void CMagInv::DumpToStream(std::ostream & str, std::string comment, char * close_str, char * delimit,
						   std::string offset) {

   str << comment.c_str()  << "Inversion" << std::endl;
   str << comment.c_str()  << "INCL="             << I                        << delimit 
	 	          << "DECL="             << D                        << delimit 
                  << "LINEARDISCR="      << m_lfLinearDiscrepancy    << delimit
				  << "LINEARMAXDIFF="    << m_lfLinearMaxDiff        << delimit 
				  << "NONLINEARDISCR="   << m_lfNonLinearDiscrepancy << delimit 
				  << "NONLINEARMAXDIFF=" << m_lfNonLinearMaxDiff     << std::endl;

   if(m_pCollection) {
	   std::string comment1 = comment;
	   comment1.insert(0, offset.c_str());
		m_pCollection->DumpToStream(str, comment1, close_str, delimit, offset);
   }

   str << comment.c_str()  << close_str << "Inversion" << std::endl;
}



void CMagInv::DumpToXMLStream(std::ostream & str) {

   char delimit =' ';
   char quot = '\"';
   str << "<Inversion>" << std::endl;
   str << "<params " << "ID="               << quot << m_sInvID.c_str()         << quot << delimit
	                 << "INCL="             << quot << I                        << quot << delimit 
	 	             << "DECL="             << quot << D                        << quot << delimit 
                     << "LINEARDISCR="      << quot << m_lfLinearDiscrepancy    << quot << delimit
				     << "LINEARMAXDIFF="    << quot << m_lfLinearMaxDiff        << quot << delimit 
				     << "NONLINEARDISCR="   << quot << m_lfNonLinearDiscrepancy << quot << delimit 
				     << "NONLINEARMAXDIFF=" << quot << m_lfNonLinearMaxDiff     << quot << "/>" << std::endl;


   str << "<properties ";

   for(int i=0; i<int(m_vPropNames.size()); i++) {
	   str <<  m_vPropNames[i] << "=" << quot << m_sPropVals[i] << quot;
	   if(i==m_vPropNames.size()-1) str << "/>" << std::endl;
	   else str << delimit;
   }


   if(m_pCollection) m_pCollection->DumpToXMLStream(str);
   if(m_pData)       m_pData->DumpToXMLStream(str);

   str << "</Inversion>" << std::endl;
}


int CMagInv::GetTotalAbilities() {

 if(!m_pCollection) return 0;

 int lin, non_lin;
 m_pCollection->GetAbilities(&lin, &non_lin);

 return (lin < non_lin) ? lin : non_lin;

}

int CMagInv::ComputeField(int type, int replace) {
	
    double * T_synt = NULL;
    
    CMagDataCollection * pData = GetMagData();
    if(!pData) return 0; 
    if(!m_pCollection) return 0;

    int m, n, nlin;

     m_pCollection->ReScanAllObjects(&m, &nlin);

    if(!GetTotalParameters(&n, &m, &nlin)) return 0;

    T_synt = new double[n];
    if(!T_synt) return 0;

    ComputeField(T_synt, n, type);

    m_pData->set_synthetic_data(T_synt, replace);
	     
    if(T_synt)  delete [] T_synt;
    return 1;
}




void   CMagInv::AddDoubleProperty(char * name, double prop) {

	m_vPropNames.push_back(name);
	m_sPropVals.push_back(prop);

}



double  CMagInv::GetDoubleProperty(char * name) {

	int n = m_vPropNames.size();

	for(int i=0; i<n; i++) {
		if(name == m_vPropNames[i]) {
			return m_sPropVals[i];
		}

	}
		return NO_DATA;
}

int CMagInv::EstimateStandardDeviations() {

	int n, m, nlin;
	DNLS1_STORAGE dstorage;
	dstorage.iopt   =  1;
	dstorage.maxfev = 10;
	char message[1024];

	int std_flag = GetTheoryStdDevFlag();
	SetTheoryStdDevFlag(1);

	EstimateAllParameters(Estimator_dnls1, &dstorage, message);

	SetTheoryStdDevFlag(std_flag);

	double av_fit  = GetNonLinearDiscrepancy();

	if(!GetTotalParameters(&n, &m, &nlin)) return -1;

	int m_total = m + nlin;

	double * x = new double [m_total];
	if(!x) return -1;

	GetAllVariables(x, m_total, 1);

    for(int i =0; i<m_total; i++) {
		x[i] *=av_fit;
	}

	SetAllVariables(x, m_total, 1);

	delete [] x;

	return 0;

}