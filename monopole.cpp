//
// @(#)monopole.cpp  1.1  misha1-09oct108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  09oct108  by  misha
// Version:  09oct108  by  misha
//


#include "monopole.h"
#include "math.h"


/* -----------------------  computional functions -----------------------------*/



#define R(x,y,z)  sqrt((x)*(x)+(y)*(y)+(z)*(z))


//         Center of monopole is origin  of local coordinat system 



double Vm_z(double x, double y, double z) {
	double r;
	r=sqrt(x*x+y*y+z*z);
	return -z/(r*r*r);
}


double Vm_x(double x, double y, double z) {
	double r;
	r=sqrt(x*x+y*y+z*z);
	return -x/(r*r*r);
}


double T_mono(double y,  double x,  double z,
	      double y0, double x0, double z0,
	      double A,  double B,  double C) {

 double X, Y, Z, deltaT;

 x=x-x0; y=y-y0; z=z-z0;

 X = Vm_x(x,y,z);
 Y = Vm_x(y,x,z);
 Z = Vm_z(x,y,z);

 deltaT=(X*A+Y*B+Z*C);
 return(deltaT);

}


int dTmono_d(double y,  double x,  double z,  
         double y0, double x0, double z0, double M,
         double A,  double B,  double C,
         double *dTdx, double *dTdy, double *dTdz) {

	double Xx, Xy, Xz, Yx, Yy, Yz, Zx, Zy, Zz;
	x=x-x0; y=y-y0; z=z-z0;

	double r = sqrt(x*x+y*y+z*z);
	double r3 = r*r*r;
	double r5 = r3*r*r;
	double r1 = 1/r3;

	Xx = -(3*x*x/r5 - r1);  Xy = -(3*x*y)/r5;       Xz = -(3*x*z)/r5;
	Yx = Xy;                Yy = -(3*y*y/r5 - r1);  Yz = -(3*y*z)/r5;
	Zx = Xz;                Zy = Yz;                Zz = -(3*z*z/r5 - r1);


	*dTdx = M*(Xx*A+ Yx*B+ Zx*C);
	*dTdy = M*(Xy*A+ Yy*B+ Zy*C);
	*dTdz = M*(Xz*A+ Yz*B+ Zz*C);

	return 1;
}



CMagMonopole::CMagMonopole() : CMagObject(0, 1,3, "CMagMonopole", "monopole") { 

  m_iObjectLinearAbilities    = 2; 
  m_iObjectNonLinearAbilities = 2;

}



CMagMonopole::~CMagMonopole() {
                                                    
}





double CMagMonopole::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
					double x0,    double y0,    double z0,
				        double A,     double B,     double C,
				        int data_pos) {
    
    if(!m_bIsValid) return 0.;
    if(m_iFieldType != anyfield  && m_iFieldType != type) return 0;

    x_obs += add_params[X_OFFSET];
    y_obs += add_params[Y_OFFSET];
    z_obs += add_params[Z_OFFSET];

    switch(type) {

	case totalmagfield:
	    return ComputeTotalField(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

	    /*	case totalmaggradient:   
	    return ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

        case totalfinitemaggrad:
	    return ComputeFiniteGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

        case absmaggradient:	
    	    return ComputeAbsGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);
    	    
    	 case absfinitegradient:
            return ComputeAbsFiniteGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);   

	    */

	default:
	return CMagObject::ComputeField(type, add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, data_pos);

    }

    return 0.;
}




double CMagMonopole::ComputeTotalField(double * add_params, double x_obs, double y_obs, double z_obs,
					double x0,    double y0,    double z0,
				        double A,     double B,     double C) {


  return m_lfJ[0]*T_mono(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C);
}





int CMagMonopole::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
				     double * pos_deriv,
				     double x_obs, double y_obs, double z_obs,
				     double x0,    double y0,    double z0,
				     double A,     double B,     double C,
				     int data_pos, int *n_params, int * start) { 

    if(!m_bIsValid) return 0;
 
    if(m_iFieldType != anyfield && m_iFieldType != type) {
	deriv[0] = 0.;
	return 0;
    }

    x_obs += add_params[X_OFFSET];
    y_obs += add_params[Y_OFFSET];
    z_obs += add_params[Z_OFFSET];
    
    switch(type) {

	case totalmagfield:
	    return GetTotalFieldLinearDerivatives(add_params, deriv, size, pos_deriv, 
						  x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

	    /*	case totalmaggradient:   
	    return GetGradientLinearDerivatives(add_params, deriv, size, pos_deriv, 
						x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

        case totalfinitemaggrad:
	    return GetFiniteGradientLinearDerivatives(add_params, deriv, size, pos_deriv, 
						      x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

	case absmaggradient:
	    return GetAbsGradientLinearDerivatives(add_params, deriv, size, pos_deriv, 
						   x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);
						   
	case absfinitegradient: 
         return GetAbsFiniteGradientLinearDerivatives(add_params, deriv, size, pos_deriv, 
						   x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);					   
	    */

	default:
	    return CMagObject::GetLinearDerivatives(type, add_params, deriv, size, pos_deriv, 
						    x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, 
						    data_pos, n_params, start);
    }

    return 0;

}


int CMagMonopole::GetTotalFieldLinearDerivatives(double * add_params, double * deriv, int size, 
						 double * pos_deriv,
						 double x_obs, double y_obs, double z_obs,
						 double x0,    double y0,    double z0,
						 double A,     double B,     double C,
						 int *n_params, int * start) { 
				   


 if(start) *start =  m_iLinearStart;
 *n_params = 0;

if(size < 1) return 0;
if(m_pJfixed[0]) return 0;

 deriv[0] =   T_mono(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C);
 *n_params  = 1;
 return 1;

}




int CMagMonopole::GetNonLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
					double * pos_deriv,
				        double x_obs, double y_obs, double z_obs,
                                        double x0,    double y0,    double z0,
				        double A,     double B,     double C,
					int data_pos, int * n_params, int * start) {


    if(!m_bIsValid) return 0;
 
    if(m_iFieldType != anyfield && m_iFieldType != type) {
	deriv[0] = deriv[1] = deriv[2] = 0.;
	return 0;
    }

    x_obs += add_params[X_OFFSET];
    y_obs += add_params[Y_OFFSET];
    z_obs += add_params[Z_OFFSET];

    switch(type) {

     case totalmagfield:
	    return GetTotalFieldNonLinearDerivatives(add_params, deriv, size, pos_deriv, 
						     x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);
	/*
	case totalmaggradient:   
	    return GetGradientNonLinearDerivatives(add_params, deriv, size, pos_deriv, 
						   x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);


       case totalfinitemaggrad:
	    return GetFiniteGradientNonLinearDerivatives(add_params, deriv, size, pos_deriv, 
						      x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

	case absmaggradient:
	    return GetAbsGradientNonLinearDerivatives(add_params, deriv, size, pos_deriv, 
						   x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);


    case absfinitegradient:
         return GetAbsFiniteGradientNonLinearDerivatives(add_params, deriv, size, pos_deriv, 
						   x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);
      */

	default:
	    return CMagObject::GetNonLinearDerivatives(type, add_params, deriv, size, pos_deriv, 
						       x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, 
						       data_pos, n_params, start);
    }

    return 0;

}


int CMagMonopole::GetTotalFieldNonLinearDerivatives(double * add_params, double * deriv, int size, 
													double * pos_deriv,
													double x_obs, double y_obs, double z_obs,
													double x0,    double y0,    double z0,
													double A,     double B,     double C,
													int * n_params, int * start) {
	
	// non-linear parameters NOTE that d/dx and d/dy swapped!
	
	int i, j;
	
	double deriv_tmp[3];
	if(!m_bIsValid || size < 3) return 0;
	
	if(start) *start =  m_iNonLinearStart;
	
	*n_params = 0;
	
	dTmono_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, m_lfJ[0], A, B, C,
			&deriv_tmp[1], &deriv_tmp[0], &deriv_tmp[2]);

	
	for(i=0; i<3; i++) pos_deriv[i] = deriv_tmp[i];
	j=0;
	for(i=0; i<3; i++) {
		if(!m_pPosfixed[i]) { 
			deriv[j] = deriv_tmp[i]; 
			j++;
			(*n_params)++;
		}
	}
    
	return 1;
	
}					





void CMagMonopole::DumpParameters(FILE * dat){
	
	CMagObject::DumpParameters(dat);
	fprintf(dat,"Position:\n");
	
	if(m_pPosfixed[0]) 
		fprintf(dat," X: %.3lf FIXED\n", m_lfPos[0]);
	else 
		fprintf(dat," X: %.3lf +/- %lg\n", m_lfPos[0], m_lfStdDevPos[0]);
	
	if(m_pPosfixed[1]) 
		fprintf(dat," Y: %.3lf FIXED\n", m_lfPos[1]);
	else 
		fprintf(dat," Y: %.3lf +/- %lg\n", m_lfPos[1], m_lfStdDevPos[1]);
	
	if(m_pPosfixed[2]) 
		fprintf(dat," Z: %.3lf FIXED\n", m_lfPos[2]);
	else 
		fprintf(dat," Z: %.3lf +/- %lg\n", m_lfPos[2], m_lfStdDevPos[2]);
	
	
	fprintf(dat,"Monopole Ampitude:\n");
    if(m_pJfixed[0]) 
		fprintf(dat," M: %lg FIXED\n",  m_lfJ[0]);
    else
		fprintf(dat," M: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);
	
}



void CMagMonopole::DumpToXMLStream(std::ostream & str)  {

   char delimit = ' ';
   char quot = '\"';
   str << "<" << m_sObjectType << ">" << std::endl;

   str << "<params ";

   str << "X="     << quot << m_lfPos[0]       << quot << delimit;
   str << "Y="     << quot << m_lfPos[1]       << quot << delimit;
   str << "Z="     << quot << m_lfPos[2]       << quot << delimit;
   str << "X_std=" << quot << m_lfStdDevPos[0] << quot << delimit;
   str << "Y_std=" << quot << m_lfStdDevPos[1] << quot << delimit;
   str << "Z_std=" << quot << m_lfStdDevPos[2] << quot << delimit;

   str << "Jtotal="     << quot << m_lfJ[0]       << quot << delimit;
   str << "Jtotal_std=" << quot << m_lfStdDevJ[0] << quot << delimit;

   str << "/>" << std::endl;
   str << "</" << m_sObjectType << ">" << std::endl;

}