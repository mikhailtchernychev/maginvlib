//
// @(#)magdipole.cpp  1.1  misha-10may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  10may103  by  misha@misha.local
// Version:  10may103  by  misha@misha.local
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "magdipole.h"


// ------------------- start computational functions -------------------------------
 
// magnetic field & derivatives computational functions.


#define R(x,y,z)  sqrt((x)*(x)+(y)*(y)+(z)*(z))


/*            Center of sphere or rectangular is origin
              of local coordinat system 
*/

/* ------ for sphere -------------*/

double Vh(double x, double y, double z) {
	double r;
	r=sqrt(x*x+y*y+z*z);
	return(3.*x*y/(r*r*r*r*r));
	}

double Vv(double x, double y, double z) {
	double r;
	r=sqrt(x*x+y*y+z*z);
	return( (2.*x*x-y*y-z*z)/(r*r*r*r*r));
	}


double dT_2(double y,  double x,  double z,
           double y0, double x0, double z0,
           double A,  double B,  double C,
           double Jx, double Jy, double Jz) {

double X, Y, Z, Vxy, Vxz, Vyz, deltaT;

x=x-x0; y=y-y0; z=z-z0;

Vxy=Vh(x,y,z); Vxz=Vh(x,z,y); Vyz=Vh(y,z,x); 

X = Jx*Vv(x,y,z) + Jy*Vxy       + Jz*Vxz;
Y = Jx*Vxy       + Jy*Vv(y,x,z) + Jz*Vyz;
Z = Jx*Vxz	 + Jy*Vyz       + Jz*Vv(z,y,x);

deltaT=X*A+Y*B+Z*C;
return(deltaT);
}


int dT_inv_2(double y, double x, double z, 
             double y0, double x0, double z0, 
             double A, double B, double C,
             double *a, double *b, double *c) {
double  Vxy, Vxz, Vyz;

x=x-x0; y=y-y0; z=z-z0;

Vxy=Vh(x,y,z); Vxz=Vh(x,z,y); Vyz=Vh(y,z,x); 

*a =  A*Vv(x,y,z) + B*Vxy       + C*Vxz;
*b =  A*Vxy       + B*Vv(y,x,z) + C*Vyz;
*c =  A*Vxz	  + B*Vyz       + C*Vv(z,y,x);
return 0;
}


int dT_d(double y,  double x,  double z,  
         double y0, double x0, double z0,
         double A,  double B,  double C,
         double Jx, double Jy, double Jz,
         double *dTdx, double *dTdy, double *dTdz) {

double deltaT, R7;

double dVxy_dx0, dVxz_dx0, dVyz_dx0;
double dVxx_dx0, dVyy_dx0, dVzz_dx0; 
double dVxy_dy0, dVxz_dy0, dVyz_dy0;
double dVxx_dy0, dVyy_dy0, dVzz_dy0; 
double dVxy_dz0, dVxz_dz0, dVyz_dz0;
double dVxx_dz0, dVyy_dz0, dVzz_dz0; 

x=x-x0; y=y-y0; z=z-z0;

dVxy_dx0 = 12.*x*x*y - 3.*y*y*y - 3.*y*z*z;
dVxz_dx0 = 12.*x*x*z - 3.*z*z*z - 3.*z*y*y;
dVyz_dx0 = 15.*x*y*z;
dVxx_dx0 =  6.*x*x*x - 9.*x*y*y - 9.*x*z*z;
dVyy_dx0 = 12.*y*y*x - 3.*x*x*x - 3.*x*z*z;
dVzz_dx0 = 12.*z*z*x - 3.*x*x*x - 3.*y*y*x;

dVxy_dy0 = dVyy_dx0;
dVxz_dy0 = dVyz_dx0;
dVyz_dy0 = 12.*y*y*z - 3.*z*x*x - 3.*z*z*z;
dVxx_dy0 = dVxy_dx0;
dVyy_dy0 =  6.*y*y*y - 9.*y*x*x - 9.*y*z*z;
dVzz_dy0 = 12.*z*z*y - 3.*y*x*x - 3.*y*y*y;

dVxy_dz0 = dVyz_dx0;
dVxz_dz0 = dVzz_dx0; 
dVyz_dz0 = dVzz_dy0;
dVxx_dz0 = dVxz_dx0;
dVyy_dz0 = dVyz_dy0;
dVzz_dz0 =  6.*z*z*z - 9.*z*x*x - 9.*y*y*z;


R7= R(x,y,z);
R7 = R7*R7*R7*R7*R7*R7*R7;

deltaT = dVxy_dx0*(A*Jy + B*Jx) + 
         dVxz_dx0*(A*Jz + C*Jx) + 
         dVyz_dx0*(B*Jz + C*Jy) +
         A*Jx*dVxx_dx0 + B*Jy*dVyy_dx0 + C*Jz*dVzz_dx0;


*dTdx = (double)deltaT/R7;

deltaT = dVxy_dy0*(A*Jy + B*Jx) + 
         dVxz_dy0*(A*Jz + C*Jx) + 
         dVyz_dy0*(B*Jz + C*Jy) +
         A*Jx*dVxx_dy0 + B*Jy*dVyy_dy0 + C*Jz*dVzz_dy0;

*dTdy = (double)deltaT/R7;

deltaT = dVxy_dz0*(A*Jy + B*Jx) + 
         dVxz_dz0*(A*Jz + C*Jx) + 
         dVyz_dz0*(B*Jz + C*Jy) +
         A*Jx*dVxx_dz0 + B*Jy*dVyy_dz0 + C*Jz*dVzz_dz0;

*dTdz = (double)deltaT/R7;

return 0;
}





int dT_d_dx0(double y,  double x,  double z,  
             double y0, double x0, double z0,
             double A,  double B,  double C,
             double Jx, double Jy, double Jz,
             double *dTdx_dx0, double *dTdy_dx0, double *dTdz_dx0) {

double  R1, R5, R2, R7, R9;

double dVxy_dx_dx0, dVxz_dx_dx0, dVyz_dx_dx0;
double dVxx_dx_dx0, dVyy_dx_dx0, dVzz_dx_dx0; 
double dVxy_dy_dx0, dVxz_dy_dx0, dVyz_dy_dx0;
double dVxx_dy_dx0, dVyy_dy_dx0, dVzz_dy_dx0; 
double dVxy_dz_dx0, dVxz_dz_dx0, dVyz_dz_dx0;
double dVxx_dz_dx0, dVyy_dz_dx0, dVzz_dz_dx0; 

 double x2, y2, z2, x3, y3, z3;

x=x-x0; y=y-y0; z=z-z0;

R1 = R(x,y,z);
R2 = R1*R1; 
R5 = R2*R2*R1;
R7 = R5*R2;
R9 = R7*R2;

 x2 = x*x;  x3 = x2*x;
 y2 = y*y;  y3 = y2*y;
 z2 = z*z;  z3 = z2*z;

/*1*/  dVxy_dx_dx0 = (21.*x*y*z2 + 21.*x*y3 +-84*x3*y)/R9 + 24.*x*y/R7;  //(7.*x*(3.*y*R2-15.*x2*y))/R9 + 24.*x*y/R7; 
/*2*/  dVxz_dx_dx0 = (21.*x*z3 + 21*x*y2*z - 84.*x3*z)/R9 + 24.*x*z/R7;  //(7.*x*(3.*z*R2-15.*x2*z))/R9 + 24.*x*z/R7;
/*3*/  dVyz_dx_dx0 = 15.*y*z/R7 - 105.*x2*y*z/R9;
/*4*/  dVxx_dx_dx0 = (63.*x2*z2 + 63.*x2*y2 - 42.*x2*x2)/R9 + (-9.*z2-9.*y2+18.*x2)/R7;  //(7.*x2*(4.*R2-5.*(-z2-y2+2.*x2)))/R9 - (4.*R2-5.*(-z2-y2+2.*x2) - 12.*x2)/R7;
/*11*/ dVyy_dx_dx0 = (21.*x2*z2 - 84.*x2*y2 + 21.*x2*x2)/R9 + (-3.*z2+12.*y2 -9.*x2)/R7; //(7.*x2*(-2.*R2-5.*(-z2+2.*y2-x2)))/R9- (-2.*R2-5.*(-z2+2.*y2-x2) +  6.*x2)/R7;
/*12*/ dVzz_dx_dx0 = -(dVxx_dx_dx0+dVyy_dx_dx0); //(-84.*x2*z2 + 21.*x2*y2 + 21.*x2*x2)/R9 + (12.*z2-3.*y2-9.*x2)/R7; //(7.*x2*(-2.*R2-5.*(2*z2-y2-x2)))/R9  - (-2.*R2-5.*(2.*z2-y2-x2) +  6.*x2)/R7;

/*5*/  dVxy_dy_dx0 = dVyy_dx_dx0;
/*6*/  dVxz_dy_dx0 = dVyz_dx_dx0;
/*7*/  dVyz_dy_dx0 = (21.*x*z3-84.*x*y2*z + 21.*x3*z)/R9 - 6.*x*z/R7;  //(7.*x*(3.*z*R2 - 15.*y2*z))/R9 - 6.*x*z/R7;
/*8*/  dVxx_dy_dx0 =  dVxy_dx_dx0;
/*9*/  dVyy_dy_dx0 = (63.*x*y*z2 - 42.*x*y3 + 63.*x3*y)/R9 - 18.*x*y/R7;  //(7.*x*(4.*y*R2 - 5.*y*(-z2+2*y2-x2)))/R9 - 18.*x*y/R7;
/*10*/ dVzz_dy_dx0 = -(dVyy_dy_dx0+dVxx_dy_dx0); //(-84.*x*y*z2 + 21.*x*y3 + 21*x3*y)/R9 - 6.*x*y/R7; //(7.*x*(-5.*y*(2.*z2-y2-x2) - 2.*y*R2))/R9 - 6*x*y/R7;

 dVxy_dz_dx0 = dVyz_dx_dx0;
 dVxz_dz_dx0 = dVzz_dx_dx0; 
 dVyz_dz_dx0 = dVzz_dy_dx0;
 dVxx_dz_dx0 = dVxz_dx_dx0;
 dVyy_dz_dx0 = dVyz_dy_dx0;
 /*18*/  dVzz_dz_dx0 = (-42.*x*z3 + 63.*x*y2*z + 63.*x3*z)/R9 - 18.*x*z/R7;  //(7.*x*(4.*z*R2 - 5.*z*(2.*z2-y2-x2)))/R9 - 18*x*z/R7;

/* - this part is identical in all routines */   

*dTdx_dx0 = dVxy_dx_dx0*(A*Jy + B*Jx) + 
            dVxz_dx_dx0*(A*Jz + C*Jx) + 
            dVyz_dx_dx0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dx_dx0 + B*Jy*dVyy_dx_dx0 + C*Jz*dVzz_dx_dx0;

*dTdy_dx0 = dVxy_dy_dx0*(A*Jy + B*Jx) + 
            dVxz_dy_dx0*(A*Jz + C*Jx) + 
            dVyz_dy_dx0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dy_dx0 + B*Jy*dVyy_dy_dx0 + C*Jz*dVzz_dy_dx0;

*dTdz_dx0 = dVxy_dz_dx0*(A*Jy + B*Jx) + 
            dVxz_dz_dx0*(A*Jz + C*Jx) + 
            dVyz_dz_dx0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dz_dx0 + B*Jy*dVyy_dz_dx0 + C*Jz*dVzz_dz_dx0;

// take "-" because gradient is always bottom mag - top mag.

 *dTdx_dx0 = -(*dTdx_dx0);
 *dTdy_dx0 = -(*dTdy_dx0);
 *dTdz_dx0 = -(*dTdz_dx0);


return 0;
}




int dT_d_dy0(double y,  double x,  double z,  
             double y0, double x0, double z0,
             double A,  double B,  double C,
             double Jx, double Jy, double Jz,
             double *dTdx_dy0, double *dTdy_dy0, double *dTdz_dy0) {

double R1, R5, R2, R7, R9;

double dVxy_dx_dy0, dVxz_dx_dy0, dVyz_dx_dy0;
double dVxx_dx_dy0, dVyy_dx_dy0, dVzz_dx_dy0; 
double dVxy_dy_dy0, dVxz_dy_dy0, dVyz_dy_dy0;
double dVxx_dy_dy0, dVyy_dy_dy0, dVzz_dy_dy0; 
double dVxy_dz_dy0, dVxz_dz_dy0, dVyz_dz_dy0;
double dVxx_dz_dy0, dVyy_dz_dy0, dVzz_dz_dy0; 

 double x2, y2, z2, x3, y3, z3;

x=x-x0; y=y-y0; z=z-z0;

R1 = R(x,y,z);
R2 = R1*R1; 
R5 = R2*R2*R1;
R7 = R5*R2;
R9 = R7*R2;

 x2 = x*x; x3 = x2*x;
 y2 = y*y; y3 = y2*y;
 z2 = z*z; z3 = z2*z;

/*1*/  dVxy_dx_dy0 = (21.*y2*z2+21.*y2*y2-84.*x2*y2)/R9 +  (-3.*z2-9.*y2+12.*x2)/R7;  //(7.*y2*(3.*R2-15.*x2))/R9 - (3.*R2+6.*y2-15.*x2)/R7;
/*2*/  dVxz_dx_dy0 = (21.*y*z3+21.*y3*z-84.*x2*y*z)/R9  - 6.*y*z/R7; //(7.*y*(3.*z*R2-15.*x2*z))/R9 - 6.*y*z/R7;
/*3*/  dVyz_dx_dy0 = 15.*x*z/R7 - 105.*x*y2*z/R9;
/*4*/  dVxx_dx_dy0 = (63.*x*y*z2 + 63.*x*y3 - 42.*x3*y)/R9 - 18.*x*y/R7; //(7.*y*(4.*x*R2 - 5.*x*(-z2-y2+2.*x2)))/R9 - 18.*x*y/R7;
/*11*/ dVyy_dx_dy0 = (21.*x*y*z2-84.*x*y3+21.*x3*y)/R9  + 24.*x*y/R7;   //(7.*y*(-2.*x*R2-5.*x*(-z2+2.*y2-x2)))/R9 + 24*x*y/R7;
/*12*/ dVzz_dx_dy0 =  -(dVxx_dx_dy0+dVyy_dx_dy0); //(-84.*x*y*z2+21.*x*y3+21.*x3*y)/R9 - 6.*x*y/R7; //(7.*y*(-5.*x*(2.*z2-y2-x2) - 2.*x*R2))/R9 - 6.*x*y/R7;

/*5*/  dVxy_dy_dy0 = dVyy_dx_dy0;
/*6*/  dVxz_dy_dy0 = dVyz_dx_dy0;
/*7*/  dVyz_dy_dy0 = (21.*y*z3-84.*y3*z+21.*x2*y*z)/R9 + 24.*y*z/R7; //(7.*y*(3.*z*R2-15.*y2*z))/R9 + 24.*y*z/R7;
/*8*/  dVxx_dy_dy0 = dVxy_dx_dy0;
/*9*/  dVyy_dy_dy0 = (63.*y2*z2-42.*y2*y2+63.*x2*y2)/R9 + (-9.*z2+18.*y2-9.*x2)/R7;  //(7.*y2*(4.*R2-5.*(-z2+2.*y2-x2)))/R9 - (4.*R2-5.*(-z2+2.*y2-x2)-12.*y2)/R7;
/*10*/ dVzz_dy_dy0 = -(dVxx_dy_dy0+dVyy_dy_dy0); //(-84.*y2*z2+21.*y2*y2+21.*x2*y2)/R9 + (12.*z2-9.*y2-3.*x2)/R7; //(7.*y2*(-5.*(2.*z2-y2-x2) - 2.*R2))/R9 - (-5.*(2.*z2-y2-x2) - 2.*R2 + 6.*y2)/R7;

 dVxy_dz_dy0 = dVyz_dx_dy0;
 dVxz_dz_dy0 = dVzz_dx_dy0; 
 dVyz_dz_dy0 = dVzz_dy_dy0;
 dVxx_dz_dy0 = dVxz_dx_dy0;
 dVyy_dz_dy0 = dVyz_dy_dy0;
 /*18*/  dVzz_dz_dy0 = (-42.*y*z3+63.*y3*z+63.*x2*y*z)/R9 - 18.*y*z/R7; //(7.*y*(4.*z*R2 - 5.*z*(2.*z2-y2-x2)))/R9 - 18.*y*z/R7;


/* - this part is identical in all routines */   

*dTdx_dy0 = dVxy_dx_dy0*(A*Jy + B*Jx) + 
            dVxz_dx_dy0*(A*Jz + C*Jx) + 
            dVyz_dx_dy0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dx_dy0 + B*Jy*dVyy_dx_dy0 + C*Jz*dVzz_dx_dy0;

*dTdy_dy0 = dVxy_dy_dy0*(A*Jy + B*Jx) + 
            dVxz_dy_dy0*(A*Jz + C*Jx) + 
            dVyz_dy_dy0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dy_dy0 + B*Jy*dVyy_dy_dy0 + C*Jz*dVzz_dy_dy0;

*dTdz_dy0 = dVxy_dz_dy0*(A*Jy + B*Jx) + 
            dVxz_dz_dy0*(A*Jz + C*Jx) + 
            dVyz_dz_dy0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dz_dy0 + B*Jy*dVyy_dz_dy0 + C*Jz*dVzz_dz_dy0;

// take "-" because gradient is always bottom mag - top mag.

 *dTdx_dy0 = -(*dTdx_dy0);
 *dTdy_dy0 = -(*dTdy_dy0);
 *dTdz_dy0 = -(*dTdz_dy0);


return 0;
}


int dT_d_dz0(double y,  double x,  double z,  
             double y0, double x0, double z0,
             double A,  double B,  double C,
             double Jx, double Jy, double Jz,
             double *dTdx_dz0, double *dTdy_dz0, double *dTdz_dz0) {

double R1, R5, R2, R7, R9;

double dVxy_dx_dz0, dVxz_dx_dz0, dVyz_dx_dz0;
double dVxx_dx_dz0, dVyy_dx_dz0, dVzz_dx_dz0; 
double dVxy_dy_dz0, dVxz_dy_dz0, dVyz_dy_dz0;
double dVxx_dy_dz0, dVyy_dy_dz0, dVzz_dy_dz0; 
double dVxy_dz_dz0, dVxz_dz_dz0, dVyz_dz_dz0;
double dVxx_dz_dz0, dVyy_dz_dz0, dVzz_dz_dz0; 

 double x2, y2, z2, x3, y3,z3;

x=x-x0; y=y-y0; z=z-z0;

R1 = R(x,y,z);
R2 = R1*R1; 
R5 = R2*R2*R1;
R7 = R5*R2;
R9 = R7*R2;

 x2 = x*x; x3 = x2*x;
 y2 = y*y; y3 = y2*y;
 z2 = z*z; z3 = z2*z;

/*1*/  dVxy_dx_dz0 = (21.*y*z3+21.*y3*z-84.*x2*y*z)/R9 - 6.*y*z/R7;  //(7.*z*(3.*y*R2-15.*x2*y))/R9 - 6.*y*z/R7;
/*2*/  dVxz_dx_dz0 = (21.*z2*z2+21.*y2*z2-84.*x2*z2)/R9 + (-9*z2-3.*y2+12.*x2)/R7; //(7.*z2*(3.*R2-15.*x2))/R9 - (3.*R2+6.*z2-15.*x2)/R7;
/*3*/  dVyz_dx_dz0 = 15.*x*y/R7 - 105.*x*y*z2/R9;
/*4*/  dVxx_dx_dz0 = (63.*x*z3+63.*x*y2*z-42.*x3*z)/R9 - 18.*x*z/R7; //(7.*z*(4.*x*R2-5.*x*(-z2-y2+2.*x2)))/R9 - 18.*x*z/R7;
/*11*/ dVyy_dx_dz0 = (21.*x*z3-84.*x*y2*z+21.*x3*z)/R9  - 6.*x*z/R7; //(7.*z*(-2.*x*R2 -5.*x*(-z2+2.*y2-x2)))/R9 - 6.*x*z/R7;
/*12*/ dVzz_dx_dz0 = -(dVxx_dx_dz0+dVyy_dx_dz0); //(-84.*x*z3+21.*x*y2*z+21.*x3*z)/R9 + 24.*x*z/R7;  //(7.*z*(-5.*x*(2.*z2-x2-y2) - 2.*x*R2))/R9 + 24.*x*z/R7;

/*5*/  dVxy_dy_dz0 = dVyy_dx_dz0;
/*6*/  dVxz_dy_dz0 = dVyz_dx_dz0;
/*7*/  dVyz_dy_dz0 = (21.*z2*z2-84.*y2*z2+21*x2*z2)/R9 + (-9.*z2+12.*y2-3.*x2)/R7; //(7.*z2*(3.*R2-15.*y2))/R9 - (3.*R2+6*z2-15.*y2)/R7;
/*8*/  dVxx_dy_dz0 = dVxy_dx_dz0;
/*9*/  dVyy_dy_dz0 = (63.*y*z3-42.*y3*z+63.*x2*y*z)/R9 - 18.*y*z/R7; //(7.*z*(4.*y*R2 - 5.*y*(-z2+2.*y2-x2)))/R9 - 18.*y*z/R7;
/*10*/ dVzz_dy_dz0 = -(dVxx_dy_dz0+dVyy_dy_dz0); //(-84.*y*z3+21.*y3*z+21*x2*y*z)/R9 + 24.*y*z/R7; // (7.*z*(-5.*y*(2*z2-y2-x2) - 2.*y*R2))/R9 + 24.*y*z/R7;

 dVxy_dz_dz0 = dVyz_dx_dz0;
 dVxz_dz_dz0 = dVzz_dx_dz0; 
 dVyz_dz_dz0 = dVzz_dy_dz0;
 dVxx_dz_dz0 = dVxz_dx_dz0;
 dVyy_dz_dz0 = dVyz_dy_dz0;
 /*18*/  dVzz_dz_dz0 = (-42*z2*z2+63.*y2*z2+63.*x2*z2)/R9 + (18.*z2-9.*y2-9.*x2)/R7; //(7.*z2*(4.*R2-5.*(2.*z2-y2-x2)))/R9 - (-5.*(2.*z2-y2-x2) + 4.*R2 - 12*z2)/R7;

/* - this part is identical in all routines */   

*dTdx_dz0 = dVxy_dx_dz0*(A*Jy + B*Jx) + 
            dVxz_dx_dz0*(A*Jz + C*Jx) + 
            dVyz_dx_dz0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dx_dz0 + B*Jy*dVyy_dx_dz0 + C*Jz*dVzz_dx_dz0;

*dTdy_dz0 = dVxy_dy_dz0*(A*Jy + B*Jx) + 
            dVxz_dy_dz0*(A*Jz + C*Jx) + 
            dVyz_dy_dz0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dy_dz0 + B*Jy*dVyy_dy_dz0 + C*Jz*dVzz_dy_dz0;

*dTdz_dz0 = dVxy_dz_dz0*(A*Jy + B*Jx) + 
            dVxz_dz_dz0*(A*Jz + C*Jx) + 
            dVyz_dz_dz0*(B*Jz + C*Jy) +
            A*Jx*dVxx_dz_dz0 + B*Jy*dVyy_dz_dz0 + C*Jz*dVzz_dz_dz0;

// take "-" because gradient is always bottom mag - top mag.

 *dTdx_dz0 = -(*dTdx_dz0);
 *dTdy_dz0 = -(*dTdy_dz0);
 *dTdz_dz0 = -(*dTdz_dz0);

return 0;
}





// ---------------------- end computational functions -------------------------


CMagDipole::CMagDipole() : CMagObject(0, 3,3, "CMagDipole", "dipole") { 

  m_iObjectLinearAbilities    = 2; 
  m_iObjectNonLinearAbilities = 2;

}


CMagDipole::CMagDipole(int is_induced) : CMagObject(is_induced, 3,3, "CMagDipole", "dipole") {

  m_iObjectLinearAbilities    = 2; 
  m_iObjectNonLinearAbilities = 2;

}

CMagDipole::~CMagDipole() {
                          
                                                    
}


int CMagDipole::GetTotalLinearParams() {

	if(m_nLinearParamsReport==0) return 0;

    if(m_bInduced) return 1;
    else return m_nLinearParamsReport;

}



double CMagDipole::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
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

	case totalmaggradient:   
	    return ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

        case totalfinitemaggrad:
	    return ComputeFiniteGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

        case absmaggradient:	
    	    return ComputeAbsGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);
    	    
    	 case absfinitegradient:
            return ComputeAbsFiniteGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);   

	default:
	    return CMagObject::ComputeField(type, add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, data_pos);

    }

    return 0.;
}


double CMagDipole::ComputeTotalField(double * add_params, double x_obs, double y_obs, double z_obs,
					double x0,    double y0,    double z0,
				        double A,     double B,     double C) {
	
	
  if(m_bInduced) return m_lfJ[0]*dT_2(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);

  return dT_2(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);

}

double CMagDipole::ComputeFiniteGradient(double * add_params, double x_obs, double y_obs, double z_obs,
					 double x0,    double y0,    double z0,
					 double A,     double B,     double C) {

 double s = add_params[FINITE_SEPARATION]; 
 if(s <= 0.) return 0; 
	
 double x_obs2 = x_obs + add_params[DIRECTION_X]*s;
 double y_obs2 = y_obs + add_params[DIRECTION_Y]*s;
 double z_obs2 = z_obs + add_params[DIRECTION_Z]*s;
 double T1, T2;

  if(m_bInduced){
      T1 =  m_lfJ[0]*dT_2(x_obs,  y_obs, z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
      T2 =  m_lfJ[0]*dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
  }
  else {
      T1 =  dT_2(x_obs,  y_obs, z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);
      T2 =  dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);
  }

      return (T2-T1)/s;
}

int  CMagDipole::GetAbsFiniteGradientLinearDerivatives(double * add_params, double * deriv, int size, 
			       double * pos_deriv,
			       double x_obs, double y_obs, double z_obs,
			       double x0,    double y0,    double z0,
			       double A,     double B,     double C,
			       int *n_params, int * start) {
          
    
   *n_params = 0;
    int i, j;
    if(start) *start =  m_iLinearStart;
    
    if(m_bInduced) { // for induced it's correct
      double Jtmp = m_lfJ[0];
       m_lfJ[0] = 1.;
       deriv[0] = ComputeAbsFiniteGradient(add_params, x_obs, y_obs, z_obs, x0, y0, z0, A, B, C);     
       m_lfJ[0] = Jtmp;
       *n_params = 1;
    }
    else { // this is approximation only. J included non-linearly  
      double Jorig[3], deriv_tmp[3];
      for(i=0; i<3; i++) Jorig[i] = m_lfJ[i];
    
      m_lfJ[0] = 1;  m_lfJ[1] = 0.;  m_lfJ[2] = 0.;
      deriv_tmp[0] = ComputeAbsFiniteGradient(add_params, x_obs, y_obs, z_obs, x0, y0, z0, A, B, C);
    
      m_lfJ[0] = 0;  m_lfJ[1] = 1.;  m_lfJ[2] = 0.;
      deriv_tmp[1] = ComputeAbsFiniteGradient(add_params, x_obs, y_obs, z_obs, x0, y0, z0, A, B, C);
    
      m_lfJ[0] = 0;  m_lfJ[1] = 0.;  m_lfJ[2] = 1.;
      deriv_tmp[2] = ComputeAbsFiniteGradient(add_params, x_obs, y_obs, z_obs, x0, y0, z0, A, B, C);
    
      for(i=j=0; i<3; i++)  {
          m_lfJ[i] = Jorig[i];
          if(!m_pJfixed[i]) { 
              deriv[j] = deriv_tmp[i];  
              j++;
              (*n_params)++;
             }
        }
     }
    
                                                           
return 1;                       
}	


// this is not working. It should be mixture of linear and no-linear
// derivates calls

int  CMagDipole::GetAbsFiniteGradientNonLinearDerivatives(double * add_params, double * deriv, int size, 
			       double * pos_deriv,
			       double x_obs, double y_obs, double z_obs,
			       double x0,    double y0,    double z0,
			       double A,     double B,     double C,
			       int *n_params, int * start) {
                       
       
       double deriv2[9], T0, T2[3], deriv0[3];      
       double x_obs2, y_obs2, z_obs2;
       *n_params = 0;
       int i, j;

       double s = add_params[FINITE_SEPARATION]; 
       
        // compute field and derivatives for 0 point
        if(m_bInduced)
          T0 =  m_lfJ[0]*dT_2(x_obs,  y_obs, z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
        else 
          T0 =  dT_2(x_obs,  y_obs, z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);
 
        int ret = GetTotalFieldLinearDerivatives(add_params, deriv0, 3, pos_deriv,  x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);
        if(!ret) return 0;
            
       // now compute in 3 points
            
       for(i=0; i<3; i++) {
      
          GetAbsFinalGradSecondPoint(i, s, x_obs, y_obs, z_obs, &x_obs2, &y_obs2, &z_obs2);
       
         // field
         if(m_bInduced)
           T2[i] =  m_lfJ[0]*dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
         else
           T2[i] =  dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);
           
         // derivatives
         ret = GetTotalFieldLinearDerivatives(add_params, &deriv2[3*i], 3, pos_deriv,  x_obs2, y_obs2,  z_obs2, x0, y0, z0, A, B, C, n_params, start);
         if(!ret) return 0;
                    
   }

   // gradient and final differences
   
   double gx   = (T2[0]-T0);
   double gy   = (T2[1]-T0);
   double gz   = (T2[2]-T0);   
   double grad = sqrt(gx*gx + gy*gy + gz*gz);                

   if(grad <= 0.) { // special case when field is 0 (magnetic momentum is 0).
           gx = gy = gz = grad = 1.;
   }

   double derivtmp;
   *n_params = 0; 
   for(i=j=0; i<*n_params; i++) {
    derivtmp = pos_deriv[i] = (gx*(deriv2[i]   - deriv0[i]) +
                               gy*(deriv2[3+i] - deriv0[i]) +
                               gz*(deriv2[6+i] - deriv0[i]))/(s*grad);                               
    if(!m_pPosfixed[i]) { 
         deriv[j] = derivtmp; 
         j++;
         (*n_params)++;
      }                   
    }
    
  if(start) *start =  m_iNonLinearStart;
                                                           
return 1;                       
}	



 // presently finite gradient is taken only along main axis
 // it will be a problem with different directions of the profile

void CMagDipole::GetAbsFinalGradSecondPoint(int axis, double s, double x1, double y1, double z1,
                                                double * x2, double *y2, double *z2){
 *x2 = x1; *y2 = y1; *z2 = z1;
 
  switch(axis) {
       case 0: *x2 = x1 + s;  break;
       case 1: *y2 = y1 + s;  break;                        
       case 2: *z2 = z1 + s; break;
       }                                     
                                                                            
} 
                                
                                
double CMagDipole::ComputeAbsFiniteGradient(double * add_params, double x_obs, double y_obs, double z_obs,
					 double x0,    double y0,    double z0,
					 double A,     double B,     double C) {

 double s = add_params[FINITE_SEPARATION]; 
 if(s <= 0.) return 0; 

 double x_obs2, y_obs2, z_obs2;
 double T1, T2, grad;

 // this is for 4 sensor setup - disable now

/*  if(m_bInduced)
      T1 =  m_lfJ[0]*dT_2(x_obs,  y_obs, z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
  else 
      T1 =  dT_2(x_obs,  y_obs, z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);
  
  grad  = 0.;
  for(int i=0; i<3; i++) {
          
  GetAbsFinalGradSecondPoint(i, s, x_obs, y_obs, z_obs, &x_obs2, &y_obs2, &z_obs2);
       
  if(m_bInduced)
      T2 =  m_lfJ[0]*dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
   else 
      T2 =  dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);

    grad += (T2-T1)*(T2-T1);

  }*/
  
  
  // for 6 sensor setup
  
  double s2 = s/2;
  grad = 0.;
  //fprintf(stderr,"s=%lf x_obs=%lf y_obs=%lf z_obs=%lf x=%lf y=%lf z=%lf %lf %lf %lf\n",
  //                      s, x_obs, y_obs, z_obs, x0, y0, z0, m_lfPos[0], m_lfPos[1], m_lfPos[2]); 
  
   for(int i=0; i<3; i++) {
                    
   // positive along the axis
                    
  GetAbsFinalGradSecondPoint(i, s2, x_obs, y_obs, z_obs, &x_obs2, &y_obs2, &z_obs2);
       
  if(m_bInduced)
      T2 =  m_lfJ[0]*dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
   else 
      T2 =  dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);

  // negative along the axis
  
  GetAbsFinalGradSecondPoint(i, -s2, x_obs, y_obs, z_obs, &x_obs2, &y_obs2, &z_obs2);
       
  if(m_bInduced)
      T1 =  m_lfJ[0]*dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
   else 
      T1 =  dT_2(x_obs2,y_obs2,z_obs2, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, m_lfJ[0], m_lfJ[1], m_lfJ[2]);


    grad += (T2-T1)*(T2-T1);

  }

      return sqrt(grad)/s;
}


double CMagDipole::ComputeGradient(double * add_params, double x_obs, double y_obs, double z_obs,
				   double x0,    double y0,    double z0,
				   double A,     double B,     double C) {
	

    double deriv_tmp[3];
    if(m_bInduced) {
    dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     m_lfJ[0]*A, m_lfJ[0]*B, m_lfJ[0]*C, &deriv_tmp[1], &deriv_tmp[0], &deriv_tmp[2]);
    }
    else {	
      dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
         m_lfJ[0], m_lfJ[1], m_lfJ[2], &deriv_tmp[1], &deriv_tmp[0], &deriv_tmp[2]);
    }

    for(int i=0; i<3; i++) deriv_tmp[i] = -deriv_tmp[i]; // gradinet always bottom - top

   return add_params[DIRECTION_X]*deriv_tmp[0] + 
          add_params[DIRECTION_Y]*deriv_tmp[1] +
          add_params[DIRECTION_Z]*deriv_tmp[2];

}

double CMagDipole::ComputeAbsGradient(double * add_params, double x_obs, double y_obs, double z_obs,
				   double x0,    double y0,    double z0,
				   double A,     double B,     double C) {

    double gx, gy, gz;

    add_params[DIRECTION_X] = 1;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z]=0;  
    gx = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 1;  add_params[DIRECTION_Z]=0;  
    gy = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z]=1;  
    gz = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    return sqrt(gx*gx + gy*gy + gz*gz);

}


void CMagDipole::DumpParameters(FILE * dat){

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


  if(m_bInduced) {
    fprintf(dat,"Mag. Moment - induced only:\n");
    if(m_pJfixed[0]) 
      fprintf(dat," Jt: %lg FIXED\n",  m_lfJ[0]);
    else
      fprintf(dat," Jt: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);
   }
  else {
    fprintf(dat,"Mag. Moment:\n");

    if(m_pJfixed[0]) 
      fprintf(dat," Jx: %lg FIXED\n",  m_lfJ[1]);
    else
      fprintf(dat," Jx: %lg +/- %lg\n",  m_lfJ[1], m_lfStdDevJ[1]);

    if(m_pJfixed[1]) 
      fprintf(dat," Jy: %lg FIXED\n",  m_lfJ[0]);
    else
      fprintf(dat," Jy: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);

    if(m_pJfixed[2]) 
      fprintf(dat," Jz: %lg FIXED\n",  m_lfJ[2]);
    else
      fprintf(dat," Jz: %lg +/- %lg\n",  m_lfJ[2], m_lfStdDevJ[2]);
      
      fprintf(stderr,"Total mag. moment: %lg\n", 
            sqrt(m_lfJ[0]*m_lfJ[0] + m_lfJ[1]*m_lfJ[1] + m_lfJ[2]*m_lfJ[2]));
  }


}


void CMagDipole::DumpToStream(std::ostream & str, std::string comment, 
						      char * close_str, char * delimit, std::string offset) {

	str << comment.c_str() << m_sObjectType << std::endl;

   str << comment.c_str();

   str << "X="     << m_lfPos[0] << delimit;
   str << "Y="     << m_lfPos[1] << delimit;
   str << "Z="     << m_lfPos[2]       << delimit;
   str << "X_std=" << m_lfStdDevPos[0] << delimit;
   str << "Y_std=" << m_lfStdDevPos[1] << delimit;
   str << "Z_std=" << m_lfStdDevPos[2] << delimit;


   if(m_bInduced) {
		str << "INDUCED_ONLY=1" << delimit << "Jtotal=" << sqrt(m_lfJ[0]*m_lfJ[0] + m_lfJ[1]*m_lfJ[1] + m_lfJ[2]*m_lfJ[2]);
		str << "Jtotal_std=" << m_lfStdDevJ[0] << delimit;
   }
   else {
 	    str << "INDUCED_ONLY=0" << delimit;
		str << "Jx="            << m_lfJ[0] << delimit;
		str << "Jy="            << m_lfJ[1] << delimit;
		str << "Jz="            << m_lfJ[2] << delimit;
		str << "Jtotal="        << sqrt(m_lfJ[0]*m_lfJ[0] + m_lfJ[1]*m_lfJ[1] + m_lfJ[2]*m_lfJ[2]) << delimit;
		str << "Jx_std="        << m_lfStdDevJ[0] << delimit;
		str << "Jy_std="        << m_lfStdDevJ[1] << delimit;
		str << "Jz_std="        << m_lfStdDevJ[2] << delimit;
   }

   str << std::endl << comment.c_str()  << close_str << m_sObjectType << std::endl;


}



void CMagDipole::DumpToXMLStream(std::ostream & str) {

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


   if(m_bInduced) {
		str << "INDUCED_ONLY=\"1\"" << delimit << "Jtotal=" << quot << sqrt(m_lfJ[0]*m_lfJ[0] + m_lfJ[1]*m_lfJ[1] + m_lfJ[2]*m_lfJ[2]) << quot << delimit;
		str << "Jtotal_std=" << quot << m_lfStdDevJ[0] << quot << delimit;
   }
   else {
 	    str << "INDUCED_ONLY=\"0\"" << delimit;
		str << "Jx="            << quot << m_lfJ[0] << quot << delimit;
		str << "Jy="            << quot << m_lfJ[1] << quot << delimit;
		str << "Jz="            << quot << m_lfJ[2] << quot << delimit;
		str << "Jtotal="        << quot << sqrt(m_lfJ[0]*m_lfJ[0] + m_lfJ[1]*m_lfJ[1] + m_lfJ[2]*m_lfJ[2]) << quot << delimit;
		str << "Jx_std="        << quot << m_lfStdDevJ[0] << quot << delimit;
		str << "Jy_std="        << quot << m_lfStdDevJ[1] << quot << delimit;
		str << "Jz_std="        << quot << m_lfStdDevJ[2] << quot << delimit;
   }

   str << "/>" << std::endl;
   str << "</" << m_sObjectType << ">" << std::endl;

}


int CMagDipole::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
				     double * pos_deriv,
				     double x_obs, double y_obs, double z_obs,
				     double x0,    double y0,    double z0,
				     double A,     double B,     double C,
				     int data_pos, int *n_params, int * start) { 

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
	    return GetTotalFieldLinearDerivatives(add_params, deriv, size, pos_deriv, 
						  x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

	case totalmaggradient:   
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

	default:
	    return CMagObject::GetLinearDerivatives(type, add_params, deriv, size, pos_deriv, 
						    x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, 
						    data_pos, n_params, start);
    }

    return 0;

}


int CMagDipole::GetTotalFieldLinearDerivatives(double * add_params, double * deriv, int size, 
				               double * pos_deriv,
				               double x_obs, double y_obs, double z_obs,
				               double x0,    double y0,    double z0,
				               double A,     double B,     double C,
				               int *n_params, int * start) { 
				   

 int i, j;

 double deriv_tmp[3];

 if(start) *start =  m_iLinearStart;
 *n_params = 0;


if(m_bInduced) {
if(size < 1) return 0;
if(m_pJfixed[0]) return 0;
 deriv[0] = dT_2(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, A, B, C);
 *n_params  = 1;
 return 1;
 }

if(size < 3) return 0;
dT_inv_2(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C, 
	 &deriv_tmp[0], &deriv_tmp[1], &deriv_tmp[2]);
 *n_params  = 0;
 j=0;
 for(i=0; i<3; i++) {
   if(!m_pJfixed[i]) { 
     deriv[j] = deriv_tmp[i];  
     j++;
     (*n_params)++;
   }
 }


 return 1;

}	   



int CMagDipole::GetFiniteGradientLinearDerivatives(double * add_params, double * deriv, int size, 
						   double * pos_deriv,
						   double x_obs, double y_obs, double z_obs,
						   double x0,    double y0,    double z0,
						   double A,     double B,     double C,
						   int *n_params, int * start) { 

int i;
 *n_params = 0;
 double s = add_params[FINITE_SEPARATION];
 if(s <= 0.) return 0; 
	
 double x_obs2 = x_obs + add_params[DIRECTION_X]*s;
 double y_obs2 = y_obs + add_params[DIRECTION_Y]*s;
 double z_obs2 = z_obs + add_params[DIRECTION_Z]*s;

 double deriv1[3], deriv2[3];


 int ret1 = GetTotalFieldLinearDerivatives(add_params, deriv1, 3, pos_deriv,  x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);
 int ret2 = GetTotalFieldLinearDerivatives(add_params, deriv2, 3, pos_deriv,  x_obs2, y_obs2,  z_obs2, x0, y0, z0, A, B, C, n_params, start);

 if(!ret1 || !ret2) {
     *n_params = 0;
     return 0;
 }

 for(i=0; i<*n_params; i++) {
     deriv[i] = (deriv2[i] - deriv1[i])/s;
 }


 return 1;

}

				   

int CMagDipole::GetGradientLinearDerivatives(double * add_params, double * deriv, int size, 
				             double * pos_deriv,
				             double x_obs, double y_obs, double z_obs,
				             double x0,    double y0,    double z0,
				             double A,     double B,     double C,
				             int *n_params, int * start) { 
				   
 int i, j;


 double deriv_tmp[3];
 double deriv_tmp2[3];

 if(start) *start =  m_iLinearStart;
 *n_params = 0;

 if(m_bInduced) {
  if(size < 1) return 0;
  if(m_pJfixed[0]) return 0;
  dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
       A, B, C, &deriv_tmp[1], &deriv_tmp[0], &deriv_tmp[2]);

  for(int i=0; i<3; i++) deriv_tmp[i] = -deriv_tmp[i]; // gradinet always bottom - top

  deriv[0] = add_params[DIRECTION_X]*deriv_tmp[0] + 
             add_params[DIRECTION_Y]*deriv_tmp[1] +
             add_params[DIRECTION_Z]*deriv_tmp[2];

  *n_params  = 1;
  return 1;
  }

  // any direction;
  if(size < 3) return 0;

   dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     1., 0., 0., &deriv_tmp2[1], &deriv_tmp2[0], &deriv_tmp2[2]);

   for(i=0; i<3; i++) deriv_tmp2[i] = -deriv_tmp2[i]; // gradinet always bottom - top

   deriv_tmp[0] = add_params[DIRECTION_X]*deriv_tmp2[0] + 
                  add_params[DIRECTION_Y]*deriv_tmp2[1] +
                  add_params[DIRECTION_Z]*deriv_tmp2[2];


   dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     0., 1., 0., &deriv_tmp2[1], &deriv_tmp2[0], &deriv_tmp2[2]);

   for(i=0; i<3; i++) deriv_tmp2[i] = -deriv_tmp2[i]; // gradinet always bottom - top


   deriv_tmp[1] = add_params[DIRECTION_X]*deriv_tmp2[0] + 
                  add_params[DIRECTION_Y]*deriv_tmp2[1] +
                  add_params[DIRECTION_Z]*deriv_tmp2[2];

   dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     0., 0., 1., &deriv_tmp2[1], &deriv_tmp2[0], &deriv_tmp2[2]);

   for(i=0; i<3; i++) deriv_tmp2[i] = -deriv_tmp2[i]; // gradinet always bottom - top

   deriv_tmp[2] = add_params[DIRECTION_X]*deriv_tmp2[0] + 
                  add_params[DIRECTION_Y]*deriv_tmp2[1] +
 	          add_params[DIRECTION_Z]*deriv_tmp2[2];

   *n_params  = 0;
   j=0;
   for(i=0; i<3; i++) {
    if(!m_pJfixed[i]) { 
      deriv[j] = deriv_tmp[i];
     j++;
     (*n_params)++;
   }
 }

  return 1;

}


int CMagDipole::GetAbsGradientLinearDerivatives(double * add_params, double * deriv, int size, 
				             double * pos_deriv,
				             double x_obs, double y_obs, double z_obs,
				             double x0,    double y0,    double z0,
				             double A,     double B,     double C,
				             int *n_params, int * start) { 


    double gx, gy, gz;
    int i;

    *n_params = 0;


    add_params[DIRECTION_X] = 1;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z]=0;  
    gx = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 1;  add_params[DIRECTION_Z]=0;  
    gy = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z]=1;  
    gz = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    double grad = sqrt(gx*gx + gy*gy + gz*gz);

    if(grad <= 0.) grad = gx = gy = gz = 1.;

    double derivx[3], derivy[3], derivz[3];

    add_params[DIRECTION_X] = 1;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z] = 0;  
    GetGradientLinearDerivatives(add_params, derivx, size, pos_deriv, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 1;  add_params[DIRECTION_Z] = 0;  
    GetGradientLinearDerivatives(add_params, derivy, size, pos_deriv, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);
    
    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z] = 1;  
    GetGradientLinearDerivatives(add_params, derivz, size, pos_deriv, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

    for(i=0; i<*n_params; i++) {
	deriv[i] = (gx*derivx[i] + gy*derivy[i] + gz*derivz[i])/grad;
    }
	
    return 1;

}


int CMagDipole::GetNonLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
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

	default:
	    return CMagObject::GetNonLinearDerivatives(type, add_params, deriv, size, pos_deriv, 
						       x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, 
						       data_pos, n_params, start);
    }

    return 0;

}


int CMagDipole::GetTotalFieldNonLinearDerivatives(double * add_params, double * deriv, int size, 
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

if(m_bInduced) {
dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
     m_lfJ[0]*A, m_lfJ[0]*B, m_lfJ[0]*C, &deriv_tmp[1], &deriv_tmp[0], &deriv_tmp[2]);
}
else 	
  dT_d(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
     m_lfJ[0], m_lfJ[1], m_lfJ[2], &deriv_tmp[1], &deriv_tmp[0], &deriv_tmp[2]); 

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

int CMagDipole::GetGradientNonLinearDerivatives(double * add_params, double * deriv, int size, 
				 	        double * pos_deriv,
				                double x_obs, double y_obs, double z_obs,
                                                double x0,    double y0,    double z0,
				                double A,     double B,     double C,
					        int * n_params, int * start) {

// non-linear parameters NOTE that d/dx and d/dy swapped!


 int i, j;

if(!m_bIsValid || size < 3) return 0;

 if(start) *start =  m_iNonLinearStart;

 *n_params = 0;

   double deriv_tmp_x[3];
   double deriv_tmp_y[3];
   double deriv_tmp_z[3];

   if(m_bInduced) {
     	dT_d_dx0(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     m_lfJ[0]*A, m_lfJ[0]*B, m_lfJ[0]*C, &deriv_tmp_x[1], &deriv_tmp_x[0], &deriv_tmp_x[2]);
 
    	dT_d_dy0(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     m_lfJ[0]*A, m_lfJ[0]*B, m_lfJ[0]*C, &deriv_tmp_y[1], &deriv_tmp_y[0], &deriv_tmp_y[2]);

    	dT_d_dz0(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     m_lfJ[0]*A, m_lfJ[0]*B, m_lfJ[0]*C, &deriv_tmp_z[1], &deriv_tmp_z[0], &deriv_tmp_z[2]);
   }
   else {
     	dT_d_dx0(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     m_lfJ[0], m_lfJ[1], m_lfJ[2], &deriv_tmp_x[1], &deriv_tmp_x[0], &deriv_tmp_x[2]);
 
    	dT_d_dy0(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     m_lfJ[0], m_lfJ[1], m_lfJ[2], &deriv_tmp_y[1], &deriv_tmp_y[0], &deriv_tmp_y[2]);

    	dT_d_dz0(x_obs,y_obs,z_obs, m_lfPos[0]-x0, m_lfPos[1]-y0, m_lfPos[2]-z0, A, B, C,
	     m_lfJ[0], m_lfJ[1], m_lfJ[2], &deriv_tmp_z[1], &deriv_tmp_z[0], &deriv_tmp_z[2]);

   }


	 pos_deriv[1] =  -add_params[DIRECTION_X]*deriv_tmp_x[0] + 
	                 -add_params[DIRECTION_Y]*deriv_tmp_x[1] +
	                 -add_params[DIRECTION_Z]*deriv_tmp_x[2];

	 pos_deriv[0] =  -add_params[DIRECTION_X]*deriv_tmp_y[0] + 
	                 -add_params[DIRECTION_Y]*deriv_tmp_y[1] +
	                 -add_params[DIRECTION_Z]*deriv_tmp_y[2];

	 pos_deriv[2] =  -add_params[DIRECTION_X]*deriv_tmp_z[0] + 
	                 -add_params[DIRECTION_Y]*deriv_tmp_z[1] +
	                 -add_params[DIRECTION_Z]*deriv_tmp_z[2];


	 j=0;
	 for(i=0; i<3; i++) {
	  if(!m_pPosfixed[i]) { 
	    deriv[j] = pos_deriv[i];
	    j++;
	    (*n_params)++;
	  }
	}
	

   return 1;   
 }


int CMagDipole::GetAbsGradientNonLinearDerivatives(double * add_params, double * deriv, int size, 
				             double * pos_deriv,
				             double x_obs, double y_obs, double z_obs,
				             double x0,    double y0,    double z0,
				             double A,     double B,     double C,
				             int *n_params, int * start) { 


    double gx, gy, gz;
    int i;

    *n_params = 0;

    add_params[DIRECTION_X] = 1;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z]=0;  
    gx = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 1;  add_params[DIRECTION_Z]=0;  
    gy = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z]=1;  
    gz = ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

    double grad = sqrt(gx*gx + gy*gy + gz*gz);

    if(grad <= 0.) grad = gx = gy = gz = 1.;

    double derivx[3], derivy[3], derivz[3];
    double pos_derivx[3], pos_derivy[3], pos_derivz[3];

    add_params[DIRECTION_X] = 1;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z] = 0;  
    GetGradientNonLinearDerivatives(add_params, derivx, size, pos_derivx, x_obs, y_obs,  z_obs, 
				    x0, y0, z0, A, B, C, n_params, start);

    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 1;  add_params[DIRECTION_Z] = 0;  
    GetGradientNonLinearDerivatives(add_params, derivy, size, pos_derivy, x_obs, y_obs,  z_obs, 
				    x0, y0, z0, A, B, C, n_params, start);
    
    add_params[DIRECTION_X] = 0;  add_params[DIRECTION_Y] = 0;  add_params[DIRECTION_Z] = 1;  
    GetGradientNonLinearDerivatives(add_params, derivz, size, pos_derivz, x_obs, y_obs,  z_obs, 
				    x0, y0, z0, A, B, C, n_params, start);

    for(i=0; i<*n_params; i++) {
	deriv[i] = (gx*derivx[i] + gy*derivy[i] + gz*derivz[i])/grad;
    }


   for(i=0; i<3; i++) {
	pos_deriv[i] = (gx*pos_derivx[i] + gy*pos_derivy[i] + gz*pos_derivz[i])/grad;
    }

    return 1;

}




int CMagDipole::GetFiniteGradientNonLinearDerivatives(double * add_params, double * deriv, int size, 
						   double * pos_deriv,
						   double x_obs, double y_obs, double z_obs,
						   double x0,    double y0,    double z0,
						   double A,     double B,     double C,
						   int *n_params, int * start) { 

 int i;
 *n_params = 0;
 double s = add_params[FINITE_SEPARATION];
 if(s <= 0.) return 0; 
	
 double x_obs2 = x_obs + add_params[DIRECTION_X]*s;
 double y_obs2 = y_obs + add_params[DIRECTION_Y]*s;
 double z_obs2 = z_obs + add_params[DIRECTION_Z]*s;

 double deriv1[3],     deriv2[3];
 double pos_deriv1[3], pos_deriv2[3];

 int ret1 = GetTotalFieldNonLinearDerivatives(add_params, deriv1, size, pos_deriv1, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);
 int ret2 = GetTotalFieldNonLinearDerivatives(add_params, deriv2, size, pos_deriv2, x_obs2, y_obs2,  z_obs2, x0, y0, z0, A, B, C, 
					      n_params, start);

 if(!ret1 || !ret2) {
     *n_params = 0;
     return 0;
 }

 for(i=0; i<*n_params; i++) 
     deriv[i] = (deriv2[i] - deriv1[i])/s;

 for(i=0; i<3; i++) pos_deriv[i] = (pos_deriv2[i] - pos_deriv1[i])/s;

 return 1;

}




void CMagDipole::GetDipoleMoment(double A, double B, double C, double Az,
                       double * Jtotal, double * incl, double  *decl,
                       double * angle) {
                           
     double rad=M_PI/180.;
                                                   
      if(m_bInduced) {
          *Jtotal = m_lfJ[0];
          *angle  = 0.;                                                        
 
           *incl = asin(C)/rad;
           *decl = atan2(B,A)/rad; //acos(A/(sqrt(1-C*C)))/rad;
           return;
       }   
       
   *Jtotal = sqrt(m_lfJ[0]* m_lfJ[0]+m_lfJ[1]* m_lfJ[1] + m_lfJ[2]* m_lfJ[2]);
    *angle = acos((m_lfJ[0]*A+m_lfJ[1]*B+m_lfJ[2]*C)/(*Jtotal))/rad;
    
    *incl = asin(m_lfJ[2]/(*Jtotal))/rad;  
    *decl = (atan2(m_lfJ[1], m_lfJ[0])/rad) - 90. + Az;
                       
}    



int CMagDipole::GetInfoHeader(char * str, int size) {



	return 1;
}

int CMagDipole::GetDataHeader(char * str, int size) {


	return 1;
}












