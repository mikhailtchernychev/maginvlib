//
// @(#)magpipe.cpp  1.1  misha-04nov103
//
// Copyright (c) 2003 Geometrics, Inc.
//
// Created:  04nov103  by  misha@misha2.geometrics.com
// Version:  04nov103  by  misha@misha2.geometrics.com
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "magpipe.h"

#define R(x,y,z)  sqrt((x)*(x)+(y)*(y)+(z)*(z))

#define get_der(op, x, z, y1, y2) (op)((x),(y2), (z)) -  (op)((x),(y1), (z))

int mult_vect(double * a, double * x, double * b);

double Vz_sub(double x, double y, double z) {
  return y/(R(x,y,z));
}



double Vz_line(double x, double y, double z, double y_1, double y_2) {
double mult;

mult = z/(x*x+z*z);

return mult*((y_2/(R(x,y_2,z))) - (y_1/(R(x,y_1,z)))); 

}


double Vz_line1(double x, double y, double z, double y_1, double y_2) {

double V1, V2;

V2 = z/((y_2+R(x,y_2,z))*R(x,y_2,z));
V1 = z/((y_1+R(x,y_1,z))*R(x,y_1,z));

return V1 - V2; 

}


double Vzz_2(double x, double y, double z) {
double r = R(x,y,z);
return ((2.*z*z-x*x)*y/r) - (z*z*y*y*y/(r*r*r));
}


double Vzz_line(double x, double y, double z, double y_1, double y_2) {

double mult = ((x*x+z*z)*(x*x+z*z));

if(mult < 1.e-14)  return (1./(2.*y_2*y_2)) -  (1./(2.*y_1*y_1));


return (Vzz_2(x,y_2,z) - Vzz_2(x,y_1,z))/mult;

}


double Vxz_2(double x, double y, double z) {
double r = R(x,y,z);

return y*((3.*(x*x+z*z)+2.*y*y)/(r*r*r));

}



double Vxz_line(double x, double y, double z, double y_1, double y_2) {

double mult = ((x*x+z*z)*(x*x+z*z));

if(mult < 1.e-14) return 0.; 

return x*z*(Vxz_2(x,y_2,z) - Vxz_2(x,y_1,z))/mult;

}


double Vyy_line(double x, double y, double z, double y_1, double y_2) {

double r1 = R(x,y_1,z);
double r2 = R(x,y_2,z);

return y_1/(r1*r1*r1) - y_2/(r2*r2*r2);
}


double Vyx_line(double x, double y, double z, double y_1, double y_2) {

double r1 = R(x,y_1,z);
double r2 = R(x,y_2,z);

return x/(r1*r1*r1) - x/(r2*r2*r2);

}



// new functions 11/07/03

double Vzz_half(double x, double y, double z) {

  double r = R(x,y,z);
  double r2 = r*r;

  return (r2*(r+y)-z*z*(r+y)-z*z*r)/(r2*r*(r+y)*(r+y));

}


double Vxx_half(double x, double y, double z) {

  double r = R(x,y,z);
  double r2 = r*r;

  return (r2*(r+y)-x*x*(r+y)-x*x*r)/(r2*r*(r+y)*(r+y));

}


double Vxy_half(double x, double y, double z) {
 
  double r = R(x,y,z);
  double r2 = r*r;
  
  return (-x*r-x*y)/(r2*r*(r+y));

}


double Vxz_half(double x, double y, double z) {

  double r = R(x,y,z);
  double r2 = r*r;

  return (-x*z*(r+y)-x*z*r)/(r2*r*(r+y)*(r+y));

}

double Vyz_half(double x, double y, double z) {

  double r = R(x,y,z);
  double r2 = r*r;

  return (-z*r-y*z)/(r2*r*(r+y));

}


// derivatives functions, 11/10/02


/* ----------------- for X grad component ------------------*/

double Vzz_x_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double z2  = z*z;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(3.*x*yr*z2 + 2.*r*x*z2 - r2*x*yr) + 3.*x*yr2*z2 - r2*x*yr2) / (r2*r2*r*yr2*yr);
}


double Vxx_x_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double x3  = x*x*x;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(3.*x3*yr - 3.*r2*x*yr + 2.*r*x3) + 3.*x3*yr2 -3.* r2*x*yr2) / (r2*r2*r*yr2*yr);
}


double Vxy_x_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double x2  = x*x;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(x2*yr-r2*yr+2.*x2*y+2.*r*x2) + 3.*x2*y*yr - r2*y*yr) / (r2*r2*r*yr2);
}


double Vxz_x_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double x2  = x*x;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(3.*x2*yr*z - r2*yr*z +2.*r*x2*z) + 3.*x2*yr2*z - r2*yr2*z) / (r2*r2*r*yr2*yr);
}



double Vyz_x_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(x*yr*z+2.*x*y*z + 2.*r*x*z) + 3.*x*y*yr*z) / (r2*r2*r*yr2); 
}




/* ----------------- for Y grad component ------------------*/

double Vzz_y_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double z2  = z*z;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(yr*z2+2.*y*z2+2.*r*z2 - r2*yr) + 3.*y*yr*z2 - r2*y*yr) / (r2*r2*r*yr2);
}


double Vxx_y_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double x2  = x*x;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(x2*yr-r2*yr+2.*x2*y+2.*r*x2) + 3.*x2*y*yr-r2*y*yr) /  (r2*r2*r*yr2);
}


double Vxy_y_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double x2  = x*x;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(-x*(x2+z*z) + 2.*x*y*yr + 2.*r*x*yr) + 3.*x*y*y*yr - r2*x*yr) / (r2*r2*r*yr2); 
}


double Vxz_y_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double x2  = x*x;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(x*yr*z+2.*x*y*z+2.*r*x*z) + 3.*x*y*yr*z) / (r2*r2*r*yr2); 
}



double Vyz_y_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double yr  = y+r;
  double yr2 = yr*yr;
  double x2 = x*x;
  double y2 = y*y;
  double z2 = z*z;

  return (r*(-z*(x2+z2) + 2.*y*yr*z + 2*r*yr*z) + 3.*y2*yr*z - r2*yr*z) / (r2*r2*r*yr2); 
}





/* ----------------- for Z grad component ------------------*/

double Vzz_z_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double z3  = z*z*z;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(3.*yr*z3+2.*r*z3-3.*r2*yr*z) + 3.*yr2*z3 - 3.*r2*yr2*z) / (r2*r2*r*yr2*yr);
}


double Vxx_z_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double yr  = y+r;
  double yr2 = yr*yr;
  double x2  = x*x;

  return (r*(3.*x2*yr*z - r2*yr*z + 2*r*x2*z) + 3.*x2*yr2*z - r2*yr2*z) / (r2*r2*r*yr2*yr);
}


double Vxy_z_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double z2  = z*z;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(x*yr*z+2.*x*y*z+2.*r*x*z) + 3.*x*y*yr*z) / (r2*r2*r*yr2);
}


double Vxz_z_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double z2  = z*z;
  double yr  = y+r;
  double yr2 = yr*yr;

  return (r*(3.*x*yr*z2+2.*r*x*z2 - r2*x*yr) + 3.*x*yr2*z2 - r2*x*yr2) / (r2*r2*r*yr2*yr);
}



double Vyz_z_half(double x, double y, double z) {

  double r   = R(x,y,z);
  double r2  = r*r;
  double yr  = y+r;
  double yr2 = yr*yr;
  double z2 = z*z;

  return (r*(yr*z2+2.*y*z2+2.*r*z2-r2*yr) +  3.*y*yr*z2 - r2*y*yr) / (r2*r2*r*yr2);
}






/** Calculation of total magnetic field from line segment

r  - observation position
r0 - center of the segment position
n  - unit vector along main magnetic field
J  - magnetisation vector
a  - half of the size
A_tr - transformation matrix

*/


double dT_line(double * r, double * r0, double *n, double *J, double a, double * A_tr) {

double X, Y, Z, Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, yy=0.;
double rt[3], Jt[3], ri[3], nt[3];
int i;

/* transform coordinates */

for(i=0; i<3; i++) ri[i] = r[i] - r0[i];  
mult_vect(A_tr, ri, rt);
mult_vect(A_tr, J,  Jt); 
mult_vect(A_tr, n,  nt); 

Vyy = Vyy_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vzz = Vzz_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vxx = -Vyy - Vzz;

Vxy = Vyx_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vxz = Vxz_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vyz = Vyx_line(rt[2],yy,rt[0], rt[1]-a, rt[1]+a);

X = Jt[0]*Vxx + Jt[1]*Vxy  +Jt[2]*Vxz;
Y = Jt[0]*Vxy + Jt[1]*Vyy  +Jt[2]*Vyz;
Z = Jt[0]*Vxz + Jt[1]*Vyz  +Jt[2]*Vzz;

return(X*nt[0] + Y*nt[1] + Z*nt[2]);
}


double dT_line2(double * r, double * r0, double *n, double *J, double a, double * A_tr) {

double X, Y, Z, Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, yy=0.;
double rt[3], Jt[3], ri[3], nt[3];
int i;

/* transform coordinates */

for(i=0; i<3; i++) ri[i] = r[i] - r0[i];  
mult_vect(A_tr, ri, rt);
mult_vect(A_tr, J,  Jt); 
mult_vect(A_tr, n,  nt); 

Vxx = get_der(Vxx_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);

X = Jt[0]*Vxx + Jt[1]*Vxy  +Jt[2]*Vxz;
Y = Jt[0]*Vxy + Jt[1]*Vyy  +Jt[2]*Vyz;
Z = Jt[0]*Vxz + Jt[1]*Vyz  +Jt[2]*Vzz;

return(X*nt[0] + Y*nt[1] + Z*nt[2]);
}



int  grad_line(double * r, double * r0, double *n, double *J, double a, double * A_tr,
	       double *grad_x, double *grad_y,  double * grad_z ) {

double X, Y, Z, Vxx, Vyy, Vzz, Vxy, Vxz, Vyz,  yy=0.;
double rt[3], Jt[3], ri[3], nt[3];
int i;

/* transform coordinates */

for(i=0; i<3; i++) ri[i] = r[i] - r0[i];  
mult_vect(A_tr, ri, rt);
mult_vect(A_tr, J,  Jt); 
mult_vect(A_tr, n,  nt); 

Vxx = get_der(Vxx_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);

X = Jt[0]*Vxx + Jt[1]*Vxy  +Jt[2]*Vxz;
Y = Jt[0]*Vxy + Jt[1]*Vyy  +Jt[2]*Vyz;
Z = Jt[0]*Vxz + Jt[1]*Vyz  +Jt[2]*Vzz;

*grad_x = X*nt[0] + Y*nt[1] + Z*nt[2];


Vxx = get_der(Vxx_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);

X = Jt[0]*Vxx + Jt[1]*Vxy  +Jt[2]*Vxz;
Y = Jt[0]*Vxy + Jt[1]*Vyy  +Jt[2]*Vyz;
Z = Jt[0]*Vxz + Jt[1]*Vyz  +Jt[2]*Vzz;

*grad_y = X*nt[0] + Y*nt[1] + Z*nt[2];


Vxx = get_der(Vxx_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);

X = Jt[0]*Vxx + Jt[1]*Vxy  +Jt[2]*Vxz;
Y = Jt[0]*Vxy + Jt[1]*Vyy  +Jt[2]*Vyz;
Z = Jt[0]*Vxz + Jt[1]*Vyz  +Jt[2]*Vzz;

*grad_z = X*nt[0] + Y*nt[1] + Z*nt[2];

return 1;
}




int  grad_line_matrix(double * r, double * r0, double *n,  double a, double * A_tr,
	       double * out ) {

double X, Y, Z, Vxx, Vyy, Vzz, Vxy, Vxz, Vyz,  yy=0.;
double rt[3], ri[3], nt[3], J[3], Jtx[3], Jty[3], Jtz[3];
int i;


for(i=0; i<9; i++) out[i] = 0;

/* transform coordinates */

for(i=0; i<3; i++) ri[i] = r[i] - r0[i];  
mult_vect(A_tr, ri, rt);
mult_vect(A_tr, n,  nt); 

// transform unit J's into new system

J[0] = 1.; J[1] = 0.; J[2] = 0.; mult_vect(A_tr, J,  Jtx); 
J[0] = 0.; J[1] = 1.; J[2] = 0.; mult_vect(A_tr, J,  Jty); 
J[0] = 0.; J[1] = 0.; J[2] = 1.; mult_vect(A_tr, J,  Jtz); 


/* ------------------- d/dx -------------------------- */

Vxx = get_der(Vxx_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_x_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);

X = Jtx[0]*Vxx + Jtx[1]*Vxy  +Jtx[2]*Vxz;
Y = Jtx[0]*Vxy + Jtx[1]*Vyy  +Jtx[2]*Vyz;
Z = Jtx[0]*Vxz + Jtx[1]*Vyz  +Jtx[2]*Vzz;
out[0] = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jty[0]*Vxx + Jty[1]*Vxy  +Jty[2]*Vxz;
Y = Jty[0]*Vxy + Jty[1]*Vyy  +Jty[2]*Vyz;
Z = Jty[0]*Vxz + Jty[1]*Vyz  +Jty[2]*Vzz;
out[1] = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jtz[0]*Vxx + Jtz[1]*Vxy  +Jtz[2]*Vxz;
Y = Jtz[0]*Vxy + Jtz[1]*Vyy  +Jtz[2]*Vyz;
Z = Jtz[0]*Vxz + Jtz[1]*Vyz  +Jtz[2]*Vzz;
out[2] = X*nt[0] + Y*nt[1] + Z*nt[2];

/* ------------------- d/dy -------------------------- */

Vxx = get_der(Vxx_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_y_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);

X = Jtx[0]*Vxx + Jtx[1]*Vxy  +Jtx[2]*Vxz;
Y = Jtx[0]*Vxy + Jtx[1]*Vyy  +Jtx[2]*Vyz;
Z = Jtx[0]*Vxz + Jtx[1]*Vyz  +Jtx[2]*Vzz;
out[3] = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jty[0]*Vxx + Jty[1]*Vxy  +Jty[2]*Vxz;
Y = Jty[0]*Vxy + Jty[1]*Vyy  +Jty[2]*Vyz;
Z = Jty[0]*Vxz + Jty[1]*Vyz  +Jty[2]*Vzz;
out[4] = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jtz[0]*Vxx + Jtz[1]*Vxy  +Jtz[2]*Vxz;
Y = Jtz[0]*Vxy + Jtz[1]*Vyy  +Jtz[2]*Vyz;
Z = Jtz[0]*Vxz + Jtz[1]*Vyz  +Jtz[2]*Vzz;
out[5] = X*nt[0] + Y*nt[1] + Z*nt[2];


/* ------------------- d/dz -------------------------- */

Vxx = get_der(Vxx_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_z_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);

X = Jtx[0]*Vxx + Jtx[1]*Vxy  +Jtx[2]*Vxz;
Y = Jtx[0]*Vxy + Jtx[1]*Vyy  +Jtx[2]*Vyz;
Z = Jtx[0]*Vxz + Jtx[1]*Vyz  +Jtx[2]*Vzz;
out[6] = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jty[0]*Vxx + Jty[1]*Vxy  +Jty[2]*Vxz;
Y = Jty[0]*Vxy + Jty[1]*Vyy  +Jty[2]*Vyz;
Z = Jty[0]*Vxz + Jty[1]*Vyz  +Jty[2]*Vzz;
out[7] = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jtz[0]*Vxx + Jtz[1]*Vxy  +Jtz[2]*Vxz;
Y = Jtz[0]*Vxy + Jtz[1]*Vyy  +Jtz[2]*Vyz;
Z = Jtz[0]*Vxz + Jtz[1]*Vyz  +Jtz[2]*Vzz;
out[8] = X*nt[0] + Y*nt[1] + Z*nt[2];


return 1;
}







/** calculate matrix elements for J estimation */

double dT_line_matr(double * r, double * r0, double *n, double a, double * A_tr,
	      double *Tx, double *Ty, double *Tz) {

double X, Y, Z, Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, yy=0.;
double rt[3], J[3], ri[3], nt[3], Jtx[3], Jty[3], Jtz[3];
int i;

/* transform coordinates */

for(i=0; i<3; i++) ri[i] = r[i] - r0[i];  
mult_vect(A_tr, ri, rt);
mult_vect(A_tr, n,  nt); 

// transform unit J's into new system

J[0] = 1.; J[1] = 0.; J[2] = 0.; mult_vect(A_tr, J,  Jtx); 
J[0] = 0.; J[1] = 1.; J[2] = 0.; mult_vect(A_tr, J,  Jty); 
J[0] = 0.; J[1] = 0.; J[2] = 1.; mult_vect(A_tr, J,  Jtz); 

Vyy = Vyy_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vzz = Vzz_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vxx = -Vyy - Vzz;

Vxy = Vyx_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vxz = Vxz_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vyz = Vyx_line(rt[2],yy,rt[0], rt[1]-a, rt[1]+a);

X = Jtx[0]*Vxx + Jtx[1]*Vxy  +Jtx[2]*Vxz;
Y = Jtx[0]*Vxy + Jtx[1]*Vyy  +Jtx[2]*Vyz;
Z = Jtx[0]*Vxz + Jtx[1]*Vyz  +Jtx[2]*Vzz;
*Tx = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jty[0]*Vxx + Jty[1]*Vxy  +Jty[2]*Vxz;
Y = Jty[0]*Vxy + Jty[1]*Vyy  +Jty[2]*Vyz;
Z = Jty[0]*Vxz + Jty[1]*Vyz  +Jty[2]*Vzz;
*Ty = X*nt[0] + Y*nt[1] + Z*nt[2];

X = Jtz[0]*Vxx + Jtz[1]*Vxy  +Jtz[2]*Vxz;
Y = Jtz[0]*Vxy + Jtz[1]*Vyy  +Jtz[2]*Vyz;
Z = Jtz[0]*Vxz + Jtz[1]*Vyz  +Jtz[2]*Vzz;
*Tz = X*nt[0] + Y*nt[1] + Z*nt[2];

//*Tx = Vxx*nt[0] + Vxy*nt[1] + Vxz*nt[2];
//*Ty = Vxy*nt[0] + Vyy*nt[1] + Vyz*nt[2];
//*Tz = Vxz*nt[0] + Vyz*nt[1] + Vzz*nt[2];

return 0;
}


double dT_line_matr2(double * r, double * r0, double *n, double a, double * A_tr,
	      double *Tx, double *Ty, double *Tz) {

double Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, yy=0.;
double rt[3], ri[3], nt[3];
int i;

/* transform coordinates */

for(i=0; i<3; i++) ri[i] = r[i] - r0[i];  
mult_vect(A_tr, ri, rt);
mult_vect(A_tr, n,  nt); 

Vxx = get_der(Vxx_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vzz = get_der(Vzz_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyy = -Vxx  - Vzz;

Vxy = get_der(Vxy_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vxz = get_der(Vxz_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);
Vyz = get_der(Vyz_half, rt[0], rt[2],  rt[1]-a, rt[1]+a);


*Tx = Vxx*nt[0] + Vxy*nt[1] + Vxz*nt[2];
*Ty = Vxy*nt[0] + Vyy*nt[1] + Vyz*nt[2];
*Tz = Vxz*nt[0] + Vyz*nt[1] + Vzz*nt[2];

return 0;
}


/** Prepare rotation matrix for given direction

 an - vector along new Y-axis
 a  - space for matrix (3x3 array)

 */

int get_rotated(double * an, double * a, double * alen) {
double s;

*alen = s = sqrt(an[0]*an[0] + an[1]*an[1] + an[2]*an[2]);

if(s==0.) {
  fprintf(stderr,"Abnormal get_rotated call\n");
  return 0;
}


an[0] /= s; an[1] /= s; an[2] /= s;

/* fill in Y row */

a[3] = an[0]; a[4] = an[1]; a[5] = an[2];

/* X perpendicular  to an Z plane or an Y plan */

if(fabs(an[2]) < fabs(an[1])) { 
    a[0]=an[1];  a[1]=-an[0]; a[2] = 0.;
}
else {
    a[0]=-an[2];  a[1]=0.; a[2] = an[0];
}

s = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
a[0] = a[0]/s; a[1] = a[1]/s; a[2] = a[2]/s;

/* now choose Z axis as X x Y */

a[6] =   a[1]*a[5] - a[2]*a[4];
a[7] = -(a[0]*a[5] - a[2]*a[3]);
a[8]  =  a[0]*a[4] - a[1]*a[3];

s = sqrt(a[6]*a[6] + a[7]*a[7] + a[8]*a[8]);
a[6] = a[6]/s; a[7] = a[7]/s; a[8] = a[8]/s;


return 1;
}


/**  specific 3x3 matrix - vector multiplication
*/

int mult_vect(double * a, double * x, double * b) {

b[0] = a[0]*x[0] + a[1]*x[1] + a[2]*x[2];
b[1] = a[3]*x[0] + a[4]*x[1] + a[5]*x[2];
b[2] = a[6]*x[0] + a[7]*x[1] + a[8]*x[2];

return 0;
}




// ---------------------- end computational functions -------------------------

CMagPipe::CMagPipe() : CMagObject() {

  m_iObjectLinearAbilities    = 1; 
  m_iObjectNonLinearAbilities = 1;

}


CMagPipe::CMagPipe(int is_induced) : CMagObject(is_induced, 3, 6, "CMagPipe", "pipe_segment") {

  m_iObjectLinearAbilities    = 1; 
  m_iObjectNonLinearAbilities = 1;

}


int CMagPipe::GetTotalLinearParams() {
    if(m_bInduced) return 1;
    else return m_nLinearParams;

}



void CMagPipe::GetCoefficients(double * add_params, double &x_obs, double &y_obs, double &z_obs,
				 double x0,    double y0,    double z0,
			         double A,     double B,     double C,
				 double A_tr[9], double n_vect[3], double r[3],
				 double r0[3], double n_earth[3], double J[3],  double * alen) {


 x_obs += add_params[X_OFFSET];
 y_obs += add_params[Y_OFFSET];
 z_obs += add_params[Z_OFFSET];

 n_vect[1] = cos(m_lfPos[INCL_POS])*cos(m_lfPos[DECL_POS]);
 n_vect[0] = cos(m_lfPos[INCL_POS])*sin(m_lfPos[DECL_POS]);
 n_vect[2] = sin(m_lfPos[INCL_POS]);

 r[0] = x_obs;  r[1] = y_obs;  r[2] = z_obs;

 r0[0] = m_lfPos[X0_POS]-x0;
 r0[1] = m_lfPos[Y0_POS]-y0;
 r0[2] = m_lfPos[Z0_POS]-z0;

 n_earth[1] = A;
 n_earth[0] = B;
 n_earth[2] = C;

 get_rotated(n_vect, A_tr, alen);

 if(m_bInduced) {
         J[0] = m_lfJ[0]*A;
         J[1] = m_lfJ[0]*B;
         J[2] = m_lfJ[0]*C;
        }
      else {
         J[0] = m_lfJ[0];
         J[1] = m_lfJ[1];
         J[2] = m_lfJ[2];
      }

}

double CMagPipe::ComputeField(FIELDTYPE type, double * add_params, double x_obs, double y_obs, double z_obs,
					double x0,    double y0,    double z0,
				        double A,     double B,     double C, 
			                int data_pos) {
    
    if(!m_bIsValid) return 0.;
    if(m_iFieldType != anyfield  && m_iFieldType != type) return 0;

    switch(type) {

	case totalmagfield:
	    return ComputeTotalField(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

	case totalmaggradient:   
	    return ComputeGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);

	case absmaggradient:
	    return ComputeAbsGradient(add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C);	    

	default:
	    return CMagObject::ComputeField(type, add_params, x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, data_pos);

    }

    return 0.;
}



double CMagPipe::ComputeTotalField(double * add_params, double x_obs, double y_obs, double z_obs,
					double x0,    double y0,    double z0,
				        double A,     double B,     double C) {


 double A_tr[9];
 double n_vect[3], alen;
 double r[3], r0[3], n_earth[3], J[3];

 GetCoefficients(add_params, x_obs, y_obs, z_obs,  x0,  y0,  z0,  A,  B, C, A_tr, n_vect, r, r0, n_earth, J, &alen);

 return dT_line2(r, r0, n_earth, J, m_lfPos[AHALF_POS], A_tr); 

}
 


double CMagPipe::ComputeGradient(double * add_params, double x_obs, double y_obs, double z_obs,
					double x0,    double y0,    double z0,
				        double A,     double B,     double C) {


 double A_tr[9];
 double n_vect[3], alen;
 double r[3], r0[3], n_earth[3], J[3];

 GetCoefficients(add_params, x_obs, y_obs, z_obs,  x0,  y0,  z0,  A,  B, C, A_tr, n_vect, r, r0, n_earth, J, &alen);

 double grad[3], direction[3], dir_nt[3];

 grad_line(r, r0, n_earth, J, m_lfPos[AHALF_POS], A_tr, &grad[0], &grad[1], &grad[2]);

 direction[0] = add_params[DIRECTION_X];
 direction[1] = add_params[DIRECTION_Y];
 direction[2] = add_params[DIRECTION_Z];
    
 mult_vect(A_tr, direction,  dir_nt);

 return dir_nt[0]*grad[0] + dir_nt[1]*grad[1] + dir_nt[2]*grad[2];

}


double CMagPipe::ComputeAbsGradient(double * add_params, double x_obs, double y_obs, double z_obs,
				   double x0,    double y0,    double z0,
				   double A,     double B,     double C) {

 double A_tr[9];
 double n_vect[3], alen;
 double r[3], r0[3], n_earth[3], J[3];

 GetCoefficients(add_params, x_obs, y_obs, z_obs,  x0,  y0,  z0,  A,  B, C, A_tr, n_vect, r, r0, n_earth, J, &alen);

 double grad[3];

 grad_line(r, r0, n_earth, J, m_lfPos[AHALF_POS], A_tr, &grad[0], &grad[1], &grad[2]);

 return sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);

}





int CMagPipe::GetLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size, 
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

    switch(type) {

	case totalmagfield:
	    return GetTotalFieldLinearDerivatives(add_params, deriv, size, pos_deriv, 
						  x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

	case totalmaggradient:   
	    return GetGradientLinearDerivatives(add_params, deriv, size, pos_deriv, 
						x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, n_params, start);

	default:
	    return CMagObject::GetLinearDerivatives(type, add_params, deriv, size, pos_deriv, 
						    x_obs, y_obs,  z_obs, x0, y0, z0, A, B, C, data_pos,
						    n_params, start);
    }

    return 0;

}


int CMagPipe::GetTotalFieldLinearDerivatives(double * add_params, double * deriv, int size, 
					     double * pos_deriv,
					     double x_obs, double y_obs, double z_obs,
					     double x0,    double y0,    double z0,
					     double A,     double B,     double C,
					     int *n_params, int * start) { 
				    int i, j;
 double deriv_tmp[3], deriv_tmp2[3];

 if(start) *start =  m_iLinearStart;
 *n_params = 0;

 double A_tr[9];
 double n_vect[3], alen;
 double r[3], r0[3], n_earth[3], J[3];

 GetCoefficients(add_params, x_obs, y_obs, z_obs,  x0,  y0,  z0,  A,  B, C, A_tr, n_vect, r, r0, n_earth, J, &alen);

 if(m_bInduced) {
        J[0] = A;
        J[1] = B;
        J[2] = C;
	if(size < 1) return 0; 
	if(m_pJfixed[0]) return 0;
	 deriv[0] = dT_line(r, r0, n_earth, J, m_lfPos[AHALF_POS], A_tr);
	 *n_params = 1;
	 return 1;
   }

   if(size < 3) return 0;

   dT_line_matr(r, r0, n_earth, m_lfPos[AHALF_POS], A_tr,  
		  &deriv_tmp[0], &deriv_tmp[1], &deriv_tmp[2]);
   

   deriv_tmp2[0] = deriv_tmp[0];
   deriv_tmp2[1] = deriv_tmp[1];
   deriv_tmp2[2] = deriv_tmp[2];

   *n_params  = 0;
        j=0;
	for(i=0; i<3; i++) {
	    if(!m_pJfixed[i]) { 
		deriv[j] = deriv_tmp2[i];
		j++;
		(*n_params)++;
	    }
	}

      return 1;
  
}



int CMagPipe::GetGradientLinearDerivatives(double * add_params, double * deriv, int size, 
					   double * pos_deriv,
					   double x_obs, double y_obs, double z_obs,
					   double x0,    double y0,    double z0,
					   double A,     double B,     double C,
					   int *n_params, int * start) { 
				   

 int i, j;
 double deriv_tmp[3], deriv_tmp2[3];

 if(start) *start =  m_iLinearStart;
 *n_params = 0;

 double A_tr[9];
 double n_vect[3], alen;
 double r[3], r0[3], n_earth[3], J[3];

 GetCoefficients(add_params, x_obs, y_obs, z_obs,  x0,  y0,  z0,  A,  B, C, A_tr, n_vect, r, r0, n_earth, J, &alen);

 double grad[3], direction[3], dir_nt[3], grads[9];

 direction[0] = add_params[DIRECTION_X];
 direction[1] = add_params[DIRECTION_Y];
 direction[2] = add_params[DIRECTION_Z];
	     
 mult_vect(A_tr, direction,  dir_nt);

  if(m_bInduced) {
       J[0] = A;
       J[1] = B;
       J[2] = C;
       if(size < 1) return 0; 
       if(m_pJfixed[0]) return 0;

       grad_line(r, r0, n_earth, J, m_lfPos[AHALF_POS], A_tr, &grad[0], &grad[1], &grad[2]);

       deriv[0] = dir_nt[0]*grad[0] + dir_nt[1]*grad[1] + dir_nt[2]*grad[2];
       *n_params = 1;
       return 1;
   }

   if(size < 3) return 0;

   grad_line_matrix(r, r0, n_earth, m_lfPos[AHALF_POS], A_tr, grads);

   deriv_tmp[0] = dir_nt[0]*grads[0] + dir_nt[1]*grads[3] + dir_nt[2]*grads[6];
   deriv_tmp[1] = dir_nt[0]*grads[1] + dir_nt[1]*grads[4] + dir_nt[2]*grads[7];
   deriv_tmp[2] = dir_nt[0]*grads[2] + dir_nt[1]*grads[5] + dir_nt[2]*grads[8];

   deriv_tmp2[0] = deriv_tmp[0];
   deriv_tmp2[1] = deriv_tmp[1];
   deriv_tmp2[2] = deriv_tmp[2];

   *n_params  = 0;
       j=0;
       for(i=0; i<3; i++) {
	   if(!m_pJfixed[i]) { 
	       deriv[j] = deriv_tmp2[i];
	       j++;
	      (*n_params)++;
	   }
       }

      return 1;

 }



void CMagPipe::DumpParameters(FILE * dat){

  CMagObject::DumpParameters(dat);
  fprintf(dat,"Position:\n"); 
  if(m_pPosfixed[X0_POS]) 
    fprintf(dat," X: %lg FIXED\n", m_lfPos[X0_POS]);
  else 
    fprintf(dat," X: %lg +/- %lg\n", m_lfPos[X0_POS], m_lfStdDevPos[X0_POS]);

 if(m_pPosfixed[Y0_POS]) 
    fprintf(dat," Y: %lg FIXED\n", m_lfPos[Y0_POS]);
  else 
    fprintf(dat," Y: %lg +/- %lg\n", m_lfPos[Y0_POS], m_lfStdDevPos[Y0_POS]);

 if(m_pPosfixed[Z0_POS]) 
    fprintf(dat," Z: %lg FIXED\n", m_lfPos[Z0_POS]);
  else 
    fprintf(dat," Z: %lg +/- %lg\n", m_lfPos[Z0_POS], m_lfStdDevPos[Z0_POS]);


 if(m_pPosfixed[AHALF_POS]) 
    fprintf(dat," Half length: %lg FIXED\n", m_lfPos[AHALF_POS]);
  else 
    fprintf(dat," Half length: %lg +/- %lg\n", m_lfPos[AHALF_POS], m_lfStdDevPos[AHALF_POS]);


  double rad=M_PI/180.;

 if(m_pPosfixed[DECL_POS]) 
    fprintf(dat," Direction declination: %lg FIXED\n", m_lfPos[DECL_POS]/rad);
  else 
    fprintf(dat," Direction declination: %lg +/- %lg\n", m_lfPos[DECL_POS]/rad, m_lfStdDevPos[DECL_POS]/rad);


 if(m_pPosfixed[INCL_POS]) 
    fprintf(dat," Direction inclination: %lg FIXED\n", m_lfPos[INCL_POS]/rad);
  else 
    fprintf(dat," Direction inclination: %lg +/- %lg\n", m_lfPos[INCL_POS]/rad, m_lfStdDevPos[INCL_POS]/rad);



  if(m_bInduced) {
    fprintf(dat,"Mag. Moment - induced only:\n");
    if(m_pJfixed[0]) 
      fprintf(dat," Jtotal: %lg FIXED\n",  m_lfJ[0]);
    else
      fprintf(dat," Jtotal: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);
   }
  else {
    fprintf(dat,"Mag. Moment:\n");

    if(m_pJfixed[0]) 
      fprintf(dat," Jx: %lg FIXED\n",  m_lfJ[0]);
    else
      fprintf(dat," Jx: %lg +/- %lg\n",  m_lfJ[0], m_lfStdDevJ[0]);

    if(m_pJfixed[1]) 
      fprintf(dat," Jy: %lg FIXED\n",  m_lfJ[1]);
    else
      fprintf(dat," Jy: %lg +/- %lg\n",  m_lfJ[1], m_lfStdDevJ[1]);

    if(m_pJfixed[2]) 
      fprintf(dat," Jz: %lg FIXED\n",  m_lfJ[2]);
    else
      fprintf(dat," Jz: %lg +/- %lg\n",  m_lfJ[2], m_lfStdDevJ[2]);
  }


}










