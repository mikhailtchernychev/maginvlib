#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "line.h"


#define R(x,y,z)  sqrt((x)*(x)+(y)*(y)+(z)*(z))

int mult_vect(float * a, float * x, float * b);

float Vz_sub(float x, float y, float z) {
  return y/(R(x,y,z));
}



float Vz_line(float x, float y, float z, float y_1, float y_2) {
float mult;

mult = z/(x*x+z*z);

return mult*((y_2/(R(x,y_2,z))) - (y_1/(R(x,y_1,z)))); 

}


float Vz_line1(float x, float y, float z, float y_1, float y_2) {

float V1, V2;

V2 = z/((y_2+R(x,y_2,z))*R(x,y_2,z));
V1 = z/((y_1+R(x,y_1,z))*R(x,y_1,z));

return V1 - V2; 

}


float Vzz_2(float x, float y, float z) {
float r = R(x,y,z);
return float(((2.*z*z-x*x)*y/r) - (z*z*y*y*y/(r*r*r)));
}


float Vzz_line(float x, float y, float z, float y_1, float y_2) {

float mult = ((x*x+z*z)*(x*x+z*z));

if(mult < 1.e-14)  return float((1./(2.*y_2*y_2)) -  (1./(2.*y_1*y_1)));


return (Vzz_2(x,y_2,z) - Vzz_2(x,y_1,z))/mult;

}


float Vxz_2(float x, float y, float z) {
float r = R(x,y,z);

return float(y*((3.*(x*x+z*z)+2.*y*y)/(r*r*r)));

}



float Vxz_line(float x, float y, float z, float y_1, float y_2) {

float mult = ((x*x+z*z)*(x*x+z*z));

if(mult < 1.e-14) return 0.; 

return x*z*(Vxz_2(x,y_2,z) - Vxz_2(x,y_1,z))/mult;

}


float Vyy_line(float x, float y, float z, float y_1, float y_2) {

float r1 = R(x,y_1,z);
float r2 = R(x,y_2,z);

return y_1/(r1*r1*r1) - y_2/(r2*r2*r2);
}


float Vyx_line(float x, float y, float z, float y_1, float y_2) {

float r1 = R(x,y_1,z);
float r2 = R(x,y_2,z);

return x/(r1*r1*r1) - x/(r2*r2*r2);

}


/** Calculation of total magnetic field from line segment

r  - observation position
r0 - center of the segment position
n  - unit vector along main magnetic field
J  - magnetisation vector
a  - half of the size
A_tr - transformation matrix

*/


float dT_line(float * r, float * r0, float *n, float *J, float a, float * A_tr) {

float X, Y, Z, Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, yy=0.;
float rt[3], Jt[3], ri[3], nt[3];
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


/** calculate matrix elements for J estimation */

float dT_line_matr(float * r, float * r0, float *n, float a, float * A_tr,
	      float *Tx, float *Ty, float *Tz) {

float Vxx, Vyy, Vzz, Vxy, Vxz, Vyz, yy=0.;
float rt[3], ri[3], nt[3];
int i;

/* transform coordinates */

for(i=0; i<3; i++) ri[i] = r[i] - r0[i];  
mult_vect(A_tr, ri, rt);
mult_vect(A_tr, n,  nt); 

Vyy = Vyy_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vzz = Vzz_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vxx = -Vyy - Vzz;

Vxy = Vyx_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vxz = Vxz_line(rt[0],yy,rt[2], rt[1]-a, rt[1]+a);
Vyz = Vyx_line(rt[2],yy,rt[0], rt[1]-a, rt[1]+a);


*Tx = Vxx*nt[0] + Vxy*nt[1] + Vxz*nt[2];
*Ty = Vxy*nt[0] + Vyy*nt[1] + Vyz*nt[2];
*Tz = Vxz*nt[0] + Vyz*nt[1] + Vzz*nt[2];

return 0;
}


/** Prepare rotation matrix for given direction

 an - vector along new Y-axis
 a  - space for matrix (3x3 array)

 */

int get_rotated(float * an, float * a, float * alen) {
float s;


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

/* check det */

/*s1 = a[0]*(a[4]*a[8] - a[5]*a[7]) -
     a[1]*(a[3]*a[8] - a[5]*a[6]) +
     a[2]*(a[3]*a[7] - a[4]*a[6]);


fprintf(stderr,"det = %f\n",s1);

for(i=0; i<3; i++) {
  for(j=0; j<3; j++)
    fprintf(stderr,"%f ",a[3*i + j]);
  fprintf(stderr,"\n");
}
*/

return 1;
}


/**  specific 3x3 matrix - vector multiplication
*/

int mult_vect(float * a, float * x, float * b) {

b[0] = a[0]*x[0] + a[1]*x[1] + a[2]*x[2];
b[1] = a[3]*x[0] + a[4]*x[1] + a[5]*x[2];
b[2] = a[6]*x[0] + a[7]*x[1] + a[8]*x[2];

return 0;
}

