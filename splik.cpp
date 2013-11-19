#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define far

int splik(int n, double * x, double  * u,  double  dp,
          double * a, double * b, double  * c, double  * d, float *ex=NULL)

{
 int n1,n2,n3,n4,k,k1,k2,j,j1,j2;
 double a1,a2,a3,b1,b2,p1,p2,p3,h,h1,h2,h3,r,r1,r2,s,t,z;
 double far *y2, far *p; 

 if(n< 5) {
	 for(j=0; j<n; j++) {
		 d[j] = u[j];
	 }
	 return -1;
 }

 y2=(double *)calloc(n,sizeof(double));
 p=(double *)calloc(n,sizeof(double));

 n--;

  if(!ex) 
   for(j=0;j<=n;j++) p[j]=dp;
  else
   for(j=0;j<=n;j++) p[j]=ex[j];

 y2[0]=0.; y2[n]=0.;
 n1=n-1; n2=n-2; n3=n-3; n4=n-4;

 for(k=0; k<=n1; k++) d[k]=1./(x[k+1]-x[k]);

 d[n]=0.; a1=0.; a2=0.; b1=0.; a[0]=0.; b[0]=0.; b[1]=0.; j1=-1;
 p1=1./p[0]; p2=1./p[1]; h=0.;  h1=d[0]; h2=d[1]; r1=(u[1]-u[0])*h1;

 for(k=0; k<=n2; k++)
    {
     k1=k+1; k2=k+2; if(k>0) h=b1-a1*a[j1];
     h3=d[k2]; p3=1./p[k2]; s=h1+h2;
     t=2./h1+2./h2+6.*(h1*h1*p1+s*s*p2+h2*h2*p3);
     r2=(u[k2]-u[k1])*h2;
     b2=1./h2-6.*h2*(p2*s+p3*(h2+h3));
     a3=6.*p3*h2*h3;
     z=1./(t-a1*b[k]-h*a[k]);

     if(k<=n3) a[k1]=z*(b2-h*b[k1]);
     if(k<=n4) b[k2]=z*a3;
     r=6.*(r2-r1);
     if(k>=1) r=r-h*c[j1];
     if(k>1)  r=r-a1*c[j2];
     c[k]=z*r;
     j2=j1; j1=k; a1=a2; a2=a3; b1=b2; h1=h2; h2=h3; p1=p2; p2=p3; r1=r2;
    }

   y2[n1]=c[n2]; if(n3==-1) goto m4;
   y2[n2]=c[n3]-a[n2]*y2[n1];
   if(n4==-1) goto m4;

     for(j=0; j<=n4; j++)
        { k=n2-j-1; k1=k+1; k2=k+2; y2[k]=c[k-1]-a[k]*y2[k1]-b[k1]*y2[k2]; }

 m4: h1=0.;

     for(k=0; k<=n1; k++)
       {
        j2=k+1; c[k]=d[k]; h2=d[k]*(y2[j2]-y2[k]); a[k]=h2/6.;
        d[k]=u[k]-(h2-h1)/p[k]; b[k]=0.5*y2[k]; h1=h2;
       }

      d[n]=u[n]+h1/p[n];

       for(k=0; k<=n1; k++)
       { j2=k+1; h=c[k]; c[k]=(d[j2]-d[k])*h-(y2[j2]+2.*y2[k])/(6.*h); }

       c[n]=(d[n]-d[n1])*h+(2.*y2[n]+y2[n1])/(6.*h);

      free(y2); free(p);

return(0);
      }

