/*   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/*#include "f2c.h"*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define real    double
#define integer int
typedef long int logical;
typedef double doublereal;

int l1_driver_d(int m, int n, double * ct,  double *f, double *a, 
	      double prec, double eps, int *irank, int * iter, 
	      double *r, double *z, int *ind, int (*cancel_func)());


/* Subroutine */ int l1_d(integer *m, integer *n, integer *mm, integer *mmm, 
	integer *nn, real *ct, real *f, real *prec, real *eps, integer *ic, 
	integer *ir, integer *ib, real *uf, real *bp, real *xp, real *t, real 
	*alfa, real *p, real *ginv, real *vb, integer *irank, integer *iter, 
	real *r__, real *a, real *z__, integer *ind, int (*cancel_func)());



/* test function for l1_driver */


#ifdef TEST_MAIN

main() {
  int m,n, i,j, err =1;
  double *ct, *f, *a, *r;
  double prec, eps,z;
  int irank, iter, ind;

  prec = 1.e-6;
  eps = 1.e-4;

  fscanf(stdin,"%d%d", &m, &n);

  ct = f = a = r = NULL;
  fprintf(stderr,"m=%d n=%d\n", m,n);

  ct = (double *)malloc(m*n*sizeof(double)); if(!ct) goto error;
  f  = (double *)malloc(n*sizeof(double));   if(!f)  goto error;
  a  = (double *)malloc(m*sizeof(double));   if(!a)  goto error;
  r  = (double *)malloc(n*sizeof(double));   if(!r)  goto error;

  for(i=0; i<m; i++)
    for(j=0; j<n; j++) fscanf(stdin,"%f", &ct[j*m+i]);

  //for(j=0; j<m*n; j++) fscanf(stdin,"%f", &ct[j]);

  for(i=0; i<n; i++) fscanf(stdin,"%f", &f[i]);

  fprintf(stderr,"f[0]=%f f[n-1]=%f\n", f[0], f[n-1]);

  l1_driver_d(m, n, ct,  f, a, prec, eps, &irank, &iter, r, &z, &ind, NULL); 
  err =0;

  fprintf(stderr,"Solution:\n");
  for(j=0; j<m; j++) fprintf(stderr,"%f\n", a[j]);

  for(j=0; j<n; j++) fprintf(stdout,"%f %f %f\n", f[j], f[j]+r[j], r[j]);

  fprintf(stderr,"ind=%d\n", ind);



error:
  if(err) fprintf(stderr,"Allocation error in l1 test\n");
  
  if(ct) free(ct);
  if(f)  free(f);
  if(a)  free(a);
  if(r)  free(r);


}

#endif



/* driver function for L1 */ 

int l1_driver_d(int m, int n, double * ct,  double *f, double *a, 
	      double prec, double eps, int *irank, int * iter, 
	      double *r, double *z, int *ind, int (*cancel_func)()) {

  /* - parameters - 
  
   input:

    m  - number of variables
    n  - number of equitions (data)
    ct - matrix m x n
    f  - data - vector n

    prec  THE ROUND-OFF LEVEL OF THE COMPUTER. FOR THE IBM
          360/67 COMPUTER, PREC IS ABOUT 1.E-6 AND 1.E-16 FOR
          SINGLE AND DOUBLE PRECISION RESPECTIVELY.
    eps   A SPECIFIED TOLERANCE SUCH THAT A CALCULATED NUMBER
          X IS CONSIDERED ZERO IF  ABS(X) <EPS. FOR THE IBM 360/67
          COMPUTER, EPS IS USUALLY TAKEN 1.E-4 AND 1.E-11
          RESPECTIVELY.   
    
  	  
    output:
    
   irank  RANK  THE CALCULATED RANK OF MATRIX C.
   iter   THE NUMBER OF ITERATIONS WHICH THE SOLUTION NEEDED.
   a      AN MM-VECTOR. ITS FIRST M ELEMENTS ARE THE SOLUTION
          VECTOR A*.
   r      AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE 
          THE RESIDUAL VECTOR R=C*A-F.
   z      THE MINIMUM L1 NORM OF THE RESIDUAL VECTOR R.
   ind    A RETURN INDICATOR. IND=0 INDICATES THAT THE
          SOLUTION VECTOR A* IS UNIQUE. IND=1, INDICATES THAT
          A* IS MOST PROBABLY NOT UNIQUE. IND=-1 INDICATES
          PREMATURE TERMINATION OF THE CALCULATION DUE TO
          VERY ILL-CONDITIONING OF MATRIX C.
          THE MEANING OF THE OTHER PARAMETERS.

	  */

        int mm, nn, mmm, ret;

        int * ic, * ir, * ib; 
	double * uf, * bp, *xp, *t, * alfa, * p, * ginv, *vb, * r__, *z__;

        mm = m;
	nn = n;
	mmm = (mm*(mm+3))/2.;
     
	z__ = z;
	r__ = r;


        ib = ir = ic = NULL;
	uf = bp = xp = t = alfa = p = ginv = vb = NULL;

	ret = -1;

        ic = (int *)malloc(sizeof(int)*n); if(!ic) goto error;
        ir = (int *)malloc(sizeof(int)*m); if(!ir) goto error;
        ib = (int *)malloc(sizeof(int)*n); if(!ib) goto error;

	uf = (double *)malloc(sizeof(double)*mm);  if(!uf)    goto error; 
	bp = (double *)malloc(sizeof(double)*mm);  if(!bp)    goto error; 
	xp = (double *)malloc(sizeof(double)*mm);  if(!xp)    goto error; 
	t  = (double *)malloc(sizeof(double)*n);   if(!t)     goto error; 
     alfa  = (double *)malloc(sizeof(double)*n);   if(!alfa)  goto error; 
	p  = (double *)malloc(sizeof(double)*mmm); if(!p)     goto error; 
     ginv  = (double *)malloc(sizeof(double)*m*m); if(!ginv)  goto error; 
       vb  = (double *)malloc(sizeof(double)*m);   if(!vb)    goto error; 
	
   l1_d(&m, &n, &mm, &mmm, &nn, ct, f, &prec, &eps, ic, ir, ib, uf, bp,
       xp, t, alfa, p, ginv, vb, irank, iter, r__, a, z__, ind, cancel_func);
       ret = 0;


      error:
       if(ic) free(ic);
       if(ir) free(ir);
       if(ib) free(ib);

       if(uf)   free(uf);
       if(bp)   free(bp);
       if(xp)   free(xp);
       if(t)    free(t);
       if(alfa) free(alfa);     
       if(p)    free(p);
       if(ginv) free(ginv);
       if(vb)   free(vb);


       return ret;

}








/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int pslv_d(integer *id, integer *k, integer *mm, integer *mmm,
	 real *p, real *b, real *x)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j, l;
    doublereal s;
    integer kd;
    doublereal sa, sb;
    integer jj, kk, ll, km1, kkd, kkm1;

/*     THIS SUBROUTINE SOLVES THE SQUARE NON-SINGULAR SYSTEM */
/* OF LINEAR EQUATIONS */
/*                            P*X=B, */
/* OR THE SQUARE NON-SINGULAR SYSTEM OF LINEAR EQUATIONS */
/*                         P(TRANSPOSE)*X=B, */
/* WHERE P IS AN UPPER TRIANGULAR MATRIX, B IS THE RIGHT HAND */
/* SIDE VECTOR AND X IS THE SOLUTION VECTOR. */
/*     THE INPUT DATA TO THE SUBROUTINE. */
/* ID     AN INTEGER INDICATOR SPECIFIED BY THE USER. */
/*        IF ID=1 THE EQUATION P*X=B IS SOLVED. */
/*        IF ID= ANY INTEGER OTHER THAN 1, THE EQUATION */
/*        P(TRANSPOSE)*X=B IS SOLVED. */
/* K      AN INTEGER = THE NUMBER OF EQUATIONS OF THE GIVEN */
/*        SYSTEM. */
/* MM     AN INTEGER GREATER THAN OR EQUAL TO K. */
/* MMM    AN INTEGER = (MM*(MM+3))/2 */
/* P      AN MMM-VECTOR. ITS FIRST (K+1) ELEMENTS CONTAIN THE */
/*        FIRST K ELEMENTS OF ROW 1 OF THE UPPER TRIANGULAR */
/*        MATRIX PLUS AN EXTRA LOCATION TO THE RIGHT. ITS */
/*        NEXT K ELEMENTS CONTAIN THE (K-1) ELEMENTS OF ROW 2 */
/*        OF THE UPPER TRIANGULAR MATRIX PLUS AN EXTRA */
/*        LOCATION TO THE RIGHT, ..., ETC. */
/* B      AN MM-VECTOR. ITS FIRST K ELEMENTS CONTAIN THE */
/*        R.H.S. VECTOR OF THE GIVEN SYSTEM. */
/*     THE OUTPUT OF THE SUBROUTINE. */
/* X      AN MM-VECTOR. ON EXIT, ITS FIRST K ELEMENTS CONTAIN */
/*        THE SOLUTION TO THE GIVEN SYSTEM. */
    /* Parameter adjustments */
    --x;
    --b;
    --p;

    /* Function Body */
    if (*id != 1) {
	goto L3;
    }
/* SOLUTION OF THE UPPER TRIANGULAR SYSTEM. */
    l = *k - 1 + *k * (*k + 1) / 2;
    x[*k] = b[*k] / p[l];
    if (*k == 1) {
	return 0;
    }
    kd = 3;
    km1 = *k - 1;
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *k - i__;
	l -= kd;
	++kd;
	s = b[j];
	ll = l;
	jj = j;
	i__2 = i__;
	for (kk = 1; kk <= i__2; ++kk) {
	    ++ll;
	    sa = p[ll];
	    ++jj;
	    sb = x[jj];
/* L1: */
	    s -= sa * sb;
	}
	x[j] = s / p[l];
/* L2: */
    }
    return 0;
/* SOLUTION OF THE LOWER TRIANGULAR SYSTEM. */
L3:
    x[1] = b[1] / p[1];
    if (*k == 1) {
	return 0;
    }
    l = 1;
    kd = *k + 1;
    i__1 = *k;
    for (i__ = 2; i__ <= i__1; ++i__) {
	l += kd;
	--kd;
	s = b[i__];
	kk = i__;
	kkm1 = i__ - 1;
	kkd = *k;
	i__2 = kkm1;
	for (j = 1; j <= i__2; ++j) {
	    sa = p[kk];
	    sb = x[j];
	    s -= sa * sb;
	    kk += kkd;
	    --kkd;
/* L4: */
	}
/* L5: */
	x[i__] = s / p[l];
    }
    return 0;
} /* pslv_ */

/* Subroutine */ int l1_d(integer *m, integer *n, integer *mm, integer *mmm, 
	integer *nn, real *ct, real *f, real *prec, real *eps, integer *ic, 
	integer *ir, integer *ib, real *uf, real *bp, real *xp, real *t, real 
	*alfa, real *p, real *ginv, real *vb, integer *irank, integer *iter, 
	real *r__, real *a, real *z__, integer *ind, int (*cancel_func)())
{
    /* System generated locals */
    integer ct_dim1, ct_offset, ginv_dim1, ginv_offset, i__1, i__2;

    /* Local variables */
    integer iout;
    extern /* Subroutine */ int pslv_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *);
    real d__, e, g;
    integer i__, j, k, l;
    doublereal s;
    real alpha;
    integer icjin;
    real alfmn, alfmx;
    integer itest, i1;
    real tpeps;
    integer k1;
    real pivot;
    integer irank1, irnkm1;
    real gg;
    integer kd, ii;
    doublereal sa, sb;
    integer li, lj;
    real xb;
    integer kk, kl, iciout;
    real pivotn, pivoto;
    integer icj;
    real bxb;
    integer jin, nmm, ivo;
    real piv;

/*     THIS SUBROUTINE SOLVES AN OVERDETERMINED SYSTEM OF */
/* LINEAR EQUATIONS IN THE L1 NORM BY USING A DUAL SIMPLEX */
/* ALGORITHM TO THE LINEAR PROGRAMMING FORMULATION OF THE */
/* GIVEN PROBLEM. IN THIS ALGORITHM CERTAIN INTERMEDIATE */
/* SIMPLEX ITERATIONS ARE SKIPPED. */
/*     FOR PURPOSE OF NUMERICAL STABILITY, THIS SUBROUTINE */
/* USES A TRIANGULAR DECOMPOSITION TO THE BASIS MATRIX. */
/*     THE SYSTEM OF LINEAR EQUATIONS HAS THE FORM */
/*                          C*A=F, */
/* WHERE C IS A GIVEN REAL N BY M MATRIX OF RANK K.LE.M.LE.N */
/* AND F IS A GIVEN REAL N-VECTOR. */
/*     THE PROBLEM TO BE SOLVED IS TO CALCULATE THE ELEMENTS */
/* OF THE M-VECTOR  A*  WHICH GIVES THE MINIMUM NORM Z. */
/*               Z= ABS(R(1)) + ABS(R(2)) + ... + ABS(R(N)) , */
/* WHERE R(I) IS THE ITH RESIDUAL AND IS GIVEN BY */
/*   R(I)=C(I,1)*A(1)+C(I,2)*A(2)+ ... +C(I,M)*A(M)-F(I). */
/*     SUBROUTINE L1 IS COMPLETELY SELF CONTAINED (CONSISTS */
/* OF TWO SUBROUTINES L1   AND PSLV). */
/*     THE INPUT DATA TO THE SUBROUTINE. */
/* M      THE NUMBER OF COLUMNS OF MATRIX C. */
/* N      THE NUMBER OF ROWS OF MATRIX C. */
/* MM     AN INTEGER GREATER THAN OR EQUAL TO M. */
/* MMM    AN INTEGER = (MM*(MM+3))/2 . */
/* NN     AN INTEGER GREATER THAN OR EQUAL TO N. */
/* CT     A MATRIX OF DIMENSIONS MM BY NN. ON ENTRY, ITS FIRST */
/*        M ROWS AND FIRST N COLUMNS CONTAIN THE TRANSPOSE OF */
/*        MATRIX C IN THE GIVEN SYSTEM C*A=F. */
/* F      AN NN-VECTOR. ON ENTRY, ITS FIRST N ELEMENTS CONTAIN */
/*        THE R.H.S. OF THE GIVEN SYSTEM C*A=F. */
/* PREC   THE ROUND-OFF LEVEL OF THE COMPUTER. FOR THE IBM */
/*        360/67 COMPUTER, PREC IS ABOUT 1.E-6 AND 1.E-16 FOR */
/*        SINGLE AND DOUBLE PRECISION RESPECTIVELY. */
/* EPS    A SPECIFIED TOLERANCE SUCH THAT A CALCULATED NUMBER */
/*        X IS CONSIDERED ZERO IF  ABS(X) <EPS. FOR THE IBM 360/67 */
/*        COMPUTER, EPS IS USUALLY TAKEN 1.E-4 AND 1.E-11 */
/*        RESPECTIVELY. */
/*     THE RESULTS OF THE PROBLEM. */
/* IRANK  THE CALCULATED RANK OF MATRIX C. */
/* ITER   THE NUMBER OF ITERATIONS WHICH THE SOLUTION NEEDED. */
/* A      AN MM-VECTOR. ITS FIRST M ELEMENTS ARE THE SOLUTION */
/*        VECTOR A*. */
/* R      AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE */
/*        THE RESIDUAL VECTOR R=C*A-F. */
/* Z      THE MINIMUM L1 NORM OF THE RESIDUAL VECTOR R. */
/* IND    A RETURN INDICATOR. IND=0 INDICATES THAT THE */
/*        SOLUTION VECTOR A* IS UNIQUE. IND=1, INDICATES THAT */
/*        A* IS MOST PROBABLY NOT UNIQUE. IND=-1 INDICATES */
/*        PREMATURE TERMINATION OF THE CALCULATION DUE TO */
/*        VERY ILL-CONDITIONING OF MATRIX C. */
/*     THE MEANING OF THE OTHER PARAMETERS. */
/* GINV   AN MM-SQUARE MATRIX. ITS FIRST IRANK COLUMNS AND */
/*        FIRST IRANK ROWS CONTAIN THE INVERSE OF THE INITIAL */
/*        BASIS MATRIX AND ITS UPDATE. */
/* VB     AN MM-VECTOR. ITS FIRST IRANK ELEMENTS CONTAIN THE */
/*        INITIAL BASIC SOLUTION AND ITS UPDATE. */
/* T      AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE */
/*        ELEMENTS OF THE ROW IN THE SIMPLEX TABLEAU, THAT */
/*        CORRESPONDS TO THE COLUMN WHICH LEAVES THE BASIS. */
/* ALFA   AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE */
/*        RATIOS:  ALFA(J) = R(J)/T(J). */
/* IC     AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE */
/*        COLUMN INDICES OF MATRIX CT. */
/* IR     AN MM VECTOR. ITS FIRST IRANK ELEMENTS CONTAIN THE */
/*        ROW INDICES OF THE LINEARLY INDEPENDENT ROWS OF CT. */
/* IB     A SIGN NN-VECTOR. ITS FIRST N ELEMENTS HAVE THE */
/*        VALUES +1 OR -1. IB(J)=+1 INDICATES THAT COLUMN J */
/*        IN THE TABLEAU IS AT ITS LOWER BOUND 0. IB(J)=-1 */
/*        INDICATES THAT COLUMN J IS AT ITS UPPER BOUND 2. */
/* P      AN MMM-VECTOR. ITS FIRST ((IRANK*(IRANK+3))/2)-1 */
/*        ELEMENTS CONTAIN THE (IRANK*(IRANK+1))/2 ELEMENTS OF */
/*        THE UPPER TRIANGULAR MATRIX + EXTRA (IRANK-1) WORKING */
/*        LOCATIONS. SEE THE COMMENTS IN SUBROUTINE PSLV. */
/* BP     AN MM-VECTOR. ITS FIRST IRANK ELEMENTS ARE THE */
/*        R.H.S. OF THE TRIANGULAR EQUATIONS AS P*XP=BP. */
/* XP     AN MM-VECTOR . ITS FIRST IRANK ELEMENTS ARE THE */
/*        SOLUTION OF THE TRIANGULAR EQUATIONS AS P*XP=BP. */
/* UF     AN MM-WORKING VECTOR. */
    /* Parameter adjustments */
    --a;
    --vb;
    ginv_dim1 = *mm;
    ginv_offset = ginv_dim1 + 1;
    ginv -= ginv_offset;
    --xp;
    --bp;
    --uf;
    --ir;
    --p;
    --r__;
    --alfa;
    --t;
    --ib;
    --ic;
    --f;
    ct_dim1 = *mm;
    ct_offset = ct_dim1 + 1;
    ct -= ct_offset;

    /* Function Body */
    *ind = 0;
    tpeps = *eps + 2.f;
    nmm = *m * (*m + 3) / 2;
    *irank = *m;
    *iter = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	ib[j] = 1;
/* L1: */
	ic[j] = j;
    }
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	ir[j] = j;
	a[j] = 0.f;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L2: */
	    ginv[i__ + j * ginv_dim1] = 0.f;
	}
/* L3: */
	ginv[j + j * ginv_dim1] = 1.f;
    }
    iout = 0;
/* PART 1 OF THE ALGORITHM. */
L4:
    ++iout;
    if (iout > *irank) {
	goto L16;
    }
    piv = 0.f;
    i__1 = *n;
    for (j = iout; j <= i__1; ++j) {
	icj = ic[j];
	i__2 = *irank;
	for (i__ = iout; i__ <= i__2; ++i__) {
	    d__ = ct[i__ + icj * ct_dim1];
	    if (d__ < 0.f) {
		d__ = -d__;
	    }
	    if (d__ <= piv) {
		goto L5;
	    }
	    li = i__;
	    jin = icj;
	    lj = j;
	    piv = d__;
L5:
	    ;
	}
/* L6: */
    }
/* DETECTION OF RANK DEFICIENCY. */
    if (piv > *eps) {
	goto L7;
    }
    *irank = iout - 1;
    *ind = 1;
    goto L16;
L7:
    if (li == iout) {
	goto L10;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	g = ct[li + j * ct_dim1];
	ct[li + j * ct_dim1] = ct[iout + j * ct_dim1];
/* L8: */
	ct[iout + j * ct_dim1] = g;
    }
    k = ir[li];
    ir[li] = ir[iout];
    ir[iout] = k;
    if (iout == 1) {
	goto L10;
    }
    k = iout - 1;
    i__1 = k;
    for (j = 1; j <= i__1; ++j) {
	d__ = ginv[li + j * ginv_dim1];
	ginv[li + j * ginv_dim1] = ginv[iout + j * ginv_dim1];
/* L9: */
	ginv[iout + j * ginv_dim1] = d__;
    }
/* A GAUSS-JORDAN ELIMINATION STEP. */
L10:
    pivot = ct[iout + jin * ct_dim1];
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L11: */
	ct[iout + j * ct_dim1] /= pivot;
    }
    i__1 = iout;
    for (j = 1; j <= i__1; ++j) {
/* L12: */
	ginv[iout + j * ginv_dim1] /= pivot;
    }
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == iout) {
	    goto L15;
	}
	d__ = ct[i__ + jin * ct_dim1];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L13: */
	    ct[i__ + j * ct_dim1] -= d__ * ct[iout + j * ct_dim1];
	}
	i__2 = iout;
	for (j = 1; j <= i__2; ++j) {
/* L14: */
	    ginv[i__ + j * ginv_dim1] -= d__ * ginv[iout + j * ginv_dim1];
	}
L15:
	;
    }
    ++(*iter);
    k = ic[lj];
    ic[lj] = ic[iout];
    ic[iout] = k;
    goto L4;
/* PART 2 OF THE ALGORITHM. */
L16:
    irank1 = *irank + 1;
    irnkm1 = *irank - 1;
/* INITIAL RESIDUALS AND INITIAL BASIC SOLUTION. */
    i__1 = *irank;
    for (j = 1; j <= i__1; ++j) {
	icj = ic[j];
	r__[icj] = 0.f;
/* L17: */
	uf[j] = f[icj];
    }
    i__1 = *n;
    for (j = irank1; j <= i__1; ++j) {
	icj = ic[j];
	s = -f[icj];
	i__2 = *irank;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sa = uf[i__];
	    sb = ct[i__ + icj * ct_dim1];
/* L18: */
	    s += sa * sb;
	}
	r__[icj] = s;
	if (s >= 0.f) {
	    goto L19;
	}
	ib[icj] = -1;
L19:
	;
    }
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = 0.f;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    sa = ct[i__ + j * ct_dim1];
	    if (ib[j] == -1) {
		sa = -sa;
	    }
/* L20: */
	    s += sa;
	}
/* L21: */
	vb[i__] = s;
    }
/* INITIALIZING THE TRIANGULAR MATRIX. */
    i__1 = nmm;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L22: */
	p[i__] = 0.f;
    }
    k = 1;
    kd = irank1;
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[k] = 1.f;
	k += kd;
/* L23: */
	--kd;
    }
/* DETERMINE THE VECTOR WHICH LEAVES THE BASIS. */
L24:
    ivo = 0;
    pslv_d(&c__1, irank, mm, mmm, &p[1], &vb[1], &xp[1]);
    g = 1.f;
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	e = xp[i__];
	if (e < -(*eps)) {
	    goto L25;
	}
	if (e <= tpeps) {
	    goto L27;
	}
	d__ = 2.f - e;
	if (d__ >= g) {
	    goto L27;
	}
	ivo = 1;
	goto L26;
L25:
	d__ = e;
	if (d__ >= g) {
	    goto L27;
	}
	ivo = -1;
L26:
	g = d__;
	iout = i__;
	xb = e;
L27:
	;
    }
    if (ivo == 0) {
	goto L66;
    }
/* CALCULATION OF ROW (IOUT) IN THE TABLEAU. */
    iciout = ic[iout];
    t[iciout] = 1.f;
    bxb = xb;
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L28: */
	bp[i__] = 0.f;
    }
    bp[iout] = 1.f;
    pslv_d(&c__2, irank, mm, mmm, &p[1], &bp[1], &xp[1]);
    alfmx = 0.f;
    i__1 = *n;
    for (j = irank1; j <= i__1; ++j) {
	icj = ic[j];
	alfa[icj] = 0.f;
	s = 0.f;
	i__2 = *irank;
	for (i__ = iout; i__ <= i__2; ++i__) {
	    sa = xp[i__];
	    sb = ct[i__ + icj * ct_dim1];
/* L29: */
	    s += sa * sb;
	}
	e = s;
	t[icj] = e;
	if (e < *eps && e > -(*eps)) {
	    goto L31;
	}
	d__ = r__[icj];
	if (d__ != 0.f) {
	    goto L30;
	}
	d__ = *prec * *prec * (real) ib[icj];
L30:
	alfa[icj] = d__ / e;
	gg = alfa[icj];
	if (gg < 0.f) {
	    gg = -gg;
	}
	if (gg <= alfmx) {
	    goto L31;
	}
	alfmx = gg;
L31:
	;
    }
    pivoto = 1.f;
    itest = 0;
/* DETERMINE THE VECTOR WHICH ENTERS THE BASIS. */
    gg = alfmx + alfmx;
L32:
	if(cancel_func) {
		if((*cancel_func)()) {
			*ind = -1;
			return 0;
		}
	}

    alfmx = gg;
    alfmn = -gg;
    i__1 = *n;
    for (j = irank1; j <= i__1; ++j) {
	icj = ic[j];
	e = alfa[icj];
	d__ = e * (real) ivo;
	if (d__ <= 0.f) {
	    goto L35;
	}
	if (ivo == 1) {
	    goto L33;
	}
	if (e <= alfmn) {
	    goto L35;
	}
	alfmn = e;
	goto L34;
L33:
	if (e >= alfmx) {
	    goto L35;
	}
	alfmx = e;
L34:
	jin = j;
	itest = 1;
L35:
	;
    }
    if (itest == 1) {
	goto L36;
    }
/* NO FEASIBLE SOLUTION HAS BEEN FOUND. */
    *ind = -1;
    goto L66;
L36:
    icjin = ic[jin];
    pivot = t[icjin];
    alpha = alfa[icjin];
    pivotn = pivot / pivoto;
    if (xb > tpeps) {
	goto L37;
    }
    if (pivotn > 0.f) {
	goto L39;
    }
    goto L41;
L37:
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	e = ct[i__ + iciout * ct_dim1];
/* L38: */
	vb[i__] = vb[i__] - e - e;
    }
    e = t[iciout];
    bxb = bxb - e - e;
    ib[iciout] = -1;
    if (pivotn > 0.f) {
	goto L41;
    }
L39:
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	e = ct[i__ + icjin * ct_dim1];
/* L40: */
	vb[i__] = vb[i__] + e + e;
    }
    e = t[icjin];
    bxb = bxb + e + e;
    ib[icjin] = 1;
L41:
    xb = bxb / pivot;
    if (xb >= -(*eps) && xb <= tpeps) {
	goto L42;
    }
    itest = 0;
L42:
    alfa[icjin] = 0.f;
    if (itest == 1) {
	goto L43;
    }
    pivoto = pivot;
    iciout = icjin;
    goto L32;
L43:
    r__[icjin] = 0.f;
    ++(*iter);
    if (iout == *irank) {
	goto L46;
    }
/* UPDATING MATRIX (P,GINV,VB,CT). */
    i__1 = irnkm1;
    for (j = iout; j <= i__1; ++j) {
	k = j;
	k1 = k + 1;
	kd = *irank;
	l = ic[k];
	ic[k] = ic[k1];
	ic[k1] = l;
	i__2 = k1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p[k] = p[k + 1];
	    k += kd;
/* L44: */
	    --kd;
	}
/* L45: */
    }
L46:
    l = ic[*irank];
    ic[*irank] = ic[jin];
    ic[jin] = l;
    k = *irank;
    kd = *irank;
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[k] = ct[i__ + icjin * ct_dim1];
	k += kd;
/* L47: */
	--kd;
    }
    if (iout == *irank) {
	goto L58;
    }
    i__1 = irnkm1;
    for (i__ = iout; i__ <= i__1; ++i__) {
	ii = i__;
	i1 = i__ + 1;
	k = 0;
	kd = irank1;
	i__2 = ii;
	for (j = 1; j <= i__2; ++j) {
	    k += kd;
/* L48: */
	    --kd;
	}
	kk = k;
	kl = k - kd;
	l = kl;
	g = p[k];
	d__ = p[l];
	if (g < 0.f) {
	    g = -g;
	}
	if (d__ < 0.f) {
	    d__ = -d__;
	}
	if (g <= d__) {
	    goto L53;
	}
	i__2 = *irank;
	for (j = ii; j <= i__2; ++j) {
	    e = p[k];
	    p[k] = p[l];
	    p[l] = e;
	    ++k;
/* L49: */
	    ++l;
	}
	j = ir[i__];
	ir[i__] = ir[i1];
	ir[i1] = j;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    e = ct[i__ + j * ct_dim1];
	    ct[i__ + j * ct_dim1] = ct[i1 + j * ct_dim1];
/* L50: */
	    ct[i1 + j * ct_dim1] = e;
	}
	i__2 = *irank;
	for (j = 1; j <= i__2; ++j) {
	    e = ginv[i__ + j * ginv_dim1];
	    ginv[i__ + j * ginv_dim1] = ginv[i1 + j * ginv_dim1];
/* L51: */
	    ginv[i1 + j * ginv_dim1] = e;
	}
	i__2 = *irank;
	for (j = 1; j <= i__2; ++j) {
	    e = ginv[j + i__ * ginv_dim1];
	    ginv[j + i__ * ginv_dim1] = ginv[j + i1 * ginv_dim1];
/* L52: */
	    ginv[j + i1 * ginv_dim1] = e;
	}
	e = vb[i__];
	vb[i__] = vb[i1];
	vb[i1] = e;
L53:
	e = p[kk] / p[kl];
	if (e < *prec && e > -(*prec)) {
	    goto L57;
	}
	k = kk;
	l = kl;
	i__2 = *irank;
	for (j = ii; j <= i__2; ++j) {
	    p[k] -= e * p[l];
	    ++k;
/* L54: */
	    ++l;
	}
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L55: */
	    ct[i1 + j * ct_dim1] -= e * ct[i__ + j * ct_dim1];
	}
	i__2 = *irank;
	for (j = 1; j <= i__2; ++j) {
/* L56: */
	    ginv[i1 + j * ginv_dim1] -= e * ginv[i__ + j * ginv_dim1];
	}
	vb[i1] -= e * vb[i__];
L57:
	;
    }
L58:
    i__1 = *irank;
    for (j = 1; j <= i__1; ++j) {
	icj = ic[j];
/* L59: */
	uf[j] = f[icj];
    }
    if (alpha > 1.f || alpha < -1.f) {
	goto L61;
    }
    i__1 = *n;
    for (j = irank1; j <= i__1; ++j) {
	icj = ic[j];
/* L60: */
	r__[icj] -= alpha * t[icj];
    }
    goto L64;
L61:
    pslv_d(&c__2, irank, mm, mmm, &p[1], &uf[1], &xp[1]);
    i__1 = *n;
    for (j = irank1; j <= i__1; ++j) {
	icj = ic[j];
	s = -f[icj];
	i__2 = *irank;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sa = xp[i__];
	    sb = ct[i__ + icj * ct_dim1];
/* L62: */
	    s += sa * sb;
	}
/* L63: */
	r__[icj] = s;
    }
L64:
    i__1 = *n;
    for (j = irank1; j <= i__1; ++j) {
	icj = ic[j];
	d__ = r__[icj] * (real) ib[icj];
	if (d__ >= 0.f) {
	    goto L65;
	}
	r__[icj] = 0.f;
L65:
	;
    }
    goto L24;
/* CALCULATING THE ANSWER OF THE PROBLEM. */
L66:
    pslv_d(&c__2, irank, mm, mmm, &p[1], &uf[1], &vb[1]);
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = 0.f;
	i__2 = *irank;
	for (k = 1; k <= i__2; ++k) {
	    sa = vb[k];
	    sb = ginv[k + i__ * ginv_dim1];
/* L67: */
	    s += sa * sb;
	}
	k = ir[i__];
/* L68: */
	a[k] = s;
    }
    s = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sa = r__[j];
	if (sa < 0.f) {
	    sa = -sa;
	}
/* L69: */
	s += sa;
    }
    *z__ = s;
    if (*ind != 0) {
	return 0;
    }
    e = 2.f - *eps;
    i__1 = *irank;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__ = xp[i__];
	if (d__ < *eps || d__ > e) {
	    *ind = 1;
	}
/* L70: */
    }
    return 0;
} /* l1_ */

/*

 ORIGINAL FORTRAM TEXT:

C     THIS PROGRAM IS A DRIVER FOR THE SUBROUTINE  L1
C WHICH SOLVES AN OVERDETERMINED SYSTEM OF LINEAR EQUATIONS
C IN THE L1 NORM, USING A DUAL SIMPLEX METHOD.
C     THE OVERDETERMINED SYSTEM HAS THE FORM  CA=F
C   C IS A GIVEN REAL N BY M MATRIX OF RANK K, K.LE.M.LE.N
C   F IS A GIVEN REAL N-VECTOR.
C   A IS THE SOLUTION M-VECTOR.
C
C
C
      DIMENSION AA(252,27)
      DIMENSION FAY(8,250),EF(250)
      DIMENSION CT(25,250),F(250),R(250),A(25),BINV(25,25)
      DIMENSION P(350),IC(250),IB(250),Y(250),TH(250),IR(25)
      DIMENSION U(25),V(25),W(25),GINV(25,25),VB(25)
      DIMENSION BV(25),IBOUND(250),ICBAS(25),IRBAS(25),ZC(250)
      DIMENSION F1(50),F2(50),BA(6,15)
      INTEGER S(250)
C
  101 FORMAT(10H1EXAMPLE  ,I3,8H SIZE   ,I3,8H  BY    ,I3)
  102 FORMAT(40H0   R.H.S.   OF THE SYSTEM              )
  103 FORMAT(40H0TRANSPOSE OF COEFFICIENT MATRIX        )
  104 FORMAT(17H0EXECUTION TIME =,F12.5,13HSECONDS      )
  105 FORMAT(9H0L1 NORM=,F10.5,6H,RANK=,I4,6H,ITER=,I4,5H,IND=,I4)
  106 FORMAT(40H0THE ANSWER :  A(I)   OR   X(I)         )
  107 FORMAT(40H0RESIDUALS  R(J)    OR E(J)             )
  108 FORMAT(40H0IC(I)                                  )
  109 FORMAT(40H0IR(I)                                  )
  110 FORMAT(1H ,10F12.5)
  111 FORMAT(1H ,20I4)
  115 FORMAT(9H0L1 NORM=,F10.5,6H,RANK=,F3.0,6H,ITER=,F3.0,5H,IND=,F2.0)
  120 FORMAT(1H ,8F15.5)
  125 FORMAT(5F10.6)
  126 FORMAT(6F10.6)
C
      CALL SUNDER
      PREC=1.E-6
      EPS=1.E-4
      TOLER=1.0E-4
      MM=25
      MM2=MM+2
      MMM=(MM*(MM+3))/2
      NN=250
      NN2=NN+2
      IEXMPL=0
    1 IEXMPL=IEXMPL+1
      GO TO (2,4,5,6,7,8,9,13,17,19,20,21,22,26,27,28,29,32,35,100),
     *IEXMPL
C EXAMPLE 1.
    2 N=201
      DX=0.02
      DO 3 I=1,N
      X1=DX*FLOAT(I-1)
      X2=X1*X1
      X3=X2*X1
      X4=X2*X2
      X5=X3*X2
      X6=X3*X3
      X7=X3*X4
      FAY(1,I)=1.
      FAY(2,I)=X1
      FAY(3,I)=X2
      FAY(4,I)=X3
      FAY(5,I)=X4
      FAY(6,I)=X5
      FAY(7,I)=X6
    3 EF(I)= (EXP(-X1))*SIN(X1)
      M=1
      GO TO 10
    4 M=2
      GO TO 10
    5 M=3
      GO TO 10
    6 M=4
      GO TO 10
    7 M=5
      GO TO 10
    8 M=6
      GO TO 10
    9 M=7
   10 DO 12 J=1,N
      F(J)=EF(J)
        DO 11 I=1,M
   11   CT(I,J)=FAY(I,J)
   12 CONTINUE
      GO TO 90
C EXAMPLE 2.
   13 N=100
      N2=50
      M=6
      READ(5,125)(F1(I), I=1,N2)
      READ(5,125)(F2(I), I=1,N2)
      DO 14 J=1,N2
      F(J)=F1(J)
      JN2=J+N2
   14 F(JN2)=F2(J)
      DO 16 J=1,10
      D=0.1*FLOAT(J)
      DD=D*D
        DO 15 I=1,10
        E=0.1*FLOAT(I)
        EE=E*E
        K=10*(J-1)+I
        CT(1,K)=1.0
        CT(2,K)=E
        CT(3,K)=D
        CT(4,K)=E*D
        CT(5,K)=EE
        CT(6,K)=DD
   15   CONTINUE
   16 CONTINUE
      GO TO 90
C EXAMPLE 3.
   17 N=51
      M=8
      DO 18 I=1,N
      D=0.02*FLOAT(I-1)
      DD=D*D
      DDD=D*DD
      E1=D-0.1
      E2=D-0.2
      E3=D-0.4
      E4=D-0.7
      E13=E1*E1*E1
      E23=E2*E2*E2
      E33=E3*E3*E3
      E43=E4*E4*E4
      CT(1,I)=1.0
      CT(2,I)=D
      CT(3,I)=DD
      CT(4,I)=DDD
      IF(E1.GT.0.0) CT(5,I)=E13
      IF(E1.LE.0.0) CT(5,I)=0.
      IF(E2.GT.0.0) CT(6,I)=E23
      IF(E2.LE.0.0) CT(6,I)=0.
      IF(E3.GT.0.0) CT(7,I)=E33
      IF(E3.LE.0.0) CT(7,I)=0.
      IF(E4.GT.0.0) CT(8,I)=E43
      IF(E4.LE.0.0) CT(8,I)=0.
   18 F(I)= SQRT(D)
      GO TO 90
C EXAMPLE 4.
   19 M=6
      PI=4.0*ATAN(1.0)
      X01=PI/2.0
      X02=X01*X01
      X03=X02*X01
      X04=X02*X02
      N=23
      GO TO 23
   20 N=53
      GO TO 23
   21 N=103
      GO TO 23
   22 N=203
   23 N2=N-2
      N3=N-3
      DX=PI/FLOAT(N3)
      DO 24 I=1,N2
      X1=-X01+DX*FLOAT(I-1)
      X2=X1*X1
      X3=X2*X1
      X4=X2*X2
      X5=X3*X2
      CT(1,I)=1.0
      CT(2,I)=X1
      CT(3,I)=X2
      CT(4,I)=X3
      CT(5,I)=X4
      CT(6,I)=X5
   24 F(I)= SIN(X1)
      G=-1.0
      DO 25 I=1,2
      IF(I.EQ.2) G=1.0
      J=N2+I
      CT(1,J)=0.0
      CT(2,J)=1000.0
      CT(3,J)=2000.0*G*X01
      CT(4,J)=3000.0*X02
      CT(5,J)=4000.0*G*X03
      CT(6,J)=5000.0*X04
   25 F(J)=0.0
      GO TO 90
C EXAMPLE 5.
   26 M=11
      D=20.0
      DX=1.0/D
      N=21
      GO TO 30
   27 D=50.0
      DX=1.0/D
      N=51
      GO TO 30
   28 D=100.0
      DX=1.0/D
      N=101
      GO TO 30
   29 D=200.0
      DX=1.0/D
      N=201
   30 DO 31 I=1,N
      X1=DX*FLOAT(I-1)
      X2=X1+X1
      X3=X1+X2
      X4=X2+X2
      X5=X2+X3
      CT(1,I)=1.0
      CT(2,I)= SIN(X1)
      CT(3,I)= COS(X1)
      CT(4,I)= SIN(X2)
      CT(5,I)= COS(X2)
      CT(6,I)= SIN(X3)
      CT(7,I)= COS(X3)
      CT(8,I)= SIN(X4)
      CT(9,I)= COS(X4)
      CT(10,I)= SIN(X5)
      CT(11,I)= COS(X5)
      G= EXP(X1)
      G1= EXP(0.5)
      F(I)=G1
      IF(G.LT.G1) F(I)=G
   31 CONTINUE
      GO TO 90
C EXAMPLE 6.
   32 M=5
      M1=M+1
      N=15
      DO 34 J=1,N
      READ(5,126)(BA(I,J),I=1,M1)
        DO 33 I=1,M
   33   CT(I,J)=BA(I,J)
      F(J)=BA(M1,J)
   34 CONTINUE
      GO TO 90
C EXAMPLE 7.
   35 M=5
      N=51
      DX=0.02
      DO 36 I=1,N
      X1=DX*FLOAT(I-1)
      X2=X1*X1
      X3=X2*X1
      X4=X2*X2
      CT(1,I)=1.0
      CT(2,I)=X1
      CT(3,I)=X2
      CT(4,I)=X3
      CT(5,I)=X4
      GX=0.0
      IF(X1.GT.0.93) GX=5.0
      F(I)=X1*(1.0+X1*(1.0+X1*(1.0+X1)))+GX
   36 CONTINUE
C
   90 WRITE (3,101)IEXMPL,N,M
C     WRITE(3,102)
C     WRITE(3,110)(F(J),J=1,N)
C     WRITE(3,103)
C     DO 93 I=1,M
C  93 WRITE(3,110)(CT(I,J),J=1,N)
      CALL SETCLK

      WRITE(*,1001) IEXMPL, M,N
 1001 FORMAT(10X, 'EXAMPLE', I5, I5,I5)
      WRITE(*,*) ((CT(JJ,II), II=1,N), JJ=1,M)
      WRITE(*,1002)
 1002 FORMAT(10X,'F values:')
      WRITE(*,*) (F(II), II=1,N)

      CALL L1(M,N,MM,MMM,NN,CT,F,PREC,EPS,IC,IR,IB,
     *U,V,W,Y,TH,P,GINV,VB,IRANK,ITER,R,A,Z,IND)
      CALL RDCLK(TX)
      WRITE(3,108)
      WRITE(3,111)(IC(I),I=1,M)
      WRITE(3,109)
      WRITE(3,111)(IR(I),I=1,M)
      WRITE(3,107)
      WRITE(3,110)(R(J),J=1,N)
      WRITE(3,106)
      WRITE(3,120)(A(I),I=1,M)
      WRITE(3,105)Z,IRANK,ITER,IND
      WRITE(3,104)TX
      GO TO 1
  100 STOP
      END
      SUBROUTINE PSLV(ID,K,MM,MMM,P,B,X)                                PSLV0010
C     THIS SUBROUTINE SOLVES THE SQUARE NON-SINGULAR SYSTEM
C OF LINEAR EQUATIONS
C                            P*X=B,
C OR THE SQUARE NON-SINGULAR SYSTEM OF LINEAR EQUATIONS
C                         P(TRANSPOSE)*X=B,
C WHERE P IS AN UPPER TRIANGULAR MATRIX, B IS THE RIGHT HAND
C SIDE VECTOR AND X IS THE SOLUTION VECTOR.
C     THE INPUT DATA TO THE SUBROUTINE.
C ID     AN INTEGER INDICATOR SPECIFIED BY THE USER.
C        IF ID=1 THE EQUATION P*X=B IS SOLVED.
C        IF ID= ANY INTEGER OTHER THAN 1, THE EQUATION
C        P(TRANSPOSE)*X=B IS SOLVED.
C K      AN INTEGER = THE NUMBER OF EQUATIONS OF THE GIVEN
C        SYSTEM.
C MM     AN INTEGER GREATER THAN OR EQUAL TO K.
C MMM    AN INTEGER = (MM*(MM+3))/2
C P      AN MMM-VECTOR. ITS FIRST (K+1) ELEMENTS CONTAIN THE
C        FIRST K ELEMENTS OF ROW 1 OF THE UPPER TRIANGULAR
C        MATRIX PLUS AN EXTRA LOCATION TO THE RIGHT. ITS
C        NEXT K ELEMENTS CONTAIN THE (K-1) ELEMENTS OF ROW 2
C        OF THE UPPER TRIANGULAR MATRIX PLUS AN EXTRA
C        LOCATION TO THE RIGHT, ..., ETC.
C B      AN MM-VECTOR. ITS FIRST K ELEMENTS CONTAIN THE
C        R.H.S. VECTOR OF THE GIVEN SYSTEM.
C     THE OUTPUT OF THE SUBROUTINE.
C X      AN MM-VECTOR. ON EXIT, ITS FIRST K ELEMENTS CONTAIN
C        THE SOLUTION TO THE GIVEN SYSTEM.
      DOUBLE PRECISION S,SA,SB
      DIMENSION P(MMM),B(MM),X(MM)
      IF(ID.NE.1) GO TO 3
C SOLUTION OF THE UPPER TRIANGULAR SYSTEM.
      L=(K-1)+(K*(K+1))/2
      X(K)=B(K)/P(L)
      IF(K.EQ.1) RETURN
      KD=3
      KM1=K-1
      DO 2 I=1,KM1
        J=K-I
        L=L-KD
        KD=KD+1
        S=B(J)
        LL=L
        JJ=J
        DO 1 KK=1,I
          LL=LL+1
          SA=P(LL)
          JJ=JJ+1
          SB=X(JJ)
    1     S=S-SA*SB
        X(J)=S/P(L)
    2   CONTINUE
      RETURN
C SOLUTION OF THE LOWER TRIANGULAR SYSTEM.
    3 X(1)=B(1)/P(1)
      IF(K.EQ.1) RETURN
      L=1
      KD=K+1
      DO 5 I=2,K
        L=L+KD
        KD=KD-1
        S=B(I)
        KK=I
        KKM1=I-1
        KKD=K
        DO 4 J=1,KKM1
          SA=P(KK)
          SB=X(J)
          S=S-SA*SB
          KK=KK+KKD
          KKD=KKD-1
    4     CONTINUE
    5   X(I)=S/P(L)
      RETURN
      END
      SUBROUTINE L1(M,N,MM,MMM,NN,CT,F,PREC,EPS,IC,IR,IB,               L1  0010
     *UF,BP,XP,T,ALFA,P,GINV,VB,IRANK,ITER,R,A,Z,IND)
C     THIS SUBROUTINE SOLVES AN OVERDETERMINED SYSTEM OF
C LINEAR EQUATIONS IN THE L1 NORM BY USING A DUAL SIMPLEX
C ALGORITHM TO THE LINEAR PROGRAMMING FORMULATION OF THE
C GIVEN PROBLEM. IN THIS ALGORITHM CERTAIN INTERMEDIATE
C SIMPLEX ITERATIONS ARE SKIPPED.
C     FOR PURPOSE OF NUMERICAL STABILITY, THIS SUBROUTINE
C USES A TRIANGULAR DECOMPOSITION TO THE BASIS MATRIX.
C     THE SYSTEM OF LINEAR EQUATIONS HAS THE FORM
C                          C*A=F,
C WHERE C IS A GIVEN REAL N BY M MATRIX OF RANK K.LE.M.LE.N
C AND F IS A GIVEN REAL N-VECTOR.
C     THE PROBLEM TO BE SOLVED IS TO CALCULATE THE ELEMENTS
C OF THE M-VECTOR  A*  WHICH GIVES THE MINIMUM NORM Z.
C               Z= ABS(R(1)) + ABS(R(2)) + ... + ABS(R(N)) ,
C WHERE R(I) IS THE ITH RESIDUAL AND IS GIVEN BY
C   R(I)=C(I,1)*A(1)+C(I,2)*A(2)+ ... +C(I,M)*A(M)-F(I).
C     SUBROUTINE L1 IS COMPLETELY SELF CONTAINED (CONSISTS
C OF TWO SUBROUTINES L1   AND PSLV).
C     THE INPUT DATA TO THE SUBROUTINE.
C M      THE NUMBER OF COLUMNS OF MATRIX C.
C N      THE NUMBER OF ROWS OF MATRIX C.
C MM     AN INTEGER GREATER THAN OR EQUAL TO M.
C MMM    AN INTEGER = (MM*(MM+3))/2 .
C NN     AN INTEGER GREATER THAN OR EQUAL TO N.
C CT     A MATRIX OF DIMENSIONS MM BY NN. ON ENTRY, ITS FIRST
C        M ROWS AND FIRST N COLUMNS CONTAIN THE TRANSPOSE OF
C        MATRIX C IN THE GIVEN SYSTEM C*A=F.
C F      AN NN-VECTOR. ON ENTRY, ITS FIRST N ELEMENTS CONTAIN
C        THE R.H.S. OF THE GIVEN SYSTEM C*A=F.
C PREC   THE ROUND-OFF LEVEL OF THE COMPUTER. FOR THE IBM
C        360/67 COMPUTER, PREC IS ABOUT 1.E-6 AND 1.E-16 FOR
C        SINGLE AND DOUBLE PRECISION RESPECTIVELY.
C EPS    A SPECIFIED TOLERANCE SUCH THAT A CALCULATED NUMBER
C        X IS CONSIDERED ZERO IF  ABS(X) <EPS. FOR THE IBM 360/67
C        COMPUTER, EPS IS USUALLY TAKEN 1.E-4 AND 1.E-11
C        RESPECTIVELY.
C     THE RESULTS OF THE PROBLEM.
C IRANK  THE CALCULATED RANK OF MATRIX C.
C ITER   THE NUMBER OF ITERATIONS WHICH THE SOLUTION NEEDED.
C A      AN MM-VECTOR. ITS FIRST M ELEMENTS ARE THE SOLUTION
C        VECTOR A*.
C R      AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE
C        THE RESIDUAL VECTOR R=C*A-F.
C Z      THE MINIMUM L1 NORM OF THE RESIDUAL VECTOR R.
C IND    A RETURN INDICATOR. IND=0 INDICATES THAT THE
C        SOLUTION VECTOR A* IS UNIQUE. IND=1, INDICATES THAT
C        A* IS MOST PROBABLY NOT UNIQUE. IND=-1 INDICATES
C        PREMATURE TERMINATION OF THE CALCULATION DUE TO
C        VERY ILL-CONDITIONING OF MATRIX C.
C     THE MEANING OF THE OTHER PARAMETERS.
C GINV   AN MM-SQUARE MATRIX. ITS FIRST IRANK COLUMNS AND
C        FIRST IRANK ROWS CONTAIN THE INVERSE OF THE INITIAL
C        BASIS MATRIX AND ITS UPDATE.
C VB     AN MM-VECTOR. ITS FIRST IRANK ELEMENTS CONTAIN THE
C        INITIAL BASIC SOLUTION AND ITS UPDATE.
C T      AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE
C        ELEMENTS OF THE ROW IN THE SIMPLEX TABLEAU, THAT
C        CORRESPONDS TO THE COLUMN WHICH LEAVES THE BASIS.
C ALFA   AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE
C        RATIOS:  ALFA(J) = R(J)/T(J).
C IC     AN NN-VECTOR. ITS FIRST N ELEMENTS CONTAIN THE
C        COLUMN INDICES OF MATRIX CT.
C IR     AN MM VECTOR. ITS FIRST IRANK ELEMENTS CONTAIN THE
C        ROW INDICES OF THE LINEARLY INDEPENDENT ROWS OF CT.
C IB     A SIGN NN-VECTOR. ITS FIRST N ELEMENTS HAVE THE
C        VALUES +1 OR -1. IB(J)=+1 INDICATES THAT COLUMN J
C        IN THE TABLEAU IS AT ITS LOWER BOUND 0. IB(J)=-1
C        INDICATES THAT COLUMN J IS AT ITS UPPER BOUND 2.
C P      AN MMM-VECTOR. ITS FIRST ((IRANK*(IRANK+3))/2)-1
C        ELEMENTS CONTAIN THE (IRANK*(IRANK+1))/2 ELEMENTS OF
C        THE UPPER TRIANGULAR MATRIX + EXTRA (IRANK-1) WORKING
C        LOCATIONS. SEE THE COMMENTS IN SUBROUTINE PSLV.
C BP     AN MM-VECTOR. ITS FIRST IRANK ELEMENTS ARE THE
C        R.H.S. OF THE TRIANGULAR EQUATIONS AS P*XP=BP.
C XP     AN MM-VECTOR . ITS FIRST IRANK ELEMENTS ARE THE
C        SOLUTION OF THE TRIANGULAR EQUATIONS AS P*XP=BP.
C UF     AN MM-WORKING VECTOR.
      DOUBLE PRECISION S,SA,SB
      DIMENSION CT(MM,NN),F(NN),A(MM),GINV(MM,MM),P(MMM)
      DIMENSION IC(NN),IB(NN),R(NN),T(NN),ALFA(NN),IR(MM)
      DIMENSION UF(MM),BP(MM),XP(MM),VB(MM)
      IND=0
      TPEPS=2.+EPS
      NMM=(M*(M+3))/2
      IRANK=M
      ITER=0
      DO 1 J=1,N
        IB(J)=1
    1   IC(J)=J
      DO 3 J=1,M
        IR(J)=J
        A(J)=0.
        DO 2 I=1,M
    2     GINV(I,J)=0.
    3   GINV(J,J)=1.
      IOUT=0
C PART 1 OF THE ALGORITHM.
    4 IOUT=IOUT+1
      IF(IOUT.GT.IRANK) GO TO 16
      PIV=0.
      DO 6 J=IOUT,N
        ICJ=IC(J)
        DO 5 I=IOUT,IRANK
          D=CT(I,ICJ)
          IF(D.LT.0.0) D=-D
          IF(D.LE.PIV) GO TO 5
          LI=I
          JIN=ICJ
          LJ=J
          PIV=D
    5     CONTINUE
    6   CONTINUE
C DETECTION OF RANK DEFICIENCY.
      IF(PIV.GT.EPS) GO TO 7
      IRANK=IOUT-1
      IND=1
      GO TO 16
    7 IF(LI.EQ.IOUT) GO TO 10
      DO 8 J=1,N
        G=CT(LI,J)
        CT(LI,J)=CT(IOUT,J)
    8   CT(IOUT,J)=G
      K=IR(LI)
      IR(LI)=IR(IOUT)
      IR(IOUT)=K
      IF(IOUT.EQ.1) GO TO 10
      K=IOUT-1
      DO 9 J=1,K
        D=GINV(LI,J)
        GINV(LI,J)=GINV(IOUT,J)
    9   GINV(IOUT,J)=D
C A GAUSS-JORDAN ELIMINATION STEP.
   10 PIVOT=CT(IOUT,JIN)
      DO 11 J=1,N
   11   CT(IOUT,J)=CT(IOUT,J)/PIVOT
      DO 12 J=1,IOUT
   12   GINV(IOUT,J)=GINV(IOUT,J)/PIVOT
      DO 15 I=1,IRANK
        IF(I.EQ.IOUT) GO TO 15
        D=CT(I,JIN)
        DO 13 J=1,N
   13     CT(I,J)=CT(I,J)-D*CT(IOUT,J)
        DO 14 J=1,IOUT
   14     GINV(I,J)=GINV(I,J)-D*GINV(IOUT,J)
   15   CONTINUE
      ITER=ITER+1
      K=IC(LJ)
      IC(LJ)=IC(IOUT)
      IC(IOUT)=K
      GO TO 4
C PART 2 OF THE ALGORITHM.
   16 IRANK1=IRANK+1
      IRNKM1=IRANK-1
C INITIAL RESIDUALS AND INITIAL BASIC SOLUTION.
      DO 17 J=1,IRANK
        ICJ=IC(J)
        R(ICJ)=0.
   17   UF(J)=F(ICJ)
      DO 19 J=IRANK1,N
        ICJ=IC(J)
        S=-F(ICJ)
        DO 18 I=1,IRANK
          SA=UF(I)
          SB=CT(I,ICJ)
   18     S=S+SA*SB
        R(ICJ)=S
        IF(S.GE.0.0) GO TO 19
        IB(ICJ)=-1
   19   CONTINUE
      DO 21 I=1,IRANK
        S=0.
        DO 20 J=1,N
          SA=CT(I,J)
          IF(IB(J).EQ.(-1)) SA=-SA
   20     S=S+SA
   21   VB(I)=S
C INITIALIZING THE TRIANGULAR MATRIX.
      DO 22 I=1,NMM
   22   P(I)=0.
      K=1
      KD=IRANK1
      DO 23 I=1,IRANK
        P(K)=1.
        K=K+KD
   23   KD=KD-1
C DETERMINE THE VECTOR WHICH LEAVES THE BASIS.
   24 IVO=0
      CALL PSLV(1,IRANK,MM,MMM,P,VB,XP)
      G=1.
      DO 27 I=1,IRANK
        E=XP(I)
        IF(E.LT.(-EPS)) GO TO 25
        IF(E.LE.TPEPS) GO TO 27
        D=2.0-E
        IF(D.GE.G) GO TO 27
        IVO=1
        GO TO 26
   25   D=E
        IF(D.GE.G) GO TO 27
        IVO=-1
   26   G=D
        IOUT=I
        XB=E
   27   CONTINUE
      IF(IVO.EQ.0) GO TO 66
C CALCULATION OF ROW (IOUT) IN THE TABLEAU.
      ICIOUT=IC(IOUT)
      T(ICIOUT)=1.
      BXB=XB
      DO 28 I=1,IRANK
   28   BP(I)=0.
      BP(IOUT)=1.
      CALL PSLV(2,IRANK,MM,MMM,P,BP,XP)
      ALFMX=0.0
      DO 31 J=IRANK1,N
        ICJ=IC(J)
        ALFA(ICJ)=0.
        S=0.
        DO 29 I=IOUT,IRANK
          SA=XP(I)
          SB=CT(I,ICJ)
   29     S=S+SA*SB
        E=S
        T(ICJ)=E
        IF(E.LT.EPS.AND.E.GT.(-EPS)) GO TO 31
        D=R(ICJ)
        IF(D.NE.0.0) GO TO 30
        D=PREC*PREC*FLOAT(IB(ICJ))
   30   ALFA(ICJ)=D/E
        GG=ALFA(ICJ)
        IF(GG.LT.0.0) GG=-GG
        IF(GG.LE.ALFMX) GO TO 31
        ALFMX=GG
   31   CONTINUE
      PIVOTO=1.
      ITEST=0
C DETERMINE THE VECTOR WHICH ENTERS THE BASIS.
      GG=ALFMX+ALFMX
   32 ALFMX=GG
      ALFMN=-GG
      DO 35 J=IRANK1,N
        ICJ=IC(J)
        E=ALFA(ICJ)
        D=E*FLOAT(IVO)
        IF(D.LE.0.0) GO TO 35
        IF(IVO.EQ.1) GO TO 33
        IF(E.LE.ALFMN) GO TO 35
        ALFMN=E
        GO TO 34
   33   IF(E.GE.ALFMX) GO TO 35
        ALFMX=E
   34   JIN=J
        ITEST=1
   35   CONTINUE
      IF(ITEST.EQ.1) GO TO 36
C NO FEASIBLE SOLUTION HAS BEEN FOUND.
      IND=-1
      GO TO 66
   36 ICJIN=IC(JIN)
      PIVOT=T(ICJIN)
      ALPHA=ALFA(ICJIN)
      PIVOTN=PIVOT/PIVOTO
      IF(XB.GT.TPEPS) GO TO 37
      IF(PIVOTN.GT.0.) GO TO 39
      GO TO 41
   37 DO 38 I=1,IRANK
        E=CT(I,ICIOUT)
   38   VB(I)=VB(I)-E-E
      E=T(ICIOUT)
      BXB=BXB-E-E
      IB(ICIOUT)=-1
      IF(PIVOTN.GT.0.) GO TO 41
   39 DO 40 I=1,IRANK
        E=CT(I,ICJIN)
   40   VB(I)=VB(I)+E+E
      E=T(ICJIN)
      BXB=BXB+E+E
      IB(ICJIN)=1
   41 XB=BXB/PIVOT
      IF(XB.GE.(-EPS).AND.XB.LE.TPEPS) GO TO 42
      ITEST=0
   42 ALFA(ICJIN)=0.
      IF(ITEST.EQ.1) GO TO 43
      PIVOTO=PIVOT
      ICIOUT=ICJIN
      GO TO 32
   43 R(ICJIN)=0.
      ITER=ITER+1
      IF(IOUT.EQ.IRANK) GO TO 46
C UPDATING MATRIX (P,GINV,VB,CT).
      DO 45 J=IOUT,IRNKM1
        K=J
        K1=K+1
        KD=IRANK
        L=IC(K)
        IC(K)=IC(K1)
        IC(K1)=L
        DO 44 I=1,K1
          P(K)=P(K+1)
          K=K+KD
   44     KD=KD-1
   45   CONTINUE
   46 L=IC(IRANK)
      IC(IRANK)=IC(JIN)
      IC(JIN)=L
      K=IRANK
      KD=IRANK
      DO 47 I=1,IRANK
        P(K)=CT(I,ICJIN)
        K=K+KD
   47   KD=KD-1
      IF(IOUT.EQ.IRANK) GO TO 58
      DO 57 I=IOUT,IRNKM1
        II=I
        I1=I+1
        K=0
        KD=IRANK1
        DO 48 J=1,II
          K=K+KD
   48     KD=KD-1
        KK=K
        KL=K-KD
        L=KL
        G=P(K)
        D=P(L)
        IF(G.LT.0.0) G=-G
        IF(D.LT.0.0) D=-D
        IF(G.LE.D) GO TO 53
        DO 49 J=II,IRANK
          E=P(K)
          P(K)=P(L)
          P(L)=E
          K=K+1
   49     L=L+1
        J=IR(I)
        IR(I)=IR(I1)
        IR(I1)=J
        DO 50 J=1,N
          E=CT(I,J)
          CT(I,J)=CT(I1,J)
   50     CT(I1,J)=E
        DO 51 J=1,IRANK
          E=GINV(I,J)
          GINV(I,J)=GINV(I1,J)
   51     GINV(I1,J)=E
        DO 52 J=1,IRANK
          E=GINV(J,I)
          GINV(J,I)=GINV(J,I1)
   52     GINV(J,I1)=E
        E=VB(I)
        VB(I)=VB(I1)
        VB(I1)=E
   53   E=P(KK)/P(KL)
        IF(E.LT.PREC.AND.E.GT.(-PREC)) GO TO 57
        K=KK
        L=KL
        DO 54 J=II,IRANK
          P(K)=P(K)-E*P(L)
          K=K+1
   54     L=L+1
        DO 55 J=1,N
   55     CT(I1,J)=CT(I1,J)-E*CT(I,J)
        DO 56 J=1,IRANK
   56     GINV(I1,J)=GINV(I1,J)-E*GINV(I,J)
        VB(I1)=VB(I1)-E*VB(I)
   57   CONTINUE
   58 DO 59 J=1,IRANK
        ICJ=IC(J)
   59   UF(J)=F(ICJ)
      IF(ALPHA.GT.(1.0).OR.ALPHA.LT.(-1.0)) GO TO 61
      DO 60 J=IRANK1,N
        ICJ=IC(J)
   60   R(ICJ)=R(ICJ)-ALPHA*T(ICJ)
      GO TO 64
   61 CALL PSLV(2,IRANK,MM,MMM,P,UF,XP)
      DO 63 J=IRANK1,N
        ICJ=IC(J)
        S=-F(ICJ)
        DO 62 I=1,IRANK
          SA=XP(I)
          SB=CT(I,ICJ)
   62     S=S+SA*SB
   63   R(ICJ)=S
   64 DO 65 J=IRANK1,N
        ICJ=IC(J)
        D=R(ICJ)*FLOAT(IB(ICJ))
        IF(D.GE.0.0) GO TO 65
        R(ICJ)=0.0
   65   CONTINUE
      GO TO 24
C CALCULATING THE ANSWER OF THE PROBLEM.
   66 CALL PSLV(2,IRANK,MM,MMM,P,UF,VB)
      DO 68 I=1,IRANK
        S=0.
        DO 67 K=1,IRANK
          SA=VB(K)
          SB=GINV(K,I)
   67     S=S+SA*SB
        K=IR(I)
   68   A(K)=S
      S=0.
      DO 69 J=1,N
        SA=R(J)
        IF(SA.LT.0.0) SA=-SA
   69   S=S+SA
      Z=S
      IF(IND.NE.0) RETURN
      E=2.-EPS
      DO 70 I=1,IRANK
        D=XP(I)
        IF(D.LT.EPS.OR.D.GT.E) IND=1
   70   CONTINUE
      RETURN
      END
      SUBROUTINE SUNDER                                                 SUN 0010
C     THIS IS A DUMMY SUBROUTINE. ITS FUNCTION IS TO
C SUPPRESS THE UNDERFLOWS OCCURING IN THE CALCULATION.
      RETURN
      END
      SUBROUTINE SETCLK                                                 SET 0010
C     THIS IS A DUMMY SUBROUTINE. ITS FUNCTION IS TO START
C READING THE CPU TIME.
      RETURN
      END
      SUBROUTINE RDCLK(TX)                                              RDC 0010
C     THIS IS A DUMMY SUBROUTINE. ITS FUNCTION IS TO RECORD
C THE CPU TIME  TX  IN SECONDS, WHERE TX IS A REAL VARIABLE.
      TX=0.00
      RETURN
      END


-------------------------------cut here -------------------------------------


TEST DATA:

11   21
  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.  1.
  1.  1.  0.  0.0499791689  0.0998334214  0.149438143  0.198669329  0.247403964
  0.295520216  0.342897803  0.389418334  0.434965551  0.47942555  0.522687256
  0.564642489  0.605186462  0.64421767  0.681638777  0.717356086  0.751280427
  0.783326924  0.813415527  0.841470957  1.  0.998750269  0.995004177
  0.988771081  0.980066597  0.968912423  0.955336511  0.939372718  0.921060979
  0.900447071  0.87758255  0.852524519  0.825335622  0.796083748  0.764842212
  0.731688857  0.696706712  0.659983099  0.621609926  0.581683099  0.540302277
  0.  0.0998334214  0.198669329  0.295520216  0.389418334  0.47942555
  0.564642489  0.64421767  0.717356086  0.783326924  0.841470957  0.891207397
  0.932039082  0.963558197  0.985449731  0.997494996  0.999573588  0.991664827
  0.973847628  0.946300089  0.909297407  1.  0.995004177  0.980066597
  0.955336511  0.921060979  0.87758255  0.825335622  0.764842212  0.696706712
  0.621609926  0.540302277  0.453596085  0.362357706  0.267498761  0.16996716
  0.070737198 -0.0291995462 -0.128844544 -0.227202162 -0.323289543 -0.416146845
  0.  0.149438143  0.295520216  0.434965551  0.564642489  0.681638777
  0.783326924  0.867423177  0.932039082  0.975723386  0.997494996  0.996865034
  0.973847628  0.928959668  0.863209426  0.778073192  0.67546314  0.557683587
  0.427379847  0.287478089  0.141120002  1.  0.988771081  0.955336511
  0.900447071  0.825335622  0.731688857  0.621609926  0.497571081  0.362357706
  0.219006658  0.070737198 -0.0791209862 -0.227202162 -0.370180875 -0.504846036
 -0.628173649 -0.737393796 -0.830053627 -0.904072165 -0.957787216 -0.989992499
  0.  0.198669329  0.389418334  0.564642489  0.717356086  0.841470957
  0.932039082  0.985449731  0.999573588  0.973847628  0.909297407  0.808496356
  0.67546314  0.515501261  0.334988207  0.141120002 -0.0583741926 -0.255541205
 -0.442520559 -0.611857831 -0.756802499  1.  0.980066597  0.921060979
  0.825335622  0.696706712  0.540302277  0.362357706  0.16996716 -0.0291995462
 -0.227202162 -0.416146845 -0.588501155 -0.737393796 -0.856888831 -0.942222297
 -0.989992499 -0.998294771 -0.966798186 -0.896758378 -0.790967762 -0.653643608
  0.  0.247403964  0.47942555  0.681638777  0.841470957  0.948984623
  0.997494996  0.98398596  0.909297407  0.778073192  0.598472118  0.381660998
  0.141120002 -0.108195134 -0.350783229 -0.571561337 -0.756802499 -0.894989371
 -0.977530122 -0.999292791 -0.958924294  1.  0.968912423  0.87758255
  0.731688857  0.540302277  0.315322369  0.070737198 -0.178246051 -0.416146845
 -0.628173649 -0.801143587 -0.924302399 -0.989992499 -0.994129658 -0.93645668
 -0.820559382 -0.653643608 -0.44608748 -0.210795805  0.0376021527  0.2836622
   
  1.  1.05127108  1.10517097  1.16183424  1.22140276  1.28402543  1.34985888
  1.4190675  1.49182475  1.56831217  1.64872122  1.64872122  1.64872122
  1.64872122  1.64872122  1.64872122  1.64872122  1.64872122  1.64872122
  1.64872122  1.64872122

-----------------------------cut here----------------------------------------------


Part of printout from original code (for these data)

1EXAMPLE   14 SIZE    21  BY     11
0IC(I)                                  
    1   9  21  13   6  19   3  16   2  10  11
0IR(I)                                  
    1  10  11   5   8   4   7   9   2   3   6
0RESIDUALS  R(J)    OR E(J)             
      0.00000     0.00301     0.00000    -0.00307    -0.00316     0.00000     0.00428     0.00584     0.00000    -0.01789
     -0.05160    -0.01881     0.00000     0.00679     0.00524     0.00000    -0.00429    -0.00440     0.00000     0.00480
      0.00000
0THE ANSWER :  A(I)   OR   X(I)         
        -5.83387        0.00000        0.00000        2.34378       20.37653        0.00000      -26.16526        0.29254
        17.42299       -0.93676       -4.79977
0L1 NORM=   0.13317,RANK=   8,ITER=  12,IND=   1
0EXECUTION TIME =     0.00000SECONDS      

*/
