/* z29svd.f -- translated by f2c (version of 16 May 1991  13:06:06).
   You must link the resulting object file with the libraries:
	-link <S|C|M|L>f2c.lib   (in that order)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#include "f2c.h"*/

/*typedef long int integer;*/
typedef int integer;
typedef char *address;
typedef short int shortint;
//typedef float real;
typedef double real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

/* Common Block Declarations */

/*struct {
    integer kode;
} koder_;
*/

#define koder_1 koder_

/*      INTERFACE TO INTEGER*4 FUNCTION KEY_PR[C] (TLOC[REFERENCE]) */
/*      END */
/*<    >*/
/* Subroutine */ int dz29svd_(integer *nm, integer *m, integer *n, real *a,
	real *w, logical *matu, real *u, logical *matv, real *v, integer *
	ierr, real *rv1)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2,
	    i__3;
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    static real c, f, g, h;
    static integer i, j, k, l;
    static real s, scale, y, z, x, anorm;
    static integer i1, k1, l1, ii, kk, ll, mn;
    extern doublereal sign_m__(real *, real *);
    static integer its;


/* ******* èOÑèPOÉPAMMA METOÑA CàHÉìã˚PHOÉO PAáãOÜEHà˚ ********** */
/*   (SINGULAR FACTORIZATION OF MATRIX A(M,N)) */
/* CMOTPàTE : ˆOPCAâT Ñ.,MAãúKOãúM M., MOìãEP K. */
/*   " MAﬂàHH¢E METOÑ¢ MATEMATàóECKàX B¢óàCãEHàâ ",M.,"MàP",1980. */


/*# 14 "z29svd.f"*/
/*<       DIMENSION A(NM,N),W(N),U(NM,N),V(N,NM),RV1(N) >*/
/*# 15 "z29svd.f"*/
/*<       LOGICAL MATU,MATV >*/
/*# 16 "z29svd.f"*/
/*<       COMMON /KODER/ KODE >*/
/*# 17 "z29svd.f"*/
/*<       INTEGER*4 KEY_PR >*/
/*# 18 "z29svd.f"*/
/*<       INTEGER*4 TLOC >*/
/*# 19 "z29svd.f"*/
/*<       IERR=0 >*/
    /* Parameter adjustments */
    --rv1;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    u_dim1 = *nm;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --w;
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *ierr = 0;
/*# 20 "z29svd.f"*/
/*<       DO 100 I=1,M >*/
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/*# 21 "z29svd.f"*/
/*<       DO 100 J=1,N >*/
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/*# 22 "z29svd.f"*/
/*<       U(I,J)=A(I,J) >*/
	    u[i + j * u_dim1] = a[i + j * a_dim1];
/*# 23 "z29svd.f"*/
/*< 100   CONTINUE >*/
/* L100: */
	}
    }
/*# 24 "z29svd.f"*/
/*<       G=0.0 >*/
    g = 0.f;
/*# 25 "z29svd.f"*/
/*<       SCALE=0.0 >*/
    scale = 0.f;
/*# 26 "z29svd.f"*/
/*<       ANORM=0.0 >*/
    anorm = 0.f;
/*# 27 "z29svd.f"*/
/*<       DO 300 I=1,N >*/
    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
/*# 28 "z29svd.f"*/
/*<       L=I+1 >*/
	l = i + 1;
/*# 29 "z29svd.f"*/
/*<       RV1(I)=SCALE*G >*/
	rv1[i] = scale * g;
/*# 30 "z29svd.f"*/
/*<       G=0.0 >*/
	g = 0.f;
/*# 31 "z29svd.f"*/
/*<       S=0.0 >*/
	s = 0.f;
/*# 32 "z29svd.f"*/
/*<       SCALE=0.0 >*/
	scale = 0.f;
/*# 33 "z29svd.f"*/
/*<       IF(I.GT.M) GOTO 210 >*/
	if (i > *m) {
	    goto L210;
	}
/*# 34 "z29svd.f"*/
/*<       DO 120 K=I,M >*/
	i__1 = *m;
	for (k = i; k <= i__1; ++k) {
/*# 35 "z29svd.f"*/
/*< 120   SCALE=SCALE+ABS(U(K,I)) >*/
/* L120: */
	    scale += (double)(r__1 = u[k + i * u_dim1], dabs(r__1));
	}
/*# 36 "z29svd.f"*/
/*<       IF(SCALE.EQ.0.0) GOTO 210 >*/
	if (scale == 0.f) {
	    goto L210;
	}
/*# 37 "z29svd.f"*/
/*<       DO 130 K=I,M >*/
	i__1 = *m;
	for (k = i; k <= i__1; ++k) {
/*# 38 "z29svd.f"*/
/*<       U(K,I)=U(K,I)/SCALE >*/
	    u[k + i * u_dim1] /= scale;
/*# 39 "z29svd.f"*/
/*<       S=S+U(K,I)**2 >*/
/* Computing 2nd power */
	    r__1 = u[k + i * u_dim1];
	    s += r__1 * r__1;
/*# 40 "z29svd.f"*/
/*< 130   CONTINUE >*/
/* L130: */
	}
/*# 41 "z29svd.f"*/
/*<       F=U(I,I) >*/
	f = u[i + i * u_dim1];
/*# 42 "z29svd.f"*/
/*<       G=-SIGN_M(SQRT(S),F) >*/
	r__1 = (double)sqrt(s);
	g = -(double)sign_m__(&r__1, &f);
/*# 43 "z29svd.f"*/
/*<       H=F*G-S >*/
	h = f * g - s;
/*# 44 "z29svd.f"*/
/*<       U(I,I)=F-G >*/
	u[i + i * u_dim1] = f - g;
/*# 45 "z29svd.f"*/
/*<       IF(I.EQ.N) GOTO 190 >*/
	if (i == *n) {
	    goto L190;
	}
/*# 46 "z29svd.f"*/
/*<       DO 150 J=L,N >*/
	i__1 = *n;
	for (j = l; j <= i__1; ++j) {
/*# 47 "z29svd.f"*/
/*<       S=0.0 >*/
	    s = 0.f;
/*# 48 "z29svd.f"*/
/*<       DO 140 K=I,M >*/
	    i__3 = *m;
	    for (k = i; k <= i__3; ++k) {
/*# 49 "z29svd.f"*/
/*< 140   S=S+U(K,I)*U(K,J) >*/
/* L140: */
		s += u[k + i * u_dim1] * u[k + j * u_dim1];
	    }
/*# 50 "z29svd.f"*/
/*<       F=S/H >*/
	    f = s / h;
/*# 51 "z29svd.f"*/
/*<       DO 150 K=I,M >*/
	    i__3 = *m;
	    for (k = i; k <= i__3; ++k) {
/*# 52 "z29svd.f"*/
/*<       U(K,J)=U(K,J)+F*U(K,I) >*/
		u[k + j * u_dim1] += f * u[k + i * u_dim1];
/*# 53 "z29svd.f"*/
/*< 150   CONTINUE >*/
/* L150: */
	    }
	}
/*# 54 "z29svd.f"*/
/*< 190   DO 200 K=I,M >*/
L190:
	i__3 = *m;
	for (k = i; k <= i__3; ++k) {
/*# 55 "z29svd.f"*/
/*< 200   U(K,I)=SCALE*U(K,I) >*/
/* L200: */
	    u[k + i * u_dim1] = scale * u[k + i * u_dim1];
	}
/*# 56 "z29svd.f"*/
/*< 210   W(I)=SCALE*G >*/
L210:
	w[i] = scale * g;
/*# 57 "z29svd.f"*/
/*<       G=0.0 >*/
	g = 0.f;
/*# 58 "z29svd.f"*/
/*<       S=0.0 >*/
	s = 0.f;
/*# 59 "z29svd.f"*/
/*<       SCALE=0.0 >*/
	scale = 0.f;
/*# 60 "z29svd.f"*/
/*<       IF(I.GT.M.OR.I.EQ.N) GOTO 290 >*/
	if (i > *m || i == *n) {
	    goto L290;
	}
/*# 61 "z29svd.f"*/
/*<       DO 220 K=L,N >*/
	i__3 = *n;
	for (k = l; k <= i__3; ++k) {
/*# 62 "z29svd.f"*/
/*< 220   SCALE=SCALE+ABS(U(I,K)) >*/
/* L220: */
	    scale += (double)(r__1 = u[i + k * u_dim1], dabs(r__1));
	}
/*# 63 "z29svd.f"*/
/*<       IF(SCALE.EQ.0.0) GOTO 290 >*/
	if (scale == 0.f) {
	    goto L290;
	}
/*# 64 "z29svd.f"*/
/*<       DO 230 K=L,N >*/
	i__3 = *n;
	for (k = l; k <= i__3; ++k) {
/*# 65 "z29svd.f"*/
/*<       U(I,K)=U(I,K)/SCALE >*/
	    u[i + k * u_dim1] /= scale;
/*# 66 "z29svd.f"*/
/*<       S=S+U(I,K)**2 >*/
/* Computing 2nd power */
	    r__1 = u[i + k * u_dim1];
	    s += r__1 * r__1;
/*# 67 "z29svd.f"*/
/*< 230   CONTINUE >*/
/* L230: */
	}
/*# 68 "z29svd.f"*/
/*<       F=U(I,L) >*/
	f = u[i + l * u_dim1];
/*# 69 "z29svd.f"*/
/*<       G=-SIGN_M(SQRT(S),F) >*/
	r__1 = (double)sqrt(s);
	g = -(double)sign_m__(&r__1, &f);
/*# 70 "z29svd.f"*/
/*<       H=F*G-S >*/
	h = f * g - s;
/*# 71 "z29svd.f"*/
/*<       U(I,L)=F-G >*/
	u[i + l * u_dim1] = f - g;
/*# 72 "z29svd.f"*/
/*<       DO 240 K=L,N >*/
	i__3 = *n;
	for (k = l; k <= i__3; ++k) {
/*# 73 "z29svd.f"*/
/*< 240   RV1(K)=U(I,K)/H >*/
/* L240: */
	    rv1[k] = u[i + k * u_dim1] / h;
	}
/*# 74 "z29svd.f"*/
/*<       IF(I.EQ.M) GOTO 270 >*/
	if (i == *m) {
	    goto L270;
	}
/*# 75 "z29svd.f"*/
/*<       DO 260 J=L,M >*/
	i__3 = *m;
	for (j = l; j <= i__3; ++j) {
/*# 76 "z29svd.f"*/
/*<       S=0.0 >*/
	    s = 0.f;
/*# 77 "z29svd.f"*/
/*<       DO 250 K=L,N >*/
	    i__1 = *n;
	    for (k = l; k <= i__1; ++k) {
/*# 78 "z29svd.f"*/
/*< 250   S=S+U(J,K)*U(I,K) >*/
/* L250: */
		s += u[j + k * u_dim1] * u[i + k * u_dim1];
	    }
/*# 79 "z29svd.f"*/
/*<       DO 260 K=L,N >*/
	    i__1 = *n;
	    for (k = l; k <= i__1; ++k) {
/*# 80 "z29svd.f"*/
/*<       U(J,K)=U(J,K)+S*RV1(K) >*/
		u[j + k * u_dim1] += s * rv1[k];
/*# 81 "z29svd.f"*/
/*<   260 CONTINUE >*/
/* L260: */
	    }
	}
/*# 82 "z29svd.f"*/
/*< 270   DO 280 K=L,N >*/
L270:
	i__1 = *n;
	for (k = l; k <= i__1; ++k) {
/*# 83 "z29svd.f"*/
/*< 280   U(I,K)=SCALE*U(I,K) >*/
/* L280: */
	    u[i + k * u_dim1] = scale * u[i + k * u_dim1];
	}
/*# 84 "z29svd.f"*/
/*< 290   ANORM=AMAX1(ANORM,ABS(W(I))+ABS(RV1(I))) >*/
L290:
/* Computing MAX */
	r__3 = anorm, r__4 = (r__1 = w[i], (double)dabs(r__1)) + (r__2 = rv1[i], (double)dabs(
		r__2));
	anorm = (double)dmax(r__3,r__4);
/*# 85 "z29svd.f"*/
/*< 300   CONTINUE >*/
/* L300: */
    }
/*# 86 "z29svd.f"*/
/*<       IF(.NOT.MATV) GOTO 410 >*/
    if (! (*matv)) {
	goto L410;
    }
/*# 87 "z29svd.f"*/
/*<       DO 400 II=1,N >*/
    i__2 = *n;
    for (ii = 1; ii <= i__2; ++ii) {
/*# 88 "z29svd.f"*/
/*<       I=N+1-II >*/
	i = *n + 1 - ii;
/*# 89 "z29svd.f"*/
/*<       IF(I.EQ.N) GOTO 390 >*/
	if (i == *n) {
	    goto L390;
	}
/*# 90 "z29svd.f"*/
/*<       IF(G.EQ.0.0) GOTO 360 >*/
	if (g == 0.f) {
	    goto L360;
	}
/*# 91 "z29svd.f"*/
/*<       DO 320 J=L,N >*/
	i__1 = *n;
	for (j = l; j <= i__1; ++j) {
/*# 92 "z29svd.f"*/
/*< 320   V(J,I)=(U(I,J)/U(I,L))/G >*/
/* L320: */
	    v[j + i * v_dim1] = u[i + j * u_dim1] / u[i + l * u_dim1] / g;
	}
/*# 93 "z29svd.f"*/
/*<       DO 350 J=L,N >*/
	i__1 = *n;
	for (j = l; j <= i__1; ++j) {
/*# 94 "z29svd.f"*/
/*<       S=0.0 >*/
	    s = 0.f;
/*# 95 "z29svd.f"*/
/*<       DO 340 K=L,N >*/
	    i__3 = *n;
	    for (k = l; k <= i__3; ++k) {
/*# 96 "z29svd.f"*/
/*< 340   S=S+U(I,K)*V(K,J) >*/
/* L340: */
		s += u[i + k * u_dim1] * v[k + j * v_dim1];
	    }
/*# 97 "z29svd.f"*/
/*<       DO 350 K=L,N >*/
	    i__3 = *n;
	    for (k = l; k <= i__3; ++k) {
/*# 98 "z29svd.f"*/
/*<       V(K,J)=V(K,J)+S*V(K,I) >*/
		v[k + j * v_dim1] += s * v[k + i * v_dim1];
/*# 99 "z29svd.f"*/
/*< 350   CONTINUE >*/
/* L350: */
	    }
	}
/*# 100 "z29svd.f"*/
/*< 360   DO 380 J=L,N >*/
L360:
	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
/*# 101 "z29svd.f"*/
/*<       V(I,J)=0.0 >*/
	    v[i + j * v_dim1] = 0.f;
/*# 102 "z29svd.f"*/
/*<       V(J,I)=0.0 >*/
	    v[j + i * v_dim1] = 0.f;
/*# 103 "z29svd.f"*/
/*< 380   CONTINUE >*/
/* L380: */
	}
/*# 104 "z29svd.f"*/
/*< 390   V(I,I)=1.0 >*/
L390:
	v[i + i * v_dim1] = 1.f;
/*# 105 "z29svd.f"*/
/*<       G=RV1(I) >*/
	g = rv1[i];
/*# 106 "z29svd.f"*/
/*<       L=I >*/
	l = i;
/*# 107 "z29svd.f"*/
/*< 400   CONTINUE >*/
/* L400: */
    }
/*# 108 "z29svd.f"*/
/*< 410   IF(.NOT.MATU)GOTO 510 >*/
L410:
    if (! (*matu)) {
	goto L510;
    }
/*# 109 "z29svd.f"*/
/*<       MN=N >*/
    mn = *n;
/*# 110 "z29svd.f"*/
/*<       IF(M.LT.N) MN=M >*/
    if (*m < *n) {
	mn = *m;
    }
/*# 111 "z29svd.f"*/
/*<       DO 500 II=1,MN >*/
    i__2 = mn;
    for (ii = 1; ii <= i__2; ++ii) {
/*# 112 "z29svd.f"*/
/*<       I=MN+1-II >*/
	i = mn + 1 - ii;
/*# 113 "z29svd.f"*/
/*<       L=I+1 >*/
	l = i + 1;
/*# 114 "z29svd.f"*/
/*<       G=W(I) >*/
	g = w[i];
/*# 115 "z29svd.f"*/
/*<       IF(I.EQ.N) GOTO 430 >*/
	if (i == *n) {
	    goto L430;
	}
/*# 116 "z29svd.f"*/
/*<       DO 420 J=L,N >*/
	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
/*# 117 "z29svd.f"*/
/*< 420   U(I,J)=0.0 >*/
/* L420: */
	    u[i + j * u_dim1] = 0.f;
	}
/*# 118 "z29svd.f"*/
/*< 430   IF(G.EQ.0.0) GOTO 475 >*/
L430:
	if (g == 0.f) {
	    goto L475;
	}
/*# 119 "z29svd.f"*/
/*<       IF(I.EQ.MN) GOTO 460 >*/
	if (i == mn) {
	    goto L460;
	}
/*# 120 "z29svd.f"*/
/*<       DO 450 J=L,N >*/
	i__3 = *n;
	for (j = l; j <= i__3; ++j) {
/*# 121 "z29svd.f"*/
/*<       S=0.0 >*/
	    s = 0.f;
/*# 122 "z29svd.f"*/
/*<       DO 440 K=L,M >*/
	    i__1 = *m;
	    for (k = l; k <= i__1; ++k) {
/*# 123 "z29svd.f"*/
/*< 440   S=S+U(K,I)*U(K,J) >*/
/* L440: */
		s += u[k + i * u_dim1] * u[k + j * u_dim1];
	    }
/*# 124 "z29svd.f"*/
/*<       F=(S/U(I,I))/G >*/
	    f = s / u[i + i * u_dim1] / g;
/*# 125 "z29svd.f"*/
/*<       DO 450 K=I,M >*/
	    i__1 = *m;
	    for (k = i; k <= i__1; ++k) {
/*# 126 "z29svd.f"*/
/*<       U(K,J)=U(K,J)+F*U(K,I) >*/
		u[k + j * u_dim1] += f * u[k + i * u_dim1];
/*# 127 "z29svd.f"*/
/*< 450   CONTINUE >*/
/* L450: */
	    }
	}
/*# 128 "z29svd.f"*/
/*< 460   DO 470 J=I,M >*/
L460:
	i__1 = *m;
	for (j = i; j <= i__1; ++j) {
/*# 129 "z29svd.f"*/
/*< 470   U(J,I)=U(J,I)/G >*/
/* L470: */
	    u[j + i * u_dim1] /= g;
	}
/*# 130 "z29svd.f"*/
/*<       GOTO 490 >*/
	goto L490;
/*# 131 "z29svd.f"*/
/*< 475   DO 480 J=I,M >*/
L475:
	i__1 = *m;
	for (j = i; j <= i__1; ++j) {
/*# 132 "z29svd.f"*/
/*< 480   U(J,I)=0.0 >*/
/* L480: */
	    u[j + i * u_dim1] = 0.f;
	}
/*# 133 "z29svd.f"*/
/*< 490   U(I,I)=U(I,I)+1.0 >*/
L490:
	u[i + i * u_dim1] += 1.f;
/*# 134 "z29svd.f"*/
/*< 500   CONTINUE >*/
/* L500: */
    }
/*# 135 "z29svd.f"*/
/*<   510 DO 700 KK=1,N >*/
L510:
    i__2 = *n;
    for (kk = 1; kk <= i__2; ++kk) {
/*# 136 "z29svd.f"*/
/*<       K1=N-KK >*/
	k1 = *n - kk;
/*# 137 "z29svd.f"*/
/*<       K=K1+1 >*/
	k = k1 + 1;
/*# 138 "z29svd.f"*/
/*<       ITS=0 >*/
	its = 0;
/*# 139 "z29svd.f"*/
/*< 520   DO 530 LL=1,K >*/
L520:
	i__1 = k;
	for (ll = 1; ll <= i__1; ++ll) {
/*# 140 "z29svd.f"*/
/*<       L1=K-LL >*/
	    l1 = k - ll;
/*# 141 "z29svd.f"*/
/*<       L=L1+1 >*/
	    l = l1 + 1;
/*# 142 "z29svd.f"*/
/*<       IF (ABS(RV1(L))+ANORM.EQ.ANORM) GOTO 565 >*/
	    if ((r__1 = rv1[l], dabs(r__1)) + anorm == anorm) {
		goto L565;
	    }
/*# 143 "z29svd.f"*/
/*<       IF (ABS(W(L1))+ANORM.EQ.ANORM) GOTO 540 >*/
	    if ((r__1 = w[l1], dabs(r__1)) + anorm == anorm) {
		goto L540;
	    }
/*# 144 "z29svd.f"*/
/*< 530   CONTINUE >*/
/* L530: */
	}
/*# 145 "z29svd.f"*/
/*< 540   C=0.0 >*/
L540:
	c = 0.f;
/*# 146 "z29svd.f"*/
/*<       S=1.0 >*/
	s = 1.f;
/*# 147 "z29svd.f"*/
/*<       DO 560 I=L,K >*/
	i__1 = k;
	for (i = l; i <= i__1; ++i) {
/*# 148 "z29svd.f"*/
/*<       F=S*RV1(I) >*/
	    f = s * rv1[i];
/*# 149 "z29svd.f"*/
/*<       RV1(I)=C*RV1(I) >*/
	    rv1[i] = c * rv1[i];
/*# 150 "z29svd.f"*/
/*<       IF(ABS(F)+ANORM.EQ.ANORM) GOTO 565 >*/
	    if (dabs(f) + anorm == anorm) {
		goto L565;
	    }
/*# 151 "z29svd.f"*/
/*<       G=W(I) >*/
	    g = w[i];
/*# 152 "z29svd.f"*/
/*<       H=SQRT(F*F+G*G) >*/
	    h = (double)sqrt(f * f + g * g);
/*# 153 "z29svd.f"*/
/*<       W(I)=H >*/
	    w[i] = h;
/*# 154 "z29svd.f"*/
/*<       C=G/H >*/
	    c = g / h;
/*# 155 "z29svd.f"*/
/*<       S=-F/H >*/
	    s = -(double)f / h;
/*# 156 "z29svd.f"*/
/*<       IF(.NOT.MATU) GOTO 560 >*/
	    if (! (*matu)) {
		goto L560;
	    }
/*# 157 "z29svd.f"*/
/*<       DO 550 J=1,M >*/
	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
/*# 158 "z29svd.f"*/
/*<       Y=U(J,L1) >*/
		y = u[j + l1 * u_dim1];
/*# 159 "z29svd.f"*/
/*<       Z=U(J,I) >*/
		z = u[j + i * u_dim1];
/*# 160 "z29svd.f"*/
/*<       U(J,L1)=Y*C+Z*S >*/
		u[j + l1 * u_dim1] = y * c + z * s;
/*# 161 "z29svd.f"*/
/*<       U(J,I)=-Y*S+Z*C >*/
		u[j + i * u_dim1] = -(double)y * s + z * c;
/*# 162 "z29svd.f"*/
/*< 550   CONTINUE >*/
/* L550: */
	    }
/*# 163 "z29svd.f"*/
/*< 560   CONTINUE >*/
L560:
	    ;
	}
/*# 164 "z29svd.f"*/
/*< 565   Z=W(K) >*/
L565:
	z = w[k];
/*# 165 "z29svd.f"*/
/*<       IF(L.EQ.K) GOTO 650 >*/
	if (l == k) {
	    goto L650;
	}
/*# 166 "z29svd.f"*/
/*<       IF(ITS.EQ.50) GOTO 650 >*/
	if (its == 50) {
	    goto L650;
	}
/*# 167 "z29svd.f"*/
/*<       ITS=ITS+1 >*/
	++its;

/*      WRITE(*,1002) ITS */
/* 1002  FORMAT('+S',I3) */
/*      call mess('Singular factorizations  step',ITS) */
/*      KEY=KEY_PR(TLOC) */
/*      IF(KEY.EQ.27) THEN */
/* 	KODE=131 */
/*        IERR=K */
/*        return */
/*      ENDIF */
/*# 179 "z29svd.f"*/
/*<       X=W(L) >*/
	x = w[l];
/*# 180 "z29svd.f"*/
/*<       Y=W(K1) >*/
	y = w[k1];
/*# 181 "z29svd.f"*/
/*<       G=RV1(K1) >*/
	g = rv1[k1];
/*# 182 "z29svd.f"*/
/*<       H=RV1(K) >*/
	h = rv1[k];
/*# 183 "z29svd.f"*/
/*<       F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y) >*/
	f = ((y - z) * (y + z) + (g - h) * (g + h)) / (h * 2.f * y);
/*# 184 "z29svd.f"*/
/*<       G=SQRT(F*F+1.0) >*/
	g = (double)sqrt(f * f + 1.f);
/*# 185 "z29svd.f"*/
/*<       F=((X-Z)*(X+Z)+H*(Y/(F+SIGN_M(G,F))-H))/X >*/
	f = (double)((x - z) * (x + z) + h * (y / (f + sign_m__(&g, &f)) - h)) / x;
/*# 186 "z29svd.f"*/
/*<       C=1.0 >*/
	c = 1.f;
/*# 187 "z29svd.f"*/
/*<       S=1.0 >*/
	s = 1.f;
/*# 188 "z29svd.f"*/
/*<       DO 600 I1=L,K1 >*/
	i__1 = k1;
	for (i1 = l; i1 <= i__1; ++i1) {
/*# 189 "z29svd.f"*/
/*<       I=I1+1 >*/
	    i = i1 + 1;
/*# 190 "z29svd.f"*/
/*<       G=RV1(I) >*/
	    g = rv1[i];
/*# 191 "z29svd.f"*/
/*<       Y=W(I) >*/
	    y = w[i];
/*# 192 "z29svd.f"*/
/*<       H=S*G >*/
	    h = s * g;
/*# 193 "z29svd.f"*/
/*<       G=C*G >*/
	    g = c * g;
/*# 194 "z29svd.f"*/
/*<       Z=SQRT(F*F+H*H) >*/
	    z = (double)sqrt(f * f + h * h);
/*# 195 "z29svd.f"*/
/*<       RV1(I1)=Z >*/
	    rv1[i1] = z;
/*# 196 "z29svd.f"*/
/*<       C=F/Z >*/
	    c = f / z;
/*# 197 "z29svd.f"*/
/*<       S=H/Z >*/
	    s = h / z;
/*# 198 "z29svd.f"*/
/*<       F=X*C+G*S >*/
	    f = x * c + g * s;
/*# 199 "z29svd.f"*/
/*<       G=-X*S+G*C >*/
	    g = -(double)x * s + g * c;
/*# 200 "z29svd.f"*/
/*<       H=Y*S >*/
	    h = y * s;
/*# 201 "z29svd.f"*/
/*<       Y=Y*C >*/
	    y *= c;
/*# 202 "z29svd.f"*/
/*<       IF(.NOT.MATV) GOTO 575 >*/
	    if (! (*matv)) {
		goto L575;
	    }
/*# 203 "z29svd.f"*/
/*<       DO 570 J=1,N >*/
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/*# 204 "z29svd.f"*/
/*<       X=V(J,I1) >*/
		x = v[j + i1 * v_dim1];
/*# 205 "z29svd.f"*/
/*<       Z=V(J,I) >*/
		z = v[j + i * v_dim1];
/*# 206 "z29svd.f"*/
/*<       V(J,I1)=X*C+Z*S >*/
		v[j + i1 * v_dim1] = x * c + z * s;
/*# 207 "z29svd.f"*/
/*<       V(J,I)=-X*S+Z*C >*/
		v[j + i * v_dim1] = -(double)x * s + z * c;
/*# 208 "z29svd.f"*/
/*< 570   CONTINUE >*/
/* L570: */
	    }
/*# 209 "z29svd.f"*/
/*< 575   Z=SQRT(F*F+H*H) >*/
L575:
	    z = (double)sqrt(f * f + h * h);
/*# 210 "z29svd.f"*/
/*<       W(I1)=Z >*/
	    w[i1] = z;
/*# 211 "z29svd.f"*/
/*<       IF(Z.EQ.0.0) GOTO 580 >*/
	    if (z == 0.f) {
		goto L580;
	    }
/*# 212 "z29svd.f"*/
/*<       C=F/Z >*/
	    c = f / z;
/*# 213 "z29svd.f"*/
/*<       S=H/Z >*/
	    s = h / z;
/*# 214 "z29svd.f"*/
/*< 580   F=C*G+S*Y >*/
L580:
	    f = c * g + s * y;
/*# 215 "z29svd.f"*/
/*<       X=-S*G+C*Y >*/
	    x = -(double)s * g + c * y;
/*# 216 "z29svd.f"*/
/*<       IF(.NOT.MATU) GOTO 600 >*/
	    if (! (*matu)) {
		goto L600;
	    }
/*# 217 "z29svd.f"*/
/*<       DO 590 J=1,M >*/
	    i__3 = *m;
	    for (j = 1; j <= i__3; ++j) {
/*# 218 "z29svd.f"*/
/*<       Y=U(J,I1) >*/
		y = u[j + i1 * u_dim1];
/*# 219 "z29svd.f"*/
/*<       Z=U(J,I) >*/
		z = u[j + i * u_dim1];
/*# 220 "z29svd.f"*/
/*<       U(J,I1)=Y*C+Z*S >*/
		u[j + i1 * u_dim1] = y * c + z * s;
/*# 221 "z29svd.f"*/
/*<       U(J,I)=-Y*S+Z*C >*/
		u[j + i * u_dim1] = -(double)y * s + z * c;
/*# 222 "z29svd.f"*/
/*<   590 CONTINUE >*/
/* L590: */
	    }
/*# 223 "z29svd.f"*/
/*< 600   CONTINUE >*/
L600:
	    ;
	}
/*# 224 "z29svd.f"*/
/*<       RV1(L)=0.0 >*/
	rv1[l] = 0.f;
/*# 225 "z29svd.f"*/
/*<       RV1(K)=F >*/
	rv1[k] = f;
/*# 226 "z29svd.f"*/
/*<       W(K)=X >*/
	w[k] = x;
/*# 227 "z29svd.f"*/
/*<       GOTO 520 >*/
	goto L520;
/*# 228 "z29svd.f"*/
/*< 650   IF(Z.GE.0.0) GOTO 700 >*/
L650:
	if (z >= 0.f) {
	    goto L700;
	}
/*# 229 "z29svd.f"*/
/*<       W(K)=-Z >*/
	w[k] = -(double)z;
/*# 230 "z29svd.f"*/
/*<       IF(.NOT.MATV) GOTO 700 >*/
	if (! (*matv)) {
	    goto L700;
	}
/*# 231 "z29svd.f"*/
/*<       DO 690 J=1,N >*/
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/*# 232 "z29svd.f"*/
/*< 690   V(J,K)=-V(J,K) >*/
/* L690: */
	    v[j + k * v_dim1] = -(double)v[j + k * v_dim1];
	}
/*# 233 "z29svd.f"*/
/*< 700   CONTINUE >*/
L700:
	;
    }
/*# 234 "z29svd.f"*/
/*<       GOTO 1001 >*/
    goto L1001;
/*# 235 "z29svd.f"*/
/*< 1000  IERR=K >*/
/* L1000: */
    *ierr = k;
/*# 236 "z29svd.f"*/
/*< 1001  CONTINUE >*/
L1001:
/*      WRITE(*,1003) */
/* 1003  FORMAT('+‹‹‹‹') */
/*# 239 "z29svd.f"*/
/*<       RETURN >*/
    return 0;
/*# 240 "z29svd.f"*/
/*<       END >*/
} /* z29svd_ */

/*# 242 "z29svd.f"*/
/*<       REAL*4 FUNCTION SIGN_M(A, B) >*/
doublereal sign_m__(real *a, real *b)
{
    /* System generated locals */
    real ret_val;

/*# 243 "z29svd.f"*/
/*<       IF(B.GE.0.) THEN  >*/
    if (*b >= 0.f) {
/*# 244 "z29svd.f"*/
/*< 		SIGN_M=ABS(A)               >*/
	ret_val = (double)dabs(*a);
/*# 245 "z29svd.f"*/
/*<       ELSE	 >*/
    } else {
/*# 246 "z29svd.f"*/
/*<                SIGN_M=-ABS(A) >*/
	ret_val = -(double)dabs(*a);
/*# 247 "z29svd.f"*/
/*<       ENDIF >*/
    }
/*# 248 "z29svd.f"*/
/*<       END >*/
    return ret_val;
} /* sign_m__ */

