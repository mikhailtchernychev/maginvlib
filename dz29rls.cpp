/* z29rls.f -- translated by f2c (version of 16 May 1991  13:06:06).
   You must link the resulting object file with the libraries:
	-link <S|C|M|L>f2c.lib   (in that order)
*/

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

struct {
    integer kode;
} koder_;

#define koder_1 koder_

/* Subroutine */ int dz29svd_(integer *nm, integer *m, integer *n, real *a,
	real *w, logical *matu, real *u, logical *matv, real *v, integer *
	ierr, real *rv1);

/*<    >*/
/* Subroutine */ int dz29rls_(real *a, real *b, real *x, real *u, real *v,
	real *s, real *rv, real *vu, integer *mn, integer *m, integer *n,
	real *eps, real *epsm, real *dm, real *am, integer *k)
{
    /* Initialized data */

    static logical matv = TRUE_;
    static logical matu = TRUE_;

    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2,
	    i__3;
    real r__1, r__2;

    /* Local variables */
    static real d;
    static integer i, j, l;
    static real w, d1;
    static integer n1, n2;
    extern /* Subroutine */ int dz29svd_(integer *, integer *, integer *, real
	    *, real *, logical *, real *, logical *, real *, integer *, real *
	    );
    static integer ii;
    static real si, sm, sq, spe;

/****************** ≈O‰≈PO‚PAMMA T.¸. ñHOBCKOÎ **************************
***/
/**PEòAET CÍCTEMÙ ÔÍHEÎHÛX ÙPABHEHÍÎ METO‰OM CÍH‚ÙÔñPHO‚O PAÁÔOÂEHÍñ**
C*/
/* ≈POBO‰ÍT PE‚ÙÔñPÍÁA˚ÍÒ METO‰OM ÙCE˘EHÍñ C≈EKTPA CÍH‚ÙÔñPHÛX ˘ÍCEÔC */
/* Í BBE‰EHÍEM ≈APAMETPA PE‚ÙÔñPÍÁA˚ÍÍ  C */
/* ÍC≈OÔ£ÁÙET ≈O‰≈PO‚PAMMÙ ‰Â.˜OPCAÎTA Z29SVD  C */
/* BXO‰HÛE ‰AHHÛE : A - MATPÍ˚A CÍCTEMÛ (A*X=B),PAÁMEP N*M   C */
/*    B - ≈PABAñ ˘ACT£ CÍCTEMÛ,PAÁMEP - M      C */
/*    M,N - PAÁMEPHOCTÍ   C */
/*    MN - MAKCÍMAÔ£HOE ˘ÍCÔO ÍÁ M,N    C */
/*    EPS - ≈APAMETP OT¸PACÛBAHÍñ MAÔÛX CÍH‚. ˘ÍCEÔ   C */
/*    DM - ≈APAMETP,≈O KOTOPOMÙ BÛ¸ÍPAETCñ ≈APAMETP PE‚.     C */
/*    (DM ÍMEET CMÛCÔ MAKCÍMAÔ£HOÎ ‰O≈ÙCTÍMOÎ ‰ÍC≈EPCÍÍ)     C */
/* BÛXO‰HÛE ‰AHHÛE :X - PEòEHÍE CÍCTEMÛ C */
/*    U,V - OPTOHOPMAÔ£HÛE MATPÍ˚Û CÍH‚.PAÁÔOÂEHÍñ    C */
/*    S - PE‚ÙÔñPÍÁOBAHHÛE  CÍH‚.˘ÍCEÔA C */
/*    RV - HEPñ‚ÙÔñPÍÁOBAHHÛE CÍH‚.˘ÍCÔAC */
/*    VU - CÔÙÂE¸HÛÎ MACCÍB      C */
/*    AM - HAÎ‰EHHÛÎ ≈O DM ≈APAMETP PE‚ÙÔñPÍÁA˚ÍÍ     C */
/*    K - HAÎ‰EHHÛÎ (≈O EPS) PAH‚ A(KOÔ.HEHÙÔEBÛX CÍH‚.˘ÍCEÔ)C */
/************************************************************************
***/

/*# 24 "z29rls.f"*/
/*<       DIMENSION A(MN,N),U(MN,N),V(N,MN),B(M),X(N),S(N),RV(N) >*/
/*# 25 "z29rls.f"*/
/*<       DIMENSION VU(MN) >*/
/*# 26 "z29rls.f"*/
/*<       LOGICAL MATV,MATU >*/
/*# 27 "z29rls.f"*/
/*<       COMMON /KODER/ KODE >*/
/*# 28 "z29rls.f"*/
/*<       DATA  MATV/.TRUE./,MATU/.TRUE./ >*/
    /* Parameter adjustments */
    --vu;
    --rv;
    --s;
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    u_dim1 = *mn;
    u_offset = u_dim1 + 1;
    u -= u_offset;
    --x;
    --b;
    a_dim1 = *mn;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
/*# 29 "z29rls.f"*/
/*<       CALL Z29SVD(MN,M,N,A,S,MATU,U,MATV,V,II,RV) >*/
    dz29svd_(mn, m, n, &a[a_offset], &s[1], &matu, &u[u_offset], &matv, &v[
	    v_offset], &ii, &rv[1]);

/*  ‰AÔEE Ù≈OPñ‰O˘ÍBAEM MACCÍB CÍH‚.˘ÍCEÔ S Í COOTBETCTBEHHO */
/*  ≈EPECTABÔñEM CTOÔ¸˚Û U Í V */
/*     CALL PRTMT4(A,M,N) */
/*     CALL PRTMT4(U,M,N) */
/*     CALL PRTMT4(V,M,N) */
/*# 36 "z29rls.f"*/
/*<       IF(II.GT.0) GOTO 170 >*/
    if (ii > 0) {
	goto L170;
    }
/*# 37 "z29rls.f"*/
/*<       IF(N.LE.1) GOTO 4444 >*/
    if (*n <= 1) {
	goto L4444;
    }
/*# 38 "z29rls.f"*/
/*<       N1=N-1 >*/
    n1 = *n - 1;
/*# 39 "z29rls.f"*/
/*<       DO 45 I=1,N1 >*/
    i__1 = n1;
    for (i = 1; i <= i__1; ++i) {
/*# 40 "z29rls.f"*/
/*<       N2=N-I >*/
	n2 = *n - i;
/*# 41 "z29rls.f"*/
/*<       DO 40 J=1,N2 >*/
	i__2 = n2;
	for (j = 1; j <= i__2; ++j) {
/*# 42 "z29rls.f"*/
/*<       IF(ABS(S(J)).GT.ABS(S(J+1)))GOTO 40 >*/
	    if ((r__1 = s[j], dabs(r__1)) > (r__2 = s[j + 1], dabs(r__2))) {
		goto L40;
	    }
/*# 43 "z29rls.f"*/
/*<       SM=S(J) >*/
	    sm = s[j];
/*# 44 "z29rls.f"*/
/*<       S(J)=S(J+1) >*/
	    s[j] = s[j + 1];
/*# 45 "z29rls.f"*/
/*<       S(J+1)=SM >*/
	    s[j + 1] = sm;
/*# 46 "z29rls.f"*/
/*<       DO 35 K=1,N >*/
	    i__3 = *n;
	    for (*k = 1; *k <= i__3; ++(*k)) {
/*# 47 "z29rls.f"*/
/*<    35 VU(K)=V(K,J) >*/
/* L35: */
		vu[*k] = v[*k + j * v_dim1];
	    }
/*# 48 "z29rls.f"*/
/*<       DO 36 K=1,N >*/
	    i__3 = *n;
	    for (*k = 1; *k <= i__3; ++(*k)) {
/*# 49 "z29rls.f"*/
/*<       V(K,J)=V(K,J+1) >*/
		v[*k + j * v_dim1] = v[*k + (j + 1) * v_dim1];
/*# 50 "z29rls.f"*/
/*<    36 V(K,J+1)=VU(K) >*/
/* L36: */
		v[*k + (j + 1) * v_dim1] = vu[*k];
	    }
/*# 51 "z29rls.f"*/
/*<       DO 37 K=1,M >*/
	    i__3 = *m;
	    for (*k = 1; *k <= i__3; ++(*k)) {
/*# 52 "z29rls.f"*/
/*<    37 VU(K)=U(K,J) >*/
/* L37: */
		vu[*k] = u[*k + j * u_dim1];
	    }
/*# 53 "z29rls.f"*/
/*<       DO 38 K=1,M >*/
	    i__3 = *m;
	    for (*k = 1; *k <= i__3; ++(*k)) {
/*# 54 "z29rls.f"*/
/*<       U(K,J)=U(K,J+1) >*/
		u[*k + j * u_dim1] = u[*k + (j + 1) * u_dim1];
/*# 55 "z29rls.f"*/
/*<    38 U(K,J+1)=VU(K) >*/
/* L38: */
		u[*k + (j + 1) * u_dim1] = vu[*k];
	    }
/*# 56 "z29rls.f"*/
/*<    40 CONTINUE >*/
L40:
	    ;
	}
/*# 57 "z29rls.f"*/
/*<    45 CONTINUE >*/
/* L45: */
    }
/*# 58 "z29rls.f"*/
/*< 4444  CONTINUE >*/
L4444:

/*OT¸PACÛBAEM MAÔÛE CÍH‚ÙÔñPHÛE ˘ÍCÔA Í BÛ¸ÍPAEM ≈APAMETP PE‚ÙÔñPÍÁA˚ÍÍ AM
*/
/*     CALL PRTMT4(S,1,N) */
/*# 62 "z29rls.f"*/
/*<       SPE=ABS(S(1))*EPS >*/
    spe = dabs(s[1]) * *eps;
/*# 63 "z29rls.f"*/
/*<       I=0 >*/
    i = 0;
/*# 64 "z29rls.f"*/
/*<       IF(S(1).EQ.0.) GOTO 1000 >*/
    if (s[1] == 0.f) {
	goto L1000;
    }
/*# 65 "z29rls.f"*/
/*<    50 I=I+1 >*/
L50:
    ++i;
/*# 66 "z29rls.f"*/
/*<       SI=S(I) >*/
    si = s[i];
/*     SPE=ABS(S(1))*EPS */
/*# 68 "z29rls.f"*/
/*<       IF(ABS(SI).LE.SPE)GOTO 60 >*/
    if (dabs(si) <= spe) {
	goto L60;
    }
/*# 69 "z29rls.f"*/
/*<       IF(I.LT.N)GOTO 50 >*/
    if (i < *n) {
	goto L50;
    }
/*# 70 "z29rls.f"*/
/*<       K=N >*/
    *k = *n;
/*# 71 "z29rls.f"*/
/*<       GOTO 70 >*/
    goto L70;
/*# 72 "z29rls.f"*/
/*<  60   IF(I.EQ.2) continue >*/
L60:
    if (i == 2) {
    }
/* 60   IF(I.EQ.2) WRITE(6,102) */
/*# 74 "z29rls.f"*/
/*<       K=I-1 >*/
    *k = i - 1;
/*# 75 "z29rls.f"*/
/*<    70 D=0. >*/
L70:
    d = 0.f;
/*      AM=1.E-7 */
/*      AM=1.E-10 */

/* ≈POBEPKA HA MAÔOCT£ —ÔEMEHTOB MATPÍ˚Û V */

/*# 81 "z29rls.f"*/
/*<       DO 180 I=1,K >*/
    i__1 = *k;
    for (i = 1; i <= i__1; ++i) {
/*# 82 "z29rls.f"*/
/*<       DO 180 J=1,N >*/
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/*# 83 "z29rls.f"*/
/*<       IF(ABS(V(J,I)).LT.EPSM) V(J,I)=0. >*/
	    if ((r__1 = v[j + i * v_dim1], dabs(r__1)) < *epsm) {
		v[j + i * v_dim1] = 0.f;
	    }
/*# 84 "z29rls.f"*/
/*<   180 CONTINUE >*/
/* L180: */
	}
    }
/*     IF(EPS.GT.0.)AM=EPS/100. */
/*# 86 "z29rls.f"*/
/*<    71 DO 90 I=1,N >*/
L71:
    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
/*# 87 "z29rls.f"*/
/*<       D1=0. >*/
	d1 = 0.f;
/*# 88 "z29rls.f"*/
/*<       DO 80 J=1,K >*/
	i__1 = *k;
	for (j = 1; j <= i__1; ++j) {
/*# 89 "z29rls.f"*/
/*<       SQ=S(J)+AM/S(J) >*/
	    sq = s[j] + *am / s[j];
/*# 90 "z29rls.f"*/
/*<    80 D1=D1+V(I,J)**2/SQ**2 >*/
/* L80: */
/* Computing 2nd power */
	    r__1 = v[i + j * v_dim1];
/* Computing 2nd power */
	    r__2 = sq;
	    d1 += r__1 * r__1 / (r__2 * r__2);
	}
/*# 91 "z29rls.f"*/
/*<       IF(D1.LE.D)GOTO 90 >*/
	if (d1 <= d) {
	    goto L90;
	}
/*# 92 "z29rls.f"*/
/*<       D=D1 >*/
	d = d1;
/*# 93 "z29rls.f"*/
/*<    90 CONTINUE >*/
L90:
	;
    }
/*# 94 "z29rls.f"*/
/*<       IF(D.LE.DM)GOTO 110 >*/
    if (d <= *dm) {
	goto L110;
    }
/*# 95 "z29rls.f"*/
/*<       AM=2*AM >*/
    *am *= 2;
/*# 96 "z29rls.f"*/
/*<       D=0. >*/
    d = 0.f;
/*# 97 "z29rls.f"*/
/*<       GOTO 71 >*/
    goto L71;
/*# 98 "z29rls.f"*/
/*<   110 CONTINUE >*/
L110:
/*     CALL PRTMT4(S,N,1) */
/* BBO‰ÍM HOBÛE CÍH‚.˘ÍCÔA Í ÁA≈OMÍHAEM CTAPÛE B MACCÍBE RV */

/*# 102 "z29rls.f"*/
/*<       DO 103 I=1,N >*/
    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
/*# 103 "z29rls.f"*/
/*<   103 RV(I)=S(I) >*/
/* L103: */
	rv[i] = s[i];
    }
/*# 104 "z29rls.f"*/
/*<       DO 100 I=1,K >*/
    i__2 = *k;
    for (i = 1; i <= i__2; ++i) {
/*# 105 "z29rls.f"*/
/*<   100 S(I)=S(I)+AM/S(I) >*/
/* L100: */
	s[i] += *am / s[i];
    }

/* HAXO‰ÍM PEòEHÍE CÍCTEMÛ - BEKTOP X */

/*# 109 "z29rls.f"*/
/*<       DO 120 I=1,N >*/
    i__2 = *n;
    for (i = 1; i <= i__2; ++i) {
/*# 110 "z29rls.f"*/
/*<   120 X(I)=0. >*/
/* L120: */
	x[i] = 0.f;
    }
/*# 111 "z29rls.f"*/
/*<       DO 150 J=1,K >*/
    i__2 = *k;
    for (j = 1; j <= i__2; ++j) {
/*# 112 "z29rls.f"*/
/*<       W=0. >*/
	w = 0.f;
/*# 113 "z29rls.f"*/
/*<       DO 130 L=1,M >*/
	i__1 = *m;
	for (l = 1; l <= i__1; ++l) {
/*# 114 "z29rls.f"*/
/*<       W=W+U(L,J)*B(L) >*/
	    w += u[l + j * u_dim1] * b[l];
/*# 115 "z29rls.f"*/
/*<   130 CONTINUE >*/
/* L130: */
	}
/*# 116 "z29rls.f"*/
/*<       W=W/S(J) >*/
	w /= s[j];
/*# 117 "z29rls.f"*/
/*<       DO 140 I=1,N >*/
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
/*# 118 "z29rls.f"*/
/*<       X(I)=X(I)+W*V(I,J) >*/
	    x[i] += w * v[i + j * v_dim1];
/*# 119 "z29rls.f"*/
/*<   140 CONTINUE >*/
/* L140: */
	}
/*# 120 "z29rls.f"*/
/*<   150 CONTINUE >*/
/* L150: */
    }
/*# 121 "z29rls.f"*/
/*<   160 RETURN >*/
/* L160: */
    return 0;
/*# 122 "z29rls.f"*/
/*<   170 CONTINUE >*/
L170:
/*  170 WRITE(6,101) II */
/*# 124 "z29rls.f"*/
/*<    >*/
/* L101: */
/*# 126 "z29rls.f"*/
/*<   102 FORMAT(10X,'BHÍMAHÍE! PEòEHÍE ≈O O‰HOMÙ CÍH‚ÙÔñPHOMÙ ˘ÍCÔÙ.') >*/
/* L102: */
/*# 127 "z29rls.f"*/
/*<       GOTO 1002 >*/
    goto L1002;
/*# 128 "z29rls.f"*/
/*<  1000 CONTINUE >*/
L1000:
/*      WRITE(6,1001)(S(I),I=1,N) */
/*# 130 "z29rls.f"*/
/*<  1 >*/
/* L1001: */
/*# 133 "z29rls.f"*/
/*<       KODE=139 >*/
    koder_1.kode = 139;
/*# 134 "z29rls.f"*/
/*<  1002 RETURN 1 >*/
L1002:
    return 1;
/*# 135 "z29rls.f"*/
/*<       END >*/
} /* z29rls_ */

