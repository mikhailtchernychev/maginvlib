/* dcov.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

// #include "f2c.h"

#include "stdio.h"

#define real    float
#define integer int
typedef long int logical;
typedef double doublereal;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)


typedef int /* int or long int */ ftnlen;
typedef int /* int or long int */ ftnint;

#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef /* Subroutine */ int (*S_fp)(int *, int *, int *, double *, double *, double *, int *, void *);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif


#define TRUE_ (1)
#define FALSE_ (0)


/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;

/* DECK DCOV */
/*<    >*/
/* Subroutine */ int dcov_(S_fp fcn, integer *iopt, integer *m, integer *n, 
	doublereal *x, doublereal *fvec, doublereal *r, integer *ldr, integer 
	*info, doublereal *wa1, doublereal *wa2, doublereal *wa3, doublereal *
	wa4, void *pParam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;

    /* Local variables */
    integer idum;
    logical sing;
    doublereal temp;
    integer nrow, i, j, k, iflag;
    doublereal sigma;
	sigma = 0;
    extern /* Subroutine */ int dfdjc3_(S_fp, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, void *), dqrfac_(integer *, integer *, 
	    doublereal *, integer *, logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal denorm_(integer *, doublereal *);
    extern /* Subroutine */ int dwupdt_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    integer kp1, nm1;

/* ***BEGIN PROLOGUE  DCOV */
/* ***PURPOSE  Calculate the covariance matrix for a nonlinear data */
/*            fitting problem.  It is intended to be used after a */
/*            successful return from either DNLS1 or DNLS1E. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  K1B1 */
/* ***TYPE      DOUBLE PRECISION (SCOV-S, DCOV-D) */
/* ***KEYWORDS  COVARIANCE MATRIX, NONLINEAR DATA FITTING, */
/*             NONLINEAR LEAST SQUARES */
/* ***AUTHOR  Hiebert, K. L., (SNLA) */
/* ***DESCRIPTION */

/*  1. Purpose. */

/*     DCOV calculates the covariance matrix for a nonlinear data */
/*     fitting problem.  It is intended to be used after a */
/*     successful return from either DNLS1 or DNLS1E. DCOV */
/*     and DNLS1 (and DNLS1E) have compatible parameters.  The */
/*     required external subroutine, FCN, is the same */
/*     for all three codes, DCOV, DNLS1, and DNLS1E. */

/*  2. Subroutine and Type Statements. */

/*     SUBROUTINE DCOV(FCN,IOPT,M,N,X,FVEC,R,LDR,INFO, */
/*                     WA1,WA2,WA3,WA4) */
/*     INTEGER IOPT,M,N,LDR,INFO */
/*     DOUBLE PRECISION X(N),FVEC(M),R(LDR,N),WA1(N),WA2(N),WA3(N),WA4(M) 
*/
/*     EXTERNAL FCN */

/*  3. Parameters. All TYPE REAL parameters are DOUBLE PRECISION */

/*      FCN is the name of the user-supplied subroutine which calculates 
*/
/*         the functions.  If the user wants to supply the Jacobian */
/*         (IOPT=2 or 3), then FCN must be written to calculate the */
/*         Jacobian, as well as the functions.  See the explanation */
/*         of the IOPT argument below. */
/*         If the user wants the iterates printed in DNLS1 or DNLS1E, */
/*         then FCN must do the printing.  See the explanation of NPRINT 
*/
/*         in DNLS1 or DNLS1E.  FCN must be declared in an EXTERNAL */
/*         statement in the calling program and should be written as */
/*         follows. */

/*         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC) */
/*         INTEGER IFLAG,LDFJAC,M,N */
/*         DOUBLE PRECISION X(N),FVEC(M) */
/*         ---------- */
/*         FJAC and LDFJAC may be ignored       , if IOPT=1. */
/*         DOUBLE PRECISION FJAC(LDFJAC,N)      , if IOPT=2. */
/*         DOUBLE PRECISION FJAC(N)             , if IOPT=3. */
/*         ---------- */
/*           If IFLAG=0, the values in X and FVEC are available */
/*           for printing in DNLS1 or DNLS1E. */
/*           IFLAG will never be zero when FCN is called by DCOV. */
/*           The values of X and FVEC must not be changed. */
/*         RETURN */
/*         ---------- */
/*           If IFLAG=1, calculate the functions at X and return */
/*           this vector in FVEC. */
/*         RETURN */
/*         ---------- */
/*           If IFLAG=2, calculate the full Jacobian at X and return */
/*           this matrix in FJAC.  Note that IFLAG will never be 2 unless 
*/
/*           IOPT=2.  FVEC contains the function values at X and must */
/*           not be altered.  FJAC(I,J) must be set to the derivative */
/*           of FVEC(I) with respect to X(J). */
/*         RETURN */
/*         ---------- */
/*           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian */
/*           and return this vector in FJAC.  Note that IFLAG will */
/*           never be 3 unless IOPT=3.  FJAC(J) must be set to */
/*           the derivative of FVEC(LDFJAC) with respect to X(J). */
/*         RETURN */
/*         ---------- */
/*         END */


/*         The value of IFLAG should not be changed by FCN unless the */
/*         user wants to terminate execution of DCOV.  In this case, set 
*/
/*         IFLAG to a negative integer. */


/*       IOPT is an input variable which specifies how the Jacobian will 
*/
/*         be calculated.  If IOPT=2 or 3, then the user must supply the 
*/
/*         Jacobian, as well as the function values, through the */
/*         subroutine FCN.  If IOPT=2, the user supplies the full */
/*         Jacobian with one call to FCN.  If IOPT=3, the user supplies */
/*         one row of the Jacobian with each call.  (In this manner, */
/*         storage can be saved because the full Jacobian is not stored.) 
*/
/*         If IOPT=1, the code will approximate the Jacobian by forward */
/*         differencing. */

/*       M is a positive integer input variable set to the number of */
/*         functions. */

/*       N is a positive integer input variable set to the number of */
/*         variables.  N must not exceed M. */

/*       X is an array of length N.  On input X must contain the value */
/*         at which the covariance matrix is to be evaluated.  This is */
/*         usually the value for X returned from a successful run of */
/*         DNLS1 (or DNLS1E).  The value of X will not be changed. */

/*    FVEC is an output array of length M which contains the functions */
/*         evaluated at X. */

/*       R is an output array.  For IOPT=1 and 2, R is an M by N array. */
/*         For IOPT=3, R is an N by N array.  On output, if INFO=1, */
/*         the upper N by N submatrix of R contains the covariance */
/*         matrix evaluated at X. */

/*     LDR is a positive integer input variable which specifies */
/*         the leading dimension of the array R.  For IOPT=1 and 2, */
/*         LDR must not be less than M.  For IOPT=3, LDR must not */
/*         be less than N. */

/*    INFO is an integer output variable.  If the user has terminated */
/*         execution, INFO is set to the (negative) value of IFLAG.  See 
*/
/*         description of FCN. Otherwise, INFO is set as follows. */

/*         INFO = 0 Improper input parameters (M.LE.0 or N.LE.0). */

/*         INFO = 1 Successful return.  The covariance matrix has been */
/*                  calculated and stored in the upper N by N */
/*                  submatrix of R. */

/*         INFO = 2 The Jacobian matrix is singular for the input value */
/*                  of X.  The covariance matrix cannot be calculated. */
/*                  The upper N by N submatrix of R contains the QR */
/*                  factorization of the Jacobian (probably not of */
/*                  interest to the user). */

/* WA1,WA2 are work arrays of length N. */
/* and WA3 */

/*     WA4 is a work array of length M. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DENORM, DFDJC3, DQRFAC, DWUPDT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810522  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Fixed an error message.  (RWC) */
/* ***END PROLOGUE  DCOV */

/*     REVISED 850601-1100 */
/*     REVISED YYMMDD HHMM */

/*<       INTEGER I,IDUM,IFLAG,INFO,IOPT,J,K,KP1,LDR,M,N,NM1,NROW >*/
/*<    >*/
/*<       EXTERNAL FCN >*/
/*<       DOUBLE PRECISION ONE,SIGMA,TEMP,ZERO,DENORM >*/
/*<       LOGICAL SING >*/
/*<       SAVE ZERO, ONE >*/
/*<       DATA ZERO/0.D0/,ONE/1.D0/ >*/
    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    r_dim1 = *ldr;
    r_offset = r_dim1 + 1;
    r -= r_offset;
    --fvec;
    --x;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DCOV */
/*<       SING=.FALSE. >*/
    sing = FALSE_;
/*<       IFLAG=0 >*/
    iflag = 0;
/*<       IF (M.LE.0 .OR. N.LE.0) GO TO 300 >*/
    if (*m <= 0 || *n <= 0) {
	goto L300;
    }

/*     CALCULATE SIGMA = (SUM OF THE SQUARED RESIDUALS) / (M-N) */
/*<       IFLAG=1 >*/
    iflag = 1;
/*<       CALL FCN(IFLAG,M,N,X,FVEC,R,LDR) >*/
    (*fcn)(&iflag, m, n, &x[1], &fvec[1], &r[r_offset], ldr, pParam);
/*<       IF (IFLAG.LT.0) GO TO 300 >*/
    if (iflag < 0) {
	goto L300;
    }
/*<       TEMP=DENORM(M,FVEC) >*/
    temp = denorm_(m, &fvec[1]);
/*<       SIGMA=ONE >*/

/*<       IF (M.NE.N) SIGMA=TEMP*TEMP/(M-N) >*/
    if (*m != *n) {
	sigma = temp * temp / (*m - *n);
    }
	
/*     CALCULATE THE JACOBIAN */
/*<       IF (IOPT.EQ.3) GO TO 200 >*/
    if (*iopt == 3) {
	goto L200;
    }

/*     STORE THE FULL JACOBIAN USING M*N STORAGE */
/*<       IF (IOPT.EQ.1) GO TO 100 >*/
    if (*iopt == 1) {
	goto L100;
    }

/*     USER SUPPLIES THE JACOBIAN */
/*<       IFLAG=2 >*/
    iflag = 2;
/*<       CALL FCN(IFLAG,M,N,X,FVEC,R,LDR) >*/
    (*fcn)(&iflag, m, n, &x[1], &fvec[1], &r[r_offset], ldr, pParam);
/*<       GO TO 110 >*/
    goto L110;
    

/*     CODE APPROXIMATES THE JACOBIAN */
/*< 100   CALL DFDJC3(FCN,M,N,X,FVEC,R,LDR,IFLAG,ZERO,WA4) >*/
L100:
    dfdjc3_((S_fp)fcn, m, n, &x[1], &fvec[1], &r[r_offset], ldr, &iflag, &
	    zero, &wa4[1], pParam);
/*< 110   IF (IFLAG.LT.0) GO TO 300 >*/
L110:
    if (iflag < 0) {
	goto L300;
    }



/*     COMPUTE THE QR DECOMPOSITION */
/*<       CALL DQRFAC(M,N,R,LDR,.FALSE.,IDUM,1,WA1,WA1,WA1) >*/
    dqrfac_(m, n, &r[r_offset], ldr, (logical*)&c__0, &idum, &c__1, &wa1[1], &
	    wa1[1], &wa1[1]);
/*<       DO 120 I=1,N >*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/*< 120   R(I,I)=WA1(I) >*/
/* L120: */
	r[i + i * r_dim1] = wa1[i];
    }
/*<       GO TO 225 >*/
    goto L225;

/*     COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX CALCULATED ONE 
*/
/*     ROW AT A TIME AND STORED IN THE UPPER TRIANGLE OF R. */
/*     ( (Q TRANSPOSE)*FVEC IS ALSO CALCULATED BUT NOT USED.) */
/*< 200   CONTINUE >*/
L200:
/*<       DO 210 J=1,N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<       WA2(J)=ZERO >*/
	wa2[j] = zero;
/*<       DO 205 I=1,N >*/
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
/*<       R(I,J)=ZERO >*/
	    r[i + j * r_dim1] = zero;
/*< 205   CONTINUE >*/
/* L205: */
	}
/*< 210   CONTINUE >*/
/* L210: */
    }
/*<       IFLAG=3 >*/
    iflag = 3;
/*<       DO 220 I=1,M >*/
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/*<       NROW = I >*/
	nrow = i;
/*<       CALL FCN(IFLAG,M,N,X,FVEC,WA1,NROW) >*/
	(*fcn)(&iflag, m, n, &x[1], &fvec[1], &wa1[1], &nrow, pParam);
/*<       IF (IFLAG.LT.0) GO TO 300 >*/
	if (iflag < 0) {
	    goto L300;
	}
/*<       TEMP=FVEC(I) >*/
	temp = fvec[i];
/*<       CALL DWUPDT(N,R,LDR,WA1,WA2,TEMP,WA3,WA4) >*/
	dwupdt_(n, &r[r_offset], ldr, &wa1[1], &wa2[1], &temp, &wa3[1], &wa4[
		1]);
/*< 220   CONTINUE >*/
/* L220: */
    }

/*     CHECK IF R IS SINGULAR. */
/*< 225   CONTINUE >*/
L225:
/*<       DO 230 I=1,N >*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/*<       IF (R(I,I).EQ.ZERO) SING=.TRUE. >*/
	if (r[i + i * r_dim1] == zero) {
	    sing = TRUE_;
	}
/*< 230   CONTINUE >*/
/* L230: */
    }
/*<       IF (SING) GO TO 300 >*/
    if (sing) {
	goto L300;
    }

/*     R IS UPPER TRIANGULAR.  CALCULATE (R TRANSPOSE) INVERSE AND STORE 
*/
/*     IN THE UPPER TRIANGLE OF R. */
/*<       IF (N.EQ.1) GO TO 275 >*/
    if (*n == 1) {
	goto L275;
    }
/*<       NM1=N-1 >*/
    nm1 = *n - 1;
/*<       DO 270 K=1,NM1 >*/
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {

/*     INITIALIZE THE RIGHT-HAND SIDE (WA1(*)) AS THE K-TH COLUMN OF T
HE */
/*     IDENTITY MATRIX. */
/*<       DO 240 I=1,N >*/
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
/*<       WA1(I)=ZERO >*/
	    wa1[i] = zero;
/*< 240   CONTINUE >*/
/* L240: */
	}
/*<       WA1(K)=ONE >*/
	wa1[k] = one;

/*<       R(K,K)=WA1(K)/R(K,K) >*/
	r[k + k * r_dim1] = wa1[k] / r[k + k * r_dim1];
/*<       KP1=K+1 >*/
	kp1 = k + 1;
/*<       DO 260 I=KP1,N >*/
	i__2 = *n;
	for (i = kp1; i <= i__2; ++i) {

/*     SUBTRACT R(K,I-1)*R(I-1,*) FROM THE RIGHT-HAND SIDE, WA1(*)
. */
/*<       DO 250 J=I,N >*/
	    i__3 = *n;
	    for (j = i; j <= i__3; ++j) {
/*<       WA1(J)=WA1(J)-R(K,I-1)*R(I-1,J) >*/
		wa1[j] -= r[k + (i - 1) * r_dim1] * r[i - 1 + j * r_dim1];
/*< 250   CONTINUE >*/
/* L250: */
	    }
/*<       R(K,I)=WA1(I)/R(I,I) >*/
	    r[k + i * r_dim1] = wa1[i] / r[i + i * r_dim1];
/*< 260   CONTINUE >*/
/* L260: */
	}
/*< 270   CONTINUE >*/
/* L270: */
    }
/*< 275   R(N,N)=ONE/R(N,N) >*/
L275:
    r[*n + *n * r_dim1] = one / r[*n + *n * r_dim1];

/*     CALCULATE R-INVERSE * (R TRANSPOSE) INVERSE AND STORE IN THE UPPER 
*/
/*     TRIANGLE OF R. */
/*<       DO 290 I=1,N >*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/*<       DO 290 J=I,N >*/
	i__2 = *n;
	for (j = i; j <= i__2; ++j) {
/*<       TEMP=ZERO >*/
	    temp = zero;
/*<       DO 280 K=J,N >*/
	    i__3 = *n;
	    for (k = j; k <= i__3; ++k) {
/*<       TEMP=TEMP+R(I,K)*R(J,K) >*/
		temp += r[i + k * r_dim1] * r[j + k * r_dim1];
/*< 280   CONTINUE >*/
/* L280: */
	    }
/*<       R(I,J)=TEMP*SIGMA >*/
	    r[i + j * r_dim1] = temp * sigma;
/*< 290   CONTINUE >*/
/* L290: */
	}
    }
/*<       INFO=1 >*/
    *info = 1;

/*< 300   CONTINUE >*/
L300:
/*<       IF (M.LE.0 .OR. N.LE.0) INFO=0 >*/
    if (*m <= 0 || *n <= 0) {
	*info = 0;
    }
/*<       IF (IFLAG.LT.0) INFO=IFLAG >*/
    if (iflag < 0) {
	*info = iflag;
    }
/*<       IF (SING) INFO=2 >*/
    if (sing) {
	*info = 2;
    }
/*      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'DCOV', */
/*     +   'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1) 
*/
/*      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DCOV', */
/*     +   'INVALID INPUT PARAMETER.', 2, 1) */
/*      IF (INFO .EQ. 2) CALL XERMSG ('SLATEC', 'DCOV', */
/*     +   'SINGULAR JACOBIAN MATRIX, COVARIANCE MATRIX CANNOT BE ' // */
/*     +   'CALCULATED.', 1, 1) */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* dcov_ */

