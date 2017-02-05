/* tn_m.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

//#include "f2c.h"

#define MY_F2C

#ifdef MY_F2C
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define real    float
#define integer int
typedef long int logical;
typedef double doublereal;
typedef char *address;


#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)

#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif

#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)


typedef int /* int or long int */ ftnlen;
typedef int /* int or long int */ ftnint;

#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef /* Subroutine */ int (*S_fp)(int *n, double *x, double *f, double *g, void *pParam);
#endif

#define TRUE_ (1)
#define FALSE_ (0)


double pow_dd(doublereal *ap, doublereal *bp)
{
return(pow(*ap, *bp) );
}


#endif


/* Common Block Declarations */

union {
    struct {
	integer lgv, lz1, lzk, lv, lsk, lyk, ldiagb, lsr, lyr, loldg, lhg,
		lhyk, lpk, lemat, lwtest;
    } _1;
    struct {
	integer lgv, lz1, lzk, lv, lsk, lyk, ldiagb, lsr, lyr, lhyr, lhg,
		lhyk, lpk, lemat, lwtest;
    } _2;
    struct {
	integer lsub[14], lwtest;
    } _3;
} subscr_;

#define subscr_1 (subscr_._1)
#define subscr_2 (subscr_._2)
#define subscr_3 (subscr_._3)

/* Table of constant values */

static logical c_false = FALSE_;
static integer c__1 = 1;
static logical c_true = TRUE_;
static doublereal c_b135 = .6666;

/* %% TRUNCATED-NEWTON METHOD:  SUBROUTINES */
/*   FOR OTHER MACHINES, MODIFY ROUTINE MCHPR1 (MACHINE EPSILON) */
/*   WRITTEN BY:  STEPHEN G. NASH */
/*                OPERATIONS RESEARCH AND APPLIED STATISTICS DEPT. */
/*                GEORGE MASON UNIVERSITY */
/*                FAIRFAX, VA 22030 */
/* ****************************************************************** */
/*<       SUBROUTINE TN (IERROR, N, X, F, G, W, LW, SFUN) >*/
/* Subroutine */ int tn_(integer *ierror, integer *n, doublereal *x,
	doublereal *f, doublereal *g, doublereal *w, integer *lw, S_fp sfun, void *pParam)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer nmax;
    extern /* Subroutine */ int lmqn_(integer *, integer *, doublereal *,
	    doublereal *, doublereal *, doublereal *, integer *, S_fp,
	    integer *, integer *, integer *, doublereal *, doublereal *,
	    doublereal *, doublereal *, void *);
    doublereal xtol;
    integer maxit;
    extern doublereal mchpr1_(void);
    doublereal accrcy;
    integer maxfun, msglvl;
    doublereal stepmx, eta;

/*<       IMPLICIT          DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER           IERROR, N, LW >*/
/*<       DOUBLE PRECISION  X(N), G(N), F, W(LW) >*/

/* THIS ROUTINE SOLVES THE OPTIMIZATION PROBLEM */

/*            MINIMIZE F(X) */
/*               X */

/* WHERE X IS A VECTOR OF N REAL VARIABLES.  THE METHOD USED IS */
/* A TRUNCATED-NEWTON ALGORITHM (SEE "NEWTON-TYPE MINIMIZATION VIA */
/* THE LANCZOS METHOD" BY S.G. NASH (SIAM J. NUMER. ANAL. 21 (1984), */
/* PP. 770-778).  THIS ALGORITHM FINDS A LOCAL MINIMUM OF F(X).  IT DOES */
/* NOT ASSUME THAT THE FUNCTION F IS CONVEX (AND SO CANNOT GUARANTEE A */
/* GLOBAL SOLUTION), BUT DOES ASSUME THAT THE FUNCTION IS BOUNDED BELOW. */
/* IT CAN SOLVE PROBLEMS HAVING ANY NUMBER OF VARIABLES, BUT IT IS */
/* ESPECIALLY USEFUL WHEN THE NUMBER OF VARIABLES (N) IS LARGE. */

/* SUBROUTINE PARAMETERS: */

/* IERROR - (INTEGER) ERROR CODE */
/*          ( 0 => NORMAL RETURN) */
/*          ( 2 => MORE THAN MAXFUN EVALUATIONS) */
/*          ( 3 => LINE SEARCH FAILED TO FIND */
/*          (          LOWER POINT (MAY NOT BE SERIOUS) */
/*          (-1 => ERROR IN INPUT PARAMETERS) */
/* N      - (INTEGER) NUMBER OF VARIABLES */
/* X      - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON INPUT, AN INITIAL */
/*          ESTIMATE OF THE SOLUTION; ON OUTPUT, THE COMPUTED SOLUTION. */
/* G      - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON OUTPUT, THE FINAL */
/*          VALUE OF THE GRADIENT */
/* F      - (REAL*8) ON INPUT, A ROUGH ESTIMATE OF THE VALUE OF THE */
/*          OBJECTIVE FUNCTION AT THE SOLUTION; ON OUTPUT, THE VALUE */
/*          OF THE OBJECTIVE FUNCTION AT THE SOLUTION */
/* W      - (REAL*8) WORK VECTOR OF LENGTH AT LEAST 14*N */
/* LW     - (INTEGER) THE DECLARED DIMENSION OF W */
/* SFUN   - A USER-SPECIFIED SUBROUTINE THAT COMPUTES THE FUNCTION */
/*          AND GRADIENT OF THE OBJECTIVE FUNCTION.  IT MUST HAVE */
/*          THE CALLING SEQUENCE */
/*             SUBROUTINE SFUN (N, X, F, G) */
/*             INTEGER           N */
/*             DOUBLE PRECISION  X(N), G(N), F */

/* THIS IS AN EASY-TO-USE DRIVER FOR THE MAIN OPTIMIZATION ROUTINE */
/* LMQN.  MORE EXPERIENCED USERS WHO WISH TO CUSTOMIZE PERFORMANCE */
/* OF THIS ALGORITHM SHOULD CALL LMQN DIRECTLY. */

/* ---------------------------------------------------------------------- */
/* THIS ROUTINE SETS UP ALL THE PARAMETERS FOR THE TRUNCATED-NEWTON */
/* ALGORITHM.  THE PARAMETERS ARE: */

/* ETA    - SEVERITY OF THE LINESEARCH */
/* MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS */
/* XTOL   - DESIRED ACCURACY FOR THE SOLUTION X* */
/* STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH */
/* ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES */
/* MSGLVL - DETERMINES QUANTITY OF PRINTED OUTPUT */
/*          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION. */
/* MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP */

/*<       DOUBLE PRECISION ETA, ACCRCY, XTOL, STEPMX, DSQRT, MCHPR1 >*/
/*<       EXTERNAL         SFUN >*/

/* SET UP PARAMETERS FOR THE OPTIMIZATION ROUTINE */

/*<       MAXIT = N/2 >*/
    /* Parameter adjustments */
    --g;
    --x;
    --w;

    /* Function Body */
    maxit = *n / 2;
/*<       IF (MAXIT .GT. 50) MAXIT = 50 >*/
    if (maxit > 50) {
	maxit = 50;
    }
/*<       IF (MAXIT .LE. 0) MAXIT = 1 >*/
    if (maxit <= 0) {
	maxit = 1;
    }
/*<       MSGLVL = 1 >*/
    msglvl = 1;
/*<       MAXFUN = 150*N >*/
    maxfun = *n * 150;
/*<       ETA = .25D0 >*/
    eta = .25;
/*<       STEPMX = 1.D1 >*/
    stepmx = 10.;
/*<       ACCRCY = 1.D2*MCHPR1() >*/
    accrcy = mchpr1_()* 100.;
/*<       XTOL = DSQRT(ACCRCY) >*/
    xtol = sqrt(accrcy);

/* MINIMIZE THE FUNCTION */

/*<    >*/
    lmqn_(ierror, n, &x[1], f, &g[1], &w[1], lw, (S_fp)sfun, &msglvl, &maxit,
	    &maxfun, &eta, &stepmx, &accrcy, &xtol, pParam);

/* PRINT THE RESULTS */

/*      IF (IERROR .NE. 0) WRITE(*,800) IERROR */
/*      WRITE(*,810) F */
/*<       IF (MSGLVL .LT. 1) RETURN >*/
    if (msglvl < 1) {
	return 0;
    }
/*      WRITE(*,820) */
/*<       NMAX = 10 >*/
    nmax = 10;
/*<       IF (N .LT. NMAX) NMAX = N >*/
    if (*n < nmax) {
	nmax = *n;
    }
/*      WRITE(*,830) (I,X(I),I=1,NMAX) */
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(//,' ERROR CODE =', I3) */
/* 810   FORMAT(//,' OPTIMAL FUNCTION VALUE = ', 1PD22.15) */
/* 820   FORMAT(10X, 'CURRENT SOLUTION IS (AT MOST 10 COMPONENTS)', /, */
/*     *       14X, 'I', 11X, 'X(I)') */
/* 830   FORMAT(10X, I5, 2X, 1PD22.15) */
/*<       END >*/
} /* tn_ */



/*<       SUBROUTINE TNBC (IERROR, N, X, F, G, W, LW, SFUN, LOW, UP, IPIVOT) >*/
/* Subroutine */ int tnbc_(integer *ierror, integer *n, doublereal *x,
	doublereal *f, doublereal *g, doublereal *w, integer *lw, S_fp sfun,
	doublereal *low, doublereal *up, integer *ipivot, void * pParam)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer nmax;
    doublereal xtol;
    integer maxit;
    extern doublereal mchpr1_(void);
    doublereal accrcy;
    extern /* Subroutine */ int lmqnbc_(integer *, integer *, doublereal *,
	    doublereal *, doublereal *, doublereal *, integer *, S_fp,
	    doublereal *, doublereal *, integer *, integer *, integer *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *, void *)
	    ;
    integer maxfun, msglvl;
    doublereal stepmx, eta;

/*<       IMPLICIT          DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER           IERROR, N, LW, IPIVOT(N) >*/
/*<       DOUBLE PRECISION  X(N), G(N), F, W(LW), LOW(N), UP(N) >*/

/* THIS ROUTINE SOLVES THE OPTIMIZATION PROBLEM */

/*   MINIMIZE     F(X) */
/*      X */
/*   SUBJECT TO   LOW <= X <= UP */

/* WHERE X IS A VECTOR OF N REAL VARIABLES.  THE METHOD USED IS */
/* A TRUNCATED-NEWTON ALGORITHM (SEE "NEWTON-TYPE MINIMIZATION VIA */
/* THE LANCZOS ALGORITHM" BY S.G. NASH (TECHNICAL REPORT 378, MATH. */
/* THE LANCZOS METHOD" BY S.G. NASH (SIAM J. NUMER. ANAL. 21 (1984), */
/* PP. 770-778).  THIS ALGORITHM FINDS A LOCAL MINIMUM OF F(X).  IT DOES */
/* NOT ASSUME THAT THE FUNCTION F IS CONVEX (AND SO CANNOT GUARANTEE A */
/* GLOBAL SOLUTION), BUT DOES ASSUME THAT THE FUNCTION IS BOUNDED BELOW. */
/* IT CAN SOLVE PROBLEMS HAVING ANY NUMBER OF VARIABLES, BUT IT IS */
/* ESPECIALLY USEFUL WHEN THE NUMBER OF VARIABLES (N) IS LARGE. */

/* SUBROUTINE PARAMETERS: */

/* IERROR  - (INTEGER) ERROR CODE */
/*           ( 0 => NORMAL RETURN */
/*           ( 2 => MORE THAN MAXFUN EVALUATIONS */
/*           ( 3 => LINE SEARCH FAILED TO FIND LOWER */
/*           (          POINT (MAY NOT BE SERIOUS) */
/*           (-1 => ERROR IN INPUT PARAMETERS */
/* N       - (INTEGER) NUMBER OF VARIABLES */
/* X       - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON INPUT, AN INITIAL */
/*           ESTIMATE OF THE SOLUTION; ON OUTPUT, THE COMPUTED SOLUTION. */
/* G       - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON OUTPUT, THE FINAL */
/*           VALUE OF THE GRADIENT */
/* F       - (REAL*8) ON INPUT, A ROUGH ESTIMATE OF THE VALUE OF THE */
/*           OBJECTIVE FUNCTION AT THE SOLUTION; ON OUTPUT, THE VALUE */
/*           OF THE OBJECTIVE FUNCTION AT THE SOLUTION */
/* W       - (REAL*8) WORK VECTOR OF LENGTH AT LEAST 14*N */
/* LW      - (INTEGER) THE DECLARED DIMENSION OF W */
/* SFUN    - A USER-SPECIFIED SUBROUTINE THAT COMPUTES THE FUNCTION */
/*           AND GRADIENT OF THE OBJECTIVE FUNCTION.  IT MUST HAVE */
/*           THE CALLING SEQUENCE */
/*             SUBROUTINE SFUN (N, X, F, G) */
/*             INTEGER           N */
/*             DOUBLE PRECISION  X(N), G(N), F */
/* LOW, UP - (REAL*8) VECTORS OF LENGTH AT LEAST N CONTAINING */
/*           THE LOWER AND UPPER BOUNDS ON THE VARIABLES.  IF */
/*           THERE ARE NO BOUNDS ON A PARTICULAR VARIABLE, SET */
/*           THE BOUNDS TO -1.D38 AND 1.D38, RESPECTIVELY. */
/* IPIVOT  - (INTEGER) WORK VECTOR OF LENGTH AT LEAST N, USED */
/*           TO RECORD WHICH VARIABLES ARE AT THEIR BOUNDS. */

/* THIS IS AN EASY-TO-USE DRIVER FOR THE MAIN OPTIMIZATION ROUTINE */
/* LMQNBC.  MORE EXPERIENCED USERS WHO WISH TO CUSTOMIZE PERFORMANCE */
/* OF THIS ALGORITHM SHOULD CALL LMQBC DIRECTLY. */

/* ---------------------------------------------------------------------- */
/* THIS ROUTINE SETS UP ALL THE PARAMETERS FOR THE TRUNCATED-NEWTON */
/* ALGORITHM.  THE PARAMETERS ARE: */

/* ETA    - SEVERITY OF THE LINESEARCH */
/* MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS */
/* XTOL   - DESIRED ACCURACY FOR THE SOLUTION X* */
/* STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH */
/* ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES */
/* MSGLVL - CONTROLS QUANTITY OF PRINTED OUTPUT */
/*          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION. */
/* MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP */

/*<       DOUBLE PRECISION  ETA, ACCRCY, XTOL, STEPMX, DSQRT, MCHPR1 >*/
/*<       EXTERNAL          SFUN >*/

/* SET PARAMETERS FOR THE OPTIMIZATION ROUTINE */

/*<       MAXIT = N/2 >*/
    /* Parameter adjustments */
    --ipivot;
    --up;
    --low;
    --g;
    --x;
    --w;

    /* Function Body */
    maxit = *n / 2;
/*<       IF (MAXIT .GT. 50) MAXIT = 50 >*/
    if (maxit > 50) {
	maxit = 50;
    }
/*<       IF (MAXIT .LE. 0) MAXIT = 1 >*/
    if (maxit <= 0) {
	maxit = 1;
    }
/*<       MSGLVL = 1 >*/
    msglvl = 1;
/*<       MAXFUN = 150*N >*/
    maxfun = *n * 150;
/*<       ETA = .25D0 >*/
    eta = .25;
/*<       STEPMX = 1.D1 >*/
    stepmx = 10.;
/*<       ACCRCY = 1.D2*MCHPR1() >*/
    accrcy = mchpr1_() * 100.;
/*<       XTOL = DSQRT(ACCRCY) >*/
    xtol = sqrt(accrcy);

/* MINIMIZE FUNCTION */

/*<    >*/
    lmqnbc_(ierror, n, &x[1], f, &g[1], &w[1], lw, (S_fp)sfun, &low[1], &up[1]
	    , &ipivot[1], &msglvl, &maxit, &maxfun, &eta, &stepmx, &accrcy, &
	    xtol, pParam);

/* PRINT RESULTS */

/*      IF (IERROR .NE. 0) WRITE(*,800) IERROR */
/*      WRITE(*,810) F */
/*<       IF (MSGLVL .LT. 1) RETURN >*/
    if (msglvl < 1) {
	return 0;
    }
/*      WRITE(*,820) */
/*<       NMAX = 10 >*/
    nmax = 10;
/*<       IF (N .LT. NMAX) NMAX = N >*/
    if (*n < nmax) {
	nmax = *n;
    }
/*      WRITE(*,830) (I,X(I),I=1,NMAX) */
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(//,' ERROR CODE =', I3) */
/* 810   FORMAT(//,' OPTIMAL FUNCTION VALUE = ', 1PD22.15) */
/* 820   FORMAT(10X, 'CURRENT SOLUTION IS (AT MOST 10 COMPONENTS)', /, */
/*     *       14X, 'I', 11X, 'X(I)') */
/* 830   FORMAT(10X, I5, 2X, 1PD22.15) */
/*<       END >*/
} /* tnbc_ */



/*<    >*/
/* Subroutine */ int lmqn_(integer *ifail, integer *n, doublereal *x,
	doublereal *f, doublereal *g, doublereal *w, integer *lw, S_fp sfun,
	integer *msglvl, integer *maxit, integer *maxfun, doublereal *eta,
	doublereal *stepmx, doublereal *accrcy, doublereal *xtol, void * pParam)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal fold, oldf;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);
    doublereal fnew;
    integer numf;
    doublereal peps;
    integer lhyr;
    doublereal zero, rtol, yksk, tiny;
    extern /* Subroutine */ int dxpy_(integer *, doublereal *, integer *,
	    doublereal *, integer *);
    integer nwhy;
    doublereal yrsr;
    extern doublereal dnrm2_(integer *, doublereal *, integer *), step1_(
	    doublereal *, doublereal *, doublereal *, doublereal *);
    integer i__;
    doublereal alpha, fkeep;
    integer ioldg;
    doublereal small;
    integer modet;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
	    doublereal *, integer *);
    integer niter;
    doublereal gnorm, ftest, fstop, pnorm, rteps, xnorm;
    integer idiagb;
    doublereal fm, pe, difold;
    integer icycle, nlincg, nfeval;
    doublereal difnew;
    integer nmodif;
    extern /* Subroutine */ int chkucp_(integer *, integer *, integer *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *);
    doublereal epsmch;
    extern /* Subroutine */ int linder_(integer *, S_fp, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, integer *, integer *, doublereal *, integer *, void *);
    doublereal epsred, abstol, oldgtp;
    extern /* Subroutine */ int modlnp_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     logical *, doublereal *, doublereal *, doublereal *, logical *,
	    S_fp, logical *, integer *, doublereal *, doublereal *,
	    doublereal *, doublereal *, void *);
    integer ireset;
    logical lreset;
    extern /* Subroutine */ int setpar_(integer *);
    doublereal reltol, gtpnew;
    integer nftotl;
    doublereal toleps;
    extern /* Subroutine */ int setucr_(doublereal *, integer *, integer *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, S_fp, doublereal *, doublereal *, void *);
    doublereal rtleps;
    integer ipivot, nm1;
    doublereal rtolsq, tnytol, gtg, one;
    integer ipk;
    doublereal gsk, spe;
    integer isk, iyk;
    logical upd1;

/*<       IMPLICIT          DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER           MSGLVL, N, MAXFUN, IFAIL, LW >*/
/*<       DOUBLE PRECISION  X(N), G(N), W(LW), ETA, XTOL, STEPMX, F, ACCRCY >*/

/* THIS ROUTINE IS A TRUNCATED-NEWTON METHOD. */
/* THE TRUNCATED-NEWTON METHOD IS PRECONDITIONED BY A LIMITED-MEMORY */
/* QUASI-NEWTON METHOD (THIS PRECONDITIONING STRATEGY IS DEVELOPED */
/* IN THIS ROUTINE) WITH A FURTHER DIAGONAL SCALING (SEE ROUTINE NDIA3). */
/* FOR FURTHER DETAILS ON THE PARAMETERS, SEE ROUTINE TN. */

/*<    >*/
/*<    >*/
/*<       LOGICAL LRESET, UPD1 >*/

/* THE FOLLOWING IMSL AND STANDARD FUNCTIONS ARE USED */

/*<       DOUBLE PRECISION DABS, DDOT, DSQRT, STEP1, DNRM2 >*/
/*<       EXTERNAL SFUN >*/
/*<    >*/

/* INITIALIZE PARAMETERS AND CONSTANTS */

/*      IF (MSGLVL .GE. -2) WRITE(*,800) */
/*<       CALL SETPAR(N) >*/
    /* Parameter adjustments */
    --g;
    --x;
    --w;

    /* Function Body */
    setpar_(n);
/*<       UPD1 = .TRUE. >*/
    upd1 = TRUE_;
/*<       IRESET = 0 >*/
    ireset = 0;
/*<       NFEVAL = 0 >*/
    nfeval = 0;
/*<       NMODIF = 0 >*/
    nmodif = 0;
/*<       NLINCG = 0 >*/
    nlincg = 0;
/*<       FSTOP = F >*/
    fstop = *f;
/*<       ZERO = 0.D0 >*/
    zero = 0.;
/*<       ONE = 1.D0 >*/
    one = 1.;
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;

/* WITHIN THIS ROUTINE THE ARRAY W(LOLDG) IS SHARED BY W(LHYR) */

/*<       LHYR = LOLDG >*/
    lhyr = subscr_1.loldg;

/* CHECK PARAMETERS AND SET CONSTANTS */

/*<    >*/
    chkucp_(&subscr_1.lwtest, maxfun, &nwhy, n, &alpha, &epsmch, eta, &peps, &
	    rteps, &rtol, &rtolsq, stepmx, &ftest, xtol, &xnorm, &x[1], lw, &
	    small, &tiny, accrcy);
/*<       IF (NWHY .LT. 0) GO TO 120 >*/
    if (nwhy < 0) {
	goto L120;
    }
/*<    >*/
    setucr_(&small, &nftotl, &niter, n, f, &fnew, &fm, &gtg, &oldf, (S_fp)
	    sfun, &g[1], &x[1], pParam);
/*<       FOLD = FNEW >*/
    fold = fnew;
/*      IF (MSGLVL .GE. 1) WRITE(*,810) NITER,NFTOTL,NLINCG,FNEW,GTG */

/* CHECK FOR SMALL GRADIENT AT THE STARTING POINT. */

/*<       FTEST = ONE + DABS(FNEW) >*/
    ftest = one + abs(fnew);
/*<       IF (GTG .LT. 1.D-4*EPSMCH*FTEST*FTEST) GO TO 90 >*/
    if (gtg < epsmch * 1e-4 * ftest * ftest) {
	goto L90;
    }

/* SET INITIAL VALUES TO OTHER PARAMETERS */

/*<       ICYCLE = NM1 >*/
    icycle = nm1;
/*<       TOLEPS = RTOL + RTEPS >*/
    toleps = rtol + rteps;
/*<       RTLEPS = RTOLSQ + EPSMCH >*/
    rtleps = rtolsq + epsmch;
/*<       GNORM  = DSQRT(GTG) >*/
    gnorm = sqrt(gtg);
/*<       DIFNEW = ZERO >*/
    difnew = zero;
/*<       EPSRED = 5.0D-2 >*/
    epsred = .05;
/*<       FKEEP  = FNEW >*/
    fkeep = fnew;

/* SET THE DIAGONAL OF THE APPROXIMATE HESSIAN TO UNITY. */

/*<       IDIAGB = LDIAGB >*/
    idiagb = subscr_1.ldiagb;
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          W(IDIAGB) = ONE >*/
	w[idiagb] = one;
/*<          IDIAGB = IDIAGB + 1 >*/
	++idiagb;
/*< 10    CONTINUE >*/
/* L10: */
    }

/* ..................START OF MAIN ITERATIVE LOOP.......... */

/* COMPUTE THE NEW SEARCH DIRECTION */

/*<       MODET = MSGLVL - 3 >*/
    modet = *msglvl - 3;
/*<    >*/
    modlnp_(&modet, &w[subscr_1.lpk], &w[subscr_1.lgv], &w[subscr_1.lz1], &w[
	    subscr_1.lv], &w[subscr_1.ldiagb], &w[subscr_1.lemat], &x[1], &g[
	    1], &w[subscr_1.lzk], n, &w[1], lw, &niter, maxit, &nfeval, &
	    nmodif, &nlincg, &upd1, &yksk, &gsk, &yrsr, &lreset, (S_fp)sfun, &
	    c_false, &ipivot, accrcy, &gtpnew, &gnorm, &xnorm, pParam);
/*< 20    CONTINUE >*/
L20:
/*<       CALL DCOPY(N,G,1,W(LOLDG),1) >*/
    dcopy_(n, &g[1], &c__1, &w[subscr_1.loldg], &c__1);
/*<       PNORM = DNRM2(N,W(LPK),1) >*/
    pnorm = dnrm2_(n, &w[subscr_1.lpk], &c__1);
/*<       OLDF = FNEW >*/
    oldf = fnew;
/*<       OLDGTP = GTPNEW >*/
    oldgtp = gtpnew;

/* PREPARE TO COMPUTE THE STEP LENGTH */

/*<       PE = PNORM + EPSMCH >*/
    pe = pnorm + epsmch;

/* COMPUTE THE ABSOLUTE AND RELATIVE TOLERANCES FOR THE LINEAR SEARCH */

/*<       RELTOL = RTEPS*(XNORM + ONE)/PE >*/
    reltol = rteps * (xnorm + one) / pe;
/*<       ABSTOL = - EPSMCH*FTEST/(OLDGTP - EPSMCH) >*/
    abstol = -epsmch * ftest / (oldgtp - epsmch);

/* COMPUTE THE SMALLEST ALLOWABLE SPACING BETWEEN POINTS IN */
/* THE LINEAR SEARCH */

/*<       TNYTOL = EPSMCH*(XNORM + ONE)/PE >*/
    tnytol = epsmch * (xnorm + one) / pe;
/*<       SPE = STEPMX/PE >*/
    spe = *stepmx / pe;

/* SET THE INITIAL STEP LENGTH. */

/*<       ALPHA = STEP1(FNEW,FM,OLDGTP,SPE) >*/
    alpha = step1_(&fnew, &fm, &oldgtp, &spe);

/* PERFORM THE LINEAR SEARCH */

/*<    >*/
    linder_(n, (S_fp)sfun, &small, &epsmch, &reltol, &abstol, &tnytol, eta, &
	    zero, &spe, &w[subscr_1.lpk], &oldgtp, &x[1], &fnew, &alpha, &g[1]
	    , &numf, &nwhy, &w[1], lw, pParam);

/*<       FOLD = FNEW >*/
    fold = fnew;
/*<       NITER = NITER + 1 >*/
    ++niter;
/*<       NFTOTL = NFTOTL + NUMF >*/
    nftotl += numf;
/*<       GTG = DDOT(N,G,1,G,1) >*/
    gtg = ddot_(n, &g[1], &c__1, &g[1], &c__1);
/*      IF (MSGLVL .GE. 1) WRITE(*,810) NITER,NFTOTL,NLINCG,FNEW,GTG */
/*<       IF (NWHY .LT. 0) GO TO 120 >*/
    if (nwhy < 0) {
	goto L120;
    }
/*<       IF (NWHY .EQ. 0 .OR. NWHY .EQ. 2) GO TO 30 >*/
    if (nwhy == 0 || nwhy == 2) {
	goto L30;
    }

/* THE LINEAR SEARCH HAS FAILED TO FIND A LOWER POINT */

/*<       NWHY = 3 >*/
    nwhy = 3;
/*<       GO TO 100 >*/
    goto L100;
/*< 30    IF (NWHY .LE. 1) GO TO 40 >*/
L30:
    if (nwhy <= 1) {
	goto L40;
    }
/*<       CALL SFUN(N,X,FNEW,G) >*/
    (*sfun)(n, &x[1], &fnew, &g[1], pParam);
/*<       NFTOTL = NFTOTL + 1 >*/
    ++nftotl;

/* TERMINATE IF MORE THAN MAXFUN EVALUTATIONS HAVE BEEN MADE */

/*< 40    NWHY = 2 >*/
L40:
    nwhy = 2;
/*<       IF (NFTOTL .GT. MAXFUN) GO TO 110 >*/
    if (nftotl > *maxfun) {
	goto L110;
    }
/*<       NWHY = 0 >*/
    nwhy = 0;

/* SET UP PARAMETERS USED IN CONVERGENCE AND RESETTING TESTS */

/*<       DIFOLD = DIFNEW >*/
    difold = difnew;
/*<       DIFNEW = OLDF - FNEW >*/
    difnew = oldf - fnew;

/* IF THIS IS THE FIRST ITERATION OF A NEW CYCLE, COMPUTE THE */
/* PERCENTAGE REDUCTION FACTOR FOR THE RESETTING TEST. */

/*<       IF (ICYCLE .NE. 1) GO TO 50 >*/
    if (icycle != 1) {
	goto L50;
    }
/*<       IF (DIFNEW .GT. 2.0D0 *DIFOLD) EPSRED = EPSRED + EPSRED >*/
    if (difnew > difold * 2.) {
	epsred += epsred;
    }
/*<       IF (DIFNEW .LT. 5.0D-1*DIFOLD) EPSRED = 5.0D-1*EPSRED >*/
    if (difnew < difold * .5) {
	epsred *= .5;
    }
/*< 50    CONTINUE >*/
L50:
/*<       GNORM = DSQRT(GTG) >*/
    gnorm = sqrt(gtg);
/*<       FTEST = ONE + DABS(FNEW) >*/
    ftest = one + abs(fnew);
/*<       XNORM = DNRM2(N,X,1) >*/
    xnorm = dnrm2_(n, &x[1], &c__1);

/* TEST FOR CONVERGENCE */

/*<    >*/
    if (alpha * pnorm < toleps * (one + xnorm) && abs(difnew) < rtleps *
	    ftest && gtg < peps * ftest * ftest || gtg < *accrcy * 1e-4 *
	    ftest * ftest) {
	goto L90;
    }

/* COMPUTE THE CHANGE IN THE ITERATES AND THE CORRESPONDING CHANGE */
/* IN THE GRADIENTS */

/*<       ISK = LSK >*/
    isk = subscr_1.lsk;
/*<       IPK = LPK >*/
    ipk = subscr_1.lpk;
/*<       IYK = LYK >*/
    iyk = subscr_1.lyk;
/*<       IOLDG = LOLDG >*/
    ioldg = subscr_1.loldg;
/*<       DO 60 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          W(IYK) = G(I) - W(IOLDG) >*/
	w[iyk] = g[i__] - w[ioldg];
/*<          W(ISK) = ALPHA*W(IPK) >*/
	w[isk] = alpha * w[ipk];
/*<          IPK = IPK + 1 >*/
	++ipk;
/*<          ISK = ISK + 1 >*/
	++isk;
/*<          IYK = IYK + 1 >*/
	++iyk;
/*<          IOLDG = IOLDG + 1 >*/
	++ioldg;
/*< 60    CONTINUE >*/
/* L60: */
    }

/* SET UP PARAMETERS USED IN UPDATING THE DIRECTION OF SEARCH. */

/*<       YKSK = DDOT(N,W(LYK),1,W(LSK),1) >*/
    yksk = ddot_(n, &w[subscr_1.lyk], &c__1, &w[subscr_1.lsk], &c__1);
/*<       LRESET = .FALSE. >*/
    lreset = FALSE_;
/*<    >*/
    if (icycle == nm1 || difnew < epsred * (fkeep - fnew)) {
	lreset = TRUE_;
    }
/*<       IF (LRESET) GO TO 70 >*/
    if (lreset) {
	goto L70;
    }
/*<       YRSR = DDOT(N,W(LYR),1,W(LSR),1) >*/
    yrsr = ddot_(n, &w[subscr_1.lyr], &c__1, &w[subscr_1.lsr], &c__1);
/*<       IF (YRSR .LE. ZERO) LRESET = .TRUE. >*/
    if (yrsr <= zero) {
	lreset = TRUE_;
    }
/*< 70    CONTINUE >*/
L70:
/*<       UPD1 = .FALSE. >*/
    upd1 = FALSE_;

/*      COMPUTE THE NEW SEARCH DIRECTION */

/*<       MODET = MSGLVL - 3 >*/
    modet = *msglvl - 3;
/*<    >*/
    modlnp_(&modet, &w[subscr_1.lpk], &w[subscr_1.lgv], &w[subscr_1.lz1], &w[
	    subscr_1.lv], &w[subscr_1.ldiagb], &w[subscr_1.lemat], &x[1], &g[
	    1], &w[subscr_1.lzk], n, &w[1], lw, &niter, maxit, &nfeval, &
	    nmodif, &nlincg, &upd1, &yksk, &gsk, &yrsr, &lreset, (S_fp)sfun, &
	    c_false, &ipivot, accrcy, &gtpnew, &gnorm, &xnorm, pParam);
/*<       IF (LRESET) GO TO 80 >*/
    if (lreset) {
	goto L80;
    }

/*      STORE THE ACCUMULATED CHANGE IN THE POINT AND GRADIENT AS AN */
/*      "AVERAGE" DIRECTION FOR PRECONDITIONING. */

/*<       CALL DXPY(N,W(LSK),1,W(LSR),1) >*/
    dxpy_(n, &w[subscr_1.lsk], &c__1, &w[subscr_1.lsr], &c__1);
/*<       CALL DXPY(N,W(LYK),1,W(LYR),1) >*/
    dxpy_(n, &w[subscr_1.lyk], &c__1, &w[subscr_1.lyr], &c__1);
/*<       ICYCLE = ICYCLE + 1 >*/
    ++icycle;
/*<       GOTO 20 >*/
    goto L20;

/* RESET */

/*< 80    IRESET = IRESET + 1 >*/
L80:
    ++ireset;

/* INITIALIZE THE SUM OF ALL THE CHANGES IN X. */

/*<       CALL DCOPY(N,W(LSK),1,W(LSR),1) >*/
    dcopy_(n, &w[subscr_1.lsk], &c__1, &w[subscr_1.lsr], &c__1);
/*<       CALL DCOPY(N,W(LYK),1,W(LYR),1) >*/
    dcopy_(n, &w[subscr_1.lyk], &c__1, &w[subscr_1.lyr], &c__1);
/*<       FKEEP = FNEW >*/
    fkeep = fnew;
/*<       ICYCLE = 1 >*/
    icycle = 1;
/*<       GO TO 20 >*/
    goto L20;

/* ...............END OF MAIN ITERATION....................... */

/*< 90    IFAIL = 0 >*/
L90:
    *ifail = 0;
/*<       F = FNEW >*/
    *f = fnew;
/*<       RETURN >*/
    return 0;
/*< 100   OLDF = FNEW >*/
L100:
    oldf = fnew;

/* LOCAL SEARCH HERE COULD BE INSTALLED HERE */

/*< 110    F = OLDF >*/
L110:
    *f = oldf;

/* SET IFAIL */

/*< 120   IFAIL = NWHY >*/
L120:
    *ifail = nwhy;
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(//' NIT   NF   CG', 9X, 'F', 21X, 'GTG',//) */
/* 810   FORMAT(' ',I3,1X,I4,1X,I4,1X,1PD22.15,2X,1PD15.8) */
/*<       END >*/
} /* lmqn_ */



/*<    >*/
/* Subroutine */ int lmqnbc_(integer *ifail, integer *n, doublereal *x,
	doublereal *f, doublereal *g, doublereal *w, integer *lw, S_fp sfun,
	doublereal *low, doublereal *up, integer *ipivot, integer *msglvl,
	integer *maxit, integer *maxfun, doublereal *eta, doublereal *stepmx,
	doublereal *accrcy, doublereal *xtol, void *pParam)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal fold, oldf;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);
    doublereal fnew;
    integer numf;
    logical conv;
    doublereal peps;
    extern /* Subroutine */ int modz_(integer *, doublereal *, doublereal *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);
    integer lhyr;
    doublereal zero, rtol, yksk, tiny;
    extern /* Subroutine */ int dxpy_(integer *, doublereal *, integer *,
	    doublereal *, integer *);
    integer nwhy;
    doublereal yrsr;
    extern doublereal dnrm2_(integer *, doublereal *, integer *), step1_(
	    doublereal *, doublereal *, doublereal *, doublereal *);
    integer i__;
    doublereal alpha, fkeep;
    integer ioldg;
    extern /* Subroutine */ int crash_(integer *, doublereal *, integer *,
	    doublereal *, doublereal *, integer *);
    doublereal small, flast;
    integer modet;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
	    doublereal *, integer *);
    integer niter;
    doublereal gnorm, ftest;
    extern /* Subroutine */ int monit_(), ztime_(integer *, doublereal *,
	    integer *);
    doublereal fstop, pnorm, rteps, xnorm;
    integer idiagb;
    doublereal fm, pe, difold;
    integer icycle, nlincg, nfeval;
    doublereal difnew;
    integer nmodif;
    doublereal epsmch;
    extern /* Subroutine */ int chkucp_(integer *, integer *, integer *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *), linder_(integer *,
	    S_fp, doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, integer *, integer *, doublereal *,
	    integer *, void *);
    doublereal epsred, abstol, oldgtp;
    logical newcon;
    integer ireset;
    extern /* Subroutine */ int modlnp_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, integer *, integer *,
	     logical *, doublereal *, doublereal *, doublereal *, logical *,
	    S_fp, logical *, integer *, doublereal *, doublereal *,
	    doublereal *, doublereal *, void *);
    logical lreset;
    extern /* Subroutine */ int setpar_(integer *);
    doublereal reltol, gtpnew;
    integer nftotl;
    doublereal toleps;
    extern /* Subroutine */ int setucr_(doublereal *, integer *, integer *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, S_fp, doublereal *, doublereal *, void *);
    doublereal rtleps;
    integer nm1;
    extern /* Subroutine */ int stpmax_(doublereal *, doublereal *,
	    doublereal *, integer *, doublereal *, doublereal *, integer *,
	    doublereal *, doublereal *), cnvtst_(logical *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, integer *, integer *, doublereal *);
    doublereal rtolsq, tnytol;
    integer ier;
    doublereal gtg, one;
    integer ipk;
    doublereal gsk, spe;
    integer isk, iyk;
    logical upd1;

/*<       IMPLICIT         DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER          MSGLVL,N,MAXFUN,IFAIL,LW >*/
/*<       INTEGER          IPIVOT(N) >*/
/*<       DOUBLE PRECISION ETA,XTOL,STEPMX,F,ACCRCY >*/
/*<       DOUBLE PRECISION X(N),G(N),W(LW),LOW(N),UP(N) >*/

/* THIS ROUTINE IS A BOUNDS-CONSTRAINED TRUNCATED-NEWTON METHOD. */
/* THE TRUNCATED-NEWTON METHOD IS PRECONDITIONED BY A LIMITED-MEMORY */
/* QUASI-NEWTON METHOD (THIS PRECONDITIONING STRATEGY IS DEVELOPED */
/* IN THIS ROUTINE) WITH A FURTHER DIAGONAL SCALING (SEE ROUTINE NDIA3). */
/* FOR FURTHER DETAILS ON THE PARAMETERS, SEE ROUTINE TNBC. */

/*<    >*/
/*<    >*/
/*<       LOGICAL CONV, LRESET, UPD1, NEWCON >*/

/* THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE USED */

/*<       DOUBLE PRECISION DABS, DDOT, DNRM2, DSQRT, STEP1 >*/
/*<       EXTERNAL SFUN >*/
/*<    >*/

/* CHECK THAT INITIAL X IS FEASIBLE AND THAT THE BOUNDS ARE CONSISTENT */

/*<       CALL CRASH(N,X,IPIVOT,LOW,UP,IER) >*/
    /* Parameter adjustments */
    --ipivot;
    --up;
    --low;
    --g;
    --x;
    --w;

    /* Function Body */
    crash_(n, &x[1], &ipivot[1], &low[1], &up[1], &ier);
/*      IF (IER .NE. 0) WRITE(*,800) */
/*<       IF (IER .NE. 0) RETURN >*/
    if (ier != 0) {
	return 0;
    }

/*      IF (MSGLVL .GE. 1) WRITE(*,810) */

/* INITIALIZE VARIABLES */

/*<       CALL SETPAR(N) >*/
    setpar_(n);
/*<       UPD1 = .TRUE. >*/
    upd1 = TRUE_;
/*<       IRESET = 0 >*/
    ireset = 0;
/*<       NFEVAL = 0 >*/
    nfeval = 0;
/*<       NMODIF = 0 >*/
    nmodif = 0;
/*<       NLINCG = 0 >*/
    nlincg = 0;
/*<       FSTOP = F >*/
    fstop = *f;
/*<       CONV = .FALSE. >*/
    conv = FALSE_;
/*<       ZERO = 0.D0 >*/
    zero = 0.;
/*<       ONE = 1.D0 >*/
    one = 1.;
/*<       NM1 = N - 1 >*/
    nm1 = *n - 1;

/* WITHIN THIS ROUTINE THE ARRAY W(LOLDG) IS SHARED BY W(LHYR) */

/*<       LHYR = LOLDG >*/
    lhyr = subscr_1.loldg;

/* CHECK PARAMETERS AND SET CONSTANTS */

/*<    >*/
    chkucp_(&subscr_1.lwtest, maxfun, &nwhy, n, &alpha, &epsmch, eta, &peps, &
	    rteps, &rtol, &rtolsq, stepmx, &ftest, xtol, &xnorm, &x[1], lw, &
	    small, &tiny, accrcy);
/*<       IF (NWHY .LT. 0) GO TO 160 >*/
    if (nwhy < 0) {
	goto L160;
    }
/*<    >*/
    setucr_(&small, &nftotl, &niter, n, f, &fnew, &fm, &gtg, &oldf, (S_fp)
	    sfun, &g[1], &x[1], pParam);
/*<       FOLD = FNEW >*/
    fold = fnew;
/*<       FLAST = FNEW >*/
    flast = fnew;

/* TEST THE LAGRANGE MULTIPLIERS TO SEE IF THEY ARE NON-NEGATIVE. */
/* BECAUSE THE CONSTRAINTS ARE ONLY LOWER BOUNDS, THE COMPONENTS */
/* OF THE GRADIENT CORRESPONDING TO THE ACTIVE CONSTRAINTS ARE THE */
/* LAGRANGE MULTIPLIERS.  AFTERWORDS, THE PROJECTED GRADIENT IS FORMED. */

/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (IPIVOT(I) .EQ. 2) GO TO 10 >*/
	if (ipivot[i__] == 2) {
	    goto L10;
	}
/*<          IF (-IPIVOT(I)*G(I) .GE. 0.D0) GO TO 10 >*/
	if (-ipivot[i__] * g[i__] >= 0.) {
	    goto L10;
	}
/*<          IPIVOT(I) = 0 >*/
	ipivot[i__] = 0;
/*< 10    CONTINUE >*/
L10:
	;
    }
/*<       CALL ZTIME(N,G,IPIVOT) >*/
    ztime_(n, &g[1], &ipivot[1]);
/*<       GTG = DDOT(N,G,1,G,1) >*/
    gtg = ddot_(n, &g[1], &c__1, &g[1], &c__1);
/*<    >*/
    if (*msglvl >= 1) {
	monit_(n, &x[1], &fnew, &g[1], &niter, &nftotl, &nfeval, &lreset, &
		ipivot[1]);
    }

/* CHECK IF THE INITIAL POINT IS A LOCAL MINIMUM. */

/*<       FTEST = ONE + DABS(FNEW) >*/
    ftest = one + abs(fnew);
/*<       IF (GTG .LT. 1.D-4*EPSMCH*FTEST*FTEST) GO TO 130 >*/
    if (gtg < epsmch * 1e-4 * ftest * ftest) {
	goto L130;
    }

/* SET INITIAL VALUES TO OTHER PARAMETERS */

/*<       ICYCLE = NM1 >*/
    icycle = nm1;
/*<       TOLEPS = RTOL + RTEPS >*/
    toleps = rtol + rteps;
/*<       RTLEPS = RTOLSQ + EPSMCH >*/
    rtleps = rtolsq + epsmch;
/*<       GNORM  = DSQRT(GTG) >*/
    gnorm = sqrt(gtg);
/*<       DIFNEW = ZERO >*/
    difnew = zero;
/*<       EPSRED = 5.0D-2 >*/
    epsred = .05;
/*<       FKEEP  = FNEW >*/
    fkeep = fnew;

/* SET THE DIAGONAL OF THE APPROXIMATE HESSIAN TO UNITY. */

/*<       IDIAGB = LDIAGB >*/
    idiagb = subscr_1.ldiagb;
/*<       DO 15 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          W(IDIAGB) = ONE >*/
	w[idiagb] = one;
/*<          IDIAGB = IDIAGB + 1 >*/
	++idiagb;
/*< 15    CONTINUE >*/
/* L15: */
    }

/* ..................START OF MAIN ITERATIVE LOOP.......... */

/* COMPUTE THE NEW SEARCH DIRECTION */

/*<       MODET = MSGLVL - 3 >*/
    modet = *msglvl - 3;
/*<    >*/
    modlnp_(&modet, &w[subscr_1.lpk], &w[subscr_1.lgv], &w[subscr_1.lz1], &w[
	    subscr_1.lv], &w[subscr_1.ldiagb], &w[subscr_1.lemat], &x[1], &g[
	    1], &w[subscr_1.lzk], n, &w[1], lw, &niter, maxit, &nfeval, &
	    nmodif, &nlincg, &upd1, &yksk, &gsk, &yrsr, &lreset, (S_fp)sfun, &
	    c_true, &ipivot[1], accrcy, &gtpnew, &gnorm, &xnorm, pParam);
/*< 20    CONTINUE >*/
L20:
/*<       CALL DCOPY(N,G,1,W(LOLDG),1) >*/
    dcopy_(n, &g[1], &c__1, &w[subscr_1.loldg], &c__1);
/*<       PNORM = DNRM2(N,W(LPK),1) >*/
    pnorm = dnrm2_(n, &w[subscr_1.lpk], &c__1);
/*<       OLDF = FNEW >*/
    oldf = fnew;
/*<       OLDGTP = GTPNEW >*/
    oldgtp = gtpnew;

/* PREPARE TO COMPUTE THE STEP LENGTH */

/*<       PE = PNORM + EPSMCH >*/
    pe = pnorm + epsmch;

/* COMPUTE THE ABSOLUTE AND RELATIVE TOLERANCES FOR THE LINEAR SEARCH */

/*<       RELTOL = RTEPS*(XNORM + ONE)/PE >*/
    reltol = rteps * (xnorm + one) / pe;
/*<       ABSTOL = - EPSMCH*FTEST/(OLDGTP - EPSMCH) >*/
    abstol = -epsmch * ftest / (oldgtp - epsmch);

/* COMPUTE THE SMALLEST ALLOWABLE SPACING BETWEEN POINTS IN */
/* THE LINEAR SEARCH */

/*<       TNYTOL = EPSMCH*(XNORM + ONE)/PE >*/
    tnytol = epsmch * (xnorm + one) / pe;
/*<       CALL STPMAX(STEPMX,PE,SPE,N,X,W(LPK),IPIVOT,LOW,UP) >*/
    stpmax_(stepmx, &pe, &spe, n, &x[1], &w[subscr_1.lpk], &ipivot[1], &low[1]
	    , &up[1]);

/* SET THE INITIAL STEP LENGTH. */

/*<       ALPHA = STEP1(FNEW,FM,OLDGTP,SPE) >*/
    alpha = step1_(&fnew, &fm, &oldgtp, &spe);

/* PERFORM THE LINEAR SEARCH */

/*<    >*/
    linder_(n, (S_fp)sfun, &small, &epsmch, &reltol, &abstol, &tnytol, eta, &
	    zero, &spe, &w[subscr_1.lpk], &oldgtp, &x[1], &fnew, &alpha, &g[1]
	    , &numf, &nwhy, &w[1], lw, pParam);
/*<       NEWCON = .FALSE. >*/
    newcon = FALSE_;
/*<       IF (DABS(ALPHA-SPE) .GT. 1.D1*EPSMCH) GO TO 30 >*/
    if ((d__1 = alpha - spe, abs(d__1)) > epsmch * 10.) {
	goto L30;
    }
/*<       NEWCON = .TRUE. >*/
    newcon = TRUE_;
/*<       NWHY   = 0 >*/
    nwhy = 0;
/*<       CALL MODZ(N,X,W(LPK),IPIVOT,EPSMCH,LOW,UP,FLAST,FNEW) >*/
    modz_(n, &x[1], &w[subscr_1.lpk], &ipivot[1], &epsmch, &low[1], &up[1], &
	    flast, &fnew);
/*<       FLAST = FNEW >*/
    flast = fnew;

/*<  30   CONTINUE >*/
L30:
/* 30    IF (MSGLVL .GE. 3) WRITE(*,820) ALPHA,PNORM */
/*<       FOLD = FNEW >*/
    fold = fnew;
/*<       NITER = NITER + 1 >*/
    ++niter;
/*<       NFTOTL = NFTOTL + NUMF >*/
    nftotl += numf;

/* IF REQUIRED, PRINT THE DETAILS OF THIS ITERATION */

/*<    >*/
    if (*msglvl >= 1) {
	monit_(n, &x[1], &fnew, &g[1], &niter, &nftotl, &nfeval, &lreset, &
		ipivot[1]);
    }
/*<       IF (NWHY .LT. 0) GO TO 160 >*/
    if (nwhy < 0) {
	goto L160;
    }
/*<       IF (NWHY .EQ. 0 .OR. NWHY .EQ. 2) GO TO 40 >*/
    if (nwhy == 0 || nwhy == 2) {
	goto L40;
    }

/* THE LINEAR SEARCH HAS FAILED TO FIND A LOWER POINT */

/*<       NWHY = 3 >*/
    nwhy = 3;
/*<       GO TO 140 >*/
    goto L140;
/*< 40    IF (NWHY .LE. 1) GO TO 50 >*/
L40:
    if (nwhy <= 1) {
	goto L50;
    }
/*<       CALL SFUN(N,X,FNEW,G) >*/
    (*sfun)(n, &x[1], &fnew, &g[1],pParam);
/*<       NFTOTL = NFTOTL + 1 >*/
    ++nftotl;

/* TERMINATE IF MORE THAN MAXFUN EVALUATIONS HAVE BEEN MADE */

/*< 50    NWHY = 2 >*/
L50:
    nwhy = 2;
/*<       IF (NFTOTL .GT. MAXFUN) GO TO 150 >*/
    if (nftotl > *maxfun) {
	goto L150;
    }
/*<       NWHY = 0 >*/
    nwhy = 0;

/* SET UP PARAMETERS USED IN CONVERGENCE AND RESETTING TESTS */

/*<       DIFOLD = DIFNEW >*/
    difold = difnew;
/*<       DIFNEW = OLDF - FNEW >*/
    difnew = oldf - fnew;

/* IF THIS IS THE FIRST ITERATION OF A NEW CYCLE, COMPUTE THE */
/* PERCENTAGE REDUCTION FACTOR FOR THE RESETTING TEST. */

/*<       IF (ICYCLE .NE. 1) GO TO 60 >*/
    if (icycle != 1) {
	goto L60;
    }
/*<       IF (DIFNEW .GT. 2.D0*DIFOLD) EPSRED = EPSRED + EPSRED >*/
    if (difnew > difold * 2.) {
	epsred += epsred;
    }
/*<       IF (DIFNEW .LT. 5.0D-1*DIFOLD) EPSRED = 5.0D-1*EPSRED >*/
    if (difnew < difold * .5) {
	epsred *= .5;
    }
/*< 60    CALL DCOPY(N,G,1,W(LGV),1) >*/
L60:
    dcopy_(n, &g[1], &c__1, &w[subscr_1.lgv], &c__1);
/*<       CALL ZTIME(N,W(LGV),IPIVOT) >*/
    ztime_(n, &w[subscr_1.lgv], &ipivot[1]);
/*<       GTG = DDOT(N,W(LGV),1,W(LGV),1) >*/
    gtg = ddot_(n, &w[subscr_1.lgv], &c__1, &w[subscr_1.lgv], &c__1);
/*<       GNORM = DSQRT(GTG) >*/
    gnorm = sqrt(gtg);
/*<       FTEST = ONE + DABS(FNEW) >*/
    ftest = one + abs(fnew);
/*<       XNORM = DNRM2(N,X,1) >*/
    xnorm = dnrm2_(n, &x[1], &c__1);

/* TEST FOR CONVERGENCE */

/*<    >*/
    cnvtst_(&conv, &alpha, &pnorm, &toleps, &xnorm, &difnew, &rtleps, &ftest,
	    &gtg, &peps, &epsmch, &gtpnew, &fnew, &flast, &g[1], &ipivot[1],
	    n, accrcy);
/*<       IF (CONV) GO TO 130 >*/
    if (conv) {
	goto L130;
    }
/*<       CALL ZTIME(N,G,IPIVOT) >*/
    ztime_(n, &g[1], &ipivot[1]);

/* COMPUTE THE CHANGE IN THE ITERATES AND THE CORRESPONDING CHANGE */
/* IN THE GRADIENTS */

/*<       IF (NEWCON) GO TO 90 >*/
    if (newcon) {
	goto L90;
    }
/*<       ISK = LSK >*/
    isk = subscr_1.lsk;
/*<       IPK = LPK >*/
    ipk = subscr_1.lpk;
/*<       IYK = LYK >*/
    iyk = subscr_1.lyk;
/*<       IOLDG = LOLDG >*/
    ioldg = subscr_1.loldg;
/*<       DO 70 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          W(IYK) = G(I) - W(IOLDG) >*/
	w[iyk] = g[i__] - w[ioldg];
/*<          W(ISK) = ALPHA*W(IPK) >*/
	w[isk] = alpha * w[ipk];
/*<          IPK = IPK + 1 >*/
	++ipk;
/*<          ISK = ISK + 1 >*/
	++isk;
/*<          IYK = IYK + 1 >*/
	++iyk;
/*<          IOLDG = IOLDG + 1 >*/
	++ioldg;
/*< 70    CONTINUE >*/
/* L70: */
    }

/* SET UP PARAMETERS USED IN UPDATING THE PRECONDITIONING STRATEGY. */

/*<       YKSK = DDOT(N,W(LYK),1,W(LSK),1) >*/
    yksk = ddot_(n, &w[subscr_1.lyk], &c__1, &w[subscr_1.lsk], &c__1);
/*<       LRESET = .FALSE. >*/
    lreset = FALSE_;
/*<    >*/
    if (icycle == nm1 || difnew < epsred * (fkeep - fnew)) {
	lreset = TRUE_;
    }
/*<       IF (LRESET) GO TO 80 >*/
    if (lreset) {
	goto L80;
    }
/*<       YRSR = DDOT(N,W(LYR),1,W(LSR),1) >*/
    yrsr = ddot_(n, &w[subscr_1.lyr], &c__1, &w[subscr_1.lsr], &c__1);
/*<       IF (YRSR .LE. ZERO) LRESET = .TRUE. >*/
    if (yrsr <= zero) {
	lreset = TRUE_;
    }
/*< 80    CONTINUE >*/
L80:
/*<       UPD1 = .FALSE. >*/
    upd1 = FALSE_;

/*      COMPUTE THE NEW SEARCH DIRECTION */

/*<  90   CONTINUE >*/
L90:
/* 90    IF (UPD1 .AND. MSGLVL .GE. 3) WRITE(*,830) */
/*      IF (NEWCON .AND. MSGLVL .GE. 3) WRITE(*,840) */
/*<       MODET = MSGLVL - 3 >*/
    modet = *msglvl - 3;
/*<    >*/
    modlnp_(&modet, &w[subscr_1.lpk], &w[subscr_1.lgv], &w[subscr_1.lz1], &w[
	    subscr_1.lv], &w[subscr_1.ldiagb], &w[subscr_1.lemat], &x[1], &g[
	    1], &w[subscr_1.lzk], n, &w[1], lw, &niter, maxit, &nfeval, &
	    nmodif, &nlincg, &upd1, &yksk, &gsk, &yrsr, &lreset, (S_fp)sfun, &
	    c_true, &ipivot[1], accrcy, &gtpnew, &gnorm, &xnorm, pParam);
/*<       IF (NEWCON) GO TO 20 >*/
    if (newcon) {
	goto L20;
    }
/*<       IF (LRESET) GO TO 110 >*/
    if (lreset) {
	goto L110;
    }

/* COMPUTE THE ACCUMULATED STEP AND ITS CORRESPONDING */
/* GRADIENT DIFFERENCE. */

/*<       CALL DXPY(N,W(LSK),1,W(LSR),1) >*/
    dxpy_(n, &w[subscr_1.lsk], &c__1, &w[subscr_1.lsr], &c__1);
/*<       CALL DXPY(N,W(LYK),1,W(LYR),1) >*/
    dxpy_(n, &w[subscr_1.lyk], &c__1, &w[subscr_1.lyr], &c__1);
/*<       ICYCLE = ICYCLE + 1 >*/
    ++icycle;
/*<       GOTO 20 >*/
    goto L20;

/* RESET */

/*< 110   IRESET = IRESET + 1 >*/
L110:
    ++ireset;

/* INITIALIZE THE SUM OF ALL THE CHANGES IN X. */

/*<       CALL DCOPY(N,W(LSK),1,W(LSR),1) >*/
    dcopy_(n, &w[subscr_1.lsk], &c__1, &w[subscr_1.lsr], &c__1);
/*<       CALL DCOPY(N,W(LYK),1,W(LYR),1) >*/
    dcopy_(n, &w[subscr_1.lyk], &c__1, &w[subscr_1.lyr], &c__1);
/*<       FKEEP = FNEW >*/
    fkeep = fnew;
/*<       ICYCLE = 1 >*/
    icycle = 1;
/*<       GO TO 20 >*/
    goto L20;

/* ...............END OF MAIN ITERATION....................... */

/*< 130   IFAIL = 0 >*/
L130:
    *ifail = 0;
/*<       F = FNEW >*/
    *f = fnew;
/*<       RETURN >*/
    return 0;
/*< 140   OLDF = FNEW >*/
L140:
    oldf = fnew;

/* LOCAL SEARCH COULD BE INSTALLED HERE */

/*< 150   F = OLDF >*/
L150:
    *f = oldf;
/*<    >*/
    if (*msglvl >= 1) {
	monit_(n, &x[1], f, &g[1], &niter, &nftotl, &nfeval, &ireset, &ipivot[
		1]);
    }

/* SET IFAIL */

/*< 160   IFAIL = NWHY >*/
L160:
    *ifail = nwhy;
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(' THERE IS NO FEASIBLE POINT; TERMINATING ALGORITHM') */
/* 810   FORMAT(//'  NIT   NF   CG', 9X, 'F', 21X, 'GTG',//) */
/* 820   FORMAT('        LINESEARCH RESULTS:  ALPHA,PNORM',2(1PD12.4)) */
/* 830   FORMAT(' UPD1 IS TRUE - TRIVIAL PRECONDITIONING') */
/* 840   FORMAT(' NEWCON IS TRUE - CONSTRAINT ADDED IN LINESEARCH') */
/*<       END >*/
} /* lmqnbc_ */



/*<       SUBROUTINE MONIT(N,X,F,G,NITER,NFTOTL,NFEVAL,IRESET,IPIVOT) >*/
/* Subroutine */ int monit_(integer *n, doublereal *x, doublereal *f,
	doublereal *g, integer *niter, integer *nftotl, integer *nfeval,
	integer *ireset, integer *ipivot)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal gtg;


/* PRINT RESULTS OF CURRENT ITERATION */

/*<       IMPLICIT         DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION X(N),F,G(N),GTG >*/
/*<       INTEGER          IPIVOT(N) >*/

/*<       GTG = 0.D0 >*/
    /* Parameter adjustments */
    --ipivot;
    --g;
    --x;

    /* Function Body */
    gtg = 0.;
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (IPIVOT(I) .NE. 0) GO TO 10 >*/
	if (ipivot[i__] != 0) {
	    goto L10;
	}
/*<          GTG = GTG + G(I)*G(I) >*/
	gtg += g[i__] * g[i__];
/*< 10    CONTINUE >*/
L10:
	;
    }
/*      WRITE(*,800) NITER,NFTOTL,NFEVAL,F,GTG */
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(' ',I4,1X,I4,1X,I4,1X,1PD22.15,2X,1PD15.8) */
/*<       END >*/
} /* monit_ */



/*<       SUBROUTINE ZTIME(N,X,IPIVOT) >*/
/* Subroutine */ int ztime_(integer *n, doublereal *x, integer *ipivot)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/*<       IMPLICIT         DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION X(N) >*/
/*<       INTEGER          IPIVOT(N) >*/

/* THIS ROUTINE MULTIPLIES THE VECTOR X BY THE CONSTRAINT MATRIX Z */

/*<       DO 10 I = 1,N >*/
    /* Parameter adjustments */
    --ipivot;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (IPIVOT(I) .NE. 0) X(I) = 0.D0 >*/
	if (ipivot[i__] != 0) {
	    x[i__] = 0.;
	}
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* ztime_ */



/*<       SUBROUTINE STPMAX(STEPMX,PE,SPE,N,X,P,IPIVOT,LOW,UP) >*/
/* Subroutine */ int stpmax_(doublereal *stepmx, doublereal *pe, doublereal *
	spe, integer *n, doublereal *x, doublereal *p, integer *ipivot,
	doublereal *low, doublereal *up)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal t;

/*<       IMPLICIT         DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION LOW(N),UP(N),X(N),P(N),STEPMX,PE,SPE,T >*/
/*<       INTEGER          IPIVOT(N) >*/

/* COMPUTE THE MAXIMUM ALLOWABLE STEP LENGTH */

/*<       SPE = STEPMX / PE >*/
    /* Parameter adjustments */
    --up;
    --low;
    --ipivot;
    --p;
    --x;

    /* Function Body */
    *spe = *stepmx / *pe;
/* SPE IS THE STANDARD (UNCONSTRAINED) MAX STEP */
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (IPIVOT(I) .NE. 0) GO TO 10 >*/
	if (ipivot[i__] != 0) {
	    goto L10;
	}
/*<          IF (P(I) .EQ. 0.D0) GO TO 10 >*/
	if (p[i__] == 0.) {
	    goto L10;
	}
/*<          IF (P(I) .GT. 0.D0) GO TO 5 >*/
	if (p[i__] > 0.) {
	    goto L5;
	}
/*<          T = LOW(I) - X(I) >*/
	t = low[i__] - x[i__];
/*<          IF (T .GT. SPE*P(I)) SPE = T / P(I) >*/
	if (t > *spe * p[i__]) {
	    *spe = t / p[i__];
	}
/*<          GO TO 10 >*/
	goto L10;
/*< 5        T = UP(I) - X(I) >*/
L5:
	t = up[i__] - x[i__];
/*<          IF (T .LT. SPE*P(I)) SPE = T / P(I) >*/
	if (t < *spe * p[i__]) {
	    *spe = t / p[i__];
	}
/*< 10    CONTINUE >*/
L10:
	;
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* stpmax_ */



/*<       SUBROUTINE MODZ(N,X,P,IPIVOT,EPSMCH,LOW,UP,FLAST,FNEW) >*/
/* Subroutine */ int modz_(integer *n, doublereal *x, doublereal *p, integer *
	ipivot, doublereal *epsmch, doublereal *low, doublereal *up,
	doublereal *flast, doublereal *fnew)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__;
    doublereal tol;

/*<       IMPLICIT         DOUBLE PRECISION (A-H,O-Z) >*/
/*<    >*/
/*<       INTEGER          IPIVOT(N) >*/

/* UPDATE THE CONSTRAINT MATRIX IF A NEW CONSTRAINT IS ENCOUNTERED */

/*<       DO 10 I = 1,N >*/
    /* Parameter adjustments */
    --up;
    --low;
    --ipivot;
    --p;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (IPIVOT(I) .NE. 0) GO TO 10 >*/
	if (ipivot[i__] != 0) {
	    goto L10;
	}
/*<          IF (P(I) .EQ. 0.D0) GO TO 10 >*/
	if (p[i__] == 0.) {
	    goto L10;
	}
/*<          IF (P(I) .GT. 0.D0) GO TO 5 >*/
	if (p[i__] > 0.) {
	    goto L5;
	}
/*<          TOL = 1.D1 * EPSMCH * (DABS(LOW(I)) + 1.D0) >*/
	tol = *epsmch * 10. * ((d__1 = low[i__], abs(d__1)) + 1.);
/*<          IF (X(I)-LOW(I) .GT. TOL) GO TO 10 >*/
	if (x[i__] - low[i__] > tol) {
	    goto L10;
	}
/*<          FLAST = FNEW >*/
	*flast = *fnew;
/*<          IPIVOT(I) = -1 >*/
	ipivot[i__] = -1;
/*<          X(I) = LOW(I) >*/
	x[i__] = low[i__];
/*<          GO TO 10 >*/
	goto L10;
/*< 5        TOL = 1.D1 * EPSMCH * (DABS(UP(I)) + 1.D0) >*/
L5:
	tol = *epsmch * 10. * ((d__1 = up[i__], abs(d__1)) + 1.);
/*<          IF (UP(I)-X(I) .GT. TOL) GO TO 10 >*/
	if (up[i__] - x[i__] > tol) {
	    goto L10;
	}
/*<          FLAST = FNEW >*/
	*flast = *fnew;
/*<          IPIVOT(I) = 1 >*/
	ipivot[i__] = 1;
/*<          X(I) = UP(I) >*/
	x[i__] = up[i__];
/*< 10    CONTINUE >*/
L10:
	;
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* modz_ */



/*<    >*/
/* Subroutine */ int cnvtst_(logical *conv, doublereal *alpha, doublereal *
	pnorm, doublereal *toleps, doublereal *xnorm, doublereal *difnew,
	doublereal *rtleps, doublereal *ftest, doublereal *gtg, doublereal *
	peps, doublereal *epsmch, doublereal *gtpnew, doublereal *fnew,
	doublereal *flast, doublereal *g, integer *ipivot, integer *n,
	doublereal *accrcy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    doublereal cmax;
    integer imax, i__;
    doublereal t;
    logical ltest;
    doublereal one;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       LOGICAL CONV,LTEST >*/
/*<       INTEGER IPIVOT(N) >*/
/*<    >*/

/* TEST FOR CONVERGENCE */

/*<       IMAX = 0 >*/
    /* Parameter adjustments */
    --ipivot;
    --g;

    /* Function Body */
    imax = 0;
/*<       CMAX = 0.D0 >*/
    cmax = 0.;
/*<       LTEST = FLAST - FNEW .LE. -5.D-1*GTPNEW >*/
    ltest = *flast - *fnew <= *gtpnew * -.5;
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (IPIVOT(I) .EQ. 0 .OR. IPIVOT(I) .EQ. 2) GO TO 10 >*/
	if (ipivot[i__] == 0 || ipivot[i__] == 2) {
	    goto L10;
	}
/*<          T = -IPIVOT(I)*G(I) >*/
	t = -ipivot[i__] * g[i__];
/*<          IF (T .GE. 0.D0) GO TO 10 >*/
	if (t >= 0.) {
	    goto L10;
	}
/*<          CONV = .FALSE. >*/
	*conv = FALSE_;
/*<          IF (LTEST) GO TO 10 >*/
	if (ltest) {
	    goto L10;
	}
/*<          IF (CMAX .LE. T) GO TO 10 >*/
	if (cmax <= t) {
	    goto L10;
	}
/*<          CMAX = T >*/
	cmax = t;
/*<          IMAX = I >*/
	imax = i__;
/*< 10    CONTINUE >*/
L10:
	;
    }
/*<       IF (IMAX .EQ. 0) GO TO 15 >*/
    if (imax == 0) {
	goto L15;
    }
/*<       IPIVOT(IMAX) = 0 >*/
    ipivot[imax] = 0;
/*<       FLAST = FNEW >*/
    *flast = *fnew;
/*<       RETURN >*/
    return 0;
/*< 15    CONTINUE >*/
L15:
/*<       CONV = .FALSE. >*/
    *conv = FALSE_;
/*<       ONE = 1.D0 >*/
    one = 1.;
/*<    >*/
    if ((*alpha * *pnorm >= *toleps * (one + *xnorm) || abs(*difnew) >= *
	    rtleps * *ftest || *gtg >= *peps * *ftest * *ftest) && *gtg >= *
	    accrcy * 1e-4 * *ftest * *ftest) {
	return 0;
    }
/*<       CONV = .TRUE. >*/
    *conv = TRUE_;

/* FOR DETAILS, SEE GILL, MURRAY, AND WRIGHT (1981, P. 308) AND */
/* FLETCHER (1981, P. 116).  THE MULTIPLIER TESTS (HERE, TESTING */
/* THE SIGN OF THE COMPONENTS OF THE GRADIENT) MAY STILL NEED TO */
/* MODIFIED TO INCORPORATE TOLERANCES FOR ZERO. */

/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* cnvtst_ */



/*<       SUBROUTINE CRASH(N,X,IPIVOT,LOW,UP,IER) >*/
/* Subroutine */ int crash_(integer *n, doublereal *x, integer *ipivot,
	doublereal *low, doublereal *up, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION X(N),LOW(N),UP(N) >*/
/*<       INTEGER IPIVOT(N) >*/

/* THIS INITIALIZES THE CONSTRAINT INFORMATION, AND ENSURES THAT THE */
/* INITIAL POINT SATISFIES  LOW <= X <= UP. */
/* THE CONSTRAINTS ARE CHECKED FOR CONSISTENCY. */

/*<       IER = 0 >*/
    /* Parameter adjustments */
    --up;
    --low;
    --ipivot;
    --x;

    /* Function Body */
    *ier = 0;
/*<       DO 30 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (X(I) .LT. LOW(I)) X(I) = LOW(I) >*/
	if (x[i__] < low[i__]) {
	    x[i__] = low[i__];
	}
/*<          IF (X(I) .GT. UP(I)) X(I) = UP(I) >*/
	if (x[i__] > up[i__]) {
	    x[i__] = up[i__];
	}
/*<          IPIVOT(I) = 0 >*/
	ipivot[i__] = 0;
/*<          IF (X(I) .EQ. LOW(I)) IPIVOT(I) = -1 >*/
	if (x[i__] == low[i__]) {
	    ipivot[i__] = -1;
	}
/*<          IF (X(I) .EQ. UP(I)) IPIVOT(I) = 1 >*/
	if (x[i__] == up[i__]) {
	    ipivot[i__] = 1;
	}
/*<          IF (UP(I) .EQ. LOW(I)) IPIVOT(I) = 2 >*/
	if (up[i__] == low[i__]) {
	    ipivot[i__] = 2;
	}
/*<          IF (LOW(I) .GT. UP(I)) IER = -I >*/
	if (low[i__] > up[i__]) {
	    *ier = -i__;
	}
/*< 30    CONTINUE >*/
/* L30: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* crash_ */


/* THE VECTORS SK AND YK, ALTHOUGH NOT IN THE CALL, */
/* ARE USED (VIA THEIR POSITION IN W) BY THE ROUTINE MSOLVE. */

/*<    >*/
/* Subroutine */ int modlnp_(integer *modet, doublereal *zsol, doublereal *gv,
	 doublereal *r__, doublereal *v, doublereal *diagb, doublereal *emat,
	doublereal *x, doublereal *g, doublereal *zk, integer *n, doublereal *
	w, integer *lw, integer *niter, integer *maxit, integer *nfeval,
	integer *nmodif, integer *nlincg, logical *upd1, doublereal *yksk,
	doublereal *gsk, doublereal *yrsr, logical *lreset, S_fp sfun,
	logical *bounds, integer *ipivot, doublereal *accrcy, doublereal *gtp,
	 doublereal *gnorm, doublereal *xnorm, void *pParam)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    doublereal beta;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);
    doublereal qold, qnew;
    extern /* Subroutine */ int ndia3_(integer *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, integer *);
    integer i__, k;
    doublereal alpha, delta;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
	    doublereal *, integer *), gtims_(doublereal *, doublereal *,
	    integer *, doublereal *, doublereal *, doublereal *, integer *,
	    S_fp, logical *, doublereal *, doublereal *, doublereal *, void *),
	    daxpy_(integer *, doublereal *, doublereal *, integer *,
	    doublereal *, integer *);
    logical first;
    extern /* Subroutine */ int ztime_(integer *, doublereal *, integer *);
    doublereal rzold, qtest, pr;
    extern /* Subroutine */ int negvec_(integer *, doublereal *);
    doublereal rz;
    extern /* Subroutine */ int initpc_(doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, logical *, doublereal *,
	    doublereal *, doublereal *, logical *), msolve_(doublereal *,
	    doublereal *, integer *, doublereal *, integer *, logical *,
	    doublereal *, doublereal *, doublereal *, logical *, logical *);
    doublereal rhsnrm, tol, vgv;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER MODET,N,NITER,IPIVOT(1) >*/
/*<       DOUBLE PRECISION ZSOL(N),G(N),GV(N),R(N),V(N),DIAGB(N),W(LW) >*/
/*<       DOUBLE PRECISION EMAT(N),ZK(N),X(N),ACCRCY >*/
/*<    >*/
/*<       DOUBLE PRECISION GNORM,XNORM >*/
/*<       DOUBLE PRECISION DDOT,DNRM2 >*/
/*<       LOGICAL FIRST,UPD1,LRESET,BOUNDS >*/
/*<       EXTERNAL SFUN >*/

/* THIS ROUTINE PERFORMS A PRECONDITIONED CONJUGATE-GRADIENT */
/* ITERATION IN ORDER TO SOLVE THE NEWTON EQUATIONS FOR A SEARCH */
/* DIRECTION FOR A TRUNCATED-NEWTON ALGORITHM.  WHEN THE VALUE OF THE */
/* QUADRATIC MODEL IS SUFFICIENTLY REDUCED, */
/* THE ITERATION IS TERMINATED. */

/* PARAMETERS */

/* MODET       - INTEGER WHICH CONTROLS AMOUNT OF OUTPUT */
/* ZSOL        - COMPUTED SEARCH DIRECTION */
/* G           - CURRENT GRADIENT */
/* GV,GZ1,V    - SCRATCH VECTORS */
/* R           - RESIDUAL */
/* DIAGB,EMAT  - DIAGONAL PRECONDITONING MATRIX */
/* NITER       - NONLINEAR ITERATION # */
/* FEVAL       - VALUE OF QUADRATIC FUNCTION */

/* ************************************************************* */
/* INITIALIZATION */
/* ************************************************************* */

/* GENERAL INITIALIZATION */

/*      IF (MODET .GT. 0) WRITE(*,800) */
/*<       IF (MAXIT .EQ. 0) RETURN >*/
    /* Parameter adjustments */
    --zk;
    --g;
    --x;
    --emat;
    --diagb;
    --v;
    --r__;
    --gv;
    --zsol;
    --w;
    --ipivot;

    /* Function Body */
    if (*maxit == 0) {
	return 0;
    }
/*<       FIRST = .TRUE. >*/
    first = TRUE_;
/*<       RHSNRM = GNORM >*/
    rhsnrm = *gnorm;
/*<       TOL = 1.D-12 >*/
    tol = 1e-12;
/*<       QOLD = 0.D0 >*/
    qold = 0.;

/* INITIALIZATION FOR PRECONDITIONED CONJUGATE-GRADIENT ALGORITHM */

/*<    >*/
    initpc_(&diagb[1], &emat[1], n, &w[1], lw, modet, upd1, yksk, gsk, yrsr,
	    lreset);
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          R(I) = -G(I) >*/
	r__[i__] = -g[i__];
/*<          V(I) = 0.D0 >*/
	v[i__] = 0.;
/*<          ZSOL(I) = 0.D0 >*/
	zsol[i__] = 0.;
/*< 10    CONTINUE >*/
/* L10: */
    }

/* ************************************************************ */
/* MAIN ITERATION */
/* ************************************************************ */

/*<       DO 30 K = 1,MAXIT >*/
    i__1 = *maxit;
    for (k = 1; k <= i__1; ++k) {
/*<          NLINCG = NLINCG + 1 >*/
	++(*nlincg);
/*         IF (MODET .GT. 1) WRITE(*,810) K */

/* CG ITERATION TO SOLVE SYSTEM OF EQUATIONS */

/*<          IF (BOUNDS) CALL ZTIME(N,R,IPIVOT) >*/
	if (*bounds) {
	    ztime_(n, &r__[1], &ipivot[1]);
	}
/*<    >*/
	msolve_(&r__[1], &zk[1], n, &w[1], lw, upd1, yksk, gsk, yrsr, lreset,
		&first);
/*<          IF (BOUNDS) CALL ZTIME(N,ZK,IPIVOT) >*/
	if (*bounds) {
	    ztime_(n, &zk[1], &ipivot[1]);
	}
/*<          RZ = DDOT(N,R,1,ZK,1) >*/
	rz = ddot_(n, &r__[1], &c__1, &zk[1], &c__1);
/*<          IF (RZ/RHSNRM .LT. TOL) GO TO 80 >*/
	if (rz / rhsnrm < tol) {
	    goto L80;
	}
/*<          IF (K .EQ. 1) BETA = 0.D0 >*/
	if (k == 1) {
	    beta = 0.;
	}
/*<          IF (K .GT. 1) BETA = RZ/RZOLD >*/
	if (k > 1) {
	    beta = rz / rzold;
	}
/*<          DO 20 I = 1,N >*/
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/*<             V(I) = ZK(I) + BETA*V(I) >*/
	    v[i__] = zk[i__] + beta * v[i__];
/*< 20       CONTINUE >*/
/* L20: */
	}
/*<          IF (BOUNDS) CALL ZTIME(N,V,IPIVOT) >*/
	if (*bounds) {
	    ztime_(n, &v[1], &ipivot[1]);
	}
/*<          CALL GTIMS(V,GV,N,X,G,W,LW,SFUN,FIRST,DELTA,ACCRCY,XNORM) >*/
	gtims_(&v[1], &gv[1], n, &x[1], &g[1], &w[1], lw, (S_fp)sfun, &first,
		&delta, accrcy, xnorm, pParam);
/*<          IF (BOUNDS) CALL ZTIME(N,GV,IPIVOT) >*/
	if (*bounds) {
	    ztime_(n, &gv[1], &ipivot[1]);
	}
/*<          NFEVAL = NFEVAL + 1 >*/
	++(*nfeval);
/*<          VGV = DDOT(N,V,1,GV,1) >*/
	vgv = ddot_(n, &v[1], &c__1, &gv[1], &c__1);
/*<          IF (VGV/RHSNRM .LT. TOL) GO TO 50 >*/
	if (vgv / rhsnrm < tol) {
	    goto L50;
	}
/*<          CALL NDIA3(N,EMAT,V,GV,R,VGV,MODET) >*/
	ndia3_(n, &emat[1], &v[1], &gv[1], &r__[1], &vgv, modet);

/* COMPUTE LINEAR STEP LENGTH */

/*<          ALPHA = RZ / VGV >*/
	alpha = rz / vgv;
/*         IF (MODET .GE. 1) WRITE(*,820) ALPHA */

/* COMPUTE CURRENT SOLUTION AND RELATED VECTORS */

/*<          CALL DAXPY(N,ALPHA,V,1,ZSOL,1) >*/
	daxpy_(n, &alpha, &v[1], &c__1, &zsol[1], &c__1);
/*<          CALL DAXPY(N,-ALPHA,GV,1,R,1) >*/
	d__1 = -alpha;
	daxpy_(n, &d__1, &gv[1], &c__1, &r__[1], &c__1);

/* TEST FOR CONVERGENCE */

/*<          GTP = DDOT(N,ZSOL,1,G,1) >*/
	*gtp = ddot_(n, &zsol[1], &c__1, &g[1], &c__1);
/*<          PR = DDOT(N,R,1,ZSOL,1) >*/
	pr = ddot_(n, &r__[1], &c__1, &zsol[1], &c__1);
/*<          QNEW = 5.D-1 * (GTP + PR) >*/
	qnew = (*gtp + pr) * .5;
/*<          QTEST = K * (1.D0 - QOLD/QNEW) >*/
	qtest = k * (1. - qold / qnew);
/*<          IF (QTEST .LT. 0.D0) GO TO 70 >*/
	if (qtest < 0.) {
	    goto L70;
	}
/*<          QOLD = QNEW >*/
	qold = qnew;
/*<          IF (QTEST .LE. 5.D-1) GO TO 70 >*/
	if (qtest <= .5) {
	    goto L70;
	}

/* PERFORM CAUTIONARY TEST */

/*<          IF (GTP .GT. 0) GO TO 40 >*/
	if (*gtp > 0.) {
	    goto L40;
	}
/*<          RZOLD = RZ >*/
	rzold = rz;
/*< 30    CONTINUE >*/
/* L30: */
    }

/* TERMINATE ALGORITHM */

/*<       K = K-1 >*/
    --k;
/*<       GO TO 70 >*/
    goto L70;

/* TRUNCATE ALGORITHM IN CASE OF AN EMERGENCY */

/*<  40   CONTINUE >*/
L40:
/* 40    IF (MODET .GE. -1) WRITE(*,830) K */
/*<       CALL DAXPY(N,-ALPHA,V,1,ZSOL,1) >*/
    d__1 = -alpha;
    daxpy_(n, &d__1, &v[1], &c__1, &zsol[1], &c__1);
/*<       GTP = DDOT(N,ZSOL,1,G,1) >*/
    *gtp = ddot_(n, &zsol[1], &c__1, &g[1], &c__1);
/*<       GO TO 90 >*/
    goto L90;
/*< 50    CONTINUE >*/
L50:
/*      IF (MODET .GT. -2) WRITE(*,840) */
/*< 60    IF (K .GT. 1) GO TO 70 >*/
/* L60: */
    if (k > 1) {
	goto L70;
    }
/*<       CALL MSOLVE(G,ZSOL,N,W,LW,UPD1,YKSK,GSK,YRSR,LRESET,FIRST) >*/
    msolve_(&g[1], &zsol[1], n, &w[1], lw, upd1, yksk, gsk, yrsr, lreset, &
	    first);
/*<       CALL NEGVEC(N,ZSOL) >*/
    negvec_(n, &zsol[1]);
/*<       IF (BOUNDS) CALL ZTIME(N,ZSOL,IPIVOT) >*/
    if (*bounds) {
	ztime_(n, &zsol[1], &ipivot[1]);
    }
/*<       GTP = DDOT(N,ZSOL,1,G,1) >*/
    *gtp = ddot_(n, &zsol[1], &c__1, &g[1], &c__1);
/*< 70    CONTINUE >*/
L70:
/*      IF (MODET .GE. -1) WRITE(*,850) K,RNORM */
/*<       GO TO 90 >*/
    goto L90;
/*< 80    CONTINUE >*/
L80:
/*      IF (MODET .GE. -1) WRITE(*,860) */
/*<       IF (K .GT. 1) GO TO 70 >*/
    if (k > 1) {
	goto L70;
    }
/*<       CALL DCOPY(N,G,1,ZSOL,1) >*/
    dcopy_(n, &g[1], &c__1, &zsol[1], &c__1);
/*<       CALL NEGVEC(N,ZSOL) >*/
    negvec_(n, &zsol[1]);
/*<       IF (BOUNDS) CALL ZTIME(N,ZSOL,IPIVOT) >*/
    if (*bounds) {
	ztime_(n, &zsol[1], &ipivot[1]);
    }
/*<       GTP = DDOT(N,ZSOL,1,G,1) >*/
    *gtp = ddot_(n, &zsol[1], &c__1, &g[1], &c__1);
/*<       GO TO 70 >*/
    goto L70;

/* STORE (OR RESTORE) DIAGONAL PRECONDITIONING */

/*< 90    CONTINUE >*/
L90:
/*<       CALL DCOPY(N,EMAT,1,DIAGB,1) >*/
    dcopy_(n, &emat[1], &c__1, &diagb[1], &c__1);
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(' ',//,' ENTERING MODLNP') */
/* 810   FORMAT(' ',//,' ### ITERATION ',I2,' ###') */
/* 820   FORMAT(' ALPHA',1PD16.8) */
/* 830   FORMAT(' G(T)Z POSITIVE AT ITERATION ',I2, */
/*     *     ' - TRUNCATING METHOD',/) */
/* 840   FORMAT(' ',10X,'HESSIAN NOT POSITIVE-DEFINITE') */
/* 850   FORMAT(' ',/,8X,'MODLAN TRUNCATED AFTER ',I3,' ITERATIONS', */
/*     *     '  RNORM = ',1PD14.6) */
/* 860   FORMAT(' PRECONDITIONING NOT POSITIVE-DEFINITE') */
/*<        END >*/
} /* modlnp_ */



/*<       SUBROUTINE NDIA3(N,E,V,GV,R,VGV,MODET) >*/
/* Subroutine */ int ndia3_(integer *n, doublereal *e, doublereal *v,
	doublereal *gv, doublereal *r__, doublereal *vgv, integer *modet)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);
    integer i__;
    doublereal vr;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION E(N),V(N),GV(N),R(N),VGV,VR,DDOT >*/

/* UPDATE THE PRECONDITIOING MATRIX BASED ON A DIAGONAL VERSION */
/* OF THE BFGS QUASI-NEWTON UPDATE. */

/*<       VR = DDOT(N,V,1,R,1) >*/
    /* Parameter adjustments */
    --r__;
    --gv;
    --v;
    --e;

    /* Function Body */
    vr = ddot_(n, &v[1], &c__1, &r__[1], &c__1);
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          E(I) = E(I) - R(I)*R(I)/VR + GV(I)*GV(I)/VGV >*/
	e[i__] = e[i__] - r__[i__] * r__[i__] / vr + gv[i__] * gv[i__] / *vgv;
/*<          IF (E(I) .GT. 1.D-6) GO TO 10 >*/
	if (e[i__] > 1e-6) {
	    goto L10;
	}
/*         IF (MODET .GT. -2) WRITE(*,800) E(I) */
/*<          E(I) = 1.D0 >*/
	e[i__] = 1.;
/*< 10    CONTINUE >*/
L10:
	;
    }
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(' *** EMAT NEGATIVE:  ',1PD16.8) */
/*<       END >*/
} /* ndia3_ */


/*      SERVICE ROUTINES FOR OPTIMIZATION */

/*<       SUBROUTINE NEGVEC(N,V) >*/
/* Subroutine */ int negvec_(integer *n, doublereal *v)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER N >*/
/*<       DOUBLE PRECISION V(N) >*/

/* NEGATIVE OF THE VECTOR V */

/*<       INTEGER I >*/
/*<       DO 10 I = 1,N >*/
    /* Parameter adjustments */
    --v;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          V(I) = -V(I) >*/
	v[i__] = -v[i__];
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* negvec_ */



/*<    >*/
/* Subroutine */ int lsout_(integer *iloc, integer *itest, doublereal *xmin,
	doublereal *fmin, doublereal *gmin, doublereal *xw, doublereal *fw,
	doublereal *gw, doublereal *u, doublereal *a, doublereal *b,
	doublereal *tol, doublereal *eps, doublereal *scxbd, doublereal *
	xlamda)
{
    doublereal ybnd, ya, yb, yu, yw;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<    >*/

/* ERROR PRINTOUTS FOR GETPTC */

/*<       DOUBLE PRECISION YA,YB,YBND,YW,YU >*/
/*<       YU = XMIN + U >*/
    yu = *xmin + *u;
/*<       YA = A + XMIN >*/
    ya = *a + *xmin;
/*<       YB = B + XMIN >*/
    yb = *b + *xmin;
/*<       YW = XW + XMIN >*/
    yw = *xw + *xmin;
/*<       YBND = SCXBD + XMIN >*/
    ybnd = *scxbd + *xmin;
/*      WRITE(*,800) */
/*      WRITE(*,810) TOL,EPS */
/*      WRITE(*,820) YA,YB */
/*      WRITE(*,830) YBND */
/*      WRITE(*,840) YW,FW,GW */
/*      WRITE(*,850) XMIN,FMIN,GMIN */
/*      WRITE(*,860) YU */
/*      WRITE(*,870) ILOC,ITEST */
/*<       RETURN >*/
    return 0;
/* 800   FORMAT(///' OUTPUT FROM LINEAR SEARCH') */
/* 810   FORMAT('  TOL AND EPS'/2D25.14) */
/* 820   FORMAT('  CURRENT UPPER AND LOWER BOUNDS'/2D25.14) */
/* 830   FORMAT('  STRICT UPPER BOUND'/D25.14) */
/* 840   FORMAT('  XW, FW, GW'/3D25.14) */
/* 850   FORMAT('  XMIN, FMIN, GMIN'/3D25.14) */
/* 860   FORMAT('  NEW ESTIMATE'/2D25.14) */
/* 870   FORMAT('  ILOC AND ITEST'/2I3) */
/*<       END >*/
} /* lsout_ */



/*<       DOUBLE PRECISION FUNCTION STEP1(FNEW,FM,GTP,SMAX) >*/
doublereal step1_(doublereal *fnew, doublereal *fm, doublereal *gtp,
	doublereal *smax)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    doublereal d__, alpha;
    extern doublereal mchpr1_(void);
    doublereal epsmch;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION FNEW,FM,GTP,SMAX >*/

/* ******************************************************** */
/* STEP1 RETURNS THE LENGTH OF THE INITIAL STEP TO BE TAKEN ALONG THE */
/* VECTOR P IN THE NEXT LINEAR SEARCH. */
/* ******************************************************** */

/*<       DOUBLE PRECISION ALPHA,D,EPSMCH >*/
/*<       DOUBLE PRECISION DABS,MCHPR1 >*/
/*<       EPSMCH = MCHPR1() >*/
    epsmch = mchpr1_();
/*<       D = DABS(FNEW-FM) >*/
    d__ = (d__1 = *fnew - *fm, abs(d__1));
/*<       ALPHA = 1.D0 >*/
    alpha = 1.;
/*<    >*/
    if (d__ * 2. <= -(*gtp) && d__ >= epsmch) {
	alpha = d__ * -2. / *gtp;
    }
/*<       IF (ALPHA .GE. SMAX) ALPHA = SMAX >*/
    if (alpha >= *smax) {
	alpha = *smax;
    }
/*<       STEP1 = ALPHA >*/
    ret_val = alpha;
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* step1_ */



/*<       DOUBLE PRECISION FUNCTION MCHPR1() >*/
doublereal mchpr1_(void)
{
    /* System generated locals */
    doublereal ret_val;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION X >*/

/* RETURNS THE VALUE OF EPSMCH, WHERE EPSMCH IS THE SMALLEST POSSIBLE */
/* REAL NUMBER SUCH THAT 1.0 + EPSMCH .GT. 1.0 */

/* FOR VAX */

/*      MCHPR1 = 1.D-17 */
/*<        MCHPR1 = 1.D-18 >*/
    ret_val = 1e-18;

/* FOR SUN */

/*     MCHPR1 = 1.0842021724855D-19 */
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* mchpr1_ */



/*<    >*/
/* Subroutine */ int chkucp_(integer *lwtest, integer *maxfun, integer *nwhy,
	integer *n, doublereal *alpha, doublereal *epsmch, doublereal *eta,
	doublereal *peps, doublereal *rteps, doublereal *rtol, doublereal *
	rtolsq, doublereal *stepmx, doublereal *test, doublereal *xtol,
	doublereal *xnorm, doublereal *x, integer *lw, doublereal *small,
	doublereal *tiny, doublereal *accrcy)
{
    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    extern doublereal dnrm2_(integer *, doublereal *, integer *), mchpr1_(
	    void);

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER LW,LWTEST,MAXFUN,NWHY,N >*/
/*<    >*/
/*<       DOUBLE PRECISION X(N) >*/

/* CHECKS PARAMETERS AND SETS CONSTANTS WHICH ARE COMMON TO BOTH */
/* DERIVATIVE AND NON-DERIVATIVE ALGORITHMS */

/*<       DOUBLE PRECISION DABS,DSQRT,MCHPR1 >*/
/*<       EPSMCH = MCHPR1() >*/
    /* Parameter adjustments */
    --x;

    /* Function Body */
    *epsmch = mchpr1_();
/*<       SMALL = EPSMCH*EPSMCH >*/
    *small = *epsmch * *epsmch;
/*<       TINY = SMALL >*/
    *tiny = *small;
/*<       NWHY = -1 >*/
    *nwhy = -1;
/*<       RTEPS = DSQRT(EPSMCH) >*/
    *rteps = sqrt(*epsmch);
/*<       RTOL = XTOL >*/
    *rtol = *xtol;
/*<       IF (DABS(RTOL) .LT. ACCRCY) RTOL = 1.D1*RTEPS >*/
    if (abs(*rtol) < *accrcy) {
	*rtol = *rteps * 10.;
    }

/* CHECK FOR ERRORS IN THE INPUT PARAMETERS */

/*<    >*/
    if (*lw < *lwtest || *n < 1 || *rtol < 0. || *eta >= 1. || *eta < 0. || *
	    stepmx < *rtol || *maxfun < 1) {
	return 0;
    }
/*<       NWHY = 0 >*/
    *nwhy = 0;

/* SET CONSTANTS FOR LATER */

/*<       RTOLSQ = RTOL*RTOL >*/
    *rtolsq = *rtol * *rtol;
/*<       PEPS = ACCRCY**0.6666D0 >*/
    *peps = pow_dd(accrcy, &c_b135);
/*<       XNORM = DNRM2(N,X,1) >*/
    *xnorm = dnrm2_(n, &x[1], &c__1);
/*<       ALPHA = 0.D0 >*/
    *alpha = 0.;
/*<       TEST = 0.D0 >*/
    *test = 0.;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* chkucp_ */



/*<    >*/
/* Subroutine */ int setucr_(doublereal *small, integer *nftotl, integer *
	niter, integer *n, doublereal *f, doublereal *fnew, doublereal *fm,
	doublereal *gtg, doublereal *oldf, S_fp sfun, doublereal *g,
	doublereal *x, void *pParam)
{
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);

/*<       IMPLICIT         DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER          NFTOTL,NITER,N >*/
/*<       DOUBLE PRECISION F,FNEW,FM,GTG,OLDF,SMALL >*/
/*<       DOUBLE PRECISION G(N),X(N) >*/
/*<       EXTERNAL         SFUN >*/

/* CHECK INPUT PARAMETERS, COMPUTE THE INITIAL FUNCTION VALUE, SET */
/* CONSTANTS FOR THE SUBSEQUENT MINIMIZATION */

/*<       FM = F >*/
    /* Parameter adjustments */
    --x;
    --g;

    /* Function Body */
    *fm = *f;

/* COMPUTE THE INITIAL FUNCTION VALUE */

/*<       CALL SFUN(N,X,FNEW,G) >*/
    (*sfun)(n, &x[1], fnew, &g[1],pParam);
/*<       NFTOTL = 1 >*/
    *nftotl = 1;

/* SET CONSTANTS FOR LATER */

/*<       NITER = 0 >*/
    *niter = 0;
/*<       OLDF = FNEW >*/
    *oldf = *fnew;
/*<       GTG = DDOT(N,G,1,G,1) >*/
    *gtg = ddot_(n, &g[1], &c__1, &g[1], &c__1);
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* setucr_ */



/*<       SUBROUTINE GTIMS(V,GV,N,X,G,W,LW,SFUN,FIRST,DELTA,ACCRCY,XNORM) >*/
/* Subroutine */ int gtims_(doublereal *v, doublereal *gv, integer *n,
	doublereal *x, doublereal *g, doublereal *w, integer *lw, S_fp sfun,
	logical *first, doublereal *delta, doublereal *accrcy, doublereal *
	xnorm, void * pParam)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal dinv, f;
    integer i__, ihg;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION V(N),GV(N),DINV,DELTA,G(N) >*/
/*<       DOUBLE PRECISION F,X(N),W(LW),ACCRCY,DSQRT,XNORM >*/
/*<       LOGICAL FIRST >*/
/*<       EXTERNAL SFUN >*/
/*<    >*/

/* THIS ROUTINE COMPUTES THE PRODUCT OF THE MATRIX G TIMES THE VECTOR */
/* V AND STORES THE RESULT IN THE VECTOR GV (FINITE-DIFFERENCE VERSION) */

/*<       IF (.NOT. FIRST) GO TO 20 >*/
    /* Parameter adjustments */
    --g;
    --x;
    --gv;
    --v;
    --w;

    /* Function Body */
    if (! (*first)) {
	goto L20;
    }
/*<       DELTA = DSQRT(ACCRCY)*(1.D0+XNORM) >*/
    *delta = sqrt(*accrcy) * (*xnorm + 1.);
/*<       FIRST = .FALSE. >*/
    *first = FALSE_;
/*< 20    CONTINUE >*/
L20:
/*<       DINV = 1.D0/DELTA >*/
    dinv = 1. / *delta;
/*<       IHG = LHG >*/
    ihg = subscr_2.lhg;
/*<       DO 30 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          W(IHG) = X(I) + DELTA*V(I) >*/
	w[ihg] = x[i__] + *delta * v[i__];
/*<          IHG = IHG + 1 >*/
	++ihg;
/*< 30    CONTINUE >*/
/* L30: */
    }
/*<       CALL SFUN(N,W(LHG),F,GV) >*/
    (*sfun)(n, &w[subscr_2.lhg], &f, &gv[1], pParam);
/*<       DO 40 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          GV(I) = (GV(I) - G(I))*DINV >*/
	gv[i__] = (gv[i__] - g[i__]) * dinv;
/*< 40    CONTINUE >*/
/* L40: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* gtims_ */



/*<    >*/
/* Subroutine */ int msolve_(doublereal *g, doublereal *y, integer *n,
	doublereal *w, integer *lw, logical *upd1, doublereal *yksk,
	doublereal *gsk, doublereal *yrsr, logical *lreset, logical *first)
{
    extern /* Subroutine */ int mslv_(doublereal *, doublereal *, integer *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     doublereal *, doublereal *, doublereal *, logical *, logical *);

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION G(N),Y(N),W(LW),YKSK,GSK,YRSR >*/
/*<       LOGICAL UPD1,LRESET,FIRST >*/

/* THIS ROUTINE SETS UPT THE ARRAYS FOR MSLV */

/*<    >*/
/*<    >*/
    /* Parameter adjustments */
    --y;
    --g;
    --w;

    /* Function Body */
    mslv_(&g[1], &y[1], n, &w[subscr_2.lsk], &w[subscr_2.lyk], &w[
	    subscr_2.ldiagb], &w[subscr_2.lsr], &w[subscr_2.lyr], &w[
	    subscr_2.lhyr], &w[subscr_2.lhg], &w[subscr_2.lhyk], upd1, yksk,
	    gsk, yrsr, lreset, first);
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* msolve_ */

/*<    >*/
/* Subroutine */ int mslv_(doublereal *g, doublereal *y, integer *n,
	doublereal *sk, doublereal *yk, doublereal *diagb, doublereal *sr,
	doublereal *yr, doublereal *hyr, doublereal *hg, doublereal *hyk,
	logical *upd1, doublereal *yksk, doublereal *gsk, doublereal *yrsr,
	logical *lreset, logical *first)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);
    doublereal ghyk, ghyr, yksr;
    integer i__;
    doublereal ykhyk, ykhyr, yrhyr, rdiagb;
    extern /* Subroutine */ int ssbfgs_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal one, gsr;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION G(N),Y(N) >*/

/* THIS ROUTINE ACTS AS A PRECONDITIONING STEP FOR THE */
/* LINEAR CONJUGATE-GRADIENT ROUTINE.  IT IS ALSO THE */
/* METHOD OF COMPUTING THE SEARCH DIRECTION FROM THE */
/* GRADIENT FOR THE NON-LINEAR CONJUGATE-GRADIENT CODE. */
/* IT REPRESENTS A TWO-STEP SELF-SCALED BFGS FORMULA. */

/*<    >*/
/*<    >*/
/*<       LOGICAL LRESET,UPD1,FIRST >*/
/*<       IF (UPD1) GO TO 100 >*/
    /* Parameter adjustments */
    --hyk;
    --hg;
    --hyr;
    --yr;
    --sr;
    --diagb;
    --yk;
    --sk;
    --y;
    --g;

    /* Function Body */
    if (*upd1) {
	goto L100;
    }
/*<       ONE = 1.D0 >*/
    one = 1.;
/*<       GSK = DDOT(N,G,1,SK,1) >*/
    *gsk = ddot_(n, &g[1], &c__1, &sk[1], &c__1);
/*<       IF (LRESET) GO TO 60 >*/
    if (*lreset) {
	goto L60;
    }

/* COMPUTE HG AND HY WHERE H IS THE INVERSE OF THE DIAGONALS */

/*<       DO 57 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          RDIAGB = 1.0D0/DIAGB(I) >*/
	rdiagb = 1. / diagb[i__];
/*<          HG(I) = G(I)*RDIAGB >*/
	hg[i__] = g[i__] * rdiagb;
/*<          IF (FIRST) HYK(I) = YK(I)*RDIAGB >*/
	if (*first) {
	    hyk[i__] = yk[i__] * rdiagb;
	}
/*<          IF (FIRST) HYR(I) = YR(I)*RDIAGB >*/
	if (*first) {
	    hyr[i__] = yr[i__] * rdiagb;
	}
/*< 57    CONTINUE >*/
/* L57: */
    }
/*<       IF (FIRST) YKSR = DDOT(N,YK,1,SR,1) >*/
    if (*first) {
	yksr = ddot_(n, &yk[1], &c__1, &sr[1], &c__1);
    }
/*<       IF (FIRST) YKHYR = DDOT(N,YK,1,HYR,1) >*/
    if (*first) {
	ykhyr = ddot_(n, &yk[1], &c__1, &hyr[1], &c__1);
    }
/*<       GSR = DDOT(N,G,1,SR,1) >*/
    gsr = ddot_(n, &g[1], &c__1, &sr[1], &c__1);
/*<       GHYR = DDOT(N,G,1,HYR,1) >*/
    ghyr = ddot_(n, &g[1], &c__1, &hyr[1], &c__1);
/*<       IF (FIRST) YRHYR = DDOT(N,YR,1,HYR,1) >*/
    if (*first) {
	yrhyr = ddot_(n, &yr[1], &c__1, &hyr[1], &c__1);
    }
/*<    >*/
    ssbfgs_(n, &one, &sr[1], &yr[1], &hg[1], &hyr[1], yrsr, &yrhyr, &gsr, &
	    ghyr, &hg[1]);
/*<    >*/
    if (*first) {
	ssbfgs_(n, &one, &sr[1], &yr[1], &hyk[1], &hyr[1], yrsr, &yrhyr, &
		yksr, &ykhyr, &hyk[1]);
    }
/*<       YKHYK = DDOT(N,HYK,1,YK,1) >*/
    ykhyk = ddot_(n, &hyk[1], &c__1, &yk[1], &c__1);
/*<       GHYK = DDOT(N,HYK,1,G,1) >*/
    ghyk = ddot_(n, &hyk[1], &c__1, &g[1], &c__1);
/*<    >*/
    ssbfgs_(n, &one, &sk[1], &yk[1], &hg[1], &hyk[1], yksk, &ykhyk, gsk, &
	    ghyk, &y[1]);
/*<       RETURN >*/
    return 0;
/*< 60    CONTINUE >*/
L60:

/* COMPUTE GH AND HY WHERE H IS THE INVERSE OF THE DIAGONALS */

/*<       DO 65 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          RDIAGB = 1.D0/DIAGB(I) >*/
	rdiagb = 1. / diagb[i__];
/*<          HG(I) = G(I)*RDIAGB >*/
	hg[i__] = g[i__] * rdiagb;
/*<          IF (FIRST) HYK(I) = YK(I)*RDIAGB >*/
	if (*first) {
	    hyk[i__] = yk[i__] * rdiagb;
	}
/*< 65    CONTINUE >*/
/* L65: */
    }
/*<       IF (FIRST) YKHYK = DDOT(N,YK,1,HYK,1) >*/
    if (*first) {
	ykhyk = ddot_(n, &yk[1], &c__1, &hyk[1], &c__1);
    }
/*<       GHYK = DDOT(N,G,1,HYK,1) >*/
    ghyk = ddot_(n, &g[1], &c__1, &hyk[1], &c__1);
/*<    >*/
    ssbfgs_(n, &one, &sk[1], &yk[1], &hg[1], &hyk[1], yksk, &ykhyk, gsk, &
	    ghyk, &y[1]);
/*<       RETURN >*/
    return 0;
/*< 100   CONTINUE >*/
L100:
/*<       DO 110 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*< 110      Y(I) = G(I) / DIAGB(I) >*/
/* L110: */
	y[i__] = g[i__] / diagb[i__];
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* mslv_ */



/*<    >*/
/* Subroutine */ int ssbfgs_(integer *n, doublereal *gamma, doublereal *sj,
	doublereal *yj, doublereal *hjv, doublereal *hjyj, doublereal *yjsj,
	doublereal *yjhyj, doublereal *vsj, doublereal *vhyj, doublereal *
	hjp1v)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    doublereal beta;
    integer i__;
    doublereal delta;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER N >*/
/*<       DOUBLE PRECISION GAMMA,YJSJ,YJHYJ,VSJ,VHYJ >*/
/*<       DOUBLE PRECISION SJ(N),YJ(N),HJV(N),HJYJ(N),HJP1V(N) >*/

/* SELF-SCALED BFGS */

/*<       INTEGER I >*/
/*<       DOUBLE PRECISION BETA,DELTA >*/
/*<    >*/
    /* Parameter adjustments */
    --hjp1v;
    --hjyj;
    --hjv;
    --yj;
    --sj;

    /* Function Body */
    delta = (*gamma * *yjhyj / *yjsj + 1.) * *vsj / *yjsj - *gamma * *vhyj / *
	    yjsj;
/*<       BETA = -GAMMA*VSJ/YJSJ >*/
    beta = -(*gamma) * *vsj / *yjsj;
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          HJP1V(I) = GAMMA*HJV(I) + DELTA*SJ(I) + BETA*HJYJ(I) >*/
	hjp1v[i__] = *gamma * hjv[i__] + delta * sj[i__] + beta * hjyj[i__];
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* ssbfgs_ */


/* ROUTINES TO INITIALIZE PRECONDITIONER */

/*<    >*/
/* Subroutine */ int initpc_(doublereal *diagb, doublereal *emat, integer *n,
	doublereal *w, integer *lw, integer *modet, logical *upd1, doublereal
	*yksk, doublereal *gsk, doublereal *yrsr, logical *lreset)
{
    extern /* Subroutine */ int initp3_(doublereal *, doublereal *, integer *,
	     logical *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, logical *);

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DOUBLE PRECISION DIAGB(N),EMAT(N),W(LW) >*/
/*<       DOUBLE PRECISION YKSK,GSK,YRSR >*/
/*<       LOGICAL LRESET,UPD1 >*/
/*<    >*/
/*<    >*/
    /* Parameter adjustments */
    --emat;
    --diagb;
    --w;

    /* Function Body */
    initp3_(&diagb[1], &emat[1], n, lreset, yksk, yrsr, &w[subscr_2.lhyk], &w[
	    subscr_2.lsk], &w[subscr_2.lyk], &w[subscr_2.lsr], &w[
	    subscr_2.lyr], modet, upd1);
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* initpc_ */

/*<    >*/
/* Subroutine */ int initp3_(doublereal *diagb, doublereal *emat, integer *n,
	logical *lreset, doublereal *yksk, doublereal *yrsr, doublereal *bsk,
	doublereal *sk, doublereal *yk, doublereal *sr, doublereal *yr,
	integer *modet, logical *upd1)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    doublereal cond;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);
    doublereal srds, yrsk;
    integer i__;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
	    doublereal *, integer *);
    doublereal d1, dn, td, sds;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<    >*/
/*<       LOGICAL LRESET,UPD1 >*/
/*<       IF (UPD1) GO TO 90 >*/
    /* Parameter adjustments */
    --yr;
    --sr;
    --yk;
    --sk;
    --bsk;
    --emat;
    --diagb;

    /* Function Body */
    if (*upd1) {
	goto L90;
    }
/*<       IF (LRESET) GO TO 60 >*/
    if (*lreset) {
	goto L60;
    }
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          BSK(I) = DIAGB(I)*SR(I) >*/
	bsk[i__] = diagb[i__] * sr[i__];
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       SDS = DDOT(N,SR,1,BSK,1) >*/
    sds = ddot_(n, &sr[1], &c__1, &bsk[1], &c__1);
/*<       SRDS = DDOT(N,SK,1,BSK,1) >*/
    srds = ddot_(n, &sk[1], &c__1, &bsk[1], &c__1);
/*<       YRSK = DDOT(N,YR,1,SK,1) >*/
    yrsk = ddot_(n, &yr[1], &c__1, &sk[1], &c__1);
/*<       DO 20 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          TD = DIAGB(I) >*/
	td = diagb[i__];
/*<          BSK(I) = TD*SK(I) - BSK(I)*SRDS/SDS+YR(I)*YRSK/YRSR >*/
	bsk[i__] = td * sk[i__] - bsk[i__] * srds / sds + yr[i__] * yrsk / *
		yrsr;
/*<          EMAT(I) = TD-TD*TD*SR(I)*SR(I)/SDS+YR(I)*YR(I)/YRSR >*/
	emat[i__] = td - td * td * sr[i__] * sr[i__] / sds + yr[i__] * yr[i__]
		 / *yrsr;
/*< 20    CONTINUE >*/
/* L20: */
    }
/*<       SDS = DDOT(N,SK,1,BSK,1) >*/
    sds = ddot_(n, &sk[1], &c__1, &bsk[1], &c__1);
/*<       DO 30 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          EMAT(I) = EMAT(I) - BSK(I)*BSK(I)/SDS+YK(I)*YK(I)/YKSK >*/
	emat[i__] = emat[i__] - bsk[i__] * bsk[i__] / sds + yk[i__] * yk[i__]
		/ *yksk;
/*< 30    CONTINUE >*/
/* L30: */
    }
/*<       GO TO 110 >*/
    goto L110;
/*< 60    CONTINUE >*/
L60:
/*<       DO 70 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          BSK(I) = DIAGB(I)*SK(I) >*/
	bsk[i__] = diagb[i__] * sk[i__];
/*< 70    CONTINUE >*/
/* L70: */
    }
/*<       SDS = DDOT(N,SK,1,BSK,1) >*/
    sds = ddot_(n, &sk[1], &c__1, &bsk[1], &c__1);
/*<       DO 80 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          TD = DIAGB(I) >*/
	td = diagb[i__];
/*<          EMAT(I) = TD - TD*TD*SK(I)*SK(I)/SDS + YK(I)*YK(I)/YKSK >*/
	emat[i__] = td - td * td * sk[i__] * sk[i__] / sds + yk[i__] * yk[i__]
		 / *yksk;
/*< 80    CONTINUE >*/
/* L80: */
    }
/*<       GO TO 110 >*/
    goto L110;
/*< 90    CONTINUE >*/
L90:
/*<       CALL DCOPY(N,DIAGB,1,EMAT,1) >*/
    dcopy_(n, &diagb[1], &c__1, &emat[1], &c__1);
/*< 110   CONTINUE >*/
L110:
/*<       IF (MODET .LT. 1) RETURN >*/
    if (*modet < 1) {
	return 0;
    }
/*<       D1 = EMAT(1) >*/
    d1 = emat[1];
/*<       DN = EMAT(1) >*/
    dn = emat[1];
/*<       DO 120 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          IF (EMAT(I) .LT. D1) D1 = EMAT(I) >*/
	if (emat[i__] < d1) {
	    d1 = emat[i__];
	}
/*<          IF (EMAT(I) .GT. DN) DN = EMAT(I) >*/
	if (emat[i__] > dn) {
	    dn = emat[i__];
	}
/*< 120   CONTINUE >*/
/* L120: */
    }
/*<       COND = DN/D1 >*/
    cond = dn / d1;
/*      WRITE(*,800) D1,DN,COND */
/* 800   FORMAT(' ',//8X,'DMIN =',1PD12.4,'  DMAX =',1PD12.4, */
/*     *     ' COND =',1PD12.4,/) */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* initp3_ */



/*<       SUBROUTINE SETPAR(N) >*/
/* Subroutine */ int setpar_(integer *n)
{
    integer i__;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER LSUB(14) >*/
/*<       COMMON/SUBSCR/ LSUB,LWTEST >*/

/* SET UP PARAMETERS FOR THE OPTIMIZATION ROUTINE */

/*<       DO 10 I = 1,14 >*/
    for (i__ = 1; i__ <= 14; ++i__) {
/*<           LSUB(I) = (I-1)*N + 1 >*/
	subscr_3.lsub[i__ - 1] = (i__ - 1) * *n + 1;
/*< 10    CONTINUE >*/
/* L10: */
    }
/*<       LWTEST = LSUB(14) + N - 1 >*/
    subscr_3.lwtest = subscr_3.lsub[13] + *n - 1;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* setpar_ */


/*      LINE SEARCH ALGORITHMS OF GILL AND MURRAY */

/*<    >*/
/* Subroutine */ int linder_(integer *n, S_fp sfun, doublereal *small,
	doublereal *epsmch, doublereal *reltol, doublereal *abstol,
	doublereal *tnytol, doublereal *eta, doublereal *sftbnd, doublereal *
	xbnd, doublereal *p, doublereal *gtp, doublereal *x, doublereal *f,
	doublereal *alpha, doublereal *g, integer *nftotl, integer *iflag,
	doublereal *w, integer *lw, void *pParam)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal oldf, fmin, gmin;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *,
	    integer *);
    integer numf;
    doublereal step, xmin, a, b, e;
    integer i__, l;
    doublereal u;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *,
	    doublereal *, integer *);
    integer itcnt;
    doublereal b1;
    integer itest, nprnt;
    extern /* Subroutine */ int lsout_(integer *, integer *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal gtest1, gtest2;
    integer lg;
    doublereal fu, gu, fw, gw;
    integer lx;
    logical braktd;
    doublereal ualpha, factor, scxbnd, xw;
    extern /* Subroutine */ int getptc_(doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, logical *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *);
    doublereal fpresn;
    integer ientry;
    doublereal rtsmll;
    integer lsprnt;
    doublereal big, tol, rmu;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER N,NFTOTL,IFLAG,LW >*/
/*<    >*/
/*<       DOUBLE PRECISION P(N),X(N),G(N),W(LW) >*/


/*<       INTEGER I,IENTRY,ITEST,L,LG,LX,NUMF,ITCNT >*/
/*<    >*/
/*<       LOGICAL BRAKTD >*/

/*      THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE */
/*      CALLED WITHIN LINDER */

/*<       DOUBLE PRECISION DDOT,DSQRT >*/
/*<       EXTERNAL SFUN >*/

/*      ALLOCATE THE ADDRESSES FOR LOCAL WORKSPACE */

/*<       LX = 1 >*/
    /* Parameter adjustments */
    --g;
    --x;
    --p;
    --w;

    /* Function Body */
    lx = 1;
/*<       LG = LX + N >*/
    lg = lx + *n;
/*<       LSPRNT = 0 >*/
    lsprnt = 0;
/*<       NPRNT  = 10000 >*/
    nprnt = 10000;
/*<       RTSMLL = DSQRT(SMALL) >*/
    rtsmll = sqrt(*small);
/*<       BIG = 1.D0/SMALL >*/
    big = 1. / *small;
/*<       ITCNT = 0 >*/
    itcnt = 0;

/*      SET THE ESTIMATED RELATIVE PRECISION IN F(X). */

/*<       FPRESN = 10.D0*EPSMCH >*/
    fpresn = *epsmch * 10.;
/*<       NUMF = 0 >*/
    numf = 0;
/*<       U = ALPHA >*/
    u = *alpha;
/*<       FU = F >*/
    fu = *f;
/*<       FMIN = F >*/
    fmin = *f;
/*<       GU = GTP >*/
    gu = *gtp;
/*<       RMU = 1.0D-4 >*/
    rmu = 1e-4;

/*      FIRST ENTRY SETS UP THE INITIAL INTERVAL OF UNCERTAINTY. */

/*<       IENTRY = 1 >*/
    ientry = 1;
/*< 10    CONTINUE >*/
L10:

/* TEST FOR TOO MANY ITERATIONS */

/*<       ITCNT = ITCNT + 1 >*/
    ++itcnt;
/*<       IFLAG = 1 >*/
    *iflag = 1;
/*<       IF (ITCNT .GT. 20) GO TO 50 >*/
    if (itcnt > 20) {
	goto L50;
    }
/*<       IFLAG = 0 >*/
    *iflag = 0;
/*<    >*/
    getptc_(&big, small, &rtsmll, reltol, abstol, tnytol, &fpresn, eta, &rmu,
	    xbnd, &u, &fu, &gu, &xmin, &fmin, &gmin, &xw, &fw, &gw, &a, &b, &
	    oldf, &b1, &scxbnd, &e, &step, &factor, &braktd, &gtest1, &gtest2,
	     &tol, &ientry, &itest);
/* LSOUT */
/*<    >*/
    if (lsprnt >= nprnt) {
	lsout_(&ientry, &itest, &xmin, &fmin, &gmin, &xw, &fw, &gw, &u, &a, &
		b, &tol, reltol, &scxbnd, xbnd);
    }

/*      IF ITEST=1, THE ALGORITHM REQUIRES THE FUNCTION VALUE TO BE */
/*      CALCULATED. */

/*<       IF (ITEST .NE. 1) GO TO 30 >*/
    if (itest != 1) {
	goto L30;
    }
/*<       UALPHA = XMIN + U >*/
    ualpha = xmin + u;
/*<       L = LX >*/
    l = lx;
/*<       DO 20 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          W(L) = X(I) + UALPHA*P(I) >*/
	w[l] = x[i__] + ualpha * p[i__];
/*<          L = L + 1 >*/
	++l;
/*< 20    CONTINUE >*/
/* L20: */
    }
/*<       CALL SFUN(N,W(LX),FU,W(LG)) >*/
    (*sfun)(n, &w[lx], &fu, &w[lg], pParam);
/*<       NUMF = NUMF + 1 >*/
    ++numf;
/*<       GU = DDOT(N,W(LG),1,P,1) >*/
    gu = ddot_(n, &w[lg], &c__1, &p[1], &c__1);

/*      THE GRADIENT VECTOR CORRESPONDING TO THE BEST POINT IS */
/*      OVERWRITTEN IF FU IS LESS THAN FMIN AND FU IS SUFFICIENTLY */
/*      LOWER THAN F AT THE ORIGIN. */

/*<    >*/
    if (fu <= fmin && fu <= oldf - ualpha * gtest1) {
	dcopy_(n, &w[lg], &c__1, &g[1], &c__1);
    }
/*<       GOTO 10 >*/
    goto L10;

/*      IF ITEST=2 OR 3 A LOWER POINT COULD NOT BE FOUND */

/*< 30    CONTINUE >*/
L30:
/*<       NFTOTL = NUMF >*/
    *nftotl = numf;
/*<       IFLAG = 1 >*/
    *iflag = 1;
/*<       IF (ITEST .NE. 0) GO TO 50 >*/
    if (itest != 0) {
	goto L50;
    }

/*      IF ITEST=0 A SUCCESSFUL SEARCH HAS BEEN MADE */

/*<       IFLAG = 0 >*/
    *iflag = 0;
/*<       F = FMIN >*/
    *f = fmin;
/*<       ALPHA = XMIN >*/
    *alpha = xmin;
/*<       DO 40 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<          X(I) = X(I) + ALPHA*P(I) >*/
	x[i__] += *alpha * p[i__];
/*< 40    CONTINUE >*/
/* L40: */
    }
/*< 50    RETURN >*/
L50:
    return 0;
/*<       END >*/
} /* linder_ */



/*<    >*/
/* Subroutine */ int getptc_(doublereal *big, doublereal *small, doublereal *
	rtsmll, doublereal *reltol, doublereal *abstol, doublereal *tnytol,
	doublereal *fpresn, doublereal *eta, doublereal *rmu, doublereal *
	xbnd, doublereal *u, doublereal *fu, doublereal *gu, doublereal *xmin,
	 doublereal *fmin, doublereal *gmin, doublereal *xw, doublereal *fw,
	doublereal *gw, doublereal *a, doublereal *b, doublereal *oldf,
	doublereal *b1, doublereal *scxbnd, doublereal *e, doublereal *step,
	doublereal *factor, logical *braktd, doublereal *gtest1, doublereal *
	gtest2, doublereal *tol, integer *ientry, integer *itest)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal half, abgw, absr, five, zero, p, q, r__, s, scale, denom,
	    three, a1, d1, d2, sumsq, point1, abgmin, chordm, eleven, chordu;
    logical convrg;
    doublereal xmidpt, twotol, one;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       LOGICAL BRAKTD >*/
/*<       INTEGER IENTRY,ITEST >*/
/*<    >*/

/* ************************************************************ */
/* GETPTC, AN ALGORITHM FOR FINDING A STEPLENGTH, CALLED REPEATEDLY BY */
/* ROUTINES WHICH REQUIRE A STEP LENGTH TO BE COMPUTED USING CUBIC */
/* INTERPOLATION. THE PARAMETERS CONTAIN INFORMATION ABOUT THE INTERVAL */
/* IN WHICH A LOWER POINT IS TO BE FOUND AND FROM THIS GETPTC COMPUTES A */
/* POINT AT WHICH THE FUNCTION CAN BE EVALUATED BY THE CALLING PROGRAM. */
/* THE VALUE OF THE INTEGER PARAMETERS IENTRY DETERMINES THE PATH TAKEN */
/* THROUGH THE CODE. */
/* ************************************************************ */

/*<       LOGICAL CONVRG >*/
/*<    >*/
/*<       DOUBLE PRECISION ZERO, POINT1,HALF,ONE,THREE,FIVE,ELEVEN >*/

/* THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE CALLED */
/* WITHIN GETPTC */

/*<       DOUBLE PRECISION DABS, DSQRT >*/

/*<       ZERO = 0.D0 >*/
    zero = 0.;
/*<       POINT1 = 1.D-1 >*/
    point1 = .1;
/*<       HALF = 5.D-1 >*/
    half = .5;
/*<       ONE = 1.D0 >*/
    one = 1.;
/*<       THREE = 3.D0 >*/
    three = 3.;
/*<       FIVE = 5.D0 >*/
    five = 5.;
/*<       ELEVEN = 11.D0 >*/
    eleven = 11.;

/*      BRANCH TO APPROPRIATE SECTION OF CODE DEPENDING ON THE */
/*      VALUE OF IENTRY. */

/*<       GOTO (10,20), IENTRY >*/
    switch (*ientry) {
	case 1:  goto L10;
	case 2:  goto L20;
    }

/*      IENTRY=1 */
/*      CHECK INPUT PARAMETERS */

/*< 10      ITEST = 2 >*/
L10:
    *itest = 2;
/*<    >*/
    if (*u <= zero || *xbnd <= *tnytol || *gu > zero) {
	return 0;
    }
/*<       ITEST = 1 >*/
    *itest = 1;
/*<       IF (XBND .LT. ABSTOL) ABSTOL = XBND >*/
    if (*xbnd < *abstol) {
	*abstol = *xbnd;
    }
/*<       TOL = ABSTOL >*/
    *tol = *abstol;
/*<       TWOTOL = TOL + TOL >*/
    twotol = *tol + *tol;

/* A AND B DEFINE THE INTERVAL OF UNCERTAINTY, X AND XW ARE POINTS */
/* WITH LOWEST AND SECOND LOWEST FUNCTION VALUES SO FAR OBTAINED. */
/* INITIALIZE A,SMIN,XW AT ORIGIN AND CORRESPONDING VALUES OF */
/* FUNCTION AND PROJECTION OF THE GRADIENT ALONG DIRECTION OF SEARCH */
/* AT VALUES FOR LATEST ESTIMATE AT MINIMUM. */

/*<       A = ZERO >*/
    *a = zero;
/*<       XW = ZERO >*/
    *xw = zero;
/*<       XMIN = ZERO >*/
    *xmin = zero;
/*<       OLDF = FU >*/
    *oldf = *fu;
/*<       FMIN = FU >*/
    *fmin = *fu;
/*<       FW = FU >*/
    *fw = *fu;
/*<       GW = GU >*/
    *gw = *gu;
/*<       GMIN = GU >*/
    *gmin = *gu;
/*<       STEP = U >*/
    *step = *u;
/*<       FACTOR = FIVE >*/
    *factor = five;

/*      THE MINIMUM HAS NOT YET BEEN BRACKETED. */

/*<       BRAKTD = .FALSE. >*/
    *braktd = FALSE_;

/* SET UP XBND AS A BOUND ON THE STEP TO BE TAKEN. (XBND IS NOT COMPUTED */
/* EXPLICITLY BUT SCXBND IS ITS SCALED VALUE.)  SET THE UPPER BOUND */
/* ON THE INTERVAL OF UNCERTAINTY INITIALLY TO XBND + TOL(XBND). */

/*<       SCXBND = XBND >*/
    *scxbnd = *xbnd;
/*<       B = SCXBND + RELTOL*DABS(SCXBND) + ABSTOL >*/
    *b = *scxbnd + *reltol * abs(*scxbnd) + *abstol;
/*<       E = B + B >*/
    *e = *b + *b;
/*<       B1 = B >*/
    *b1 = *b;

/* COMPUTE THE CONSTANTS REQUIRED FOR THE TWO CONVERGENCE CRITERIA. */

/*<       GTEST1 = -RMU*GU >*/
    *gtest1 = -(*rmu) * *gu;
/*<       GTEST2 = -ETA*GU >*/
    *gtest2 = -(*eta) * *gu;

/* SET IENTRY TO INDICATE THAT THIS IS THE FIRST ITERATION */

/*<       IENTRY = 2 >*/
    *ientry = 2;
/*<       GO TO 210 >*/
    goto L210;

/* IENTRY = 2 */

/* UPDATE A,B,XW, AND XMIN */

/*< 20      IF (FU .GT. FMIN) GO TO 60 >*/
L20:
    if (*fu > *fmin) {
	goto L60;
    }

/* IF FUNCTION VALUE NOT INCREASED, NEW POINT BECOMES NEXT */
/* ORIGIN AND OTHER POINTS ARE SCALED ACCORDINGLY. */

/*<       CHORDU = OLDF - (XMIN + U)*GTEST1 >*/
    chordu = *oldf - (*xmin + *u) * *gtest1;
/*<       IF (FU .LE. CHORDU) GO TO 30 >*/
    if (*fu <= chordu) {
	goto L30;
    }

/* THE NEW FUNCTION VALUE DOES NOT SATISFY THE SUFFICIENT DECREASE */
/* CRITERION. PREPARE TO MOVE THE UPPER BOUND TO THIS POINT AND */
/* FORCE THE INTERPOLATION SCHEME TO EITHER BISECT THE INTERVAL OF */
/* UNCERTAINTY OR TAKE THE LINEAR INTERPOLATION STEP WHICH ESTIMATES */
/* THE ROOT OF F(ALPHA)=CHORD(ALPHA). */

/*<       CHORDM = OLDF - XMIN*GTEST1 >*/
    chordm = *oldf - *xmin * *gtest1;
/*<       GU = -GMIN >*/
    *gu = -(*gmin);
/*<       DENOM = CHORDM-FMIN >*/
    denom = chordm - *fmin;
/*<       IF (DABS(DENOM) .GE. 1.D-15) GO TO 25 >*/
    if (abs(denom) >= 1e-15) {
	goto L25;
    }
/*<           DENOM = 1.D-15 >*/
    denom = 1e-15;
/*<           IF (CHORDM-FMIN .LT. 0.D0)  DENOM = -DENOM >*/
    if (chordm - *fmin < 0.) {
	denom = -denom;
    }
/*< 25    CONTINUE >*/
L25:
/*<       IF (XMIN .NE. ZERO) GU = GMIN*(CHORDU-FU)/DENOM >*/
    if (*xmin != zero) {
	*gu = *gmin * (chordu - *fu) / denom;
    }
/*<       FU = HALF*U*(GMIN+GU) + FMIN >*/
    *fu = half * *u * (*gmin + *gu) + *fmin;
/*<       IF (FU .LT. FMIN) FU = FMIN >*/
    if (*fu < *fmin) {
	*fu = *fmin;
    }
/*<       GO TO 60 >*/
    goto L60;
/*< 30      FW = FMIN >*/
L30:
    *fw = *fmin;
/*<       FMIN = FU >*/
    *fmin = *fu;
/*<       GW = GMIN >*/
    *gw = *gmin;
/*<       GMIN = GU >*/
    *gmin = *gu;
/*<       XMIN = XMIN + U >*/
    *xmin += *u;
/*<       A = A-U >*/
    *a -= *u;
/*<       B = B-U >*/
    *b -= *u;
/*<       XW = -U >*/
    *xw = -(*u);
/*<       SCXBND = SCXBND - U >*/
    *scxbnd -= *u;
/*<       IF (GU .LE. ZERO) GO TO 40 >*/
    if (*gu <= zero) {
	goto L40;
    }
/*<       B = ZERO >*/
    *b = zero;
/*<       BRAKTD = .TRUE. >*/
    *braktd = TRUE_;
/*<       GO TO 50 >*/
    goto L50;
/*< 40    A = ZERO >*/
L40:
    *a = zero;
/*< 50    TOL = DABS(XMIN)*RELTOL + ABSTOL >*/
L50:
    *tol = abs(*xmin) * *reltol + *abstol;
/*<       GO TO 90 >*/
    goto L90;

/* IF FUNCTION VALUE INCREASED, ORIGIN REMAINS UNCHANGED */
/* BUT NEW POINT MAY NOW QUALIFY AS W. */

/*< 60    IF (U .LT. ZERO) GO TO 70 >*/
L60:
    if (*u < zero) {
	goto L70;
    }
/*<       B = U >*/
    *b = *u;
/*<       BRAKTD = .TRUE. >*/
    *braktd = TRUE_;
/*<       GO TO 80 >*/
    goto L80;
/*< 70    A = U >*/
L70:
    *a = *u;
/*< 80    XW = U >*/
L80:
    *xw = *u;
/*<       FW = FU >*/
    *fw = *fu;
/*<       GW = GU >*/
    *gw = *gu;
/*< 90    TWOTOL = TOL + TOL >*/
L90:
    twotol = *tol + *tol;
/*<       XMIDPT = HALF*(A + B) >*/
    xmidpt = half * (*a + *b);

/* CHECK TERMINATION CRITERIA */

/*<    >*/
    convrg = abs(xmidpt) <= twotol - half * (*b - *a) || abs(*gmin) <= *
	    gtest2 && *fmin < *oldf && ((d__1 = *xmin - *xbnd, abs(d__1)) > *
	    tol || ! (*braktd));
/*<       IF (.NOT. CONVRG) GO TO 100 >*/
    if (! convrg) {
	goto L100;
    }
/*<       ITEST = 0 >*/
    *itest = 0;
/*<       IF (XMIN .NE. ZERO) RETURN >*/
    if (*xmin != zero) {
	return 0;
    }

/* IF THE FUNCTION HAS NOT BEEN REDUCED, CHECK TO SEE THAT THE RELATIVE */
/* CHANGE IN F(X) IS CONSISTENT WITH THE ESTIMATE OF THE DELTA- */
/* UNIMODALITY CONSTANT, TOL.  IF THE CHANGE IN F(X) IS LARGER THAN */
/* EXPECTED, REDUCE THE VALUE OF TOL. */

/*<       ITEST = 3 >*/
    *itest = 3;
/*<       IF (DABS(OLDF-FW) .LE. FPRESN*(ONE + DABS(OLDF))) RETURN >*/
    if ((d__1 = *oldf - *fw, abs(d__1)) <= *fpresn * (one + abs(*oldf))) {
	return 0;
    }
/*<       TOL = POINT1*TOL >*/
    *tol = point1 * *tol;
/*<       IF (TOL .LT. TNYTOL) RETURN >*/
    if (*tol < *tnytol) {
	return 0;
    }
/*<       RELTOL = POINT1*RELTOL >*/
    *reltol = point1 * *reltol;
/*<       ABSTOL = POINT1*ABSTOL >*/
    *abstol = point1 * *abstol;
/*<       TWOTOL = POINT1*TWOTOL >*/
    twotol = point1 * twotol;

/* CONTINUE WITH THE COMPUTATION OF A TRIAL STEP LENGTH */

/*< 100   R = ZERO >*/
L100:
    r__ = zero;
/*<       Q = ZERO >*/
    q = zero;
/*<       S = ZERO >*/
    s = zero;
/*<       IF (DABS(E) .LE. TOL) GO TO 150 >*/
    if (abs(*e) <= *tol) {
	goto L150;
    }

/* FIT CUBIC THROUGH XMIN AND XW */

/*<       R = THREE*(FMIN-FW)/XW + GMIN + GW >*/
    r__ = three * (*fmin - *fw) / *xw + *gmin + *gw;
/*<       ABSR = DABS(R) >*/
    absr = abs(r__);
/*<       Q = ABSR >*/
    q = absr;
/*<       IF (GW .EQ. ZERO .OR. GMIN .EQ. ZERO) GO TO 140 >*/
    if (*gw == zero || *gmin == zero) {
	goto L140;
    }

/* COMPUTE THE SQUARE ROOT OF (R*R - GMIN*GW) IN A WAY */
/* WHICH AVOIDS UNDERFLOW AND OVERFLOW. */

/*<       ABGW = DABS(GW) >*/
    abgw = abs(*gw);
/*<       ABGMIN = DABS(GMIN) >*/
    abgmin = abs(*gmin);
/*<       S = DSQRT(ABGMIN)*DSQRT(ABGW) >*/
    s = sqrt(abgmin) * sqrt(abgw);
/*<       IF ((GW/ABGW)*GMIN .GT. ZERO) GO TO 130 >*/
    if (*gw / abgw * *gmin > zero) {
	goto L130;
    }

/* COMPUTE THE SQUARE ROOT OF R*R + S*S. */

/*<       SUMSQ = ONE >*/
    sumsq = one;
/*<       P = ZERO >*/
    p = zero;
/*<       IF (ABSR .GE. S) GO TO 110 >*/
    if (absr >= s) {
	goto L110;
    }

/* THERE IS A POSSIBILITY OF OVERFLOW. */

/*<       IF (S .GT. RTSMLL) P = S*RTSMLL >*/
    if (s > *rtsmll) {
	p = s * *rtsmll;
    }
/*<       IF (ABSR .GE. P) SUMSQ = ONE +(ABSR/S)**2 >*/
    if (absr >= p) {
/* Computing 2nd power */
	d__1 = absr / s;
	sumsq = one + d__1 * d__1;
    }
/*<       SCALE = S >*/
    scale = s;
/*<       GO TO 120 >*/
    goto L120;

/* THERE IS A POSSIBILITY OF UNDERFLOW. */

/*< 110   IF (ABSR .GT. RTSMLL) P = ABSR*RTSMLL >*/
L110:
    if (absr > *rtsmll) {
	p = absr * *rtsmll;
    }
/*<       IF (S .GE. P) SUMSQ = ONE + (S/ABSR)**2 >*/
    if (s >= p) {
/* Computing 2nd power */
	d__1 = s / absr;
	sumsq = one + d__1 * d__1;
    }
/*<       SCALE = ABSR >*/
    scale = absr;
/*< 120   SUMSQ = DSQRT(SUMSQ) >*/
L120:
    sumsq = sqrt(sumsq);
/*<       Q = BIG >*/
    q = *big;
/*<       IF (SCALE .LT. BIG/SUMSQ) Q = SCALE*SUMSQ >*/
    if (scale < *big / sumsq) {
	q = scale * sumsq;
    }
/*<       GO TO 140 >*/
    goto L140;

/* COMPUTE THE SQUARE ROOT OF R*R - S*S */

/*< 130   Q = DSQRT(DABS(R+S))*DSQRT(DABS(R-S)) >*/
L130:
    q = sqrt((d__1 = r__ + s, abs(d__1))) * sqrt((d__2 = r__ - s, abs(d__2)));
/*<       IF (R .GE. S .OR. R .LE. (-S)) GO TO 140 >*/
    if (r__ >= s || r__ <= -s) {
	goto L140;
    }
/*<       R = ZERO >*/
    r__ = zero;
/*<       Q = ZERO >*/
    q = zero;
/*<       GO TO 150 >*/
    goto L150;

/* COMPUTE THE MINIMUM OF FITTED CUBIC */

/*< 140   IF (XW .LT. ZERO) Q = -Q >*/
L140:
    if (*xw < zero) {
	q = -q;
    }
/*<       S = XW*(GMIN - R - Q) >*/
    s = *xw * (*gmin - r__ - q);
/*<       Q = GW - GMIN + Q + Q >*/
    q = *gw - *gmin + q + q;
/*<       IF (Q .GT. ZERO) S = -S >*/
    if (q > zero) {
	s = -s;
    }
/*<       IF (Q .LE. ZERO) Q = -Q >*/
    if (q <= zero) {
	q = -q;
    }
/*<       R = E >*/
    r__ = *e;
/*<       IF (B1 .NE. STEP .OR. BRAKTD) E = STEP >*/
    if (*b1 != *step || *braktd) {
	*e = *step;
    }

/* CONSTRUCT AN ARTIFICIAL BOUND ON THE ESTIMATED STEPLENGTH */

/*< 150   A1 = A >*/
L150:
    a1 = *a;
/*<       B1 = B >*/
    *b1 = *b;
/*<       STEP = XMIDPT >*/
    *step = xmidpt;
/*<       IF (BRAKTD) GO TO 160 >*/
    if (*braktd) {
	goto L160;
    }
/*<       STEP = -FACTOR*XW >*/
    *step = -(*factor) * *xw;
/*<       IF (STEP .GT. SCXBND) STEP = SCXBND >*/
    if (*step > *scxbnd) {
	*step = *scxbnd;
    }
/*<       IF (STEP .NE. SCXBND) FACTOR = FIVE*FACTOR >*/
    if (*step != *scxbnd) {
	*factor = five * *factor;
    }
/*<       GO TO 170 >*/
    goto L170;

/* IF THE MINIMUM IS BRACKETED BY 0 AND XW THE STEP MUST LIE */
/* WITHIN (A,B). */

/*< 16 >*/
L160:
    if ((*a != zero || *xw >= zero) && (*b != zero || *xw <= zero)) {
	goto L180;
    }

/* IF THE MINIMUM IS NOT BRACKETED BY 0 AND XW THE STEP MUST LIE */
/* WITHIN (A1,B1). */

/*<       D1 = XW >*/
    d1 = *xw;
/*<       D2 = A >*/
    d2 = *a;
/*<       IF (A .EQ. ZERO) D2 = B >*/
    if (*a == zero) {
	d2 = *b;
    }
/* THIS LINE MIGHT BE */
/*     IF (A .EQ. ZERO) D2 = E */
/*<       U = - D1/D2 >*/
    *u = -d1 / d2;
/*<       STEP = FIVE*D2*(POINT1 + ONE/U)/ELEVEN >*/
    *step = five * d2 * (point1 + one / *u) / eleven;
/*<       IF (U .LT. ONE) STEP = HALF*D2*DSQRT(U) >*/
    if (*u < one) {
	*step = half * d2 * sqrt(*u);
    }
/*< 170   IF (STEP .LE. ZERO) A1 = STEP >*/
L170:
    if (*step <= zero) {
	a1 = *step;
    }
/*<       IF (STEP .GT. ZERO) B1 = STEP >*/
    if (*step > zero) {
	*b1 = *step;
    }

/* REJECT THE STEP OBTAINED BY INTERPOLATION IF IT LIES OUTSIDE THE */
/* REQUIRED INTERVAL OR IT IS GREATER THAN HALF THE STEP OBTAINED */
/* DURING THE LAST-BUT-ONE ITERATION. */

/*< 18 >*/
L180:
    if (abs(s) <= (d__1 = half * q * r__, abs(d__1)) || s <= q * a1 || s >= q
	    * *b1) {
	goto L200;
    }

/* A CUBIC INTERPOLATION STEP */

/*<       STEP = S/Q >*/
    *step = s / q;

/* THE FUNCTION MUST NOT BE EVALUTATED TOO CLOSE TO A OR B. */

/*<       IF (STEP - A .GE. TWOTOL .AND. B - STEP .GE. TWOTOL) GO TO 210 >*/
    if (*step - *a >= twotol && *b - *step >= twotol) {
	goto L210;
    }
/*<       IF (XMIDPT .GT. ZERO) GO TO 190 >*/
    if (xmidpt > zero) {
	goto L190;
    }
/*<       STEP = -TOL >*/
    *step = -(*tol);
/*<       GO TO 210 >*/
    goto L210;
/*< 190   STEP = TOL >*/
L190:
    *step = *tol;
/*<       GO TO 210 >*/
    goto L210;
/*< 200   E = B-A >*/
L200:
    *e = *b - *a;

/* IF THE STEP IS TOO LARGE, REPLACE BY THE SCALED BOUND (SO AS TO */
/* COMPUTE THE NEW POINT ON THE BOUNDARY). */

/*< 210   IF (STEP .LT. SCXBND) GO TO 220 >*/
L210:
    if (*step < *scxbnd) {
	goto L220;
    }
/*<       STEP = SCXBND >*/
    *step = *scxbnd;

/* MOVE SXBD TO THE LEFT SO THAT SBND + TOL(XBND) = XBND. */

/*<       SCXBND = SCXBND - (RELTOL*DABS(XBND)+ABSTOL)/(ONE + RELTOL) >*/
    *scxbnd -= (*reltol * abs(*xbnd) + *abstol) / (one + *reltol);
/*< 220   U = STEP >*/
L220:
    *u = *step;
/*<       IF (DABS(STEP) .LT. TOL .AND. STEP .LT. ZERO) U = -TOL >*/
    if (abs(*step) < *tol && *step < zero) {
	*u = -(*tol);
    }
/*<       IF (DABS(STEP) .LT. TOL .AND. STEP .GE. ZERO) U = TOL >*/
    if (abs(*step) < *tol && *step >= zero) {
	*u = *tol;
    }
/*<       ITEST = 1 >*/
    *itest = 1;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* getptc_ */

/* -------------------------------- BLAS ------------------------------------ */
/* %% TRUNCATED-NEWTON METHOD: BLAS */
/*   NOTE: ALL ROUTINES HERE ARE FROM LINPACK WITH THE EXCEPTION */
/*         OF DXPY (A VERSION OF DAXPY WITH A=1.0) */
/*   WRITTEN BY:  STEPHEN G. NASH */
/*                OPERATIONS RESEARCH AND APPLIED STATISTICS DEPT. */
/*                GEORGE MASON UNIVERSITY */
/*                FAIRFAX, VA 22030 */
/* **************************************************************** */
/*<       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY) >*/
doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy,
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    integer i__, m;
    doublereal dtemp;
    integer ix, iy, mp1;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/

/*     FORMS THE DOT PRODUCT OF TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */

/*<       DOUBLE PRECISION DX(1),DY(1),DTEMP >*/
/*<       INTEGER I,INCX,INCY,IX,IY,M,MP1,N >*/

/*<       DDOT = 0.0D0 >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
/*<       DTEMP = 0.0D0 >*/
    dtemp = 0.;
/*<       IF(N.LE.0)RETURN >*/
    if (*n <= 0) {
	return ret_val;
    }
/*<       IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20 >*/
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

/*<       IX = 1 >*/
    ix = 1;
/*<       IY = 1 >*/
    iy = 1;
/*<       IF(INCX.LT.0)IX = (-N+1)*INCX + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       IF(INCY.LT.0)IY = (-N+1)*INCY + 1 >*/
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DTEMP = DTEMP + DX(IX)*DY(IY) >*/
	dtemp += dx[ix] * dy[iy];
/*<         IX = IX + INCX >*/
	ix += *incx;
/*<         IY = IY + INCY >*/
	iy += *incy;
/*<    10 CONTINUE >*/
/* L10: */
    }
/*<       DDOT = DTEMP >*/
    ret_val = dtemp;
/*<       RETURN >*/
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

/*<    20 M = MOD(N,5) >*/
L20:
    m = *n % 5;
/*<       IF( M .EQ. 0 ) GO TO 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       DO 30 I = 1,M >*/
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DTEMP = DTEMP + DX(I)*DY(I) >*/
	dtemp += dx[i__] * dy[i__];
/*<    30 CONTINUE >*/
/* L30: */
    }
/*<       IF( N .LT. 5 ) GO TO 60 >*/
    if (*n < 5) {
	goto L60;
    }
/*<    40 MP1 = M + 1 >*/
L40:
    mp1 = m + 1;
/*<       DO 50 I = MP1,N,5 >*/
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
/*<    >*/
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ +
		4] * dy[i__ + 4];
/*<    50 CONTINUE >*/
/* L50: */
    }
/*<    60 DDOT = DTEMP >*/
L60:
    ret_val = dtemp;
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* ddot_ */

/*<       SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY) >*/
/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx,
	integer *incx, doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, m, ix, iy, mp1;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/

/*     CONSTANT TIMES A VECTOR PLUS A VECTOR. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */

/*<       DOUBLE PRECISION DX(1),DY(1),DA >*/
/*<       INTEGER I,INCX,INCY,IX,IY,M,MP1,N >*/

/*<       IF(N.LE.0)RETURN >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
/*<       IF (DA .EQ. 0.0D0) RETURN >*/
    if (*da == 0.) {
	return 0;
    }
/*<       IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20 >*/
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

/*<       IX = 1 >*/
    ix = 1;
/*<       IY = 1 >*/
    iy = 1;
/*<       IF(INCX.LT.0)IX = (-N+1)*INCX + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       IF(INCY.LT.0)IY = (-N+1)*INCY + 1 >*/
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DY(IY) = DY(IY) + DA*DX(IX) >*/
	dy[iy] += *da * dx[ix];
/*<         IX = IX + INCX >*/
	ix += *incx;
/*<         IY = IY + INCY >*/
	iy += *incy;
/*<    10 CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

/*<    20 M = MOD(N,4) >*/
L20:
    m = *n % 4;
/*<       IF( M .EQ. 0 ) GO TO 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       DO 30 I = 1,M >*/
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DY(I) = DY(I) + DA*DX(I) >*/
	dy[i__] += *da * dx[i__];
/*<    30 CONTINUE >*/
/* L30: */
    }
/*<       IF( N .LT. 4 ) RETURN >*/
    if (*n < 4) {
	return 0;
    }
/*<    40 MP1 = M + 1 >*/
L40:
    mp1 = m + 1;
/*<       DO 50 I = MP1,N,4 >*/
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
/*<         DY(I) = DY(I) + DA*DX(I) >*/
	dy[i__] += *da * dx[i__];
/*<         DY(I + 1) = DY(I + 1) + DA*DX(I + 1) >*/
	dy[i__ + 1] += *da * dx[i__ + 1];
/*<         DY(I + 2) = DY(I + 2) + DA*DX(I + 2) >*/
	dy[i__ + 2] += *da * dx[i__ + 2];
/*<         DY(I + 3) = DY(I + 3) + DA*DX(I + 3) >*/
	dy[i__ + 3] += *da * dx[i__ + 3];
/*<    50 CONTINUE >*/
/* L50: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* daxpy_ */

/*<       DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX) >*/
doublereal dnrm2_(integer *n, doublereal *dx, integer *incx)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal cutlo = 8.232e-11;
    static doublereal cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal xmax;
    integer next, i__, j, nn;
    doublereal hitest, sum;

    /* Assigned format variables */
    static char *next_fmt;

/*<       IMPLICIT         DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER          NEXT >*/
/*<       DOUBLE PRECISION DX(1),CUTLO,CUTHI,HITEST,SUM,XMAX,ZERO,ONE >*/
/*<       DATA   ZERO, ONE /0.0D0, 1.0D0/ >*/
    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */

/*           C.L.LAWSON, 1978 JAN 08 */

/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */

/*     BRIEF OUTLINE OF ALGORITHM.. */

/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX. */

/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE */
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */
/*<       DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / >*/

/*<       IF(N .GT. 0) GO TO 10 >*/
    if (*n > 0) {
	goto L10;
    }
/*<          DNRM2  = ZERO >*/
    ret_val = zero;
/*<          GO TO 300 >*/
    goto L300;

/*<    10 ASSIGN 30 TO NEXT >*/
L10:
    next = 0;
    next_fmt = fmt_30;
/*<       SUM = ZERO >*/
    sum = zero;
/*<       NN = N * INCX >*/
    nn = *n * *incx;
/*                                                 BEGIN MAIN LOOP */
/*<       I = 1 >*/
    i__ = 1;
/*<    20    GO TO NEXT,(30, 50, 70, 110) >*/
L20:
    switch (next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
/*<    30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85 >*/
L30:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }
/*<       ASSIGN 50 TO NEXT >*/
    next = 1;
    next_fmt = fmt_50;
/*<       XMAX = ZERO >*/
    xmax = zero;

/*                        PHASE 1.  SUM IS ZERO */

/*<    50 IF( DX(I) .EQ. ZERO) GO TO 200 >*/
L50:
    if (dx[i__] == zero) {
	goto L200;
    }
/*<       IF( DABS(DX(I)) .GT. CUTLO) GO TO 85 >*/
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }

/*                                PREPARE FOR PHASE 2. */
/*<       ASSIGN 70 TO NEXT >*/
    next = 2;
    next_fmt = fmt_70;
/*<       GO TO 105 >*/
    goto L105;

/*                                PREPARE FOR PHASE 4. */

/*<   100 I = J >*/
L100:
    i__ = j;
/*<       ASSIGN 110 TO NEXT >*/
    next = 3;
    next_fmt = fmt_110;
/*<       SUM = (SUM / DX(I)) / DX(I) >*/
    sum = sum / dx[i__] / dx[i__];
/*<   105 XMAX = DABS(DX(I)) >*/
L105:
    xmax = (d__1 = dx[i__], abs(d__1));
/*<       GO TO 115 >*/
    goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

/*<    70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75 >*/
L70:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L75;
    }

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */

/*<   110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115 >*/
L110:
    if ((d__1 = dx[i__], abs(d__1)) <= xmax) {
	goto L115;
    }
/*<          SUM = ONE + SUM * (XMAX / DX(I))**2 >*/
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
/*<          XMAX = DABS(DX(I)) >*/
    xmax = (d__1 = dx[i__], abs(d__1));
/*<          GO TO 200 >*/
    goto L200;

/*<   115 SUM = SUM + (DX(I)/XMAX)**2 >*/
L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
/*<       GO TO 200 >*/
    goto L200;


/*                  PREPARE FOR PHASE 3. */

/*<    75 SUM = (SUM * XMAX) * XMAX >*/
L75:
    sum = sum * xmax * xmax;


/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

/*<    85 HITEST = CUTHI/FLOAT( N ) >*/
L85:
    hitest = cuthi / (real) (*n);

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

/*<       DO 95 J =I,NN,INCX >*/
    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/*<       IF(DABS(DX(J)) .GE. HITEST) GO TO 100 >*/
	if ((d__1 = dx[j], abs(d__1)) >= hitest) {
	    goto L100;
	}
/*<    95    SUM = SUM + DX(J)**2 >*/
/* L95: */
/* Computing 2nd power */
	d__1 = dx[j];
	sum += d__1 * d__1;
    }
/*<       DNRM2 = DSQRT( SUM ) >*/
    ret_val = sqrt(sum);
/*<       GO TO 300 >*/
    goto L300;

/*<   200 CONTINUE >*/
L200:
/*<       I = I + INCX >*/
    i__ += *incx;
/*<       IF ( I .LE. NN ) GO TO 20 >*/
    if (i__ <= nn) {
	goto L20;
    }

/*              END OF MAIN LOOP. */

/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

/*<       DNRM2 = XMAX * DSQRT(SUM) >*/
    ret_val = xmax * sqrt(sum);
/*<   300 CONTINUE >*/
L300:
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* dnrm2_ */

/*<       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY) >*/
/* Subroutine */ int dcopy_(integer *n, doublereal *dx, integer *incx,
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, m, ix, iy, mp1;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/

/*     COPIES A VECTOR, X, TO A VECTOR, Y. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */

/*<       DOUBLE PRECISION DX(1),DY(1) >*/
/*<       INTEGER I,INCX,INCY,IX,IY,M,MP1,N >*/

/*<       IF(N.LE.0)RETURN >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
/*<       IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20 >*/
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

/*<       IX = 1 >*/
    ix = 1;
/*<       IY = 1 >*/
    iy = 1;
/*<       IF(INCX.LT.0)IX = (-N+1)*INCX + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       IF(INCY.LT.0)IY = (-N+1)*INCY + 1 >*/
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DY(IY) = DX(IX) >*/
	dy[iy] = dx[ix];
/*<         IX = IX + INCX >*/
	ix += *incx;
/*<         IY = IY + INCY >*/
	iy += *incy;
/*<    10 CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

/*<    20 M = MOD(N,7) >*/
L20:
    m = *n % 7;
/*<       IF( M .EQ. 0 ) GO TO 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       DO 30 I = 1,M >*/
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DY(I) = DX(I) >*/
	dy[i__] = dx[i__];
/*<    30 CONTINUE >*/
/* L30: */
    }
/*<       IF( N .LT. 7 ) RETURN >*/
    if (*n < 7) {
	return 0;
    }
/*<    40 MP1 = M + 1 >*/
L40:
    mp1 = m + 1;
/*<       DO 50 I = MP1,N,7 >*/
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
/*<         DY(I) = DX(I) >*/
	dy[i__] = dx[i__];
/*<         DY(I + 1) = DX(I + 1) >*/
	dy[i__ + 1] = dx[i__ + 1];
/*<         DY(I + 2) = DX(I + 2) >*/
	dy[i__ + 2] = dx[i__ + 2];
/*<         DY(I + 3) = DX(I + 3) >*/
	dy[i__ + 3] = dx[i__ + 3];
/*<         DY(I + 4) = DX(I + 4) >*/
	dy[i__ + 4] = dx[i__ + 4];
/*<         DY(I + 5) = DX(I + 5) >*/
	dy[i__ + 5] = dx[i__ + 5];
/*<         DY(I + 6) = DX(I + 6) >*/
	dy[i__ + 6] = dx[i__ + 6];
/*<    50 CONTINUE >*/
/* L50: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* dcopy_ */

/* ****************************************************************** */
/* SPECIAL BLAS FOR Y = X+Y */
/* ****************************************************************** */
/*<       SUBROUTINE DXPY(N,DX,INCX,DY,INCY) >*/
/* Subroutine */ int dxpy_(integer *n, doublereal *dx, integer *incx,
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, m, ix, iy, mp1;

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/

/*     VECTOR PLUS A VECTOR. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     STEPHEN G. NASH 5/30/89. */

/*<       DOUBLE PRECISION DX(1),DY(1) >*/
/*<       INTEGER I,INCX,INCY,IX,IY,M,MP1,N >*/

/*<       IF(N.LE.0)RETURN >*/
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
/*<       IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20 >*/
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

/*<       IX = 1 >*/
    ix = 1;
/*<       IY = 1 >*/
    iy = 1;
/*<       IF(INCX.LT.0)IX = (-N+1)*INCX + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       IF(INCY.LT.0)IY = (-N+1)*INCY + 1 >*/
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
/*<       DO 10 I = 1,N >*/
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DY(IY) = DY(IY) + DX(IX) >*/
	dy[iy] += dx[ix];
/*<         IX = IX + INCX >*/
	ix += *incx;
/*<         IY = IY + INCY >*/
	iy += *incy;
/*<    10 CONTINUE >*/
/* L10: */
    }
/*<       RETURN >*/
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

/*<    20 M = MOD(N,4) >*/
L20:
    m = *n % 4;
/*<       IF( M .EQ. 0 ) GO TO 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       DO 30 I = 1,M >*/
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<         DY(I) = DY(I) + DX(I) >*/
	dy[i__] += dx[i__];
/*<    30 CONTINUE >*/
/* L30: */
    }
/*<       IF( N .LT. 4 ) RETURN >*/
    if (*n < 4) {
	return 0;
    }
/*<    40 MP1 = M + 1 >*/
L40:
    mp1 = m + 1;
/*<       DO 50 I = MP1,N,4 >*/
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
/*<         DY(I) = DY(I) + DX(I) >*/
	dy[i__] += dx[i__];
/*<         DY(I + 1) = DY(I + 1) + DX(I + 1) >*/
	dy[i__ + 1] += dx[i__ + 1];
/*<         DY(I + 2) = DY(I + 2) + DX(I + 2) >*/
	dy[i__ + 2] += dx[i__ + 2];
/*<         DY(I + 3) = DY(I + 3) + DX(I + 3) >*/
	dy[i__ + 3] += dx[i__ + 3];
/*<    50 CONTINUE >*/
/* L50: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* dxpy_ */

