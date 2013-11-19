/* dnls1.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

//#include "f2c.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define real    double
#define integer int
typedef long int logical;
typedef double doublereal;
typedef char *address;

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

/* ----------------- Start buildin functions -------------------*/

int s_stop(char *s, ftnlen n) {return 0;}

#define log10e 0.43429448190325182765

double d_lg10(real *x)
{
return( log10e * log(*x) );
}

integer i_len(char *s, ftnlen n)
{
return(n);
}

int s_copy(char *a, char *b, ftnlen la, ftnlen lb)
{
	register char *aend, *bend;

	aend = a + la;

	if(la <= lb)
#ifndef NO_OVERWRITE
		if (a <= b || a >= b + la)
#endif
			while(a < aend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else
			for(b += la; a < aend; )
				*--aend = *--b;
#endif

	else {
		bend = b + lb;
#ifndef NO_OVERWRITE
		if (a <= b || a >= bend)
#endif
			while(b < bend)
				*a++ = *b++;
#ifndef NO_OVERWRITE
		else {
			a += lb;
			while(b < bend)
				*--a = *--bend;
			a += lb;
			}
#endif
		while(a < aend)
			*a++ = ' ';
		}
return 0;
	}



integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
register unsigned char *a, *aend, *b, *bend;
a = (unsigned char *)a0;
b = (unsigned char *)b0;
aend = a + la;
bend = b + lb;

if(la <= lb)
	{
	while(a < aend)
		if(*a != *b)
			return( *a - *b );
		else
			{ ++a; ++b; }

	while(b < bend)
		if(*b != ' ')
			return( ' ' - *b );
		else	++b;
	}

else
	{
	while(b < bend)
		if(*a == *b)
			{ ++a; ++b; }
		else
			return( *a - *b );
	while(a < aend)
		if(*a != ' ')
			return(*a - ' ');
		else	++a;
	}
return(0);
}


#ifdef KR_headers
double pow_ri(ap, bp) real *ap; integer *bp;
#else
double pow_ri(real *ap, integer *bp)
#endif
{
double pow, x;
integer n;
unsigned long u;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for(u = n; ; )
		{
		if(u & 01)
			pow *= x;
		if(u >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}


#ifdef KR_headers
void s_cat(lp, rpp, rnp, np, ll) char *lp, *rpp[]; ftnlen rnp[], *np, ll;
#else
void s_cat(char *lp, char *rpp[], ftnlen rnp[], ftnlen *np, ftnlen ll)
#endif
{
ftnlen i, n, nc;
char *f__rp;

n = (int)*np;
for(i = 0 ; i < n ; ++i)
	{
	nc = ll;
	if(rnp[i] < nc)
		nc = rnp[i];
	ll -= nc;
	f__rp = rpp[i];
	while(--nc >= 0)
		*lp++ = *f__rp++;
	}
while(--ll >= 0)
	*lp++ = ' ';
}


#ifdef KR_headers
integer i_indx(a, b, la, lb) char *a, *b; ftnlen la, lb;
#else
integer i_indx(char *a, char *b, ftnlen la, ftnlen lb)
#endif
{
ftnlen i, n;
char *s, *t, *bend;

n = la - lb + 1;
bend = b + lb;

for(i = 0 ; i < n ; ++i)
	{
	s = a + i;
	t = b;
	while(t < bend)
		if(*s++ != *t++)
			goto no;
	return(i+1);
	no: ;
	}
return(0);
}




/* ----------------- end buildin functions -------------------*/



/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__9 = 9;
static integer c__3 = 3;
static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__72 = 72;
static integer c__5 = 5;


#ifdef OLD_MACH

/*<       DOUBLE PRECISION FUNCTION D1MACH (I) >*/
doublereal d1mach_(integer *i)
{
    /* System generated locals */
    doublereal ret_val;
    static doublereal equiv_4[6];

    /* Local variables */
#define log10 ((integer *)equiv_4 + 8)
#define dmach (equiv_4)
#define large ((integer *)equiv_4 + 2)
#define small ((integer *)equiv_4)
#define diver ((integer *)equiv_4 + 6)
#define right ((integer *)equiv_4 + 4)
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D1MACH */
/* ***PURPOSE  Return floating point machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   D1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument, and can be referenced as follows: */

/*        D = D1MACH(I) */

/*   where I=1,...,5.  The (output) value of D above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   D1MACH( 3) = B**(-T), the smallest relative spacing. */
/*   D1MACH( 4) = B**(1-T), the largest relative spacing. */
/*   D1MACH( 5) = LOG10(B) */

/*   Assume double precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(14) = T, the number of base-B digits. */
/*   I1MACH(15) = EMIN, the smallest exponent E. */
/*   I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for 
*/
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890213  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   900911  Added SUN 386i constants.  (WRB) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added CONVEX -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/* ***END PROLOGUE  D1MACH */

/*<       INTEGER SMALL(4) >*/
/*<       INTEGER LARGE(4) >*/
/*<       INTEGER RIGHT(4) >*/
/*<       INTEGER DIVER(4) >*/
/*<       INTEGER LOG10(4) >*/

/*<       DOUBLE PRECISION DMACH(5) >*/
/*<       SAVE DMACH >*/

/*<       EQUIVALENCE (DMACH(1),SMALL(1)) >*/
/*<       EQUIVALENCE (DMACH(2),LARGE(1)) >*/
/*<       EQUIVALENCE (DMACH(3),RIGHT(1)) >*/
/*<       EQUIVALENCE (DMACH(4),DIVER(1)) >*/
/*<       EQUIVALENCE (DMACH(5),LOG10(1)) >*/

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 / */
/*     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 / */
/*     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 / */
/*     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA SMALL(1) / ZC00800000 / */
/*     DATA SMALL(2) / Z000000000 / */
/*     DATA LARGE(1) / ZDFFFFFFFF / */
/*     DATA LARGE(2) / ZFFFFFFFFF / */
/*     DATA RIGHT(1) / ZCC5800000 / */
/*     DATA RIGHT(2) / Z000000000 / */
/*     DATA DIVER(1) / ZCC6800000 / */
/*     DATA DIVER(2) / Z000000000 / */
/*     DATA LOG10(1) / ZD00E730E7 / */
/*     DATA LOG10(2) / ZC77800DC0 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O0000000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O0007777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O7770000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O7777777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA SMALL(1) / Z"3001800000000000" / */
/*     DATA SMALL(2) / Z"3001000000000000" / */
/*     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA LARGE(2) / Z"4FFE000000000000" / */
/*     DATA RIGHT(1) / Z"3FD2800000000000" / */
/*     DATA RIGHT(2) / Z"3FD2000000000000" / */
/*     DATA DIVER(1) / Z"3FD3800000000000" / */
/*     DATA DIVER(2) / Z"3FD3000000000000" / */
/*     DATA LOG10(1) / Z"3FFF9A209A84FBCF" / */
/*     DATA LOG10(2) / Z"3FFFF7988F8959AC" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA SMALL(1) / 00564000000000000000B / */
/*     DATA SMALL(2) / 00000000000000000000B / */
/*     DATA LARGE(1) / 37757777777777777777B / */
/*     DATA LARGE(2) / 37157777777777777777B / */
/*     DATA RIGHT(1) / 15624000000000000000B / */
/*     DATA RIGHT(2) / 00000000000000000000B / */
/*     DATA DIVER(1) / 15634000000000000000B / */
/*     DATA DIVER(2) / 00000000000000000000B / */
/*     DATA LOG10(1) / 17164642023241175717B / */
/*     DATA LOG10(2) / 16367571421742254654B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn OR -pd8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CC0000000000000' / */
/*     DATA DMACH(4) / Z'3CD0000000000000' / */
/*     DATA DMACH(5) / Z'3FF34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F900000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F910000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE CRAY */

/*     DATA SMALL(1) / 201354000000000000000B / */
/*     DATA SMALL(2) / 000000000000000000000B / */
/*     DATA LARGE(1) / 577767777777777777777B / */
/*     DATA LARGE(2) / 000007777777777777774B / */
/*     DATA RIGHT(1) / 376434000000000000000B / */
/*     DATA RIGHT(2) / 000000000000000000000B / */
/*     DATA DIVER(1) / 376444000000000000000B / */
/*     DATA DIVER(2) / 000000000000000000000B / */
/*     DATA LOG10(1) / 377774642023241175717B / */
/*     DATA LOG10(2) / 000007571421742254654B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */
/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC DMACH(5) */

/*     DATA SMALL /    20K, 3*0 / */
/*     DATA LARGE / 77777K, 3*177777K / */
/*     DATA RIGHT / 31420K, 3*0 / */
/*     DATA DIVER / 32020K, 3*0 / */
/*     DATA LOG10 / 40423K, 42023K, 50237K, 74776K / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA DMACH(1) / '0000000000000010'X / */
/*     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X / */
/*     DATA DMACH(3) / '0000000000003CC0'X / */
/*     DATA DMACH(4) / '0000000000003CD0'X / */
/*     DATA DMACH(5) / '79FF509F44133FF3'X / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FORMAT */

/*     DATA DMACH(1) / '0010000000000000'X / */
/*     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X / */
/*     DATA DMACH(3) / '3CA0000000000000'X / */
/*     DATA DMACH(4) / '3CB0000000000000'X / */
/*     DATA DMACH(5) / '3FD34413509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/ */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/ */
/*     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/ */
/*     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/ */
/*     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/ */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /        128,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /       9344,           0 / */
/*     DATA DIVER(1), DIVER(2) /       9472,           0 / */
/*     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 / */

/*     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /         16,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /      15552,           0 / */
/*     DATA DIVER(1), DIVER(2) /      15568,           0 / */
/*     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 / */

/*     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION) */

/*     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X / */
/*     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X / */
/*     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X / */
/*     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X / */
/*     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '37777577 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000334 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '10237777 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     THREE WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 / */
/*     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B / */
/*     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B / */
/*     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) /  40000B,       0 / */
/*     DATA SMALL(3), SMALL(4) /       0,       1 / */
/*     DATA LARGE(1), LARGE(2) /  77777B, 177777B / */
/*     DATA LARGE(3), LARGE(4) / 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) /  40000B,       0 / */
/*     DATA RIGHT(3), RIGHT(4) /       0,    225B / */
/*     DATA DIVER(1), DIVER(2) /  40000B,       0 / */
/*     DATA DIVER(3), DIVER(4) /       0,    227B / */
/*     DATA LOG10(1), LOG10(2) /  46420B,  46502B / */
/*     DATA LOG10(3), LOG10(4) /  76747B, 176377B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B / */
/*     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B / */
/*     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B / */
/*     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B / */
/*     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF / */

/*     MACHINE CONSTANTS FOR THE IBM PC */
/*     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION */
/*     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087. */

/*     DATA SMALL(1) / 2.23D-308  / */
/*     DATA LARGE(1) / 1.79D+308  / */
/*     DATA RIGHT(1) / 1.11D-16   / */
/*     DATA DIVER(1) / 2.22D-16   / */
/*     DATA LOG10(1) / 0.301029995663981195D0 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    8388608,           0 / */
/*     DATA LARGE(1), LARGE(2) / 2147483647,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /  612368384,           0 / */
/*     DATA DIVER(1), DIVER(2) /  620756992,           0 / */
/*     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 / */

/*     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 / */
/*     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 / */
/*     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 / */
/*     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    128,      0 / */
/*     DATA SMALL(3), SMALL(4) /      0,      0 / */
/*     DATA LARGE(1), LARGE(2) /  32767,     -1 / */
/*     DATA LARGE(3), LARGE(4) /     -1,     -1 / */
/*     DATA RIGHT(1), RIGHT(2) /   9344,      0 / */
/*     DATA RIGHT(3), RIGHT(4) /      0,      0 / */
/*     DATA DIVER(1), DIVER(2) /   9472,      0 / */
/*     DATA DIVER(3), DIVER(4) /      0,      0 / */
/*     DATA LOG10(1), LOG10(2) /  16282,   8346 / */
/*     DATA LOG10(3), LOG10(4) / -31493, -12296 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA SMALL(3), SMALL(4) / O000000, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA LARGE(3), LARGE(4) / O177777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O022200, O000000 / */
/*     DATA RIGHT(3), RIGHT(4) / O000000, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O022400, O000000 / */
/*     DATA DIVER(3), DIVER(4) / O000000, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020232 / */
/*     DATA LOG10(3), LOG10(4) / O102373, O147770 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE SUN 386i */

/*     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' / */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' / */
/*     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF' */
/*     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 / */

/* ***FIRST EXECUTABLE STATEMENT  D1MACH */
/*<    >*/
    if (*i < 1 || *i > 5) {
	xermsg_("SLATEC", "D1MACH", "I OUT OF BOUNDS", &c__1, &c__2, 6L, 6L, 
		15L);
    }

/*<       D1MACH = DMACH(I) >*/
    ret_val = dmach[*i - 1];
/*<       RETURN >*/
    return ret_val;

/*<       END >*/
} /* d1mach_ */

#undef right
#undef diver
#undef small
#undef large
#undef dmach
#undef log10


#endif






/*<       DOUBLE PRECISION FUNCTION D1MACH (I) >*/
doublereal d1mach_(integer *i)
{
    /* Initialized data */

    static struct {
	integer e_1[10];
	doublereal fill_2[1];
	doublereal e_3;
	} equiv_4 = { 0, 1048576, -1, 2146435071, 0, 1017118720, 0, 
		1018167296, 1352628735, 1070810131, {0}, 0. };


    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
#define log10 ((integer *)&equiv_4 + 8)
#define dmach ((doublereal *)&equiv_4)
#define large ((integer *)&equiv_4 + 2)
#define small ((integer *)&equiv_4)
#define diver ((integer *)&equiv_4 + 6)
#define right ((integer *)&equiv_4 + 4)
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D1MACH */
/* ***PURPOSE  Return floating point machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   D1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument, and can be referenced as follows: */

/*        D = D1MACH(I) */

/*   where I=1,...,5.  The (output) value of D above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   D1MACH( 3) = B**(-T), the smallest relative spacing. */
/*   D1MACH( 4) = B**(1-T), the largest relative spacing. */
/*   D1MACH( 5) = LOG10(B) */

/*   Assume double precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(14) = T, the number of base-B digits. */
/*   I1MACH(15) = EMIN, the smallest exponent E. */
/*   I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for 
*/
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890213  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   900911  Added SUN 386i constants.  (WRB) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added CONVEX -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/* ***END PROLOGUE  D1MACH */

/*<       INTEGER SMALL(4) >*/
/*<       INTEGER LARGE(4) >*/
/*<       INTEGER RIGHT(4) >*/
/*<       INTEGER DIVER(4) >*/
/*<       INTEGER LOG10(4) >*/

/*<       DOUBLE PRECISION DMACH(5) >*/
/*<       SAVE DMACH >*/

/*<       EQUIVALENCE (DMACH(1),SMALL(1)) >*/
/*<       EQUIVALENCE (DMACH(2),LARGE(1)) >*/
/*<       EQUIVALENCE (DMACH(3),RIGHT(1)) >*/
/*<       EQUIVALENCE (DMACH(4),DIVER(1)) >*/
/*<       EQUIVALENCE (DMACH(5),LOG10(1)) >*/

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 / */
/*     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 / */
/*     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 / */
/*     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA SMALL(1) / ZC00800000 / */
/*     DATA SMALL(2) / Z000000000 / */
/*     DATA LARGE(1) / ZDFFFFFFFF / */
/*     DATA LARGE(2) / ZFFFFFFFFF / */
/*     DATA RIGHT(1) / ZCC5800000 / */
/*     DATA RIGHT(2) / Z000000000 / */
/*     DATA DIVER(1) / ZCC6800000 / */
/*     DATA DIVER(2) / Z000000000 / */
/*     DATA LOG10(1) / ZD00E730E7 / */
/*     DATA LOG10(2) / ZC77800DC0 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O0000000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O0007777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O7770000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O7777777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA SMALL(1) / Z"3001800000000000" / */
/*     DATA SMALL(2) / Z"3001000000000000" / */
/*     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA LARGE(2) / Z"4FFE000000000000" / */
/*     DATA RIGHT(1) / Z"3FD2800000000000" / */
/*     DATA RIGHT(2) / Z"3FD2000000000000" / */
/*     DATA DIVER(1) / Z"3FD3800000000000" / */
/*     DATA DIVER(2) / Z"3FD3000000000000" / */
/*     DATA LOG10(1) / Z"3FFF9A209A84FBCF" / */
/*     DATA LOG10(2) / Z"3FFFF7988F8959AC" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA SMALL(1) / 00564000000000000000B / */
/*     DATA SMALL(2) / 00000000000000000000B / */
/*     DATA LARGE(1) / 37757777777777777777B / */
/*     DATA LARGE(2) / 37157777777777777777B / */
/*     DATA RIGHT(1) / 15624000000000000000B / */
/*     DATA RIGHT(2) / 00000000000000000000B / */
/*     DATA DIVER(1) / 15634000000000000000B / */
/*     DATA DIVER(2) / 00000000000000000000B / */
/*     DATA LOG10(1) / 17164642023241175717B / */
/*     DATA LOG10(2) / 16367571421742254654B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fn OR -pd8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CC0000000000000' / */
/*     DATA DMACH(4) / Z'3CD0000000000000' / */
/*     DATA DMACH(5) / Z'3FF34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -fi COMPILER OPTION */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX */
/*     USING THE -p8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F900000000000000000000000000000' / */
/*     DATA DMACH(4) / Z'3F910000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE CRAY */

/*     DATA SMALL(1) / 201354000000000000000B / */
/*     DATA SMALL(2) / 000000000000000000000B / */
/*     DATA LARGE(1) / 577767777777777777777B / */
/*     DATA LARGE(2) / 000007777777777777774B / */
/*     DATA RIGHT(1) / 376434000000000000000B / */
/*     DATA RIGHT(2) / 000000000000000000000B / */
/*     DATA DIVER(1) / 376444000000000000000B / */
/*     DATA DIVER(2) / 000000000000000000000B / */
/*     DATA LOG10(1) / 377774642023241175717B / */
/*     DATA LOG10(2) / 000007571421742254654B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */
/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC DMACH(5) */

/*     DATA SMALL /    20K, 3*0 / */
/*     DATA LARGE / 77777K, 3*177777K / */
/*     DATA RIGHT / 31420K, 3*0 / */
/*     DATA DIVER / 32020K, 3*0 / */
/*     DATA LOG10 / 40423K, 42023K, 50237K, 74776K / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA DMACH(1) / '0000000000000010'X / */
/*     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X / */
/*     DATA DMACH(3) / '0000000000003CC0'X / */
/*     DATA DMACH(4) / '0000000000003CD0'X / */
/*     DATA DMACH(5) / '79FF509F44133FF3'X / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FORMAT */

/*     DATA DMACH(1) / '0010000000000000'X / */
/*     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X / */
/*     DATA DMACH(3) / '3CA0000000000000'X / */
/*     DATA DMACH(4) / '3CB0000000000000'X / */
/*     DATA DMACH(5) / '3FD34413509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/ */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/ */
/*     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/ */
/*     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/ */
/*     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/ */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /        128,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /       9344,           0 / */
/*     DATA DIVER(1), DIVER(2) /       9472,           0 / */
/*     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 / */

/*     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /         16,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /      15552,           0 / */
/*     DATA DIVER(1), DIVER(2) /      15568,           0 / */
/*     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 / */

/*     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION) */

/*     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X / */
/*     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X / */
/*     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X / */
/*     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X / */
/*     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '37777577 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000334 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '10237777 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     THREE WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 / */
/*     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B / */
/*     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B / */
/*     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) /  40000B,       0 / */
/*     DATA SMALL(3), SMALL(4) /       0,       1 / */
/*     DATA LARGE(1), LARGE(2) /  77777B, 177777B / */
/*     DATA LARGE(3), LARGE(4) / 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) /  40000B,       0 / */
/*     DATA RIGHT(3), RIGHT(4) /       0,    225B / */
/*     DATA DIVER(1), DIVER(2) /  40000B,       0 / */
/*     DATA DIVER(3), DIVER(4) /       0,    227B / */
/*     DATA LOG10(1), LOG10(2) /  46420B,  46502B / */
/*     DATA LOG10(3), LOG10(4) /  76747B, 176377B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B / */
/*     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B / */
/*     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B / */
/*     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B / */
/*     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF / */

/*     MACHINE CONSTANTS FOR THE IBM PC */
/*     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION */
/*     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087. */

/*<       DATA SMALL(2), SMALL(1) / Z'00100000', Z'00000000' / >*/
/*<       DATA LARGE(2), LARGE(1) / Z'7FEFFFFF', Z'FFFFFFFF' / >*/
/*<       DATA RIGHT(2), RIGHT(1) / Z'3CA00000', Z'00000000' / >*/
/*<       DATA DIVER(2), DIVER(1) / Z'3CB00000', Z'00000000' / >*/
/*<       DATA LOG10(2), LOG10(1) / Z'3FD34413', Z'509F79FF' / >*/
/*      DATA SMALL(1) / 2.23D-308  / */
/*      DATA LARGE(1) / 1.79D+308  / */
/*      DATA RIGHT(1) / 1.11D-16   / */
/*      DATA DIVER(1) / 2.22D-16   / */
/*      DATA LOG10(1) / 0.301029995663981195D0 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    8388608,           0 / */
/*     DATA LARGE(1), LARGE(2) / 2147483647,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /  612368384,           0 / */
/*     DATA DIVER(1), DIVER(2) /  620756992,           0 / */
/*     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 / */

/*     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 / */
/*     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 / */
/*     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 / */
/*     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    128,      0 / */
/*     DATA SMALL(3), SMALL(4) /      0,      0 / */
/*     DATA LARGE(1), LARGE(2) /  32767,     -1 / */
/*     DATA LARGE(3), LARGE(4) /     -1,     -1 / */
/*     DATA RIGHT(1), RIGHT(2) /   9344,      0 / */
/*     DATA RIGHT(3), RIGHT(4) /      0,      0 / */
/*     DATA DIVER(1), DIVER(2) /   9472,      0 / */
/*     DATA DIVER(3), DIVER(4) /      0,      0 / */
/*     DATA LOG10(1), LOG10(2) /  16282,   8346 / */
/*     DATA LOG10(3), LOG10(4) / -31493, -12296 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA SMALL(3), SMALL(4) / O000000, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA LARGE(3), LARGE(4) / O177777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O022200, O000000 / */
/*     DATA RIGHT(3), RIGHT(4) / O000000, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O022400, O000000 / */
/*     DATA DIVER(3), DIVER(4) / O000000, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020232 / */
/*     DATA LOG10(3), LOG10(4) / O102373, O147770 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA DMACH(1) / Z'0010000000000000' / */
/*     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3CA0000000000000' / */
/*     DATA DMACH(4) / Z'3CB0000000000000' / */
/*     DATA DMACH(5) / Z'3FD34413509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA DMACH(1) / Z'00010000000000000000000000000000' / */
/*     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' / */
/*     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' / */
/* 4413509F79FF' / */

/*     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' / */
/*     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' / */

/*     MACHINE CONSTANTS FOR THE SUN 386i */

/*     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' / */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' / */
/*     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF' */
/*     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 / */

/* ***FIRST EXECUTABLE STATEMENT  D1MACH */
/*<    >*/
    if (*i < 1 || *i > 5) {
	xermsg_("SLATEC", "D1MACH", "I OUT OF BOUNDS", &c__1, &c__2, 6L, 6L, 
		15L);
    }

/*<       D1MACH = DMACH(I) >*/
    ret_val = dmach[*i - 1];
/*<       RETURN >*/
    return ret_val;

/*<       END >*/
} /* d1mach_ */

#undef right
#undef diver
#undef small
#undef large
#undef dmach
#undef log10













/*<    >*/
/* Subroutine */ int dckder_(integer *m, integer *n, doublereal *x, 
	doublereal *fvec, doublereal *fjac, integer *ldfjac, doublereal *xp, 
	doublereal *fvecp, integer *mode, doublereal *err)
{
    /* Initialized data */

    static doublereal factor = 100.;
    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    /*double sqrt(doublereal), d_lg10(doublereal *);*/

    /* Local variables */
    doublereal epsf, temp;
    integer i, j;
    extern doublereal d1mach_(integer *);
    doublereal epsmch, epslog, eps;

/* ***BEGIN PROLOGUE  DCKDER */
/* ***END PROLOGUE  DCKDER */
/*<       INTEGER I, J, LDFJAC, M, MODE, N >*/
/*<    >*/
/*<       SAVE FACTOR, ONE, ZERO >*/
/*<       DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/ >*/
    /* Parameter adjustments */
    --err;
    --fvecp;
    --xp;
    fjac_dim1 = *ldfjac;
    fjac_offset = fjac_dim1 + 1;
    fjac -= fjac_offset;
    --fvec;
    --x;

    /* Function Body */

/*     EPSMCH IS THE MACHINE PRECISION. */

/* ***FIRST EXECUTABLE STATEMENT  DCKDER */
/*<       EPSMCH = D1MACH(4) >*/
    epsmch = d1mach_(&c__4);

/*<       EPS = SQRT(EPSMCH) >*/
    eps = sqrt(epsmch);

/*<       IF (MODE .EQ. 2) GO TO 20 >*/
    if (*mode == 2) {
	goto L20;
    }

/*        MODE = 1. */

/*<          DO 10 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             TEMP = EPS*ABS(X(J)) >*/
	temp = eps * (d__1 = x[j], abs(d__1));
/*<             IF (TEMP .EQ. ZERO) TEMP = EPS >*/
	if (temp == zero) {
	    temp = eps;
	}
/*<             XP(J) = X(J) + TEMP >*/
	xp[j] = x[j] + temp;
/*<    10       CONTINUE >*/
/* L10: */
    }
/*<          GO TO 70 >*/
    goto L70;
/*<    20 CONTINUE >*/
L20:

/*        MODE = 2. */

/*<          EPSF = FACTOR*EPSMCH >*/
    epsf = factor * epsmch;
/*<          EPSLOG = LOG10(EPS) >*/
    epslog = d_lg10(&eps);
/*<          DO 30 I = 1, M >*/
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/*<             ERR(I) = ZERO >*/
	err[i] = zero;
/*<    30       CONTINUE >*/
/* L30: */
    }
/*<          DO 50 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             TEMP = ABS(X(J)) >*/
	temp = (d__1 = x[j], abs(d__1));
/*<             IF (TEMP .EQ. ZERO) TEMP = ONE >*/
	if (temp == zero) {
	    temp = one;
	}
/*<             DO 40 I = 1, M >*/
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
/*<                ERR(I) = ERR(I) + TEMP*FJAC(I,J) >*/
	    err[i] += temp * fjac[i + j * fjac_dim1];
/*<    40          CONTINUE >*/
/* L40: */
	}
/*<    50       CONTINUE >*/
/* L50: */
    }
/*<          DO 60 I = 1, M >*/
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/*<             TEMP = ONE >*/
	temp = one;
/*<    >*/
	if (fvec[i] != zero && fvecp[i] != zero && (d__1 = fvecp[i] - fvec[i],
		 abs(d__1)) >= epsf * (d__2 = fvec[i], abs(d__2))) {
	    temp = eps * (d__3 = (fvecp[i] - fvec[i]) / eps - err[i], abs(
		    d__3)) / ((d__4 = fvec[i], abs(d__4)) + (d__5 = fvecp[i], 
		    abs(d__5)));
	}
/*<             ERR(I) = ONE >*/
	err[i] = one;
/*<    >*/
	if (temp > epsmch && temp < eps) {
	    err[i] = (d_lg10(&temp) - epslog) / epslog;
	}
/*<             IF (TEMP .GE. EPS) ERR(I) = ZERO >*/
	if (temp >= eps) {
	    err[i] = zero;
	}
/*<    60       CONTINUE >*/
/* L60: */
    }
/*<    70 CONTINUE >*/
L70:

/*<       RETURN >*/
    return 0;

/*     LAST CARD OF SUBROUTINE DCKDER. */

/*<       END >*/
} /* dckder_ */

/*<       DOUBLE PRECISION FUNCTION DENORM (N, X) >*/
doublereal denorm_(integer *n, doublereal *x)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;
    static doublereal rdwarf = 3.834e-20;
    static doublereal rgiant = 1.304e19;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    doublereal xabs, x1max, x3max;
    integer i;
    doublereal s1, s2, s3, agiant, floatn;

/* ***BEGIN PROLOGUE  DENORM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNSQ and DNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (ENORM-S, DENORM-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an N-vector X, this function calculates the */
/*     Euclidean norm of X. */

/*     The Euclidean norm is computed by accumulating the sum of */
/*     squares in three different sums. The sums of squares for the */
/*     small and large components are scaled so that no overflows */
/*     occur. Non-destructive underflows are permitted. Underflows */
/*     and overflows do not occur in the computation of the unscaled */
/*     sum of squares for the intermediate components. */
/*     The definitions of small, intermediate and large components */
/*     depend on two constants, RDWARF and RGIANT. The main */
/*     restrictions on these constants are that RDWARF**2 not */
/*     underflow and RGIANT**2 not overflow. The constants */
/*     given here are suitable for every known computer. */

/*     The function statement is */

/*       DOUBLE PRECISION FUNCTION DENORM(N,X) */

/*     where */

/*       N is a positive integer input variable. */

/*       X is an input array of length N. */

/* ***SEE ALSO  DNSQ, DNSQE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DENORM */
/*<       INTEGER I, N >*/
/*<    >*/
/*<       SAVE ONE, ZERO, RDWARF, RGIANT >*/
/*<       DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/ >*/
    /* Parameter adjustments */
    --x;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DENORM */
/*<       S1 = ZERO >*/
    s1 = zero;
/*<       S2 = ZERO >*/
    s2 = zero;
/*<       S3 = ZERO >*/
    s3 = zero;
/*<       X1MAX = ZERO >*/
    x1max = zero;
/*<       X3MAX = ZERO >*/
    x3max = zero;
/*<       FLOATN = N >*/
    floatn = (doublereal) (*n);
/*<       AGIANT = RGIANT/FLOATN >*/
    agiant = rgiant / floatn;
/*<       DO 90 I = 1, N >*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/*<          XABS = ABS(X(I)) >*/
	xabs = (d__1 = x[i], abs(d__1));
/*<          IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70 >*/
	if (xabs > rdwarf && xabs < agiant) {
	    goto L70;
	}
/*<             IF (XABS .LE. RDWARF) GO TO 30 >*/
	if (xabs <= rdwarf) {
	    goto L30;
	}

/*              SUM FOR LARGE COMPONENTS. */

/*<                IF (XABS .LE. X1MAX) GO TO 10 >*/
	if (xabs <= x1max) {
	    goto L10;
	}
/*<                   S1 = ONE + S1*(X1MAX/XABS)**2 >*/
/* Computing 2nd power */
	d__1 = x1max / xabs;
	s1 = one + s1 * (d__1 * d__1);
/*<                   X1MAX = XABS >*/
	x1max = xabs;
/*<                   GO TO 20 >*/
	goto L20;
/*<    10          CONTINUE >*/
L10:
/*<                   S1 = S1 + (XABS/X1MAX)**2 >*/
/* Computing 2nd power */
	d__1 = xabs / x1max;
	s1 += d__1 * d__1;
/*<    20          CONTINUE >*/
L20:
/*<                GO TO 60 >*/
	goto L60;
/*<    30       CONTINUE >*/
L30:

/*              SUM FOR SMALL COMPONENTS. */

/*<                IF (XABS .LE. X3MAX) GO TO 40 >*/
	if (xabs <= x3max) {
	    goto L40;
	}
/*<                   S3 = ONE + S3*(X3MAX/XABS)**2 >*/
/* Computing 2nd power */
	d__1 = x3max / xabs;
	s3 = one + s3 * (d__1 * d__1);
/*<                   X3MAX = XABS >*/
	x3max = xabs;
/*<                   GO TO 50 >*/
	goto L50;
/*<    40          CONTINUE >*/
L40:
/*<                   IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2 >*/
	if (xabs != zero) {
/* Computing 2nd power */
	    d__1 = xabs / x3max;
	    s3 += d__1 * d__1;
	}
/*<    50          CONTINUE >*/
L50:
/*<    60       CONTINUE >*/
L60:
/*<             GO TO 80 >*/
	goto L80;
/*<    70    CONTINUE >*/
L70:

/*           SUM FOR INTERMEDIATE COMPONENTS. */

/*<             S2 = S2 + XABS**2 >*/
/* Computing 2nd power */
	d__1 = xabs;
	s2 += d__1 * d__1;
/*<    80    CONTINUE >*/
L80:
/*<    90    CONTINUE >*/
/* L90: */
	;
    }

/*     CALCULATION OF NORM. */

/*<       IF (S1 .EQ. ZERO) GO TO 100 >*/
    if (s1 == zero) {
	goto L100;
    }
/*<          DENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX) >*/
    ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
/*<          GO TO 130 >*/
    goto L130;
/*<   100 CONTINUE >*/
L100:
/*<          IF (S2 .EQ. ZERO) GO TO 110 >*/
    if (s2 == zero) {
	goto L110;
    }
/*<    >*/
    if (s2 >= x3max) {
	ret_val = sqrt(s2 * (one + x3max / s2 * (x3max * s3)));
    }
/*<    >*/
    if (s2 < x3max) {
	ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
    }
/*<             GO TO 120 >*/
    goto L120;
/*<   110    CONTINUE >*/
L110:
/*<             DENORM = X3MAX*SQRT(S3) >*/
    ret_val = x3max * sqrt(s3);
/*<   120    CONTINUE >*/
L120:
/*<   130 CONTINUE >*/
L130:
/*<       RETURN >*/
    return ret_val;

/*     LAST CARD OF FUNCTION DENORM. */

/*<       END >*/
} /* denorm_ */

/*<    >*/
/* Subroutine */ int dfdjc3_(S_fp fcn, integer *m, integer *n, doublereal *x, 
	doublereal *fvec, doublereal *fjac, integer *ldfjac, integer *iflag, 
	doublereal *epsfcn, doublereal *wa, void *pParam)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    doublereal temp, h;
    integer i, j;
    extern doublereal d1mach_(integer *);
    doublereal epsmch, eps;

/* ***BEGIN PROLOGUE  DFDJC3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNLS1 and DNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (FDJAC3-S, DFDJC3-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  **** Double Precision version of FDJAC3 **** */

/*     This subroutine computes a forward-difference approximation */
/*     to the M by N Jacobian matrix associated with a specified */
/*     problem of M functions in N variables. */

/*     The subroutine statement is */

/*       SUBROUTINE DFDJC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA) */

/*     where */

/*       FCN is the name of the user-supplied subroutine which */
/*         calculates the functions. FCN must be declared */
/*         in an external statement in the user calling */
/*         program, and should be written as follows. */

/*         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC) */
/*         INTEGER LDFJAC,M,N,IFLAG */
/*         DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N) */
/*         ---------- */
/*         When IFLAG.EQ.1 calculate the functions at X and */
/*         return this vector in FVEC. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by FCN unless */
/*         the user wants to terminate execution of DFDJC3. */
/*         In this case set IFLAG to a negative integer. */

/*       M is a positive integer input variable set to the number */
/*         of functions. */

/*       N is a positive integer input variable set to the number */
/*         of variables. N must not exceed M. */

/*       X is an input array of length N. */

/*       FVEC is an input array of length M which must contain the */
/*         functions evaluated at X. */

/*       FJAC is an output M by N array which contains the */
/*         approximation to the Jacobian matrix evaluated at X. */

/*       LDFJAC is a positive integer input variable not less than M */
/*         which specifies the leading dimension of the array FJAC. */

/*       IFLAG is an integer variable which can be used to terminate */
/*         THE EXECUTION OF DFDJC3. See description of FCN. */

/*       EPSFCN is an input variable used in determining a suitable */
/*         step length for the forward-difference approximation. This */
/*         approximation assumes that the relative errors in the */
/*         functions are of the order of EPSFCN. If EPSFCN is less */
/*         than the machine precision, it is assumed that the relative */
/*         errors in the functions are of the order of the machine */
/*         precision. */

/*       WA is a work array of length M. */

/* ***SEE ALSO  DNLS1, DNLS1E */
/* ***ROUTINES CALLED  D1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DFDJC3 */
/*<       INTEGER M,N,LDFJAC,IFLAG >*/
/*<       DOUBLE PRECISION EPSFCN >*/
/*<       DOUBLE PRECISION X(*),FVEC(*),FJAC(LDFJAC,*),WA(*) >*/
/*<       INTEGER I,J >*/
/*<       DOUBLE PRECISION EPS,EPSMCH,H,TEMP,ZERO >*/
/*<       DOUBLE PRECISION D1MACH >*/
/*<       SAVE ZERO >*/
/*<       DATA ZERO /0.0D0/ >*/
    /* Parameter adjustments */
    --wa;
    fjac_dim1 = *ldfjac;
    fjac_offset = fjac_dim1 + 1;
    fjac -= fjac_offset;
    --fvec;
    --x;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DFDJC3 */
/*<       EPSMCH = D1MACH(4) >*/
    epsmch = d1mach_(&c__4);

/*<       EPS = SQRT(MAX(EPSFCN,EPSMCH)) >*/
    eps = sqrt((max(*epsfcn,epsmch)));
/*      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES */
/*           ARE TO BE RETURNED BY FCN. */
/*<       IFLAG = 1 >*/
    *iflag = 1;
/*<       DO 20 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          TEMP = X(J) >*/
	temp = x[j];
/*<          H = EPS*ABS(TEMP) >*/
	h = eps * abs(temp);
/*<          IF (H .EQ. ZERO) H = EPS >*/
	if (h == zero) {
	    h = eps;
	}
/*<          X(J) = TEMP + H >*/
	x[j] = temp + h;
/*<          CALL FCN(IFLAG,M,N,X,WA,FJAC,LDFJAC) >*/
	(*fcn)(iflag, m, n, &x[1], &wa[1], &fjac[fjac_offset], ldfjac, pParam);
/*<          IF (IFLAG .LT. 0) GO TO 30 >*/
	if (*iflag < 0) {
	    goto L30;
	}
/*<          X(J) = TEMP >*/
	x[j] = temp;
/*<          DO 10 I = 1, M >*/
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
/*<             FJAC(I,J) = (WA(I) - FVEC(I))/H >*/
	    fjac[i + j * fjac_dim1] = (wa[i] - fvec[i]) / h;
/*<    10       CONTINUE >*/
/* L10: */
	}

/*<    20    CONTINUE >*/
/* L20: */
    }
/*<    30 CONTINUE >*/
L30:
/*<       RETURN >*/
    return 0;

/*     LAST CARD OF SUBROUTINE DFDJC3. */

/*<       END >*/
} /* dfdjc3_ */

/*<    >*/
/* Subroutine */ int dmpar_(integer *n, doublereal *r, integer *ldr, integer *
	ipvt, doublereal *diag, doublereal *qtb, doublereal *delta, 
	doublereal *par, doublereal *x, doublereal *sigma, doublereal *wa1, 
	doublereal *wa2)
{
    /* Initialized data */

    static doublereal p1 = .1;
    static doublereal p001 = .001;
    static doublereal zero = 0.;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    doublereal parc, parl;
    integer iter;
    doublereal temp, paru;
    integer i, j, k, l;
    doublereal dwarf;
    integer nsing;
    doublereal gnorm;
    extern doublereal d1mach_(integer *);
    doublereal fp;
    extern doublereal denorm_(integer *, doublereal *);
    doublereal dxnorm;
    integer jm1, jp1;
    extern /* Subroutine */ int dqrslv_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);
    doublereal sum;

/* ***BEGIN PROLOGUE  DMPAR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNLS1 and DNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (LMPAR-S, DMPAR-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   **** Double Precision version of LMPAR **** */

/*     Given an M by N matrix A, an N by N nonsingular DIAGONAL */
/*     matrix D, an M-vector B, and a positive number DELTA, */
/*     the problem is to determine a value for the parameter */
/*     PAR such that if X solves the system */

/*           A*X = B ,     SQRT(PAR)*D*X = 0 , */

/*     in the least squares sense, and DXNORM is the Euclidean */
/*     norm of D*X, then either PAR is zero and */

/*           (DXNORM-DELTA) .LE. 0.1*DELTA , */

/*     or PAR is positive and */

/*           ABS(DXNORM-DELTA) .LE. 0.1*DELTA . */

/*     This subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     QR factorization, with column pivoting, of A. That is, if */
/*     A*P = Q*R, where P is a permutation matrix, Q has orthogonal */
/*     columns, and R is an upper triangular matrix with diagonal */
/*     elements of nonincreasing magnitude, then DMPAR expects */
/*     the full upper triangle of R, the permutation matrix P, */
/*     and the first N components of (Q TRANSPOSE)*B. On output */
/*     DMPAR also provides an upper triangular matrix S such that */

/*            T   T                   T */
/*           P *(A *A + PAR*D*D)*P = S *S . */

/*     S is employed within DMPAR and may be of separate interest. */

/*     Only a few iterations are generally needed for convergence */
/*     of the algorithm. If, however, the limit of 10 iterations */
/*     is reached, then the output PAR will contain the best */
/*     value obtained so far. */

/*     The subroutine statement is */

/*       SUBROUTINE DMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA, */
/*                        WA1,WA2) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an N by N array. On input the full upper triangle */
/*         must contain the full upper triangle of the matrix R. */
/*         On output the full upper triangle is unaltered, and the */
/*         strict lower triangle contains the strict upper triangle */
/*         (transposed) of the upper triangular matrix S. */

/*       LDR is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array R. */

/*       IPVT is an integer input array of length N which defines the */
/*         permutation matrix P such that A*P = Q*R. Column J of P */
/*         is column IPVT(J) of the identity matrix. */

/*       DIAG is an input array of length N which must contain the */
/*         diagonal elements of the matrix D. */

/*       QTB is an input array of length N which must contain the first */
/*         N elements of the vector (Q TRANSPOSE)*B. */

/*       DELTA is a positive input variable which specifies an upper */
/*         bound on the Euclidean norm of D*X. */

/*       PAR is a nonnegative variable. On input PAR contains an */
/*         initial estimate of the Levenberg-Marquardt parameter. */
/*         On output PAR contains the final estimate. */

/*       X is an output array of length N which contains the least */
/*         squares solution of the system A*X = B, SQRT(PAR)*D*X = 0, */
/*         for the output PAR. */

/*       SIGMA is an output array of length N which contains the */
/*         diagonal elements of the upper triangular matrix S. */

/*       WA1 and WA2 are work arrays of length N. */

/* ***SEE ALSO  DNLS1, DNLS1E */
/* ***ROUTINES CALLED  D1MACH, DENORM, DQRSLV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DMPAR */
/*<       INTEGER N,LDR >*/
/*<       INTEGER IPVT(*) >*/
/*<       DOUBLE PRECISION DELTA,PAR >*/
/*<    >*/
/*<       INTEGER I,ITER,J,JM1,JP1,K,L,NSING >*/
/*<    >*/
/*<       DOUBLE PRECISION D1MACH,DENORM >*/
/*<       SAVE P1, P001, ZERO >*/
/*<       DATA P1,P001,ZERO /1.0D-1,1.0D-3,0.0D0/ >*/
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --sigma;
    --x;
    --qtb;
    --diag;
    --ipvt;
    r_dim1 = *ldr;
    r_offset = r_dim1 + 1;
    r -= r_offset;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DMPAR */
/*<       DWARF = D1MACH(1) >*/
    dwarf = d1mach_(&c__1);

/*     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE */
/*     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION. */

/*<       NSING = N >*/
    nsing = *n;
/*<       DO 10 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          WA1(J) = QTB(J) >*/
	wa1[j] = qtb[j];
/*<          IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1 >*/
	if (r[j + j * r_dim1] == zero && nsing == *n) {
	    nsing = j - 1;
	}
/*<          IF (NSING .LT. N) WA1(J) = ZERO >*/
	if (nsing < *n) {
	    wa1[j] = zero;
	}
/*<    10    CONTINUE >*/
/* L10: */
    }
/*<       IF (NSING .LT. 1) GO TO 50 >*/
    if (nsing < 1) {
	goto L50;
    }
/*<       DO 40 K = 1, NSING >*/
    i__1 = nsing;
    for (k = 1; k <= i__1; ++k) {
/*<          J = NSING - K + 1 >*/
	j = nsing - k + 1;
/*<          WA1(J) = WA1(J)/R(J,J) >*/
	wa1[j] /= r[j + j * r_dim1];
/*<          TEMP = WA1(J) >*/
	temp = wa1[j];
/*<          JM1 = J - 1 >*/
	jm1 = j - 1;
/*<          IF (JM1 .LT. 1) GO TO 30 >*/
	if (jm1 < 1) {
	    goto L30;
	}
/*<          DO 20 I = 1, JM1 >*/
	i__2 = jm1;
	for (i = 1; i <= i__2; ++i) {
/*<             WA1(I) = WA1(I) - R(I,J)*TEMP >*/
	    wa1[i] -= r[i + j * r_dim1] * temp;
/*<    20       CONTINUE >*/
/* L20: */
	}
/*<    30    CONTINUE >*/
L30:
/*<    40    CONTINUE >*/
/* L40: */
	;
    }
/*<    50 CONTINUE >*/
L50:
/*<       DO 60 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          L = IPVT(J) >*/
	l = ipvt[j];
/*<          X(L) = WA1(J) >*/
	x[l] = wa1[j];
/*<    60    CONTINUE >*/
/* L60: */
    }

/*     INITIALIZE THE ITERATION COUNTER. */
/*     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST */
/*     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION. */

/*<       ITER = 0 >*/
    iter = 0;
/*<       DO 70 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          WA2(J) = DIAG(J)*X(J) >*/
	wa2[j] = diag[j] * x[j];
/*<    70    CONTINUE >*/
/* L70: */
    }
/*<       DXNORM = DENORM(N,WA2) >*/
    dxnorm = denorm_(n, &wa2[1]);
/*<       FP = DXNORM - DELTA >*/
    fp = dxnorm - *delta;
/*<       IF (FP .LE. P1*DELTA) GO TO 220 >*/
    if (fp <= p1 * *delta) {
	goto L220;
    }

/*     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON */
/*     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF */
/*     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO. */

/*<       PARL = ZERO >*/
    parl = zero;
/*<       IF (NSING .LT. N) GO TO 120 >*/
    if (nsing < *n) {
	goto L120;
    }
/*<       DO 80 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          L = IPVT(J) >*/
	l = ipvt[j];
/*<          WA1(J) = DIAG(L)*(WA2(L)/DXNORM) >*/
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
/*<    80    CONTINUE >*/
/* L80: */
    }
/*<       DO 110 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM = ZERO >*/
	sum = zero;
/*<          JM1 = J - 1 >*/
	jm1 = j - 1;
/*<          IF (JM1 .LT. 1) GO TO 100 >*/
	if (jm1 < 1) {
	    goto L100;
	}
/*<          DO 90 I = 1, JM1 >*/
	i__2 = jm1;
	for (i = 1; i <= i__2; ++i) {
/*<             SUM = SUM + R(I,J)*WA1(I) >*/
	    sum += r[i + j * r_dim1] * wa1[i];
/*<    90       CONTINUE >*/
/* L90: */
	}
/*<   100    CONTINUE >*/
L100:
/*<          WA1(J) = (WA1(J) - SUM)/R(J,J) >*/
	wa1[j] = (wa1[j] - sum) / r[j + j * r_dim1];
/*<   110    CONTINUE >*/
/* L110: */
    }
/*<       TEMP = DENORM(N,WA1) >*/
    temp = denorm_(n, &wa1[1]);
/*<       PARL = ((FP/DELTA)/TEMP)/TEMP >*/
    parl = fp / *delta / temp / temp;
/*<   120 CONTINUE >*/
L120:

/*     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION. */

/*<       DO 140 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          SUM = ZERO >*/
	sum = zero;
/*<          DO 130 I = 1, J >*/
	i__2 = j;
	for (i = 1; i <= i__2; ++i) {
/*<             SUM = SUM + R(I,J)*QTB(I) >*/
	    sum += r[i + j * r_dim1] * qtb[i];
/*<   130       CONTINUE >*/
/* L130: */
	}
/*<          L = IPVT(J) >*/
	l = ipvt[j];
/*<          WA1(J) = SUM/DIAG(L) >*/
	wa1[j] = sum / diag[l];
/*<   140    CONTINUE >*/
/* L140: */
    }
/*<       GNORM = DENORM(N,WA1) >*/
    gnorm = denorm_(n, &wa1[1]);
/*<       PARU = GNORM/DELTA >*/
    paru = gnorm / *delta;
/*<       IF (PARU .EQ. ZERO) PARU = DWARF/MIN(DELTA,P1) >*/
    if (paru == zero) {
	paru = dwarf / min(*delta,p1);
    }

/*     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU), */
/*     SET PAR TO THE CLOSER ENDPOINT. */

/*<       PAR = MAX(PAR,PARL) >*/
    *par = max(*par,parl);
/*<       PAR = MIN(PAR,PARU) >*/
    *par = min(*par,paru);
/*<       IF (PAR .EQ. ZERO) PAR = GNORM/DXNORM >*/
    if (*par == zero) {
	*par = gnorm / dxnorm;
    }

/*     BEGINNING OF AN ITERATION. */

/*<   150 CONTINUE >*/
L150:
/*<          ITER = ITER + 1 >*/
    ++iter;

/*        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR. */

/*<          IF (PAR .EQ. ZERO) PAR = MAX(DWARF,P001*PARU) >*/
    if (*par == zero) {
/* Computing MAX */
	d__1 = dwarf, d__2 = p001 * paru;
	*par = max(d__1,d__2);
    }
/*<          TEMP = SQRT(PAR) >*/
    temp = sqrt(*par);
/*<          DO 160 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             WA1(J) = TEMP*DIAG(J) >*/
	wa1[j] = temp * diag[j];
/*<   160       CONTINUE >*/
/* L160: */
    }
/*<          CALL DQRSLV(N,R,LDR,IPVT,WA1,QTB,X,SIGMA,WA2) >*/
    dqrslv_(n, &r[r_offset], ldr, &ipvt[1], &wa1[1], &qtb[1], &x[1], &sigma[1]
	    , &wa2[1]);
/*<          DO 170 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             WA2(J) = DIAG(J)*X(J) >*/
	wa2[j] = diag[j] * x[j];
/*<   170       CONTINUE >*/
/* L170: */
    }
/*<          DXNORM = DENORM(N,WA2) >*/
    dxnorm = denorm_(n, &wa2[1]);
/*<          TEMP = FP >*/
    temp = fp;
/*<          FP = DXNORM - DELTA >*/
    fp = dxnorm - *delta;

/*        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE */
/*        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL */
/*        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10. */

/*<    >*/
    if (abs(fp) <= p1 * *delta || parl == zero && fp <= temp && temp < zero ||
	     iter == 10) {
	goto L220;
    }

/*        COMPUTE THE NEWTON CORRECTION. */

/*<          DO 180 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             L = IPVT(J) >*/
	l = ipvt[j];
/*<             WA1(J) = DIAG(L)*(WA2(L)/DXNORM) >*/
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
/*<   180       CONTINUE >*/
/* L180: */
    }
/*<          DO 210 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             WA1(J) = WA1(J)/SIGMA(J) >*/
	wa1[j] /= sigma[j];
/*<             TEMP = WA1(J) >*/
	temp = wa1[j];
/*<             JP1 = J + 1 >*/
	jp1 = j + 1;
/*<             IF (N .LT. JP1) GO TO 200 >*/
	if (*n < jp1) {
	    goto L200;
	}
/*<             DO 190 I = JP1, N >*/
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
/*<                WA1(I) = WA1(I) - R(I,J)*TEMP >*/
	    wa1[i] -= r[i + j * r_dim1] * temp;
/*<   190          CONTINUE >*/
/* L190: */
	}
/*<   200       CONTINUE >*/
L200:
/*<   210       CONTINUE >*/
/* L210: */
	;
    }
/*<          TEMP = DENORM(N,WA1) >*/
    temp = denorm_(n, &wa1[1]);
/*<          PARC = ((FP/DELTA)/TEMP)/TEMP >*/
    parc = fp / *delta / temp / temp;

/*        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU. */

/*<          IF (FP .GT. ZERO) PARL = MAX(PARL,PAR) >*/
    if (fp > zero) {
	parl = max(parl,*par);
    }
/*<          IF (FP .LT. ZERO) PARU = MIN(PARU,PAR) >*/
    if (fp < zero) {
	paru = min(paru,*par);
    }

/*        COMPUTE AN IMPROVED ESTIMATE FOR PAR. */

/*<          PAR = MAX(PARL,PAR+PARC) >*/
/* Computing MAX */
    d__1 = parl, d__2 = *par + parc;
    *par = max(d__1,d__2);

/*        END OF AN ITERATION. */

/*<          GO TO 150 >*/
    goto L150;
/*<   220 CONTINUE >*/
L220:

/*     TERMINATION. */

/*<       IF (ITER .EQ. 0) PAR = ZERO >*/
    if (iter == 0) {
	*par = zero;
    }
/*<       RETURN >*/
    return 0;

/*     LAST CARD OF SUBROUTINE DMPAR. */

/*<       END >*/
} /* dmpar_ */

/*<    >*/
/* Subroutine */ int dnls1_(S_fp fcn, integer *iopt, integer *m, integer *n, 
	doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, 
	doublereal *ftol, doublereal *xtol, doublereal *gtol, integer *maxfev,
	 doublereal *epsfcn, doublereal *diag, integer *mode, doublereal *
	factor, integer *nprint, integer *info, integer *nfev, integer *njev, 
	integer *ipvt, doublereal *qtf, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3, doublereal *wa4, void * pParam)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal p1 = .1;
    static doublereal p5 = .5;
    static doublereal p25 = .25;
    static doublereal p75 = .75;
    static doublereal p0001 = 1e-4;
    static doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    logical sing;
    integer iter;
    doublereal temp;
    integer nrow;
    doublereal temp1, temp2;
    integer i, j, l, iflag;
    doublereal delta;
    extern /* Subroutine */ int dmpar_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    doublereal ratio;
    integer ijunk;
    doublereal fnorm, gnorm, pnorm;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int dfdjc3_(S_fp, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, void * pParam);
    doublereal xnorm, fnorm1;
    extern /* Subroutine */ int dckder_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    integer modech;
    extern /* Subroutine */ int dqrfac_(integer *, integer *, doublereal *, 
	    integer *, logical *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    doublereal actred, dirder, epsmch, prered;
    extern doublereal denorm_(integer *, doublereal *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dwupdt_(integer *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    doublereal par, err, sum;

/* ***BEGIN PROLOGUE  DNLS1 */
/* ***PURPOSE  Minimize the sum of the squares of M nonlinear functions */
/*       IOPT is an input variable which specifies how the Jacobian will 
*/
/*       FACTOR is a positive input variable used in determining the ini- 
*/
/*         tial step bound.  This bound is set to the product of FACTOR */
/*         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR 
*/
/*         itself.  In most cases FACTOR should lie in the interval */
/*         (.1,100.).  100. is a generally recommended value. */

/*       NPRINT is an integer input variable that enables controlled */
/*         printing of iterates if it is positive.  In this case, FCN is 
*/
/*         called with IFLAG = 0 at the beginning of the first iteration 
*/
/*         and every NPRINT iterations thereafter and immediately prior */
/*         to return, with X and FVEC available for printing. Appropriate 
*/
/*         print statements must be added to FCN (see example) and */
/*         FVEC should not be altered.  If NPRINT is not positive, no */
/*         special calls to FCN with IFLAG = 0 are made. */

/*       INFO is an integer output variable.  If the user has terminated 
*/
/*        execution, INFO is set to the (negative) value of IFLAG.  See */
/*        description of FCN and JAC. Otherwise, INFO is set as follows */

/*         INFO = 0  improper input parameters. */

/*         INFO = 1  both actual and predicted relative reductions in the 
*/
/*                   sum of squares are at most FTOL. */

/*         INFO = 2  relative error between two consecutive iterates is */
/*                   at most XTOL. */

/*         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold. */

/*         INFO = 4  the cosine of the angle between FVEC and any column 
*/
/*                   of the Jacobian is at most GTOL in absolute value. */

/*         INFO = 5  number of calls to FCN for function evaluation */
/*                   has reached MAXFEV. */

/*         INFO = 6  FTOL is too small.  No further reduction in the sum 
*/
/*                   of squares is possible. */

/*         INFO = 7  XTOL is too small.  No further improvement in the */
/*                   approximate solution X is possible. */

/*         INFO = 8  GTOL is too small.  FVEC is orthogonal to the */
/*                   columns of the Jacobian to machine precision. */

/*         Sections 4 and 5 contain more details about INFO. */

/*       NFEV is an integer output variable set to the number of calls to 
*/
/*         FCN for function evaluation. */

/*       NJEV is an integer output variable set to the number of */
/*         evaluations of the full Jacobian.  If IOPT=2, only one call to 
*/
/*         FCN is required for each evaluation of the full Jacobian. */
/*         If IOPT=3, the M calls to FCN are required. */
/*         If IOPT=1, then NJEV is set to zero. */

/*       IPVT is an integer output array of length N.  IPVT defines a */
/*         permutation matrix P such that JAC*P = Q*R, where JAC is the */
/*         final calculated Jacobian, Q is orthogonal (not stored), and R 
*/
/*         is upper triangular with diagonal elements of nonincreasing */
/*         magnitude.  Column J of P is column IPVT(J) of the identity */
/*         matrix. */

/*       QTF is an output array of length N which contains the first N */
/*         elements of the vector (Q transpose)*FVEC. */

/*       WA1, WA2, and WA3 are work arrays of length N. */

/*       WA4 is a work array of length M. */


/* 4. Successful Completion. */

/*       The accuracy of DNLS1 is controlled by the convergence parame- */
/*       ters FTOL, XTOL, and GTOL.  These parameters are used in tests */
/*       which make three types of comparisons between the approximation 
*/
/*       X and a solution XSOL.  DNLS1 terminates when any of the tests */
/*       is satisfied.  If any of the convergence parameters is less than 
*/
/*       the machine precision (as defined by the function R1MACH(4)), */
/*       then DNLS1 only attempts to satisfy the test defined by the */
/*       machine precision.  Further progress is not usually possible. */

/*       The tests assume that the functions are reasonably well behaved, 
*/
/*       and, if the Jacobian is supplied by the user, that the functions 
*/
/*       and the Jacobian are coded consistently.  If these conditions */
/*       are not satisfied, then DNLS1 may incorrectly indicate conver- */
/*       gence.  If the Jacobian is coded correctly or IOPT=1, */
/*       then the validity of the answer can be checked, for example, by 
*/
/*       rerunning DNLS1 with tighter tolerances. */

/*       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm 
*/
/*         of a vector Z, then this test attempts to guarantee that */

/*               ENORM(FVEC) .LE. (1+FTOL)*ENORM(FVECS), */

/*         where FVECS denotes the functions evaluated at XSOL.  If this 
*/
/*         condition is satisfied with FTOL = 10**(-K), then the final */
/*         residual norm ENORM(FVEC) has K significant decimal digits and 
*/
/*         INFO is set to 1 (or to 3 if the second test is also satis- */
/*         fied).  Unless high precision solutions are required, the */
/*         recommended value for FTOL is the square root of the machine */
/*         precision. */

/*       Second Convergence Test.  If D is the diagonal matrix whose */
/*         entries are defined by the array DIAG, then this test attempts 
*/
/*         to guarantee that */

/*               ENORM(D*(X-XSOL)) .LE. XTOL*ENORM(D*XSOL). */

/*         If this condition is satisfied with XTOL = 10**(-K), then the 
*/
/*         larger components of D*X have K significant decimal digits and 
*/
/*         INFO is set to 2 (or to 3 if the first test is also satis- */
/*         fied).  There is a danger that the smaller components of D*X */
/*         may have large relative errors, but if MODE = 1, then the */
/*         accuracy of the components of X is usually related to their */
/*         sensitivity.  Unless high precision solutions are required, */
/*         the recommended value for XTOL is the square root of the */
/*         machine precision. */

/*       Third Convergence Test.  This test is satisfied when the cosine 
*/
/*         of the angle between FVEC and any column of the Jacobian at X 
*/
/*         is at most GTOL in absolute value.  There is no clear rela- */
/*         tionship between this test and the accuracy of DNLS1, and */
/*         furthermore, the test is equally well satisfied at other crit- 
*/
/*         ical points, namely maximizers and saddle points.  Therefore, 
*/
/*         termination caused by this test (INFO = 4) should be examined 
*/
/*         carefully.  The recommended value for GTOL is zero. */


/* 5. Unsuccessful Completion. */

/*       Unsuccessful termination of DNLS1 can be due to improper input */
/*       parameters, arithmetic interrupts, or an excessive number of */
/*       function evaluations. */

/*       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1 */
/*         or IOPT .GT. 3, or N .LE. 0, or M .LT. N, or for IOPT=1 or 2 */
/*         LDFJAC .LT. M, or for IOPT=3 LDFJAC .LT. N, or FTOL .LT. 0.E0, 
*/
/*         or XTOL .LT. 0.E0, or GTOL .LT. 0.E0, or MAXFEV .LE. 0, or */
/*         FACTOR .LE. 0.E0. */

/*       Arithmetic Interrupts.  If these interrupts occur in the FCN */
/*         subroutine during an early stage of the computation, they may 
*/
/*         be caused by an unacceptable choice of X by DNLS1.  In this */
/*         case, it may be possible to remedy the situation by rerunning 
*/
/*         DNLS1 with a smaller value of FACTOR. */

/*       Excessive Number of Function Evaluations.  A reasonable value */
/*         for MAXFEV is 100*(N+1) for IOPT=2 or 3 and 200*(N+1) for */
/*         IOPT=1.  If the number of calls to FCN reaches MAXFEV, then */
/*         this indicates that the routine is converging very slowly */
/*         as measured by the progress of FVEC, and INFO is set to 5. */
/*         In this case, it may be helpful to restart DNLS1 with MODE */
/*         set to 1. */


/* 6. Characteristics of the Algorithm. */

/*       DNLS1 is a modification of the Levenberg-Marquardt algorithm. */
/*       Two of its main characteristics involve the proper use of */
/*       implicitly scaled variables (if MODE = 1) and an optimal choice 
*/
/*       for the correction.  The use of implicitly scaled variables */
/*       achieves scale invariance of DNLS1 and limits the size of the */
/*       correction in any direction where the functions are changing */
/*       rapidly.  The optimal choice of the correction guarantees (under 
*/
/*       reasonable conditions) global convergence from starting points */
/*       far from the solution and a fast rate of convergence for */
/*       problems with small residuals. */

/*       Timing.  The time required by DNLS1 to solve a given problem */
/*         depends on M and N, the behavior of the functions, the accu- */
/*         racy requested, and the starting point.  The number of arith- 
*/
/*         metic operations needed by DNLS1 is about N**3 to process each 
*/
/*         evaluation of the functions (call to FCN) and to process each 
*/
/*         evaluation of the Jacobian it takes M*N**2 for IOPT=2 (one */
/*         call to FCN), M*N**2 for IOPT=1 (N calls to FCN) and */
/*         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN */
/*         can be evaluated quickly, the timing of DNLS1 will be */
/*         strongly influenced by the time spent in FCN. */

/*       Storage.  DNLS1 requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and */
/*         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage */
/*         locations and N integer storage locations, in addition to */
/*         the storage required by the program.  There are no internally 
*/
/*         declared storage arrays. */

/* *Long Description: */

/* 7. Example. */

/*       The problem is to determine the values of X(1), X(2), and X(3) */
/*       which provide the best fit (in the least squares sense) of */

/*             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15 */

/*       to the data */

/*             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39, */
/*                  0.37,0.58,0.73,0.96,1.34,2.10,4.39), */

/*       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The */
/*       I-th component of FVEC is thus defined by */

/*             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))). */

/*       ********** */

/*       PROGRAM TEST */
/* C */
/* C     Driver for DNLS1 example. */
/* C */
/*       INTEGER J,IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV, */
/*      *        NWRITE */
/*       INTEGER IPVT(3) */
/*       DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,FNORM,EPSFCN */
/*       DOUBLE PRECISION X(3),FVEC(15),FJAC(15,3),DIAG(3),QTF(3), */
/*      *     WA1(3),WA2(3),WA3(3),WA4(15) */
/*       DOUBLE PRECISION DENORM,D1MACH */
/*       EXTERNAL FCN */
/*       DATA NWRITE /6/ */
/* C */
/*       IOPT = 1 */
/*       M = 15 */
/*       N = 3 */
/* C */
/* C     The following starting values provide a rough fit. */
/* C */
/*       X(1) = 1.E0 */
/*       X(2) = 1.E0 */
/*       X(3) = 1.E0 */
/* C */
/*       LDFJAC = 15 */
/* C */
/* C     Set FTOL and XTOL to the square root of the machine precision */
/* C     and GTOL to zero.  Unless high precision solutions are */
/* C     required, these are the recommended settings. */
/* C */
/*       FTOL = SQRT(R1MACH(4)) */
/*       XTOL = SQRT(R1MACH(4)) */
/*       GTOL = 0.E0 */
/* C */
/*       MAXFEV = 400 */
/*       EPSFCN = 0.0 */
/*       MODE = 1 */
/*       FACTOR = 1.E2 */
/*       NPRINT = 0 */
/* C */
/*       CALL DNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL, */
/*      *           GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT, */
/*      *           INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4) */
/*       FNORM = ENORM(M,FVEC) */
/*       WRITE (NWRITE,1000) FNORM,NFEV,NJEV,INFO,(X(J),J=1,N) */
/*       STOP */
/*  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 // */
/*      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 // */
/*      *        5X,' NUMBER OF JACOBIAN EVALUATIONS',I10 // */
/*      *        5X,' EXIT PARAMETER',16X,I10 // */
/*      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7) */
/*       END */
/*       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM) */
/* C     This is the form of the FCN routine if IOPT=1, */
/* C     that is, if the user does not calculate the Jacobian. */
/*       INTEGER I,M,N,IFLAG */
/*       DOUBLE PRECISION X(N),FVEC(M),Y(15) */
/*       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4 */
/*       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8), */
/*      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15) */
/*      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1, */
/*      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/ */
/* C */
/*       IF (IFLAG .NE. 0) GO TO 5 */
/* C */
/* C     Insert print statements here when NPRINT is positive. */
/* C */
/*       RETURN */
/*     5 CONTINUE */
/*       DO 10 I = 1, M */
/*          TMP1 = I */
/*          TMP2 = 16 - I */
/*          TMP3 = TMP1 */
/*          IF (I .GT. 8) TMP3 = TMP2 */
/*          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3)) */
/*    10    CONTINUE */
/*       RETURN */
/*       END */


/*       Results obtained with different compilers or machines */
/*       may be slightly different. */

/*       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01 */

/*       NUMBER OF FUNCTION EVALUATIONS        25 */

/*       NUMBER OF JACOBIAN EVALUATIONS         0 */

/*       EXIT PARAMETER                         1 */

/*       FINAL APPROXIMATE SOLUTION */

/*        0.8241058E-01  0.1133037E+01  0.2343695E+01 */


/*       For IOPT=2, FCN would be modified as follows to also */
/*       calculate the full Jacobian when IFLAG=2. */

/*       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC) */
/* C */
/* C     This is the form of the FCN routine if IOPT=2, */
/* C     that is, if the user calculates the full Jacobian. */
/* C */
/*       INTEGER I,LDFJAC,M,N,IFLAG */
/*       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),Y(15) */
/*       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4 */
/*       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8), */
/*      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15) */
/*      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1, */
/*      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/ */
/* C */
/*       IF (IFLAG .NE. 0) GO TO 5 */
/* C */
/* C     Insert print statements here when NPRINT is positive. */
/* C */
/*       RETURN */
/*     5 CONTINUE */
/*       IF(IFLAG.NE.1) GO TO 20 */
/*       DO 10 I = 1, M */
/*          TMP1 = I */
/*          TMP2 = 16 - I */
/*          TMP3 = TMP1 */
/*          IF (I .GT. 8) TMP3 = TMP2 */
/*          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3)) */
/*    10    CONTINUE */
/*       RETURN */
/* C */
/* C     Below, calculate the full Jacobian. */
/* C */
/*    20    CONTINUE */
/* C */
/*       DO 30 I = 1, M */
/*          TMP1 = I */
/*          TMP2 = 16 - I */
/*          TMP3 = TMP1 */
/*          IF (I .GT. 8) TMP3 = TMP2 */
/*          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2 */
/*          FJAC(I,1) = -1.E0 */
/*          FJAC(I,2) = TMP1*TMP2/TMP4 */
/*          FJAC(I,3) = TMP1*TMP3/TMP4 */
/*    30    CONTINUE */
/*       RETURN */
/*       END */


/*       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3), */
/*         LDFJAC would be set to 3, and FCN would be written as */
/*         follows to calculate a row of the Jacobian when IFLAG=3. */

/*       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC) */
/* C     This is the form of the FCN routine if IOPT=3, */
/* C     that is, if the user calculates the Jacobian row by row. */
/*       INTEGER I,M,N,IFLAG */
/*       DOUBLE PRECISION X(N),FVEC(M),FJAC(N),Y(15) */
/*       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4 */
/*       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8), */
/*      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15) */
/*      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1, */
/*      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/ */
/* C */
/*       IF (IFLAG .NE. 0) GO TO 5 */
/* C */
/* C     Insert print statements here when NPRINT is positive. */
/* C */
/*       RETURN */
/*     5 CONTINUE */
/*       IF( IFLAG.NE.1) GO TO 20 */
/*       DO 10 I = 1, M */
/*          TMP1 = I */
/*          TMP2 = 16 - I */
/*          TMP3 = TMP1 */
/*          IF (I .GT. 8) TMP3 = TMP2 */
/*          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3)) */
/*    10    CONTINUE */
/*       RETURN */
/* C */
/* C     Below, calculate the LDFJAC-th row of the Jacobian. */
/* C */
/*    20 CONTINUE */

/*       I = LDFJAC */
/*          TMP1 = I */
/*          TMP2 = 16 - I */
/*          TMP3 = TMP1 */
/*          IF (I .GT. 8) TMP3 = TMP2 */
/*          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2 */
/*          FJAC(1) = -1.E0 */
/*          FJAC(2) = TMP1*TMP2/TMP4 */
/*          FJAC(3) = TMP1*TMP3/TMP4 */
/*       RETURN */
/*       END */

/* ***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm: */
/*                 implementation and theory.  In Numerical Analysis */
/*                 Proceedings (Dundee, June 28 - July 1, 1977, G. A. */
/*                 Watson, Editor), Lecture Notes in Mathematics 630, */
/*                 Springer-Verlag, 1978. */
/* ***ROUTINES CALLED  D1MACH, DCKDER, DENORM, DFDJC3, DMPAR, DQRFAC, */
/*                    DWUPDT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891006  Cosmetic changes to prologue.  (WRB) */
/*   891006  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900510  Convert XERRWV calls to XERMSG calls.  (RWC) */
/*   920205  Corrected XERN1 declaration.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DNLS1 */
/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV >*/
/*<       INTEGER IJUNK,NROW,IPVT(*) >*/
/*<       DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,EPSFCN >*/
/*<    >*/
/*<       LOGICAL SING >*/
/*<       EXTERNAL FCN >*/
/*<       INTEGER I,IFLAG,ITER,J,L,MODECH >*/
/*<    >*/
/*<       DOUBLE PRECISION D1MACH,DENORM,ERR,CHKLIM >*/
/*<       CHARACTER*8 XERN1 >*/
/*<       CHARACTER*16 XERN3 >*/
/*<       SAVE CHKLIM, ONE, P1, P5, P25, P75, P0001, ZERO >*/

/*<       DATA CHKLIM/.1D0/ >*/
    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --ipvt;
    --diag;
    fjac_dim1 = *ldfjac;
    fjac_offset = fjac_dim1 + 1;
    fjac -= fjac_offset;
    --fvec;
    --x;

    /* Function Body */
/*<    >*/
/* ***FIRST EXECUTABLE STATEMENT  DNLS1 */
/*<       EPSMCH = D1MACH(4) >*/
    epsmch = d1mach_(&c__4);

/*<       INFO = 0 >*/
    *info = 0;
/*<       IFLAG = 0 >*/
    iflag = 0;
/*<       NFEV = 0 >*/
    *nfev = 0;
/*<       NJEV = 0 >*/
    *njev = 0;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

/*<    >*/
    if (*iopt < 1 || *iopt > 3 || *n <= 0 || *m < *n || *ldfjac < *n || *ftol 
	    < zero || *xtol < zero || *gtol < zero || *maxfev <= 0 || *factor 
	    <= zero) {
	goto L300;
    }
/*<       IF (IOPT .LT. 3 .AND. LDFJAC .LT. M) GO TO 300 >*/
    if (*iopt < 3 && *ldfjac < *m) {
	goto L300;
    }
/*<       IF (MODE .NE. 2) GO TO 20 >*/
    if (*mode != 2) {
	goto L20;
    }
/*<       DO 10 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          IF (DIAG(J) .LE. ZERO) GO TO 300 >*/
	if (diag[j] <= zero) {
	    goto L300;
	}
/*<    10    CONTINUE >*/
/* L10: */
    }
/*<    20 CONTINUE >*/
L20:

/*     EVALUATE THE FUNCTION AT THE STARTING POINT */
/*     AND CALCULATE ITS NORM. */

/*<       IFLAG = 1 >*/
    iflag = 1;
/*<       IJUNK = 1 >*/
    ijunk = 1;
/*<       CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK) >*/
    (*fcn)(&iflag, m, n, &x[1], &fvec[1], &fjac[fjac_offset], &ijunk, pParam);
/*<       NFEV = 1 >*/
    *nfev = 1;
/*<       IF (IFLAG .LT. 0) GO TO 300 >*/
    if (iflag < 0) {
	goto L300;
    }
/*<       FNORM = DENORM(M,FVEC) >*/
    fnorm = denorm_(m, &fvec[1]);

/*     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER. */

/*<       PAR = ZERO >*/
    par = zero;
/*<       ITER = 1 >*/
    iter = 1;

/*     BEGINNING OF THE OUTER LOOP. */

/*<    30 CONTINUE >*/
L30:

/*        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES. */

/*<          IF (NPRINT .LE. 0) GO TO 40 >*/
    if (*nprint <= 0) {
	goto L40;
    }
/*<          IFLAG = 0 >*/
    iflag = 0;
/*<    >*/
    if ((iter - 1) % *nprint == 0) {
	(*fcn)(&iflag, m, n, &x[1], &fvec[1], &fjac[fjac_offset], &ijunk, pParam);
    }
/*<          IF (IFLAG .LT. 0) GO TO 300 >*/
    if (iflag < 0) {
	goto L300;
    }
/*<    40    CONTINUE >*/
L40:

/*        CALCULATE THE JACOBIAN MATRIX. */

/*<       IF (IOPT .EQ. 3) GO TO 475 >*/
    if (*iopt == 3) {
	goto L475;
    }

/*     STORE THE FULL JACOBIAN USING M*N STORAGE */

/*<       IF (IOPT .EQ. 1) GO TO 410 >*/
    if (*iopt == 1) {
	goto L410;
    }

/*     THE USER SUPPLIES THE JACOBIAN */

/*<          IFLAG = 2 >*/
    iflag = 2;
/*<          CALL FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC) >*/
    (*fcn)(&iflag, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, pParam);
/*<          NJEV = NJEV + 1 >*/
    ++(*njev);

/*             ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN */

/*<          IF (ITER .LE. 1) THEN >*/
    if (iter <= 1) {
/*<             IF (IFLAG .LT. 0) GO TO 300 >*/
	if (iflag < 0) {
	    goto L300;
	}

/*           GET THE INCREMENTED X-VALUES INTO WA1(*). */

/*<             MODECH = 1 >*/
	modech = 1;
/*<             CALL DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR) >*/
	dckder_(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &wa1[1], &
		wa4[1], &modech, &err);

/*           EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT IN WA4(*).
 */

/*<             IFLAG = 1 >*/
	iflag = 1;
/*<             CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,LDFJAC) >*/
	(*fcn)(&iflag, m, n, &wa1[1], &wa4[1], &fjac[fjac_offset], ldfjac, pParam);
/*<             NFEV = NFEV + 1 >*/
	++(*nfev);
/*<             IF(IFLAG .LT. 0) GO TO 300 >*/
	if (iflag < 0) {
	    goto L300;
	}
/*<             DO 350 I = 1, M >*/
	i__1 = *m;
	for (i = 1; i <= i__1; ++i) {
/*<                MODECH = 2 >*/
	    modech = 2;
/*<    >*/
	    dckder_(&c__1, n, &x[1], &fvec[i], &fjac[i + fjac_dim1], ldfjac, &
		    wa1[1], &wa4[i], &modech, &err);
/*               IF (ERR .LT. CHKLIM) THEN */
/*                  WRITE (XERN1, '(I8)') I */
/*                  WRITE (XERN3, '(1PE15.6)') ERR */
/*                  CALL XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF
 ' // */
/*     *               'FUNCTION ' // XERN1 // ' MAY BE WRONG, ERR
 = ' // */
/*     *               XERN3 // ' TOO CLOSE TO 0.', 7, 0) */
/*               ENDIF */
/*<   350       CONTINUE >*/
/* L350: */
	}
/*<          ENDIF >*/
    }

/*<          GO TO 420 >*/
    goto L420;

/*     THE CODE APPROXIMATES THE JACOBIAN */

/*< 410      IFLAG = 1 >*/
L410:
    iflag = 1;
/*<          CALL DFDJC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4) >*/
    dfdjc3_((S_fp)fcn, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &
	    iflag, epsfcn, &wa4[1], pParam);
/*<          NFEV = NFEV + N >*/
    *nfev += *n;
/*<   420    IF (IFLAG .LT. 0) GO TO 300 >*/
L420:
    if (iflag < 0) {
	goto L300;
    }

/*        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN. */

/*<          CALL DQRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3) >*/
    dqrfac_(m, n, &fjac[fjac_offset], ldfjac, (logical*)&c__1, &ipvt[1], n, &
	    wa1[1], &wa2[1], &wa3[1]);

/*        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN */
/*        QTF. */

/*<          DO 430 I = 1, M >*/
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/*<             WA4(I) = FVEC(I) >*/
	wa4[i] = fvec[i];
/*<   430         CONTINUE >*/
/* L430: */
    }
/*<          DO 470 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             IF (FJAC(J,J) .EQ. ZERO) GO TO 460 >*/
	if (fjac[j + j * fjac_dim1] == zero) {
	    goto L460;
	}
/*<             SUM = ZERO >*/
	sum = zero;
/*<             DO 440 I = J, M >*/
	i__2 = *m;
	for (i = j; i <= i__2; ++i) {
/*<                SUM = SUM + FJAC(I,J)*WA4(I) >*/
	    sum += fjac[i + j * fjac_dim1] * wa4[i];
/*<   440          CONTINUE >*/
/* L440: */
	}
/*<             TEMP = -SUM/FJAC(J,J) >*/
	temp = -sum / fjac[j + j * fjac_dim1];
/*<             DO 450 I = J, M >*/
	i__2 = *m;
	for (i = j; i <= i__2; ++i) {
/*<                WA4(I) = WA4(I) + FJAC(I,J)*TEMP >*/
	    wa4[i] += fjac[i + j * fjac_dim1] * temp;
/*<   450          CONTINUE >*/
/* L450: */
	}
/*<   460       CONTINUE >*/
L460:
/*<             FJAC(J,J) = WA1(J) >*/
	fjac[j + j * fjac_dim1] = wa1[j];
/*<             QTF(J) = WA4(J) >*/
	qtf[j] = wa4[j];
/*<   470       CONTINUE >*/
/* L470: */
    }
/*<          GO TO 560 >*/
    goto L560;

/*        ACCUMULATE THE JACOBIAN BY ROWS IN ORDER TO SAVE STORAGE. */
/*        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX */
/*        CALCULATED ONE ROW AT A TIME, WHILE SIMULTANEOUSLY */
/*        FORMING (Q TRANSPOSE)*FVEC AND STORING THE FIRST */
/*        N COMPONENTS IN QTF. */

/*<   475    DO 490 J = 1, N >*/
L475:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             QTF(J) = ZERO >*/
	qtf[j] = zero;
/*<             DO 480 I = 1, N >*/
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
/*<                FJAC(I,J) = ZERO >*/
	    fjac[i + j * fjac_dim1] = zero;
/*<   480          CONTINUE >*/
/* L480: */
	}
/*<   490        CONTINUE >*/
/* L490: */
    }
/*<          DO 500 I = 1, M >*/
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/*<             NROW = I >*/
	nrow = i;
/*<             IFLAG = 3 >*/
	iflag = 3;
/*<             CALL FCN(IFLAG,M,N,X,FVEC,WA3,NROW) >*/
	(*fcn)(&iflag, m, n, &x[1], &fvec[1], &wa3[1], &nrow, pParam);
/*<             IF (IFLAG .LT. 0) GO TO 300 >*/
	if (iflag < 0) {
	    goto L300;
	}

/*            ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN
. */

/*<             IF(ITER .GT. 1) GO TO 498 >*/
	if (iter > 1) {
	    goto L498;
	}

/*            GET THE INCREMENTED X-VALUES INTO WA1(*). */

/*<             MODECH = 1 >*/
	modech = 1;
/*<             CALL DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR) >*/
	dckder_(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &wa1[1], &
		wa4[1], &modech, &err);

/*            EVALUATE AT INCREMENTED VALUES, IF NOT ALREADY EVALUATED
. */

/*<             IF(I .NE. 1) GO TO 495 >*/
	if (i != 1) {
	    goto L495;
	}

/*            EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT INTO WA4(
*). */

/*<             IFLAG = 1 >*/
	iflag = 1;
/*<             CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,NROW) >*/
	(*fcn)(&iflag, m, n, &wa1[1], &wa4[1], &fjac[fjac_offset], &nrow, pParam);
/*<             NFEV = NFEV + 1 >*/
	++(*nfev);
/*<             IF(IFLAG .LT. 0) GO TO 300 >*/
	if (iflag < 0) {
	    goto L300;
	}
/*< 495         CONTINUE >*/
L495:
/*<             MODECH = 2 >*/
	modech = 2;
/*<             CALL DCKDER(1,N,X,FVEC(I),WA3,1,WA1,WA4(I),MODECH,ERR) >*/
	dckder_(&c__1, n, &x[1], &fvec[i], &wa3[1], &c__1, &wa1[1], &wa4[i], &
		modech, &err);
/*            IF (ERR .LT. CHKLIM) THEN */
/*               WRITE (XERN1, '(I8)') I */
/*               WRITE (XERN3, '(1PE15.6)') ERR */
/*               CALL XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF FUNCTI
ON ' */
/*     *            // XERN1 // ' MAY BE WRONG, ERR = ' // XERN3 // */
/*     *            ' TOO CLOSE TO 0.', 7, 0) */
/*            ENDIF */
/*< 498         CONTINUE >*/
L498:

/*<             TEMP = FVEC(I) >*/
	temp = fvec[i];
/*<             CALL DWUPDT(N,FJAC,LDFJAC,WA3,QTF,TEMP,WA1,WA2) >*/
	dwupdt_(n, &fjac[fjac_offset], ldfjac, &wa3[1], &qtf[1], &temp, &wa1[
		1], &wa2[1]);
/*<   500       CONTINUE >*/
/* L500: */
    }
/*<          NJEV = NJEV + 1 >*/
    ++(*njev);

/*        IF THE JACOBIAN IS RANK DEFICIENT, CALL DQRFAC TO */
/*        REORDER ITS COLUMNS AND UPDATE THE COMPONENTS OF QTF. */

/*<          SING = .FALSE. >*/
    sing = FALSE_;
/*<          DO 510 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             IF (FJAC(J,J) .EQ. ZERO) SING = .TRUE. >*/
	if (fjac[j + j * fjac_dim1] == zero) {
	    sing = TRUE_;
	}
/*<             IPVT(J) = J >*/
	ipvt[j] = j;
/*<             WA2(J) = DENORM(J,FJAC(1,J)) >*/
	wa2[j] = denorm_(&j, &fjac[j * fjac_dim1 + 1]);
/*<   510       CONTINUE >*/
/* L510: */
    }
/*<          IF (.NOT.SING) GO TO 560 >*/
    if (! sing) {
	goto L560;
    }
/*<          CALL DQRFAC(N,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3) >*/
    dqrfac_(n, n, &fjac[fjac_offset], ldfjac, (logical*)&c__1, &ipvt[1], n, &
	    wa1[1], &wa2[1], &wa3[1]);
/*<          DO 550 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             IF (FJAC(J,J) .EQ. ZERO) GO TO 540 >*/
	if (fjac[j + j * fjac_dim1] == zero) {
	    goto L540;
	}
/*<             SUM = ZERO >*/
	sum = zero;
/*<             DO 520 I = J, N >*/
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
/*<                SUM = SUM + FJAC(I,J)*QTF(I) >*/
	    sum += fjac[i + j * fjac_dim1] * qtf[i];
/*<   520         CONTINUE >*/
/* L520: */
	}
/*<             TEMP = -SUM/FJAC(J,J) >*/
	temp = -sum / fjac[j + j * fjac_dim1];
/*<             DO 530 I = J, N >*/
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
/*<                QTF(I) = QTF(I) + FJAC(I,J)*TEMP >*/
	    qtf[i] += fjac[i + j * fjac_dim1] * temp;
/*<   530          CONTINUE >*/
/* L530: */
	}
/*<   540       CONTINUE >*/
L540:
/*<             FJAC(J,J) = WA1(J) >*/
	fjac[j + j * fjac_dim1] = wa1[j];
/*<   550       CONTINUE >*/
/* L550: */
    }
/*<   560    CONTINUE >*/
L560:

/*        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING */
/*        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN. */

/*<          IF (ITER .NE. 1) GO TO 80 >*/
    if (iter != 1) {
	goto L80;
    }
/*<          IF (MODE .EQ. 2) GO TO 60 >*/
    if (*mode == 2) {
	goto L60;
    }
/*<          DO 50 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             DIAG(J) = WA2(J) >*/
	diag[j] = wa2[j];
/*<             IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE >*/
	if (wa2[j] == zero) {
	    diag[j] = one;
	}
/*<    50       CONTINUE >*/
/* L50: */
    }
/*<    60    CONTINUE >*/
L60:

/*        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X */
/*        AND INITIALIZE THE STEP BOUND DELTA. */

/*<          DO 70 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             WA3(J) = DIAG(J)*X(J) >*/
	wa3[j] = diag[j] * x[j];
/*<    70       CONTINUE >*/
/* L70: */
    }
/*<          XNORM = DENORM(N,WA3) >*/
    xnorm = denorm_(n, &wa3[1]);
/*<          DELTA = FACTOR*XNORM >*/
    delta = *factor * xnorm;
/*<          IF (DELTA .EQ. ZERO) DELTA = FACTOR >*/
    if (delta == zero) {
	delta = *factor;
    }
/*<    80    CONTINUE >*/
L80:

/*        COMPUTE THE NORM OF THE SCALED GRADIENT. */

/*<          GNORM = ZERO >*/
    gnorm = zero;
/*<          IF (FNORM .EQ. ZERO) GO TO 170 >*/
    if (fnorm == zero) {
	goto L170;
    }
/*<          DO 160 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             L = IPVT(J) >*/
	l = ipvt[j];
/*<             IF (WA2(L) .EQ. ZERO) GO TO 150 >*/
	if (wa2[l] == zero) {
	    goto L150;
	}
/*<             SUM = ZERO >*/
	sum = zero;
/*<             DO 140 I = 1, J >*/
	i__2 = j;
	for (i = 1; i <= i__2; ++i) {
/*<                SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM) >*/
	    sum += fjac[i + j * fjac_dim1] * (qtf[i] / fnorm);
/*<   140          CONTINUE >*/
/* L140: */
	}
/*<             GNORM = MAX(GNORM,ABS(SUM/WA2(L))) >*/
/* Computing MAX */
	d__2 = gnorm, d__3 = (d__1 = sum / wa2[l], abs(d__1));
	gnorm = max(d__2,d__3);
/*<   150       CONTINUE >*/
L150:
/*<   160       CONTINUE >*/
/* L160: */
	;
    }
/*<   170    CONTINUE >*/
L170:

/*        TEST FOR CONVERGENCE OF THE GRADIENT NORM. */

/*<          IF (GNORM .LE. GTOL) INFO = 4 >*/
    if (gnorm <= *gtol) {
	*info = 4;
    }
/*<          IF (INFO .NE. 0) GO TO 300 >*/
    if (*info != 0) {
	goto L300;
    }

/*        RESCALE IF NECESSARY. */

/*<          IF (MODE .EQ. 2) GO TO 190 >*/
    if (*mode == 2) {
	goto L190;
    }
/*<          DO 180 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<             DIAG(J) = MAX(DIAG(J),WA2(J)) >*/
/* Computing MAX */
	d__1 = diag[j], d__2 = wa2[j];
	diag[j] = max(d__1,d__2);
/*<   180       CONTINUE >*/
/* L180: */
    }
/*<   190    CONTINUE >*/
L190:

/*        BEGINNING OF THE INNER LOOP. */

/*<   200    CONTINUE >*/
L200:

/*           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER. */

/*<    >*/
    dmpar_(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], &delta,
	     &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

/*           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P. */

/*<             DO 210 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<                WA1(J) = -WA1(J) >*/
	wa1[j] = -wa1[j];
/*<                WA2(J) = X(J) + WA1(J) >*/
	wa2[j] = x[j] + wa1[j];
/*<                WA3(J) = DIAG(J)*WA1(J) >*/
	wa3[j] = diag[j] * wa1[j];
/*<   210          CONTINUE >*/
/* L210: */
    }
/*<             PNORM = DENORM(N,WA3) >*/
    pnorm = denorm_(n, &wa3[1]);

/*           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND. */

/*<             IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM) >*/
    if (iter == 1) {
	delta = min(delta,pnorm);
    }

/*           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM. */

/*<             IFLAG = 1 >*/
    iflag = 1;
/*<             CALL FCN(IFLAG,M,N,WA2,WA4,FJAC,IJUNK) >*/
    (*fcn)(&iflag, m, n, &wa2[1], &wa4[1], &fjac[fjac_offset], &ijunk, pParam);
/*<             NFEV = NFEV + 1 >*/
    ++(*nfev);
/*<             IF (IFLAG .LT. 0) GO TO 300 >*/
    if (iflag < 0) {
	goto L300;
    }
/*<             FNORM1 = DENORM(M,WA4) >*/
    fnorm1 = denorm_(m, &wa4[1]);

/*           COMPUTE THE SCALED ACTUAL REDUCTION. */

/*<             ACTRED = -ONE >*/
    actred = -one;
/*<             IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2 >*/
    if (p1 * fnorm1 < fnorm) {
/* Computing 2nd power */
	d__1 = fnorm1 / fnorm;
	actred = one - d__1 * d__1;
    }

/*           COMPUTE THE SCALED PREDICTED REDUCTION AND */
/*           THE SCALED DIRECTIONAL DERIVATIVE. */

/*<             DO 230 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<                WA3(J) = ZERO >*/
	wa3[j] = zero;
/*<                L = IPVT(J) >*/
	l = ipvt[j];
/*<                TEMP = WA1(L) >*/
	temp = wa1[l];
/*<                DO 220 I = 1, J >*/
	i__2 = j;
	for (i = 1; i <= i__2; ++i) {
/*<                   WA3(I) = WA3(I) + FJAC(I,J)*TEMP >*/
	    wa3[i] += fjac[i + j * fjac_dim1] * temp;
/*<   220             CONTINUE >*/
/* L220: */
	}
/*<   230          CONTINUE >*/
/* L230: */
    }
/*<             TEMP1 = DENORM(N,WA3)/FNORM >*/
    temp1 = denorm_(n, &wa3[1]) / fnorm;
/*<             TEMP2 = (SQRT(PAR)*PNORM)/FNORM >*/
    temp2 = sqrt(par) * pnorm / fnorm;
/*<             PRERED = TEMP1**2 + TEMP2**2/P5 >*/
/* Computing 2nd power */
    d__1 = temp1;
/* Computing 2nd power */
    d__2 = temp2;
    prered = d__1 * d__1 + d__2 * d__2 / p5;
/*<             DIRDER = -(TEMP1**2 + TEMP2**2) >*/
/* Computing 2nd power */
    d__1 = temp1;
/* Computing 2nd power */
    d__2 = temp2;
    dirder = -(d__1 * d__1 + d__2 * d__2);

/*           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED */
/*           REDUCTION. */

/*<             RATIO = ZERO >*/
    ratio = zero;
/*<             IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED >*/
    if (prered != zero) {
	ratio = actred / prered;
    }

/*           UPDATE THE STEP BOUND. */

/*<             IF (RATIO .GT. P25) GO TO 240 >*/
    if (ratio > p25) {
	goto L240;
    }
/*<                IF (ACTRED .GE. ZERO) TEMP = P5 >*/
    if (actred >= zero) {
	temp = p5;
    }
/*<    >*/
    if (actred < zero) {
	temp = p5 * dirder / (dirder + p5 * actred);
    }
/*<                IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1 >*/
    if (p1 * fnorm1 >= fnorm || temp < p1) {
	temp = p1;
    }
/*<                DELTA = TEMP*MIN(DELTA,PNORM/P1) >*/
/* Computing MIN */
    d__1 = delta, d__2 = pnorm / p1;
    delta = temp * min(d__1,d__2);
/*<                PAR = PAR/TEMP >*/
    par /= temp;
/*<                GO TO 260 >*/
    goto L260;
/*<   240       CONTINUE >*/
L240:
/*<                IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250 >*/
    if (par != zero && ratio < p75) {
	goto L250;
    }
/*<                DELTA = PNORM/P5 >*/
    delta = pnorm / p5;
/*<                PAR = P5*PAR >*/
    par = p5 * par;
/*<   250          CONTINUE >*/
L250:
/*<   260       CONTINUE >*/
L260:

/*           TEST FOR SUCCESSFUL ITERATION. */

/*<             IF (RATIO .LT. P0001) GO TO 290 >*/
    if (ratio < p0001) {
	goto L290;
    }

/*           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS. */

/*<             DO 270 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<                X(J) = WA2(J) >*/
	x[j] = wa2[j];
/*<                WA2(J) = DIAG(J)*X(J) >*/
	wa2[j] = diag[j] * x[j];
/*<   270          CONTINUE >*/
/* L270: */
    }
/*<             DO 280 I = 1, M >*/
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/*<                FVEC(I) = WA4(I) >*/
	fvec[i] = wa4[i];
/*<   280          CONTINUE >*/
/* L280: */
    }
/*<             XNORM = DENORM(N,WA2) >*/
    xnorm = denorm_(n, &wa2[1]);
/*<             FNORM = FNORM1 >*/
    fnorm = fnorm1;
/*<             ITER = ITER + 1 >*/
    ++iter;
/*<   290       CONTINUE >*/
L290:

/*           TESTS FOR CONVERGENCE. */

/*<    >*/
    if (abs(actred) <= *ftol && prered <= *ftol && p5 * ratio <= one) {
	*info = 1;
    }
/*<             IF (DELTA .LE. XTOL*XNORM) INFO = 2 >*/
    if (delta <= *xtol * xnorm) {
      *info = 2;
    }
/*<    >*/
    if (abs(actred) <= *ftol && prered <= *ftol && p5 * ratio <= one && *info 
	    == 2) {
	*info = 3;
    }
/*<             IF (INFO .NE. 0) GO TO 300 >*/
    if (*info != 0) {
      goto L300;
    }

/*           TESTS FOR TERMINATION AND STRINGENT TOLERANCES. */

/*<             IF (NFEV .GE. MAXFEV) INFO = 5 >*/
    if (*nfev >= *maxfev) {
	*info = 5;
    }
/*<    >*/
    if (abs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= one) {
	*info = 6;
    }
/*<             IF (DELTA .LE. EPSMCH*XNORM) INFO = 7 >*/
    if (delta <= epsmch * xnorm) {
	*info = 7;
    }
/*<             IF (GNORM .LE. EPSMCH) INFO = 8 >*/
    if (gnorm <= epsmch) {
	*info = 8;
    }
/*<             IF (INFO .NE. 0) GO TO 300 >*/
    if (*info != 0) {
      goto L300;
    }

/*           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL. */

/*<             IF (RATIO .LT. P0001) GO TO 200 >*/
    if (ratio < p0001) {
	goto L200;
    }

/*        END OF THE OUTER LOOP. */

/*<          GO TO 30 >*/
    goto L30;
/*<   300 CONTINUE >*/
L300:

/*     TERMINATION, EITHER NORMAL OR USER IMPOSED. */

/*<       IF (IFLAG .LT. 0) INFO = IFLAG >*/
    if (iflag < 0) {
	*info = iflag;
    }
/*<       IFLAG = 0 >*/
    iflag = 0;
/*<       IF (NPRINT .GT. 0) CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK) >*/
    if (*nprint > 0) {
	(*fcn)(&iflag, m, n, &x[1], &fvec[1], &fjac[fjac_offset], &ijunk, pParam);
    }
/*<    >*/
    if (*info < 0) {
	xermsg_("SLATEC", "DNLS1", "EXECUTION TERMINATED BECAUSE USER SET IF"
		"LAG NEGATIVE.", &c__1, &c__1, 6L, 5L, 53L);
    }
/*<    >*/
    if (*info == 0) {
	xermsg_("SLATEC", "DNLS1", "INVALID INPUT PARAMETER.", &c__2, &c__1, 
		6L, 5L, 24L);
    }
/*<    >*/
    if (*info == 4) {
	xermsg_("SLATEC", "DNLS1", "THIRD CONVERGENCE CONDITION, CHECK RESUL"
		"TS BEFORE ACCEPTING.", &c__1, &c__1, 6L, 5L, 60L);
    }
/*<    >*/
    if (*info == 5) {
	xermsg_("SLATEC", "DNLS1", "TOO MANY FUNCTION EVALUATIONS.", &c__9, &
		c__1, 6L, 5L, 30L);
    }
/*<    >*/
    if (*info >= 6) {
	xermsg_("SLATEC", "DNLS1", "TOLERANCES TOO SMALL, NO FURTHER IMPROVE"
		"MENT POSSIBLE.", &c__3, &c__1, 6L, 5L, 54L);
    }
/*<       RETURN >*/
    return 0;

/*     LAST CARD OF SUBROUTINE DNLS1. */

/*<       END >*/
} /* dnls1_ */

/*<    >*/
/* Subroutine */ int dqrfac_(integer *m, integer *n, doublereal *a, integer *
	lda, logical *pivot, integer *ipvt, integer *lipvt, doublereal *sigma,
	 doublereal *acnorm, doublereal *wa)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal p05 = .05;
    static doublereal zero = 0.;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    integer kmax;
    doublereal temp;
    integer i, j, k, minmn;
    extern doublereal d1mach_(integer *);
    doublereal epsmch;
    extern doublereal denorm_(integer *, doublereal *);
    doublereal ajnorm;
    integer jp1;
    doublereal sum;

/* ***BEGIN PROLOGUE  DQRFAC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNLS1, DNLS1E, DNSQ and DNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (QRFAC-S, DQRFAC-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   **** Double Precision version of QRFAC **** */

/*     This subroutine uses Householder transformations with column */
/*     pivoting (optional) to compute a QR factorization of the */
/*     M by N matrix A. That is, DQRFAC determines an orthogonal */
/*     matrix Q, a permutation matrix P, and an upper trapezoidal */
/*     matrix R with diagonal elements of nonincreasing magnitude, */
/*     such that A*P = Q*R. The Householder transformation for */
/*     column K, K = 1,2,...,MIN(M,N), is of the form */

/*                           T */
/*           I - (1/U(K))*U*U */

/*     where U has zeros in the first K-1 positions. The form of */
/*     this transformation and the method of pivoting first */
/*     appeared in the corresponding LINPACK subroutine. */

/*     The subroutine statement is */

/*       SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA) */

/*     where */

/*       M is a positive integer input variable set to the number */
/*         of rows of A. */

/*       N is a positive integer input variable set to the number */
/*         of columns of A. */

/*       A is an M by N array. On input A contains the matrix for */
/*         which the QR factorization is to be computed. On output */
/*         the strict upper trapezoidal part of A contains the strict */
/*         upper trapezoidal part of R, and the lower trapezoidal */
/*         part of A contains a factored form of Q (the non-trivial */
/*         elements of the U vectors described above). */

/*       LDA is a positive integer input variable not less than M */
/*         which specifies the leading dimension of the array A. */

/*       PIVOT is a logical input variable. If pivot is set .TRUE., */
/*         then column pivoting is enforced. If pivot is set .FALSE., */
/*         then no column pivoting is done. */

/*       IPVT is an integer output array of length LIPVT. IPVT */
/*         defines the permutation matrix P such that A*P = Q*R. */
/*         Column J of P is column IPVT(J) of the identity matrix. */
/*         If pivot is .FALSE., IPVT is not referenced. */

/*       LIPVT is a positive integer input variable. If PIVOT is */
/*             .FALSE., then LIPVT may be as small as 1. If PIVOT is */
/*             .TRUE., then LIPVT must be at least N. */

/*       SIGMA is an output array of length N which contains the */
/*         diagonal elements of R. */

/*       ACNORM is an output array of length N which contains the */
/*         norms of the corresponding columns of the input matrix A. */
/*         If this information is not needed, then ACNORM can coincide */
/*         with SIGMA. */

/*       WA is a work array of length N. If pivot is .FALSE., then WA */
/*         can coincide with SIGMA. */

/* ***SEE ALSO  DNLS1, DNLS1E, DNSQ, DNSQE */
/* ***ROUTINES CALLED  D1MACH, DENORM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DQRFAC */
/*<       INTEGER M,N,LDA,LIPVT >*/
/*<       INTEGER IPVT(*) >*/
/*<       LOGICAL PIVOT >*/
/*<       SAVE ONE, P05, ZERO >*/
/*<       DOUBLE PRECISION A(LDA,*),SIGMA(*),ACNORM(*),WA(*) >*/
/*<       INTEGER I,J,JP1,K,KMAX,MINMN >*/
/*<       DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO >*/
/*<       DOUBLE PRECISION D1MACH,DENORM >*/
/*<       DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/ >*/
    /* Parameter adjustments */
    --wa;
    --acnorm;
    --sigma;
    --ipvt;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DQRFAC */
/*<       EPSMCH = D1MACH(4) >*/
    epsmch = d1mach_(&c__4);

/*     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS. */

/*<       DO 10 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          ACNORM(J) = DENORM(M,A(1,J)) >*/
	acnorm[j] = denorm_(m, &a[j * a_dim1 + 1]);
/*<          SIGMA(J) = ACNORM(J) >*/
	sigma[j] = acnorm[j];
/*<          WA(J) = SIGMA(J) >*/
	wa[j] = sigma[j];
/*<          IF (PIVOT) IPVT(J) = J >*/
	if (*pivot) {
	    ipvt[j] = j;
	}
/*<    10    CONTINUE >*/
/* L10: */
    }

/*     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS. */

/*<       MINMN = MIN(M,N) >*/
    minmn = min(*m,*n);
/*<       DO 110 J = 1, MINMN >*/
    i__1 = minmn;
    for (j = 1; j <= i__1; ++j) {
/*<          IF (.NOT.PIVOT) GO TO 40 >*/
	if (! (*pivot)) {
	    goto L40;
	}

/*        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION. */

/*<          KMAX = J >*/
	kmax = j;
/*<          DO 20 K = J, N >*/
	i__2 = *n;
	for (k = j; k <= i__2; ++k) {
/*<             IF (SIGMA(K) .GT. SIGMA(KMAX)) KMAX = K >*/
	    if (sigma[k] > sigma[kmax]) {
		kmax = k;
	    }
/*<    20       CONTINUE >*/
/* L20: */
	}
/*<          IF (KMAX .EQ. J) GO TO 40 >*/
	if (kmax == j) {
	    goto L40;
	}
/*<          DO 30 I = 1, M >*/
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
/*<             TEMP = A(I,J) >*/
	    temp = a[i + j * a_dim1];
/*<             A(I,J) = A(I,KMAX) >*/
	    a[i + j * a_dim1] = a[i + kmax * a_dim1];
/*<             A(I,KMAX) = TEMP >*/
	    a[i + kmax * a_dim1] = temp;
/*<    30       CONTINUE >*/
/* L30: */
	}
/*<          SIGMA(KMAX) = SIGMA(J) >*/
	sigma[kmax] = sigma[j];
/*<          WA(KMAX) = WA(J) >*/
	wa[kmax] = wa[j];
/*<          K = IPVT(J) >*/
	k = ipvt[j];
/*<          IPVT(J) = IPVT(KMAX) >*/
	ipvt[j] = ipvt[kmax];
/*<          IPVT(KMAX) = K >*/
	ipvt[kmax] = k;
/*<    40    CONTINUE >*/
L40:

/*        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE */
/*        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR. */

/*<          AJNORM = DENORM(M-J+1,A(J,J)) >*/
	i__2 = *m - j + 1;
	ajnorm = denorm_(&i__2, &a[j + j * a_dim1]);
/*<          IF (AJNORM .EQ. ZERO) GO TO 100 >*/
	if (ajnorm == zero) {
	    goto L100;
	}
/*<          IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM >*/
	if (a[j + j * a_dim1] < zero) {
	    ajnorm = -ajnorm;
	}
/*<          DO 50 I = J, M >*/
	i__2 = *m;
	for (i = j; i <= i__2; ++i) {
/*<             A(I,J) = A(I,J)/AJNORM >*/
	    a[i + j * a_dim1] /= ajnorm;
/*<    50       CONTINUE >*/
/* L50: */
	}
/*<          A(J,J) = A(J,J) + ONE >*/
	a[j + j * a_dim1] += one;

/*        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS */
/*        AND UPDATE THE NORMS. */

/*<          JP1 = J + 1 >*/
	jp1 = j + 1;
/*<          IF (N .LT. JP1) GO TO 100 >*/
	if (*n < jp1) {
	    goto L100;
	}
/*<          DO 90 K = JP1, N >*/
	i__2 = *n;
	for (k = jp1; k <= i__2; ++k) {
/*<             SUM = ZERO >*/
	    sum = zero;
/*<             DO 60 I = J, M >*/
	    i__3 = *m;
	    for (i = j; i <= i__3; ++i) {
/*<                SUM = SUM + A(I,J)*A(I,K) >*/
		sum += a[i + j * a_dim1] * a[i + k * a_dim1];
/*<    60          CONTINUE >*/
/* L60: */
	    }
/*<             TEMP = SUM/A(J,J) >*/
	    temp = sum / a[j + j * a_dim1];
/*<             DO 70 I = J, M >*/
	    i__3 = *m;
	    for (i = j; i <= i__3; ++i) {
/*<                A(I,K) = A(I,K) - TEMP*A(I,J) >*/
		a[i + k * a_dim1] -= temp * a[i + j * a_dim1];
/*<    70          CONTINUE >*/
/* L70: */
	    }
/*<             IF (.NOT.PIVOT .OR. SIGMA(K) .EQ. ZERO) GO TO 80 >*/
	    if (! (*pivot) || sigma[k] == zero) {
		goto L80;
	    }
/*<             TEMP = A(J,K)/SIGMA(K) >*/
	    temp = a[j + k * a_dim1] / sigma[k];
/*<             SIGMA(K) = SIGMA(K)*SQRT(MAX(ZERO,ONE-TEMP**2)) >*/
/* Computing MAX */
/* Computing 2nd power */
	    d__3 = temp;
	    d__1 = zero, d__2 = one - d__3 * d__3;
	    sigma[k] *= sqrt((max(d__1,d__2)));
/*<             IF (P05*(SIGMA(K)/WA(K))**2 .GT. EPSMCH) GO TO 80 >*/
/* Computing 2nd power */
	    d__1 = sigma[k] / wa[k];
	    if (p05 * (d__1 * d__1) > epsmch) {
		goto L80;
	    }
/*<             SIGMA(K) = DENORM(M-J,A(JP1,K)) >*/
	    i__3 = *m - j;
	    sigma[k] = denorm_(&i__3, &a[jp1 + k * a_dim1]);
/*<             WA(K) = SIGMA(K) >*/
	    wa[k] = sigma[k];
/*<    80       CONTINUE >*/
L80:
/*<    90       CONTINUE >*/
/* L90: */
	    ;
	}
/*<   100    CONTINUE >*/
L100:
/*<          SIGMA(J) = -AJNORM >*/
	sigma[j] = -ajnorm;
/*<   110    CONTINUE >*/
/* L110: */
    }
/*<       RETURN >*/
    return 0;

/*     LAST CARD OF SUBROUTINE DQRFAC. */

/*<       END >*/
} /* dqrfac_ */

/*<       SUBROUTINE DQRSLV (N, R, LDR, IPVT, DIAG, QTB, X, SIGMA, WA) >*/
/* Subroutine */ int dqrslv_(integer *n, doublereal *r, integer *ldr, integer 
	*ipvt, doublereal *diag, doublereal *qtb, doublereal *x, doublereal *
	sigma, doublereal *wa)
{
    /* Initialized data */

    static doublereal p5 = .5;
    static doublereal p25 = .25;
    static doublereal zero = 0.;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    doublereal temp;
    integer i, j, k, l;
    doublereal cotan;
    integer nsing;
    doublereal qtbpj;
    integer jp1, kp1;
    doublereal tan_, cos_, sin_, sum;

/* ***BEGIN PROLOGUE  DQRSLV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNLS1 and DNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (QRSOLV-S, DQRSLV-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  **** Double Precision version of QRSOLV **** */

/*     Given an M by N matrix A, an N by N diagonal matrix D, */
/*     and an M-vector B, the problem is to determine an X which */
/*     solves the system */

/*           A*X = B ,     D*X = 0 , */

/*     in the least squares sense. */

/*     This subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     QR factorization, with column pivoting, of A. That is, if */
/*     A*P = Q*R, where P is a permutation matrix, Q has orthogonal */
/*     columns, and R is an upper triangular matrix with diagonal */
/*     elements of nonincreasing magnitude, then DQRSLV expects */
/*     the full upper triangle of R, the permutation matrix P, */
/*     and the first N components of (Q TRANSPOSE)*B. The system */
/*     A*X = B, D*X = 0, is then equivalent to */

/*                  T       T */
/*           R*Z = Q *B ,  P *D*P*Z = 0 , */

/*     where X = P*Z. If this system does not have full rank, */
/*     then a least squares solution is obtained. On output DQRSLV */
/*     also provides an upper triangular matrix S such that */

/*            T   T               T */
/*           P *(A *A + D*D)*P = S *S . */

/*     S is computed within DQRSLV and may be of separate interest. */

/*     The subroutine statement is */

/*       SUBROUTINE DQRSLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an N by N array. On input the full upper triangle */
/*         must contain the full upper triangle of the matrix R. */
/*         On output the full upper triangle is unaltered, and the */
/*         strict lower triangle contains the strict upper triangle */
/*         (transposed) of the upper triangular matrix S. */

/*       LDR is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array R. */

/*       IPVT is an integer input array of length N which defines the */
/*         permutation matrix P such that A*P = Q*R. Column J of P */
/*         is column IPVT(J) of the identity matrix. */

/*       DIAG is an input array of length N which must contain the */
/*         diagonal elements of the matrix D. */

/*       QTB is an input array of length N which must contain the first */
/*         N elements of the vector (Q TRANSPOSE)*B. */

/*       X is an output array of length N which contains the least */
/*         squares solution of the system A*X = B, D*X = 0. */

/*       SIGMA is an output array of length N which contains the */
/*         diagonal elements of the upper triangular matrix S. */

/*       WA is a work array of length N. */

/* ***SEE ALSO  DNLS1, DNLS1E */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DQRSLV */
/*<       INTEGER N,LDR >*/
/*<       INTEGER IPVT(*) >*/
/*<       DOUBLE PRECISION R(LDR,*),DIAG(*),QTB(*),X(*),SIGMA(*),WA(*) >*/
/*<       INTEGER I,J,JP1,K,KP1,L,NSING >*/
/*<       DOUBLE PRECISION COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO >*/
/*<       SAVE P5, P25, ZERO >*/
/*<       DATA P5,P25,ZERO /5.0D-1,2.5D-1,0.0D0/ >*/
    /* Parameter adjustments */
    --wa;
    --sigma;
    --x;
    --qtb;
    --diag;
    --ipvt;
    r_dim1 = *ldr;
    r_offset = r_dim1 + 1;
    r -= r_offset;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DQRSLV */
/*<       DO 20 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          DO 10 I = J, N >*/
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
/*<             R(I,J) = R(J,I) >*/
	    r[i + j * r_dim1] = r[j + i * r_dim1];
/*<    10       CONTINUE >*/
/* L10: */
	}
/*<          X(J) = R(J,J) >*/
	x[j] = r[j + j * r_dim1];
/*<          WA(J) = QTB(J) >*/
	wa[j] = qtb[j];
/*<    20    CONTINUE >*/
/* L20: */
    }

/*     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION. */

/*<       DO 100 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {

/*        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE */
/*        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION. */

/*<          L = IPVT(J) >*/
	l = ipvt[j];
/*<          IF (DIAG(L) .EQ. ZERO) GO TO 90 >*/
	if (diag[l] == zero) {
	    goto L90;
	}
/*<          DO 30 K = J, N >*/
	i__2 = *n;
	for (k = j; k <= i__2; ++k) {
/*<             SIGMA(K) = ZERO >*/
	    sigma[k] = zero;
/*<    30       CONTINUE >*/
/* L30: */
	}
/*<          SIGMA(J) = DIAG(L) >*/
	sigma[j] = diag[l];

/*        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D */
/*        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B */
/*        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO. */

/*<          QTBPJ = ZERO >*/
	qtbpj = zero;
/*<          DO 80 K = J, N >*/
	i__2 = *n;
	for (k = j; k <= i__2; ++k) {

/*           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE */
/*           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D. */

/*<             IF (SIGMA(K) .EQ. ZERO) GO TO 70 >*/
	    if (sigma[k] == zero) {
		goto L70;
	    }
/*<             IF (ABS(R(K,K)) .GE. ABS(SIGMA(K))) GO TO 40 >*/
	    if ((d__1 = r[k + k * r_dim1], abs(d__1)) >= (d__2 = sigma[k], 
		    abs(d__2))) {
		goto L40;
	    }
/*<                COTAN = R(K,K)/SIGMA(K) >*/
	    cotan = r[k + k * r_dim1] / sigma[k];
/*<                SIN = P5/SQRT(P25+P25*COTAN**2) >*/
/* Computing 2nd power */
	    d__1 = cotan;
	    sin_ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<                COS = SIN*COTAN >*/
	    cos_ = sin_ * cotan;
/*<                GO TO 50 >*/
	    goto L50;
/*<    40       CONTINUE >*/
L40:
/*<                TAN = SIGMA(K)/R(K,K) >*/
	    tan_ = sigma[k] / r[k + k * r_dim1];
/*<                COS = P5/SQRT(P25+P25*TAN**2) >*/
/* Computing 2nd power */
	    d__1 = tan_;
	    cos_ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<                SIN = COS*TAN >*/
	    sin_ = cos_ * tan_;
/*<    50       CONTINUE >*/
L50:

/*           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND */
/*           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0). */

/*<             R(K,K) = COS*R(K,K) + SIN*SIGMA(K) >*/
	    r[k + k * r_dim1] = cos_ * r[k + k * r_dim1] + sin_ * sigma[k];
/*<             TEMP = COS*WA(K) + SIN*QTBPJ >*/
	    temp = cos_ * wa[k] + sin_ * qtbpj;
/*<             QTBPJ = -SIN*WA(K) + COS*QTBPJ >*/
	    qtbpj = -sin_ * wa[k] + cos_ * qtbpj;
/*<             WA(K) = TEMP >*/
	    wa[k] = temp;

/*           ACCUMULATE THE TRANSFORMATION IN THE ROW OF S. */

/*<             KP1 = K + 1 >*/
	    kp1 = k + 1;
/*<             IF (N .LT. KP1) GO TO 70 >*/
	    if (*n < kp1) {
		goto L70;
	    }
/*<             DO 60 I = KP1, N >*/
	    i__3 = *n;
	    for (i = kp1; i <= i__3; ++i) {
/*<                TEMP = COS*R(I,K) + SIN*SIGMA(I) >*/
		temp = cos_ * r[i + k * r_dim1] + sin_ * sigma[i];
/*<                SIGMA(I) = -SIN*R(I,K) + COS*SIGMA(I) >*/
		sigma[i] = -sin_ * r[i + k * r_dim1] + cos_ * sigma[i];
/*<                R(I,K) = TEMP >*/
		r[i + k * r_dim1] = temp;
/*<    60          CONTINUE >*/
/* L60: */
	    }
/*<    70       CONTINUE >*/
L70:
/*<    80       CONTINUE >*/
/* L80: */
	    ;
	}
/*<    90    CONTINUE >*/
L90:

/*        STORE THE DIAGONAL ELEMENT OF S AND RESTORE */
/*        THE CORRESPONDING DIAGONAL ELEMENT OF R. */

/*<          SIGMA(J) = R(J,J) >*/
	sigma[j] = r[j + j * r_dim1];
/*<          R(J,J) = X(J) >*/
	r[j + j * r_dim1] = x[j];
/*<   100    CONTINUE >*/
/* L100: */
    }

/*     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS */
/*     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION. */

/*<       NSING = N >*/
    nsing = *n;
/*<       DO 110 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          IF (SIGMA(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1 >*/
	if (sigma[j] == zero && nsing == *n) {
	    nsing = j - 1;
	}
/*<          IF (NSING .LT. N) WA(J) = ZERO >*/
	if (nsing < *n) {
	    wa[j] = zero;
	}
/*<   110    CONTINUE >*/
/* L110: */
    }
/*<       IF (NSING .LT. 1) GO TO 150 >*/
    if (nsing < 1) {
	goto L150;
    }
/*<       DO 140 K = 1, NSING >*/
    i__1 = nsing;
    for (k = 1; k <= i__1; ++k) {
/*<          J = NSING - K + 1 >*/
	j = nsing - k + 1;
/*<          SUM = ZERO >*/
	sum = zero;
/*<          JP1 = J + 1 >*/
	jp1 = j + 1;
/*<          IF (NSING .LT. JP1) GO TO 130 >*/
	if (nsing < jp1) {
	    goto L130;
	}
/*<          DO 120 I = JP1, NSING >*/
	i__2 = nsing;
	for (i = jp1; i <= i__2; ++i) {
/*<             SUM = SUM + R(I,J)*WA(I) >*/
	    sum += r[i + j * r_dim1] * wa[i];
/*<   120       CONTINUE >*/
/* L120: */
	}
/*<   130    CONTINUE >*/
L130:
/*<          WA(J) = (WA(J) - SUM)/SIGMA(J) >*/
	wa[j] = (wa[j] - sum) / sigma[j];
/*<   140    CONTINUE >*/
/* L140: */
    }
/*<   150 CONTINUE >*/
L150:

/*     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X. */

/*<       DO 160 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          L = IPVT(J) >*/
	l = ipvt[j];
/*<          X(L) = WA(J) >*/
	x[l] = wa[j];
/*<   160    CONTINUE >*/
/* L160: */
    }
/*<       RETURN >*/
    return 0;

/*     LAST CARD OF SUBROUTINE DQRSLV. */

/*<       END >*/
} /* dqrslv_ */

/*<       SUBROUTINE DWUPDT (N, R, LDR, W, B, ALPHA, COS, SIN) >*/
/* Subroutine */ int dwupdt_(integer *n, doublereal *r, integer *ldr, 
	doublereal *w, doublereal *b, doublereal *alpha, doublereal *cos_, 
	doublereal *sin_)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal p5 = .5;
    static doublereal p25 = .25;
    static doublereal zero = 0.;

    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    /*double sqrt(doublereal);*/

    /* Local variables */
    doublereal temp, rowj;
    integer i, j;
    doublereal cotan;
    integer jm1;
    doublereal tan_;

/* ***BEGIN PROLOGUE  DWUPDT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNLS1 and DNLS1E */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (RWUPDT-S, DWUPDT-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an N by N upper triangular matrix R, this subroutine */
/*     computes the QR decomposition of the matrix formed when a row */
/*     is added to R. If the row is specified by the vector W, then */
/*     DWUPDT determines an orthogonal matrix Q such that when the */
/*     N+1 by N matrix composed of R augmented by W is premultiplied */
/*     by (Q TRANSPOSE), the resulting matrix is upper trapezoidal. */
/*     The orthogonal matrix Q is the product of N transformations */

/*           G(1)*G(2)* ... *G(N) */

/*     where G(I) is a Givens rotation in the (I,N+1) plane which */
/*     eliminates elements in the I-th plane. DWUPDT also */
/*     computes the product (Q TRANSPOSE)*C where C is the */
/*     (N+1)-vector (b,alpha). Q itself is not accumulated, rather */
/*     the information to recover the G rotations is supplied. */

/*     The subroutine statement is */

/*       SUBROUTINE DWUPDT(N,R,LDR,W,B,ALPHA,COS,SIN) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an N by N array. On input the upper triangular part of */
/*         R must contain the matrix to be updated. On output R */
/*         contains the updated triangular matrix. */

/*       LDR is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array R. */

/*       W is an input array of length N which must contain the row */
/*         vector to be added to R. */

/*       B is an array of length N. On input B must contain the */
/*         first N elements of the vector C. On output B contains */
/*         the first N elements of the vector (Q TRANSPOSE)*C. */

/*       ALPHA is a variable. On input ALPHA must contain the */
/*         (N+1)-st element of the vector C. On output ALPHA contains */
/*         the (N+1)-st element of the vector (Q TRANSPOSE)*C. */

/*       COS is an output array of length N which contains the */
/*         cosines of the transforming Givens rotations. */

/*       SIN is an output array of length N which contains the */
/*         sines of the transforming Givens rotations. */

/*     ********** */

/* ***SEE ALSO  DNLS1, DNLS1E */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DWUPDT */
/*<       INTEGER N,LDR >*/
/*<       DOUBLE PRECISION ALPHA >*/
/*<       DOUBLE PRECISION R(LDR,*),W(*),B(*),COS(*),SIN(*) >*/
/*<       INTEGER I,J,JM1 >*/
/*<       DOUBLE PRECISION COTAN,ONE,P5,P25,ROWJ,TAN,TEMP,ZERO >*/
/*<       SAVE ONE, P5, P25, ZERO >*/
/*<       DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/ >*/
    /* Parameter adjustments */
    --sin_;
    --cos_;
    --b;
    --w;
    r_dim1 = *ldr;
    r_offset = r_dim1 + 1;
    r -= r_offset;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DWUPDT */
/*<       DO 60 J = 1, N >*/
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*<          ROWJ = W(J) >*/
	rowj = w[j];
/*<          JM1 = J - 1 >*/
	jm1 = j - 1;

/*        APPLY THE PREVIOUS TRANSFORMATIONS TO */
/*        R(I,J), I=1,2,...,J-1, AND TO W(J). */

/*<          IF (JM1 .LT. 1) GO TO 20 >*/
	if (jm1 < 1) {
	    goto L20;
	}
/*<          DO 10 I = 1, JM1 >*/
	i__2 = jm1;
	for (i = 1; i <= i__2; ++i) {
/*<             TEMP = COS(I)*R(I,J) + SIN(I)*ROWJ >*/
	    temp = cos_[i] * r[i + j * r_dim1] + sin_[i] * rowj;
/*<             ROWJ = -SIN(I)*R(I,J) + COS(I)*ROWJ >*/
	    rowj = -sin_[i] * r[i + j * r_dim1] + cos_[i] * rowj;
/*<             R(I,J) = TEMP >*/
	    r[i + j * r_dim1] = temp;
/*<    10       CONTINUE >*/
/* L10: */
	}
/*<    20    CONTINUE >*/
L20:

/*        DETERMINE A GIVENS ROTATION WHICH ELIMINATES W(J). */

/*<          COS(J) = ONE >*/
	cos_[j] = one;
/*<          SIN(J) = ZERO >*/
	sin_[j] = zero;
/*<          IF (ROWJ .EQ. ZERO) GO TO 50 >*/
	if (rowj == zero) {
	    goto L50;
	}
/*<          IF (ABS(R(J,J)) .GE. ABS(ROWJ)) GO TO 30 >*/
	if ((d__1 = r[j + j * r_dim1], abs(d__1)) >= abs(rowj)) {
	    goto L30;
	}
/*<             COTAN = R(J,J)/ROWJ >*/
	cotan = r[j + j * r_dim1] / rowj;
/*<             SIN(J) = P5/SQRT(P25+P25*COTAN**2) >*/
/* Computing 2nd power */
	d__1 = cotan;
	sin_[j] = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<             COS(J) = SIN(J)*COTAN >*/
	cos_[j] = sin_[j] * cotan;
/*<             GO TO 40 >*/
	goto L40;
/*<    30    CONTINUE >*/
L30:
/*<             TAN = ROWJ/R(J,J) >*/
	tan_ = rowj / r[j + j * r_dim1];
/*<             COS(J) = P5/SQRT(P25+P25*TAN**2) >*/
/* Computing 2nd power */
	d__1 = tan_;
	cos_[j] = p5 / sqrt(p25 + p25 * (d__1 * d__1));
/*<             SIN(J) = COS(J)*TAN >*/
	sin_[j] = cos_[j] * tan_;
/*<    40    CONTINUE >*/
L40:

/*        APPLY THE CURRENT TRANSFORMATION TO R(J,J), B(J), AND ALPHA.
 */

/*<          R(J,J) = COS(J)*R(J,J) + SIN(J)*ROWJ >*/
	r[j + j * r_dim1] = cos_[j] * r[j + j * r_dim1] + sin_[j] * rowj;
/*<          TEMP = COS(J)*B(J) + SIN(J)*ALPHA >*/
	temp = cos_[j] * b[j] + sin_[j] * *alpha;
/*<          ALPHA = -SIN(J)*B(J) + COS(J)*ALPHA >*/
	*alpha = -sin_[j] * b[j] + cos_[j] * *alpha;
/*<          B(J) = TEMP >*/
	b[j] = temp;
/*<    50    CONTINUE >*/
L50:
/*<    60    CONTINUE >*/
/* L60: */
	;
    }
/*<       RETURN >*/
    return 0;

/*     LAST CARD OF SUBROUTINE DWUPDT. */

/*<       END >*/
} /* dwupdt_ */

/*<       SUBROUTINE FDUMP >*/
/* Subroutine */ int fdump_(void)
{
/* ***BEGIN PROLOGUE  FDUMP */
/* ***PURPOSE  Symbolic dump (should be locally written). */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3 */
/* ***TYPE      ALL (FDUMP-A) */
/* ***KEYWORDS  ERROR, XERMSG */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*        ***Note*** Machine Dependent Routine */
/*        FDUMP is intended to be replaced by a locally written */
/*        version which produces a symbolic dump.  Failing this, */
/*        it should be replaced by a version which prints the */
/*        subprogram nesting list.  Note that this dump must be */
/*        printed on each of up to five files, as indicated by the */
/*        XGETUA routine.  See XSETUA and XGETUA for details. */

/*     Written by Ron Jones, with SLATEC Common Math Library Subcommittee 
*/

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  FDUMP */
/* ***FIRST EXECUTABLE STATEMENT  FDUMP */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* fdump_ */

/*<       INTEGER FUNCTION I1MACH (I) >*/
integer i1mach_(integer *i)
{
    /* System generated locals */
    integer ret_val = 0;
    static integer equiv_0[16];

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
#define imach (equiv_0)
#define output (equiv_0 + 3)

/* ***BEGIN PROLOGUE  I1MACH */
/* ***PURPOSE  Return integer machine dependent constants. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  R1 */
/* ***TYPE      INTEGER (I1MACH-I) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Fox, P. A., (Bell Labs) */
/*           Hall, A. D., (Bell Labs) */
/*           Schryer, N. L., (Bell Labs) */
/* ***DESCRIPTION */

/*   I1MACH can be used to obtain machine-dependent parameters for the */
/*   local machine environment.  It is a function subprogram with one */
/*   (input) argument and can be referenced as follows: */

/*        K = I1MACH(I) */

/*   where I=1,...,16.  The (output) value of K above is determined by */
/*   the (input) value of I.  The results for various values of I are */
/*   discussed below. */

/*   I/O unit numbers: */
/*     I1MACH( 1) = the standard input unit. */
/*     I1MACH( 2) = the standard output unit. */
/*     I1MACH( 3) = the standard punch unit. */
/*     I1MACH( 4) = the standard error message unit. */

/*   Words: */
/*     I1MACH( 5) = the number of bits per integer storage unit. */
/*     I1MACH( 6) = the number of characters per integer storage unit. */

/*   Integers: */
/*     assume integers are represented in the S-digit, base-A form */

/*                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*                where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*     I1MACH( 7) = A, the base. */
/*     I1MACH( 8) = S, the number of base-A digits. */
/*     I1MACH( 9) = A**S - 1, the largest magnitude. */

/*   Floating-Point Numbers: */
/*     Assume floating-point numbers are represented in the T-digit, */
/*     base-B form */
/*                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*                where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*                0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*     I1MACH(10) = B, the base. */

/*   Single-Precision: */
/*     I1MACH(11) = T, the number of base-B digits. */
/*     I1MACH(12) = EMIN, the smallest exponent E. */
/*     I1MACH(13) = EMAX, the largest exponent E. */

/*   Double-Precision: */
/*     I1MACH(14) = T, the number of base-B digits. */
/*     I1MACH(15) = EMIN, the smallest exponent E. */
/*     I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, the desired */
/*   set of DATA statements should be activated by removing the C from */
/*   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be */
/*   checked for consistency with the local operating system. */

/* ***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for 
*/
/*                 a portable library, ACM Transactions on Mathematical */
/*                 Software 4, 2 (June 1978), pp. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891012  Added VAX G-floating constants.  (WRB) */
/*   891012  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900618  Added DEC RISC constants.  (WRB) */
/*   900723  Added IBM RS 6000 constants.  (WRB) */
/*   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16. */
/*           (RWC) */
/*   910710  Added HP 730 constants.  (SMR) */
/*   911114  Added Convex IEEE constants.  (WRB) */
/*   920121  Added SUN -r8 compiler option constants.  (WRB) */
/*   920229  Added Touchstone Delta i860 constants.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   920625  Added Convex -p8 and -pd8 compiler option constants. */
/*           (BKS, WRB) */
/*   930201  Added DEC Alpha and SGI constants.  (RWC and WRB) */
/*   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler */
/*           options.  (DWL, RWC and WRB). */
/* ***END PROLOGUE  I1MACH */

/*<       INTEGER IMACH(16),OUTPUT >*/
/*<       SAVE IMACH >*/
/*<       EQUIVALENCE (IMACH(4),OUTPUT) >*/

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        129 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1025 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA IMACH( 1) /          7 / */
/*     DATA IMACH( 2) /          2 / */
/*     DATA IMACH( 3) /          2 / */
/*     DATA IMACH( 4) /          2 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         33 / */
/*     DATA IMACH( 9) / Z1FFFFFFFF / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -256 / */
/*     DATA IMACH(13) /        255 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /       -256 / */
/*     DATA IMACH(16) /        255 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /        -50 / */
/*     DATA IMACH(16) /         76 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         48 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /          8 / */
/*     DATA IMACH(11) /         13 / */
/*     DATA IMACH(12) /        -50 / */
/*     DATA IMACH(13) /         76 / */
/*     DATA IMACH(14) /         26 / */
/*     DATA IMACH(15) /     -32754 / */
/*     DATA IMACH(16) /      32780 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 00007777777777777777B / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE CRAY */
/*     USING THE 64 BIT INTEGER COMPILER OPTION */

/*     DATA IMACH( 1) /        100 / */
/*     DATA IMACH( 2) /        101 / */
/*     DATA IMACH( 3) /        102 / */
/*     DATA IMACH( 4) /        101 / */
/*     DATA IMACH( 5) /         64 / */
/*     DATA IMACH( 6) /          8 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         63 / */
/*     DATA IMACH( 9) / 777777777777777777777B / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         47 / */
/*     DATA IMACH(12) /      -8189 / */
/*     DATA IMACH(13) /       8190 / */
/*     DATA IMACH(14) /         94 / */
/*     DATA IMACH(15) /      -8099 / */
/*     DATA IMACH(16) /       8190 / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     DATA IMACH( 1) /         11 / */
/*     DATA IMACH( 2) /         12 / */
/*     DATA IMACH( 3) /          8 / */
/*     DATA IMACH( 4) /         10 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING G_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE DEC ALPHA */
/*     USING IEEE_FLOAT */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC RISC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING D_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE DEC VAX */
/*     USING G_FLOATING */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1023 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1022 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         24 / */
/*     DATA IMACH( 6) /          3 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         23 / */
/*     DATA IMACH( 9) /    8388607 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         38 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /         43 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          6 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         63 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 730 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         39 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          4 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         23 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         55 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */
/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          7 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         32 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -126 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1015 / */
/*     DATA IMACH(16) /       1017 / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          7 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) /  Z7FFFFFFF / */
/*     DATA IMACH(10) /         16 / */
/*     DATA IMACH(11) /          6 / */
/*     DATA IMACH(12) /        -64 / */
/*     DATA IMACH(13) /         63 / */
/*     DATA IMACH(14) /         14 / */
/*     DATA IMACH(15) /        -64 / */
/*     DATA IMACH(16) /         63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1023 / */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE INTEL i860 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         54 / */
/*     DATA IMACH(15) /       -101 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          5 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/* ent. */

/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         62 / */
/*     DATA IMACH(15) /       -128 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          5 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE SUN */
/*     USING THE -r8 COMPILER OPTION */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         53 / */
/*     DATA IMACH(12) /      -1021 / */
/*     DATA IMACH(13) /       1024 / */
/*     DATA IMACH(14) /        113 / */
/*     DATA IMACH(15) /     -16381 / */
/*     DATA IMACH(16) /      16384 / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          1 / */
/*     DATA IMACH( 4) /          6 / */
/*     DATA IMACH( 5) /         36 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         27 / */
/*     DATA IMACH(12) /       -128 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         60 / */
/*     DATA IMACH(15) /      -1024 / */
/*     DATA IMACH(16) /       1023 / */

/* ected XERN1 declaration.  (WRB) */
/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA IMACH( 1) /          1 / */
/*     DATA IMACH( 2) /          1 / */
/*     DATA IMACH( 3) /          0 / */
/*     DATA IMACH( 4) /          1 / */
/*     DATA IMACH( 5) /         16 / */
/*     DATA IMACH( 6) /          2 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         15 / */
/*     DATA IMACH( 9) /      32767 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -127 / */
/*     DATA IMACH(13) /        127 / */
/*     DATA IMACH(14) /         56 / */
/*     DATA IMACH(15) /       -127 / */
/*     DATA IMACH(16) /        127 / */

/* ***FIRST EXECUTABLE STATEMENT  I1MACH */
/*<       IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10 >*/
    if (*i < 1 || *i > 16) {
	goto L10;
    }

/*<       I1MACH = IMACH(I) >*/
    ret_val = imach[*i - 1];
/*<       RETURN >*/
    return ret_val;

/*<    10 CONTINUE >*/
L10:
/*      WRITE (UNIT = OUTPUT, FMT = 9000) */
/* 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS') */

/*     CALL FDUMP */

/*<       STOP >*/
    s_stop("", 0L);
/*<       END >*/
    return ret_val;
} /* i1mach_ */

#undef output
#undef imach


/*<       FUNCTION J4SAVE (IWHICH, IVALUE, ISET) >*/
integer j4save_(integer *iwhich, integer *ivalue, logical *iset)
{
    /* Initialized data */

    static integer iparam[9] = { 0,2,0,10,1,0,0,0,0 };

    /* System generated locals */
    integer ret_val;

/* ***BEGIN PROLOGUE  J4SAVE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Save or recall global variables needed by error */
/*            handling routines. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***TYPE      INTEGER (J4SAVE-I) */
/* ***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        J4SAVE saves and recalls several global variables needed */
/*        by the library error handling routines. */

/*     Description of Parameters */
/*      --Input-- */
/*        IWHICH - Index of item desired. */
/*                = 1 Refers to current error number. */
/*                = 2 Refers to current error control flag. */
/*                = 3 Refers to current unit number to which error */
/*                    messages are to be sent.  (0 means use standard.) */
/*                = 4 Refers to the maximum number of times any */
/*                     message is to be printed (as set by XERMAX). */
/*                = 5 Refers to the total number of units to which */
/*                     each error message is to be written. */
/*                = 6 Refers to the 2nd unit for error messages */
/*                = 7 Refers to the 3rd unit for error messages */
/*                = 8 Refers to the 4th unit for error messages */
/*                = 9 Refers to the 5th unit for error messages */
/*        IVALUE - The value to be set for the IWHICH-th parameter, */
/*                 if ISET is .TRUE. . */
/*        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE */
/*                 given the value, IVALUE.  If ISET=.FALSE., the */
/*                 IWHICH-th parameter will be unchanged, and IVALUE */
/*                 is a dummy parameter. */
/*      --Output-- */
/*        The (old) value of the IWHICH-th parameter will be returned */
/*        in the function value, J4SAVE. */

/* ***SEE ALSO  XERMSG */
/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900205  Minor modifications to prologue.  (WRB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910411  Added KEYWORDS section.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  J4SAVE */
/*<       LOGICAL ISET >*/
/*<       INTEGER IPARAM(9) >*/
/*<       SAVE IPARAM >*/
/*<       DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/ >*/
/*<       DATA IPARAM(5)/1/ >*/
/*<       DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/ >*/
/* ***FIRST EXECUTABLE STATEMENT  J4SAVE */
/*<       J4SAVE = IPARAM(IWHICH) >*/
    ret_val = iparam[*iwhich - 1];
/*<       IF (ISET) IPARAM(IWHICH) = IVALUE >*/
    if (*iset) {
	iparam[*iwhich - 1] = *ivalue;
    }
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* j4save_ */

/*<       SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL) >*/
/* Subroutine */ int xercnt_(char *librar, char *subrou, char *messg, integer 
	*nerr, integer *level, integer *kontrl, ftnlen librar_len, ftnlen 
	subrou_len, ftnlen messg_len)
{
/* ***BEGIN PROLOGUE  XERCNT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Allow user control over handling of errors. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERCNT-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        Allows user control over handling of individual errors. */
/*        Just after each message is recorded, but before it is */
/*        processed any further (i.e., before it is printed or */
/*        a decision to abort is made), a call is made to XERCNT. */
/*        If the user has provided his own version of XERCNT, he */
/*        can then override the value of KONTROL used in processing */
/*        this message by redefining its value. */
/*        KONTRL may be set to any value from -2 to 2. */
/*        The meanings for KONTRL are the same as in XSETF, except */
/*        that the value of KONTRL changes only for this message. */
/*        If KONTRL is set to a value outside the range from -2 to 2, */
/*        it will be moved back into that range. */

/*     Description of Parameters */

/*      --Input-- */
/*        LIBRAR - the library that the routine is in. */
/*        SUBROU - the subroutine that XERMSG is being called from */
/*        MESSG  - the first 20 characters of the error message. */
/*        NERR   - same as in the call to XERMSG. */
/*        LEVEL  - same as in the call to XERMSG. */
/*        KONTRL - the current value of the control flag as set */
/*                 by a call to XSETF. */

/*      --Output-- */
/*        KONTRL - the new value of KONTRL.  If KONTRL is not */
/*                 defined, it will remain at its original value. */
/*                 This changed value of control affects only */
/*                 the current occurrence of the current message. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900206  Routine changed from user-callable to subsidiary.  (WRB) */
/*   900510  Changed calling sequence to include LIBRARY and SUBROUTINE */
/*           names, changed routine name from XERCTL to XERCNT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERCNT */
/*<       CHARACTER*(*) LIBRAR, SUBROU, MESSG >*/
/* ***FIRST EXECUTABLE STATEMENT  XERCNT */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* xercnt_ */

/*<       SUBROUTINE XERHLT (MESSG) >*/
/* Subroutine */ int xerhlt_(char *messg, ftnlen messg_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

/* ***BEGIN PROLOGUE  XERHLT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Abort program execution and print error message. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERHLT-A) */
/* ***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        ***Note*** machine dependent routine */
/*        XERHLT aborts the execution of the program. */
/*        The error message causing the abort is given in the calling */
/*        sequence, in case one needs it for printing on a dayfile, */
/*        for example. */

/*     Description of Parameters */
/*        MESSG is as in XERMSG. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900206  Routine changed from user-callable to subsidiary.  (WRB) */
/*   900510  Changed calling sequence to delete length of character */
/*           and changed routine name from XERABT to XERHLT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERHLT */
/*<       CHARACTER*(*) MESSG >*/
/* ***FIRST EXECUTABLE STATEMENT  XERHLT */
/*<       STOP >*/
    s_stop("", 0L);
/*<       END >*/
    return 0;
} /* xerhlt_ */

/*<       SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL) >*/
/* Subroutine */ int xermsg_(char *librar, char *subrou, char *messg, integer 
	*nerr, integer *level, ftnlen librar_len, ftnlen subrou_len, ftnlen 
	messg_len)
{
    /* System generated locals */
    address a__1[2];
    integer i__1, i__2, i__3[2];
    char ch__1[87];

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_len(char *, ftnlen);
    /* Subroutine */ void s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    integer lerr;
    char temp[72];
    integer i;
    extern /* Subroutine */ int fdump_(void);
    char xlibr[8];
    integer ltemp, kount=0;
    char xsubr[8];
    extern integer j4save_(integer *, integer *, logical *);
    integer llevel, maxmes;
    char lfirst[20];
    extern /* Subroutine */ int xercnt_(char *, char *, char *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);
    integer lkntrl;
    extern /* Subroutine */ int xerhlt_(char *, ftnlen);
    integer mkntrl;
    extern /* Subroutine */ int xersve_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, ftnlen, ftnlen, ftnlen), xerprn_(
	    char *, integer *, char *, integer *, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  XERMSG */
/* ***PURPOSE  Process error messages for SLATEC and other libraries. */
/*             to ignore trailing blanks.  Another feature is that */
/*    must be given a chance to call NUMXER or XERCLR to retrieve or */
/*    clear the error number. */
/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   880101  DATE WRITTEN */
/*   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988. 
*/
/*           THERE ARE TWO BASIC CHANGES. */
/*           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO */
/*               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES */
/*               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS */
/*               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE */
/*               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER */
/*               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY */
/*               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE */
/*               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76. */
/*           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE 
*/
/*               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE */
/*               OF LOWER CASE. */
/*   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30. */
/*           THE PRINCIPAL CHANGES ARE */
/*           1.  CLARIFY COMMENTS IN THE PROLOGUES */
/*           2.  RENAME XRPRNT TO XERPRN */
/*           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES 
*/
/*               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE / */
/*               CHARACTER FOR NEW RECORDS. */
/*   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO */
/*           CLEAN UP THE CODING. */
/*   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN 
*/
/*           PREFIX. */
/*   891013  REVISED TO CORRECT COMMENTS. */
/*   891214  Prologue converted to Version 4.0 format.  (WRB) */
/*   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but */
/*           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added */
/*           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and */
/*           XERCTL to XERCNT.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERMSG */
/*<       CHARACTER*(*) LIBRAR, SUBROU, MESSG >*/
/*<       CHARACTER*8 XLIBR, XSUBR >*/
/*<       CHARACTER*72  TEMP >*/
/*<       CHARACTER*20  LFIRST >*/
/* ***FIRST EXECUTABLE STATEMENT  XERMSG */
/*<       LKNTRL = J4SAVE (2, 0, .FALSE.) >*/
    lkntrl = j4save_(&c__2, &c__0, (logical*)&c__0);
/*<       MAXMES = J4SAVE (4, 0, .FALSE.) >*/
    maxmes = j4save_(&c__4, &c__0, (logical*)&c__0);

/*       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL. */
/*       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE */
/*          SHOULD BE PRINTED. */

/*       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN */
/*          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE, */
/*          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2. */

/*<    >*/
    if (*nerr < -9999999 || *nerr > 99999999 || *nerr == 0 || *level < -1 || *
	    level > 2) {
/*<    >*/
	xerprn_(" ***", &c_n1, "FATAL ERROR IN...$$ XERMSG -- INVALID ERROR "
		"NUMBER OR LEVEL$$ JOB ABORT DUE TO FATAL ERROR.", &c__72, 4L, 
		91L);
/*<          CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY) >*/
	//	xersve_(" ", " ", " ", &c__0, &c__0, &c__0, &kdummy, 1L, 1L, 1L);
/*<          CALL XERHLT (' ***XERMSG -- INVALID INPUT') >*/
	xerhlt_(" ***XERMSG -- INVALID INPUT", 27L);
/*<          RETURN >*/
	return 0;
/*<       ENDIF >*/
    }

/*       RECORD THE MESSAGE. */

/*<       I = J4SAVE (1, NERR, .TRUE.) >*/
    i = j4save_(&c__1, nerr, (logical*)&c__1);
/*<       CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT) >*/
    //    xersve_(librar, subrou, messg, &c__1, nerr, level, &kount, librar_len, 
    //	    subrou_len, messg_len);

/*       HANDLE PRINT-ONCE WARNING MESSAGES. */

/*<       IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN >*/
    if (*level == -1 && kount > 1) {
	return 0;
    }

/*       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG. */

/*<       XLIBR  = LIBRAR >*/
    s_copy(xlibr, librar, 8L, librar_len);
/*<       XSUBR  = SUBROU >*/
    s_copy(xsubr, subrou, 8L, subrou_len);
/*<       LFIRST = MESSG >*/
    s_copy(lfirst, messg, 20L, messg_len);
/*<       LERR   = NERR >*/
    lerr = *nerr;
/*<       LLEVEL = LEVEL >*/
    llevel = *level;
/*<       CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL) >*/
    xercnt_(xlibr, xsubr, lfirst, &lerr, &llevel, &lkntrl, 8L, 8L, 20L);

/*<       LKNTRL = MAX(-2, MIN(2,LKNTRL)) >*/
/* Computing MAX */
    i__1 = -2, i__2 = min(2,lkntrl);
    lkntrl = max(i__1,i__2);
/*<       MKNTRL = ABS(LKNTRL) >*/
    mkntrl = abs(lkntrl);

/*       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS */
/*       ZERO AND THE ERROR IS NOT FATAL. */

/*<       IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30 >*/
    if (*level < 2 && lkntrl == 0) {
	goto L30;
    }
/*<       IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30 >*/
    if (*level == 0 && kount > maxmes) {
	goto L30;
    }
/*<       IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30 >*/
    if (*level == 1 && kount > maxmes && mkntrl == 1) {
	goto L30;
    }
/*<       IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30 >*/
    if (*level == 2 && kount > max(1,maxmes)) {
	goto L30;
    }

/*       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A */
/*       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS) 
*/
/*       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG */
/*       IS NOT ZERO. */

/*<       IF (LKNTRL .NE. 0) THEN >*/
    if (lkntrl != 0) {
/*<          TEMP(1:21) = 'MESSAGE FROM ROUTINE ' >*/
	s_copy(temp, "MESSAGE FROM ROUTINE ", 21L, 21L);
/*<          I = MIN(LEN(SUBROU), 16) >*/
/* Computing MIN */
	i__1 = i_len(subrou, subrou_len);
	i = min(i__1,16);
/*<          TEMP(22:21+I) = SUBROU(1:I) >*/
	s_copy(temp + 21, subrou, i, i);
/*<          TEMP(22+I:33+I) = ' IN LIBRARY ' >*/
	i__1 = i + 21;
	s_copy(temp + i__1, " IN LIBRARY ", i + 33 - i__1, 12L);
/*<          LTEMP = 33 + I >*/
	ltemp = i + 33;
/*<          I = MIN(LEN(LIBRAR), 16) >*/
/* Computing MIN */
	i__1 = i_len(librar, librar_len);
	i = min(i__1,16);
/*<          TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I) >*/
	i__1 = ltemp;
	s_copy(temp + i__1, librar, ltemp + i - i__1, i);
/*<          TEMP(LTEMP+I+1:LTEMP+I+1) = '.' >*/
	i__1 = ltemp + i;
	s_copy(temp + i__1, ".", ltemp + i + 1 - i__1, 1L);
/*<          LTEMP = LTEMP + I + 1 >*/
	ltemp = ltemp + i + 1;
/*<          CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72) >*/
	xerprn_(" ***", &c_n1, temp, &c__72, 4L, ltemp);
/*<       ENDIF >*/
    }

/*       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE */
/*       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE */
/*       FROM EACH OF THE FOLLOWING THREE OPTIONS. */
/*       1.  LEVEL OF THE MESSAGE */
/*              'INFORMATIVE MESSAGE' */
/*              'POTENTIALLY RECOVERABLE ERROR' */
/*              'FATAL ERROR' */
/*       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE */
/*              'PROG CONTINUES' */
/*              'PROG ABORTED' */
/*       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK */
/*           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS */
/*           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.) */
/*              'TRACEBACK REQUESTED' */
/*              'TRACEBACK NOT REQUESTED' */
/*       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT */
/*       EXCEED 74 CHARACTERS. */
/*       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED. */

/*<       IF (LKNTRL .GT. 0) THEN >*/
    if (lkntrl > 0) {

/*       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL. */

/*<          IF (LEVEL .LE. 0) THEN >*/
	if (*level <= 0) {
/*<             TEMP(1:20) = 'INFORMATIVE MESSAGE,' >*/
	    s_copy(temp, "INFORMATIVE MESSAGE,", 20L, 20L);
/*<             LTEMP = 20 >*/
	    ltemp = 20;
/*<          ELSEIF (LEVEL .EQ. 1) THEN >*/
	} else if (*level == 1) {
/*<             TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,' >*/
	    s_copy(temp, "POTENTIALLY RECOVERABLE ERROR,", 30L, 30L);
/*<             LTEMP = 30 >*/
	    ltemp = 30;
/*<          ELSE >*/
	} else {
/*<             TEMP(1:12) = 'FATAL ERROR,' >*/
	    s_copy(temp, "FATAL ERROR,", 12L, 12L);
/*<             LTEMP = 12 >*/
	    ltemp = 12;
/*<          ENDIF >*/
	}

/*       THEN WHETHER THE PROGRAM WILL CONTINUE. */

/*<    >*/
	if (mkntrl == 2 && *level >= 1 || mkntrl == 1 && *level == 2) {
/*<             TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,' >*/
	    i__1 = ltemp;
	    s_copy(temp + i__1, " PROG ABORTED,", ltemp + 14 - i__1, 14L);
/*<             LTEMP = LTEMP + 14 >*/
	    ltemp += 14;
/*<          ELSE >*/
	} else {
/*<             TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,' >*/
	    i__1 = ltemp;
	    s_copy(temp + i__1, " PROG CONTINUES,", ltemp + 16 - i__1, 16L);
/*<             LTEMP = LTEMP + 16 >*/
	    ltemp += 16;
/*<          ENDIF >*/
	}

/*       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK. */

/*<          IF (LKNTRL .GT. 0) THEN >*/
	if (lkntrl > 0) {
/*<             TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED' >*/
	    i__1 = ltemp;
	    s_copy(temp + i__1, " TRACEBACK REQUESTED", ltemp + 20 - i__1, 
		    20L);
/*<             LTEMP = LTEMP + 20 >*/
	    ltemp += 20;
/*<          ELSE >*/
	} else {
/*<             TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED' >*/
	    i__1 = ltemp;
	    s_copy(temp + i__1, " TRACEBACK NOT REQUESTED", ltemp + 24 - i__1,
		     24L);
/*<             LTEMP = LTEMP + 24 >*/
	    ltemp += 24;
/*<          ENDIF >*/
	}
/*<          CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72) >*/
	xerprn_(" ***", &c_n1, temp, &c__72, 4L, ltemp);
/*<       ENDIF >*/
    }

/*       NOW SEND OUT THE MESSAGE. */

/*<       CALL XERPRN (' *  ', -1, MESSG, 72) >*/
    xerprn_(" *  ", &c_n1, messg, &c__72, 4L, messg_len);

/*       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A */
/*          TRACEBACK. */

/*<       IF (LKNTRL .GT. 0) THEN >*/
    if (lkntrl > 0) {
/*         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR */
/*<          DO 10 I=16,22 >*/
	for (i = 16; i <= 22; ++i) {
/*<             IF (TEMP(I:I) .NE. ' ') GO TO 20 >*/
	    if (temp[i - 1] != ' ') {
		goto L20;
	    }
/*<    10    CONTINUE >*/
/* L10: */
	}

/*<    20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72) >*/
L20:
/* Writing concatenation */
	i__3[0] = 15, a__1[0] = temp;
	i__3[1] = 23 - (i - 1), a__1[1] = temp + (i - 1);
	s_cat(ch__1, a__1, i__3, &c__2, 87L);
	xerprn_(" *  ", &c_n1, ch__1, &c__72, 4L, 23 - (i - 1) + 15);
/*<          CALL FDUMP >*/
	fdump_();
/*<       ENDIF >*/
    }

/*       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE. 
*/

/*<       IF (LKNTRL .NE. 0) THEN >*/
    if (lkntrl != 0) {
/*<          CALL XERPRN (' *  ', -1, ' ', 72) >*/
	xerprn_(" *  ", &c_n1, " ", &c__72, 4L, 1L);
/*<          CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72) >*/
	xerprn_(" ***", &c_n1, "END OF MESSAGE", &c__72, 4L, 14L);
/*<          CALL XERPRN ('    ',  0, ' ', 72) >*/
	xerprn_("    ", &c__0, " ", &c__72, 4L, 1L);
/*<       ENDIF >*/
    }

/*       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE */
/*       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN. */

/*<    30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN >*/
L30:
    if (*level <= 0 || *level == 1 && mkntrl <= 1) {
	return 0;
    }

/*       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A */
/*       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR */
/*       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT. 
*/

/*<       IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN >*/
    if (lkntrl > 0 && kount < max(1,maxmes)) {
/*<          IF (LEVEL .EQ. 1) THEN >*/
	if (*level == 1) {
/*<    >*/
	    xerprn_(" ***", &c_n1, "JOB ABORT DUE TO UNRECOVERED ERROR.", &
		    c__72, 4L, 35L);
/*<          ELSE >*/
	} else {
/*<             CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72) >*/
	    xerprn_(" ***", &c_n1, "JOB ABORT DUE TO FATAL ERROR.", &c__72, 
		    4L, 29L);
/*<          ENDIF >*/
	}
/*<          CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY) >*/
	//	xersve_(" ", " ", " ", &c_n1, &c__0, &c__0, &kdummy, 1L, 1L, 1L);
/*<          CALL XERHLT (' ') >*/
	xerhlt_(" ", 1L);
/*<       ELSE >*/
    } else {
/*<          CALL XERHLT (MESSG) >*/
	xerhlt_(messg, messg_len);
/*<       ENDIF >*/
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* xermsg_ */

/*<       SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP) >*/
/* Subroutine */ int xerprn_(char *prefix, integer *npref, char *messg, 
	integer *nwrap, ftnlen prefix_len, ftnlen messg_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer i_len(char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer i_indx(char *, char *, ftnlen, ftnlen), s_cmp(char *, char *, 
	    ftnlen, ftnlen);

    /* Local variables */
    integer i, n;
    char cbuff[148];
    integer lpref, nextc, lwrap, nunit;
    extern integer i1mach_(integer *);
    integer iu[5], lpiece, idelta, lenmsg;
    extern /* Subroutine */ int xgetua_(integer *, integer *);

/* ***BEGIN PROLOGUE  XERPRN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Print error messages processed by XERMSG. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XERPRN-A) */
/* ***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR */
/* ***AUTHOR  Fong, Kirby, (NMFECC at LLNL) */
/* ***DESCRIPTION */

/* This routine sends one or more lines to each of the (up to five) */
/* logical units to which error messages are to be sent.  This routine */
/* is called several times by XERMSG, sometimes with a single line to */
/* print and sometimes with a (potentially very long) message that may */
/* wrap around into multiple lines. */

/* PREFIX  Input argument of type CHARACTER.  This argument contains */
/*         characters to be put at the beginning of each line before */
/*         the body of the message.  No more than 16 characters of */
/*         PREFIX will be used. */

/* NPREF   Input argument of type INTEGER.  This argument is the number */
/*         of characters to use from PREFIX.  If it is negative, the */
/*         intrinsic function LEN is used to determine its length.  If */
/*         it is zero, PREFIX is not used.  If it exceeds 16 or if */
/*         LEN(PREFIX) exceeds 16, only the first 16 characters will be */
/*         used.  If NPREF is positive and the length of PREFIX is less */
/*         than NPREF, a copy of PREFIX extended with blanks to length */
/*         NPREF will be used. */

/* MESSG   Input argument of type CHARACTER.  This is the text of a */
/*         message to be printed.  If it is a long message, it will be */
/*         broken into pieces for printing on multiple lines.  Each line 
*/
/*         will start with the appropriate prefix and be followed by a */
/*         piece of the message.  NWRAP is the number of characters per */
/*         piece; that is, after each NWRAP characters, we break and */
/*         start a new line.  In addition the characters '$$' embedded */
/*         in MESSG are a sentinel for a new line.  The counting of */
/*         characters up to NWRAP starts over for each new line.  The */
/*         value of NWRAP typically used by XERMSG is 72 since many */
/*         older error messages in the SLATEC Library are laid out to */
/*         rely on wrap-around every 72 characters. */

/* NWRAP   Input argument of type INTEGER.  This gives the maximum size */
/*         piece into which to break MESSG for printing on multiple */
/*         lines.  An embedded '$$' ends a line, and the count restarts */
/*         at the following character.  If a line break does not occur */
/*         on a blank (it would split a word) that word is moved to the */
/*         next line.  Values of NWRAP less than 16 will be treated as */
/*         16.  Values of NWRAP greater than 132 will be treated as 132. 
*/
/*         The actual line length will be NPREF + NWRAP after NPREF has */
/*         been adjusted to fall between 0 and 16 and NWRAP has been */
/*         adjusted to fall between 16 and 132. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  I1MACH, XGETUA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   880621  DATE WRITTEN */
/*   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF */
/*           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK */
/*           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE */
/*           SLASH CHARACTER IN FORMAT STATEMENTS. */
/*   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO */
/*           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK 
*/
/*           LINES TO BE PRINTED. */
/*   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF */
/*           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH. */
/*   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH. */
/*   891214  Prologue converted to Version 4.0 format.  (WRB) */
/*   900510  Added code to break messages between words.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERPRN */
/*<       CHARACTER*(*) PREFIX, MESSG >*/
/*<       INTEGER NPREF, NWRAP >*/
/*<       CHARACTER*148 CBUFF >*/
/*<       INTEGER IU(5), NUNIT >*/
/*<       CHARACTER*2 NEWLIN >*/
/*<       PARAMETER (NEWLIN = '$$') >*/
/* ***FIRST EXECUTABLE STATEMENT  XERPRN */
/*<       CALL XGETUA(IU,NUNIT) >*/
    xgetua_(iu, &nunit);

/*       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD 
*/
/*       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD */
/*       ERROR MESSAGE UNIT. */

/*<       N = I1MACH(4) >*/
    n = i1mach_(&c__4);
/*<       DO 10 I=1,NUNIT >*/
    i__1 = nunit;
    for (i = 1; i <= i__1; ++i) {
/*<          IF (IU(I) .EQ. 0) IU(I) = N >*/
	if (iu[i - 1] == 0) {
	    iu[i - 1] = n;
	}
/*<    10 CONTINUE >*/
/* L10: */
    }

/*       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE 
*/
/*       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING 
*/
/*       THE REST OF THIS ROUTINE. */

/*<       IF ( NPREF .LT. 0 ) THEN >*/
    if (*npref < 0) {
/*<          LPREF = LEN(PREFIX) >*/
	lpref = i_len(prefix, prefix_len);
/*<       ELSE >*/
    } else {
/*<          LPREF = NPREF >*/
	lpref = *npref;
/*<       ENDIF >*/
    }
/*<       LPREF = MIN(16, LPREF) >*/
    lpref = min(16,lpref);
/*<       IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX >*/
    if (lpref != 0) {
	s_copy(cbuff, prefix, lpref, prefix_len);
    }

/*       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE 
*/
/*       TIME FROM MESSG TO PRINT ON ONE LINE. */

/*<       LWRAP = MAX(16, MIN(132, NWRAP)) >*/
/* Computing MAX */
    i__1 = 16, i__2 = min(132,*nwrap);
    lwrap = max(i__1,i__2);

/*       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS. */

/*<       LENMSG = LEN(MESSG) >*/
    lenmsg = i_len(messg, messg_len);
/*<       N = LENMSG >*/
    n = lenmsg;
/*<       DO 20 I=1,N >*/
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
/*<          IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30 >*/
	if (messg[lenmsg - 1] != ' ') {
	    goto L30;
	}
/*<          LENMSG = LENMSG - 1 >*/
	--lenmsg;
/*<    20 CONTINUE >*/
/* L20: */
    }
/*<    30 CONTINUE >*/
L30:

/*       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE. */

/*<       IF (LENMSG .EQ. 0) THEN >*/
    if (lenmsg == 0) {
/*<          CBUFF(LPREF+1:LPREF+1) = ' ' >*/
	i__1 = lpref;
	s_copy(cbuff + i__1, " ", lpref + 1 - i__1, 1L);
/*         DO 40 I=1,NUNIT */
/*            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1) */
/*   40    CONTINUE */
/*<          RETURN >*/
	return 0;
/*<       ENDIF >*/
    }

/*       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING */
/*       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL. */
/*       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT. */
/*       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED. */

/*       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE */
/*       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE 
*/
/*       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH */
/*       OF THE SECOND ARGUMENT. */

/*       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE */
/*       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER 
*/
/*       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT */
/*       POSITION NEXTC. */

/*       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE */
/*                       REMAINDER OF THE CHARACTER STRING.  LPIECE */
/*                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC, */
/*                       WHICHEVER IS LESS. */

/*       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC: */
/*                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE */
/*                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY */
/*                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION */
/*                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF */
/*                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE */
/*                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC */
/*                       SHOULD BE INCREMENTED BY 2. */

/*       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP. */

/*       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1 
*/
/*                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS */
/*                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ. 
*/
/*                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY */
/*                       AT THE END OF A LINE. */

/*<       NEXTC = 1 >*/
    nextc = 1;
/*<    50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN) >*/
L50:
    lpiece = i_indx(messg + (nextc - 1), "$$", lenmsg - (nextc - 1), 2L);
/*<       IF (LPIECE .EQ. 0) THEN >*/
    if (lpiece == 0) {

/*       THERE WAS NO NEW LINE SENTINEL FOUND. */

/*<          IDELTA = 0 >*/
	idelta = 0;
/*<          LPIECE = MIN(LWRAP, LENMSG+1-NEXTC) >*/
/* Computing MIN */
	i__1 = lwrap, i__2 = lenmsg + 1 - nextc;
	lpiece = min(i__1,i__2);
/*<          IF (LPIECE .LT. LENMSG+1-NEXTC) THEN >*/
	if (lpiece < lenmsg + 1 - nextc) {
/*<             DO 52 I=LPIECE+1,2,-1 >*/
	    for (i = lpiece + 1; i >= 2; --i) {
/*<                IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN >*/
		i__1 = nextc + i - 2;
		if (s_cmp(messg + i__1, " ", nextc + i - 1 - i__1, 1L) == 0) {
/*<                   LPIECE = I-1 >*/
		    lpiece = i - 1;
/*<                   IDELTA = 1 >*/
		    idelta = 1;
/*<                   GOTO 54 >*/
		    goto L54;
/*<                ENDIF >*/
		}
/*<    52       CONTINUE >*/
/* L52: */
	    }
/*<          ENDIF >*/
	}
/*<    54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1) >*/
L54:
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
/*<          NEXTC = NEXTC + LPIECE + IDELTA >*/
	nextc = nextc + lpiece + idelta;
/*<       ELSEIF (LPIECE .EQ. 1) THEN >*/
    } else if (lpiece == 1) {

/*       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1). */
/*       DON'T PRINT A BLANK LINE. */

/*<          NEXTC = NEXTC + 2 >*/
	nextc += 2;
/*<          GO TO 50 >*/
	goto L50;
/*<       ELSEIF (LPIECE .GT. LWRAP+1) THEN >*/
    } else if (lpiece > lwrap + 1) {

/*       LPIECE SHOULD BE SET DOWN TO LWRAP. */

/*<          IDELTA = 0 >*/
	idelta = 0;
/*<          LPIECE = LWRAP >*/
	lpiece = lwrap;
/*<          DO 56 I=LPIECE+1,2,-1 >*/
	for (i = lpiece + 1; i >= 2; --i) {
/*<             IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN >*/
	    i__1 = nextc + i - 2;
	    if (s_cmp(messg + i__1, " ", nextc + i - 1 - i__1, 1L) == 0) {
/*<                LPIECE = I-1 >*/
		lpiece = i - 1;
/*<                IDELTA = 1 >*/
		idelta = 1;
/*<                GOTO 58 >*/
		goto L58;
/*<             ENDIF >*/
	    }
/*<    56    CONTINUE >*/
/* L56: */
	}
/*<    58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1) >*/
L58:
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
/*<          NEXTC = NEXTC + LPIECE + IDELTA >*/
	nextc = nextc + lpiece + idelta;
/*<       ELSE >*/
    } else {

/*       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1. */
/*       WE SHOULD DECREMENT LPIECE BY ONE. */

/*<          LPIECE = LPIECE - 1 >*/
	--lpiece;
/*<          CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1) >*/
	i__1 = lpref;
	s_copy(cbuff + i__1, messg + (nextc - 1), lpref + lpiece - i__1, 
		nextc + lpiece - 1 - (nextc - 1));
/*<          NEXTC  = NEXTC + LPIECE + 2 >*/
	nextc = nextc + lpiece + 2;
/*<       ENDIF >*/
    }

/*       PRINT */

/*      DO 60 I=1,NUNIT */
/*         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE) */
/*   60 CONTINUE */

/*<       IF (NEXTC .LE. LENMSG) GO TO 50 >*/
    if (nextc <= lenmsg) {
	goto L50;
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* xerprn_ */

/*<       SUBROUTINE XGETUA (IUNITA, N) >*/
/* Subroutine */ int xgetua_(integer *iunita, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i, index;
    extern integer j4save_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XGETUA */
/* ***PURPOSE  Return unit number(s) to which error messages are being */
/*            sent. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (XGETUA-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        XGETUA may be called to determine the unit number or numbers */
/*        to which error messages are being sent. */
/*        These unit numbers may have been set by a call to XSETUN, */
/*        or a call to XSETUA, or may be a default value. */

/*     Description of Parameters */
/*      --Output-- */
/*        IUNIT - an array of one to five unit numbers, depending */
/*                on the value of N.  A value of zero refers to the */
/*                default unit, as defined by the I1MACH machine */
/*                constant routine.  Only IUNIT(1),...,IUNIT(N) are */
/*                defined by XGETUA.  The values of IUNIT(N+1),..., */
/*                IUNIT(5) are not defined (for N .LT. 5) or altered */
/*                in any way by XGETUA. */
/*        N     - the number of units to which copies of the */
/*                error messages are being sent.  N will be in the */
/*                range from 1 to 5. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  J4SAVE */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XGETUA */
/*<       DIMENSION IUNITA(5) >*/
/* ***FIRST EXECUTABLE STATEMENT  XGETUA */
/*<       N = J4SAVE(5,0,.FALSE.) >*/
    /* Parameter adjustments */
    --iunita;

    /* Function Body */
    *n = j4save_(&c__5, &c__0, (logical*)&c__0);
/*<       DO 30 I=1,N >*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/*<          INDEX = I+4 >*/
	index = i + 4;
/*<          IF (I.EQ.1) INDEX = 3 >*/
	if (i == 1) {
	    index = 3;
	}
/*<          IUNITA(I) = J4SAVE(INDEX,0,.FALSE.) >*/
	iunita[i] = j4save_(&index, &c__0, (logical*)&c__0);
/*<    30 CONTINUE >*/
/* L30: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* xgetua_ */






/* ------------------------------ ORIGINAL CODE FROM SLATEC LIBRARY -----------------------------*/

/*DECK DNLS1
      SUBROUTINE DNLS1 (FCN, IOPT, M, N, X, FVEC, FJAC, LDFJAC, FTOL,
     +   XTOL, GTOL, MAXFEV, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO,
     +   NFEV, NJEV, IPVT, QTF, WA1, WA2, WA3, WA4)
C***BEGIN PROLOGUE  DNLS1
C***PURPOSE  Minimize the sum of the squares of M nonlinear functions
C            in N variables by a modification of the Levenberg-Marquardt
C            algorithm.
C***LIBRARY   SLATEC
C***CATEGORY  K1B1A1, K1B1A2
C***TYPE      DOUBLE PRECISION (SNLS1-S, DNLS1-D)
C***KEYWORDS  LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
C             NONLINEAR LEAST SQUARES
C***AUTHOR  Hiebert, K. L., (SNLA)
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of DNLS1 is to minimize the sum of the squares of M
C       nonlinear functions in N variables by a modification of the
C       Levenberg-Marquardt algorithm.  The user must provide a subrou-
C       tine which calculates the functions.  The user has the option
C       of how the Jacobian will be supplied.  The user can supply the
C       full Jacobian, or the rows of the Jacobian (to avoid storing
C       the full Jacobian), or let the code approximate the Jacobian by
C       forward-differencing.   This code is the combination of the
C       MINPACK codes (Argonne) LMDER, LMDIF, and LMSTR.
C
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE DNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
C      *                 GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO
C      *                 ,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
C       INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
C       INTEGER IPVT(N)
C       DOUBLE PRECISION FTOL,XTOL,GTOL,EPSFCN,FACTOR
C       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),DIAG(N),QTF(N),
C      *     WA1(N),WA2(N),WA3(N),WA4(M)
C
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to DNLS1 and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from DNLS1.
C
C      FCN is the name of the user-supplied subroutine which calculate
C         the functions.  If the user wants to supply the Jacobian
C         (IOPT=2 or 3), then FCN must be written to calculate the
C         Jacobian, as well as the functions.  See the explanation
C         of the IOPT argument below.
C         If the user wants the iterates printed (NPRINT positive), then
C         FCN must do the printing.  See the explanation of NPRINT
C         below.  FCN must be declared in an EXTERNAL statement in the
C         calling program and should be written as follows.
C
C
C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C         INTEGER IFLAG,LDFJAC,M,N
C         DOUBLE PRECISION X(N),FVEC(M)
C         ----------
C         FJAC and LDFJAC may be ignored       , if IOPT=1.
C         DOUBLE PRECISION FJAC(LDFJAC,N)      , if IOPT=2.
C         DOUBLE PRECISION FJAC(N)             , if IOPT=3.
C         ----------
C           If IFLAG=0, the values in X and FVEC are available
C           for printing.  See the explanation of NPRINT below.
C           IFLAG will never be zero unless NPRINT is positive.
C           The values of X and FVEC must not be changed.
C         RETURN
C         ----------
C           If IFLAG=1, calculate the functions at X and return
C           this vector in FVEC.
C         RETURN
C         ----------
C           If IFLAG=2, calculate the full Jacobian at X and return
C           this matrix in FJAC.  Note that IFLAG will never be 2 unless
C           IOPT=2.  FVEC contains the function values at X and must
C           not be altered.  FJAC(I,J) must be set to the derivative
C           of FVEC(I) with respect to X(J).
C         RETURN
C         ----------
C           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
C           and return this vector in FJAC.  Note that IFLAG will
C           never be 3 unless IOPT=3.  FVEC contains the function
C           values at X and must not be altered.  FJAC(J) must be
C           set to the derivative of FVEC(LDFJAC) with respect to X(J).
C         RETURN
C         ----------
C         END
C
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of DNLS1.  In this case, set
C         IFLAG to a negative integer.
C
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=2 or 3, then the user must supply the
C         Jacobian, as well as the function values, through the
C         subroutine FCN.  If IOPT=2, the user supplies the full
C         Jacobian with one call to FCN.  If IOPT=3, the user supplies
C         one row of the Jacobian with each call.  (In this manner,
C         storage can be saved because the full Jacobian is not stored.)
C         If IOPT=1, the code will approximate the Jacobian by forward
C         differencing.
C
C       M is a positive integer input variable set to the number of
C         functions.
C
C       N is a positive integer input variable set to the number of
C         variables.  N must not exceed M.
C
C       X is an array of length N.  On input, X must contain an initial
C         estimate of the solution vector.  On output, X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length M which contains the functions
C         evaluated at the output X.
C
C       FJAC is an output array.  For IOPT=1 and 2, FJAC is an M by N
C         array.  For IOPT=3, FJAC is an N by N array.  The upper N by N
C         submatrix of FJAC contains an upper triangular matrix R with
C         diagonal elements of nonincreasing magnitude such that
C
C                T     T           T
C               P *(JAC *JAC)*P = R *R,
C
C         where P is a permutation matrix and JAC is the final calcu-
C         lated Jacobian.  Column J of P is column IPVT(J) (see below)
C         of the identity matrix.  The lower part of FJAC contains
C         information generated during the computation of R.
C
C       LDFJAC is a positive integer input variable which specifies
C         the leading dimension of the array FJAC.  For IOPT=1 and 2,
C         LDFJAC must not be less than M.  For IOPT=3, LDFJAC must not
C         be less than N.
C
C       FTOL is a non-negative input variable.  Termination occurs when
C         both the actual and predicted relative reductions in the sum
C         of squares are at most FTOL.  Therefore, FTOL measures the
C         relative error desired in the sum of squares.  Section 4 con-
C         tains more details about FTOL.
C
C       XTOL is a non-negative input variable.  Termination occurs when
C         the relative error between two consecutive iterates is at most
C         XTOL.  Therefore, XTOL measures the relative error desired in
C         the approximate solution.  Section 4 contains more details
C         about XTOL.
C
C       GTOL is a non-negative input variable.  Termination occurs when
C         the cosine of the angle between FVEC and any column of the
C         Jacobian is at most GTOL in absolute value.  Therefore, GTOL
C         measures the orthogonality desired between the function vector
C         and the columns of the Jacobian.  Section 4 contains more
C         details about GTOL.
C
C       MAXFEV is a positive integer input variable.  Termination occurs
C         when the number of calls to FCN to evaluate the functions
C         has reached MAXFEV.
C
C       EPSFCN is an input variable used in determining a suitable step
C         for the forward-difference approximation.  This approximation
C         assumes that the relative errors in the functions are of the
C         order of EPSFCN.  If EPSFCN is less than the machine preci-
C         sion, it is assumed that the relative errors in the functions
C         are of the order of the machine precision.  If IOPT=2 or 3,
C         then EPSFCN can be ignored (treat it as a dummy argument).
C
C       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
C         internally set.  If MODE = 2, DIAG must contain positive
C         entries that serve as implicit (multiplicative) scale factors
C         for the variables.
C
C       MODE is an integer input variable.  If MODE = 1, the variables
C         will be scaled internally.  If MODE = 2, the scaling is speci-
C         fied by the input DIAG.  Other values of MODE are equivalent
C         to MODE = 1.
C
C       FACTOR is a positive input variable used in determining the ini-
C         tial step bound.  This bound is set to the product of FACTOR
C         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR
C         itself.  In most cases FACTOR should lie in the interval
C         (.1,100.).  100. is a generally recommended value.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. Appropriate
C         print statements must be added to FCN (see example) and
C         FVEC should not be altered.  If NPRINT is not positive, no
C         special calls to FCN with IFLAG = 0 are made.
C
C       INFO is an integer output variable.  If the user has terminated
C        execution, INFO is set to the (negative) value of IFLAG.  See
C        description of FCN and JAC. Otherwise, INFO is set as follows
C
C         INFO = 0  improper input parameters.
C
C         INFO = 1  both actual and predicted relative reductions in the
C                   sum of squares are at most FTOL.
C
C         INFO = 2  relative error between two consecutive iterates is
C                   at most XTOL.
C
C         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
C
C         INFO = 4  the cosine of the angle between FVEC and any column
C                   of the Jacobian is at most GTOL in absolute value.
C
C         INFO = 5  number of calls to FCN for function evaluation
C                   has reached MAXFEV.
C
C         INFO = 6  FTOL is too small.  No further reduction in the sum
C                   of squares is possible.
C
C         INFO = 7  XTOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 8  GTOL is too small.  FVEC is orthogonal to the
C                   columns of the Jacobian to machine precision.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       NFEV is an integer output variable set to the number of calls to
C         FCN for function evaluation.
C
C       NJEV is an integer output variable set to the number of
C         evaluations of the full Jacobian.  If IOPT=2, only one call to
C         FCN is required for each evaluation of the full Jacobian.
C         If IOPT=3, the M calls to FCN are required.
C         If IOPT=1, then NJEV is set to zero.
C
C       IPVT is an integer output array of length N.  IPVT defines a
C         permutation matrix P such that JAC*P = Q*R, where JAC is the
C         final calculated Jacobian, Q is orthogonal (not stored), and R
C         is upper triangular with diagonal elements of nonincreasing
C         magnitude.  Column J of P is column IPVT(J) of the identity
C         matrix.
C
C       QTF is an output array of length N which contains the first N
C         elements of the vector (Q transpose)*FVEC.
C
C       WA1, WA2, and WA3 are work arrays of length N.
C
C       WA4 is a work array of length M.
C
C
C 4. Successful Completion.
C
C       The accuracy of DNLS1 is controlled by the convergence parame-
C       ters FTOL, XTOL, and GTOL.  These parameters are used in tests
C       which make three types of comparisons between the approximation
C       X and a solution XSOL.  DNLS1 terminates when any of the tests
C       is satisfied.  If any of the convergence parameters is less than
C       the machine precision (as defined by the function R1MACH(4)),
C       then DNLS1 only attempts to satisfy the test defined by the
C       machine precision.  Further progress is not usually possible.
C
C       The tests assume that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently.  If these conditions
C       are not satisfied, then DNLS1 may incorrectly indicate conver-
C       gence.  If the Jacobian is coded correctly or IOPT=1,
C       then the validity of the answer can be checked, for example, by
C       rerunning DNLS1 with tighter tolerances.
C
C       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
C         of a vector Z, then this test attempts to guarantee that
C
C               ENORM(FVEC) .LE. (1+FTOL)*ENORM(FVECS),
C
C         where FVECS denotes the functions evaluated at XSOL.  If this
C         condition is satisfied with FTOL = 10**(-K), then the final
C         residual norm ENORM(FVEC) has K significant decimal digits and
C         INFO is set to 1 (or to 3 if the second test is also satis-
C         fied).  Unless high precision solutions are required, the
C         recommended value for FTOL is the square root of the machine
C         precision.
C
C       Second Convergence Test.  If D is the diagonal matrix whose
C         entries are defined by the array DIAG, then this test attempts
C         to guarantee that
C
C               ENORM(D*(X-XSOL)) .LE. XTOL*ENORM(D*XSOL).
C
C         If this condition is satisfied with XTOL = 10**(-K), then the
C         larger components of D*X have K significant decimal digits and
C         INFO is set to 2 (or to 3 if the first test is also satis-
C         fied).  There is a danger that the smaller components of D*X
C         may have large relative errors, but if MODE = 1, then the
C         accuracy of the components of X is usually related to their
C         sensitivity.  Unless high precision solutions are required,
C         the recommended value for XTOL is the square root of the
C         machine precision.
C
C       Third Convergence Test.  This test is satisfied when the cosine
C         of the angle between FVEC and any column of the Jacobian at X
C         is at most GTOL in absolute value.  There is no clear rela-
C         tionship between this test and the accuracy of DNLS1, and
C         furthermore, the test is equally well satisfied at other crit-
C         ical points, namely maximizers and saddle points.  Therefore,
C         termination caused by this test (INFO = 4) should be examined
C         carefully.  The recommended value for GTOL is zero.
C
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of DNLS1 can be due to improper input
C       parameters, arithmetic interrupts, or an excessive number of
C       function evaluations.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1
C         or IOPT .GT. 3, or N .LE. 0, or M .LT. N, or for IOPT=1 or 2
C         LDFJAC .LT. M, or for IOPT=3 LDFJAC .LT. N, or FTOL .LT. 0.E0,
C         or XTOL .LT. 0.E0, or GTOL .LT. 0.E0, or MAXFEV .LE. 0, or
C         FACTOR .LE. 0.E0.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by DNLS1.  In this
C         case, it may be possible to remedy the situation by rerunning
C         DNLS1 with a smaller value of FACTOR.
C
C       Excessive Number of Function Evaluations.  A reasonable value
C         for MAXFEV is 100*(N+1) for IOPT=2 or 3 and 200*(N+1) for
C         IOPT=1.  If the number of calls to FCN reaches MAXFEV, then
C         this indicates that the routine is converging very slowly
C         as measured by the progress of FVEC, and INFO is set to 5.
C         In this case, it may be helpful to restart DNLS1 with MODE
C         set to 1.
C
C
C 6. Characteristics of the Algorithm.
C
C       DNLS1 is a modification of the Levenberg-Marquardt algorithm.
C       Two of its main characteristics involve the proper use of
C       implicitly scaled variables (if MODE = 1) and an optimal choice
C       for the correction.  The use of implicitly scaled variables
C       achieves scale invariance of DNLS1 and limits the size of the
C       correction in any direction where the functions are changing
C       rapidly.  The optimal choice of the correction guarantees (under
C       reasonable conditions) global convergence from starting points
C       far from the solution and a fast rate of convergence for
C       problems with small residuals.
C
C       Timing.  The time required by DNLS1 to solve a given problem
C         depends on M and N, the behavior of the functions, the accu-
C         racy requested, and the starting point.  The number of arith-
C         metic operations needed by DNLS1 is about N**3 to process each
C         evaluation of the functions (call to FCN) and to process each
C         evaluation of the Jacobian it takes M*N**2 for IOPT=2 (one
C         call to FCN), M*N**2 for IOPT=1 (N calls to FCN) and
C         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
C         can be evaluated quickly, the timing of DNLS1 will be
C         strongly influenced by the time spent in FCN.
C
C       Storage.  DNLS1 requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
C         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
C         locations and N integer storage locations, in addition to
C         the storage required by the program.  There are no internally
C         declared storage arrays.
C
C *Long Description:
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), and X(3)
C       which provide the best fit (in the least squares sense) of
C
C             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15
C
C       to the data
C
C             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
C                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
C
C       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The
C       I-th component of FVEC is thus defined by
C
C             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))).
C
C       **********
C
C       PROGRAM TEST
C C
C C     Driver for DNLS1 example.
C C
C       INTEGER J,IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,
C      *        NWRITE
C       INTEGER IPVT(3)
C       DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,FNORM,EPSFCN
C       DOUBLE PRECISION X(3),FVEC(15),FJAC(15,3),DIAG(3),QTF(3),
C      *     WA1(3),WA2(3),WA3(3),WA4(15)
C       DOUBLE PRECISION DENORM,D1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 1
C       M = 15
C       N = 3
C C
C C     The following starting values provide a rough fit.
C C
C       X(1) = 1.E0
C       X(2) = 1.E0
C       X(3) = 1.E0
C C
C       LDFJAC = 15
C C
C C     Set FTOL and XTOL to the square root of the machine precision
C C     and GTOL to zero.  Unless high precision solutions are
C C     required, these are the recommended settings.
C C
C       FTOL = SQRT(R1MACH(4))
C       XTOL = SQRT(R1MACH(4))
C       GTOL = 0.E0
C C
C       MAXFEV = 400
C       EPSFCN = 0.0
C       MODE = 1
C       FACTOR = 1.E2
C       NPRINT = 0
C C
C       CALL DNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
C      *           GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,
C      *           INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
C       FNORM = ENORM(M,FVEC)
C       WRITE (NWRITE,1000) FNORM,NFEV,NJEV,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
C      *        5X,' NUMBER OF JACOBIAN EVALUATIONS',I10 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7)
C       END
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM)
C C     This is the form of the FCN routine if IOPT=1,
C C     that is, if the user does not calculate the Jacobian.
C       INTEGER I,M,N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(M),Y(15)
C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C       END
C
C
C       Results obtained with different compilers or machines
C       may be slightly different.
C
C       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01
C
C       NUMBER OF FUNCTION EVALUATIONS        25
C
C       NUMBER OF JACOBIAN EVALUATIONS         0
C
C       EXIT PARAMETER                         1
C
C       FINAL APPROXIMATE SOLUTION
C
C        0.8241058E-01  0.1133037E+01  0.2343695E+01
C
C
C       For IOPT=2, FCN would be modified as follows to also
C       calculate the full Jacobian when IFLAG=2.
C
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C C
C C     This is the form of the FCN routine if IOPT=2,
C C     that is, if the user calculates the full Jacobian.
C C
C       INTEGER I,LDFJAC,M,N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),Y(15)
C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       IF(IFLAG.NE.1) GO TO 20
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C C
C C     Below, calculate the full Jacobian.
C C
C    20    CONTINUE
C C
C       DO 30 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
C          FJAC(I,1) = -1.E0
C          FJAC(I,2) = TMP1*TMP2/TMP4
C          FJAC(I,3) = TMP1*TMP3/TMP4
C    30    CONTINUE
C       RETURN
C       END
C
C
C       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3),
C         LDFJAC would be set to 3, and FCN would be written as
C         follows to calculate a row of the Jacobian when IFLAG=3.
C
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C C     This is the form of the FCN routine if IOPT=3,
C C     that is, if the user calculates the Jacobian row by row.
C       INTEGER I,M,N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(M),FJAC(N),Y(15)
C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       IF( IFLAG.NE.1) GO TO 20
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C C
C C     Below, calculate the LDFJAC-th row of the Jacobian.
C C
C    20 CONTINUE
C
C       I = LDFJAC
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
C          FJAC(1) = -1.E0
C          FJAC(2) = TMP1*TMP2/TMP4
C          FJAC(3) = TMP1*TMP3/TMP4
C       RETURN
C       END
C
C***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm:
C                 implementation and theory.  In Numerical Analysis
C                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
C                 Watson, Editor), Lecture Notes in Mathematics 630,
C                 Springer-Verlag, 1978.
C***ROUTINES CALLED  D1MACH, DCKDER, DENORM, DFDJC3, DMPAR, DQRFAC,
C                    DWUPDT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920205  Corrected XERN1 declaration.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNLS1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
      INTEGER IJUNK,NROW,IPVT(*)
      DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,EPSFCN
      DOUBLE PRECISION X(*),FVEC(*),FJAC(LDFJAC,*),DIAG(*),QTF(*),
     1     WA1(*),WA2(*),WA3(*),WA4(*)
      LOGICAL SING
      EXTERNAL FCN
      INTEGER I,IFLAG,ITER,J,L,MODECH
      DOUBLE PRECISION ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM,
     1     ONE,PAR,PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO,SUM,TEMP,
     2     TEMP1,TEMP2,XNORM,ZERO
      DOUBLE PRECISION D1MACH,DENORM,ERR,CHKLIM
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3
      SAVE CHKLIM, ONE, P1, P5, P25, P75, P0001, ZERO
C
      DATA CHKLIM/.1D0/
      DATA ONE,P1,P5,P25,P75,P0001,ZERO
     1     /1.0D0,1.0D-1,5.0D-1,2.5D-1,7.5D-1,1.0D-4,0.0D0/
C***FIRST EXECUTABLE STATEMENT  DNLS1
      EPSMCH = D1MACH(4)
C
      INFO = 0
      IFLAG = 0
      NFEV = 0
      NJEV = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3 .OR. N .LE. 0 .OR.
     1    M .LT. N .OR. LDFJAC .LT. N .OR. FTOL .LT. ZERO
     2    .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO
     3    .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 300
      IF (IOPT .LT. 3 .AND. LDFJAC .LT. M) GO TO 300
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 300
   10    CONTINUE
   20 CONTINUE
C
C     EVALUATE THE FUNCTION AT THE STARTING POINT
C     AND CALCULATE ITS NORM.
C
      IFLAG = 1
      IJUNK = 1
      CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      NFEV = 1
      IF (IFLAG .LT. 0) GO TO 300
      FNORM = DENORM(M,FVEC)
C
C     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
C
      PAR = ZERO
      ITER = 1
C
C     BEGINNING OF THE OUTER LOOP.
C
   30 CONTINUE
C
C        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
C
         IF (NPRINT .LE. 0) GO TO 40
         IFLAG = 0
         IF (MOD(ITER-1,NPRINT) .EQ. 0)
     1      CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
         IF (IFLAG .LT. 0) GO TO 300
   40    CONTINUE
C
C        CALCULATE THE JACOBIAN MATRIX.
C
      IF (IOPT .EQ. 3) GO TO 475
C
C     STORE THE FULL JACOBIAN USING M*N STORAGE
C
      IF (IOPT .EQ. 1) GO TO 410
C
C     THE USER SUPPLIES THE JACOBIAN
C
         IFLAG = 2
         CALL FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
         NJEV = NJEV + 1
C
C             ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN
C
         IF (ITER .LE. 1) THEN
            IF (IFLAG .LT. 0) GO TO 300
C
C           GET THE INCREMENTED X-VALUES INTO WA1(*).
C
            MODECH = 1
            CALL DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
C
C           EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT IN WA4(*).
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,LDFJAC)
            NFEV = NFEV + 1
            IF(IFLAG .LT. 0) GO TO 300
            DO 350 I = 1, M
               MODECH = 2
               CALL DCKDER(1,N,X,FVEC(I),FJAC(I,1),LDFJAC,WA1,
     1              WA4(I),MODECH,ERR)
               IF (ERR .LT. CHKLIM) THEN
                  WRITE (XERN1, '(I8)') I
                  WRITE (XERN3, '(1PE15.6)') ERR
                  CALL XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF ' //
     *               'FUNCTION ' // XERN1 // ' MAY BE WRONG, ERR = ' //
     *               XERN3 // ' TOO CLOSE TO 0.', 7, 0)
               ENDIF
  350       CONTINUE
         ENDIF
C
         GO TO 420
C
C     THE CODE APPROXIMATES THE JACOBIAN
C
410      IFLAG = 1
         CALL DFDJC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)
         NFEV = NFEV + N
  420    IF (IFLAG .LT. 0) GO TO 300
C
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
         CALL DQRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
C
C        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
C        QTF.
C
         DO 430 I = 1, M
            WA4(I) = FVEC(I)
  430         CONTINUE
         DO 470 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 460
            SUM = ZERO
            DO 440 I = J, M
               SUM = SUM + FJAC(I,J)*WA4(I)
  440          CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 450 I = J, M
               WA4(I) = WA4(I) + FJAC(I,J)*TEMP
  450          CONTINUE
  460       CONTINUE
            FJAC(J,J) = WA1(J)
            QTF(J) = WA4(J)
  470       CONTINUE
         GO TO 560
C
C        ACCUMULATE THE JACOBIAN BY ROWS IN ORDER TO SAVE STORAGE.
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX
C        CALCULATED ONE ROW AT A TIME, WHILE SIMULTANEOUSLY
C        FORMING (Q TRANSPOSE)*FVEC AND STORING THE FIRST
C        N COMPONENTS IN QTF.
C
  475    DO 490 J = 1, N
            QTF(J) = ZERO
            DO 480 I = 1, N
               FJAC(I,J) = ZERO
  480          CONTINUE
  490        CONTINUE
         DO 500 I = 1, M
            NROW = I
            IFLAG = 3
            CALL FCN(IFLAG,M,N,X,FVEC,WA3,NROW)
            IF (IFLAG .LT. 0) GO TO 300
C
C            ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN.
C
            IF(ITER .GT. 1) GO TO 498
C
C            GET THE INCREMENTED X-VALUES INTO WA1(*).
C
            MODECH = 1
            CALL DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
C
C            EVALUATE AT INCREMENTED VALUES, IF NOT ALREADY EVALUATED.
C
            IF(I .NE. 1) GO TO 495
C
C            EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT INTO WA4(*).
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,NROW)
            NFEV = NFEV + 1
            IF(IFLAG .LT. 0) GO TO 300
495         CONTINUE
            MODECH = 2
            CALL DCKDER(1,N,X,FVEC(I),WA3,1,WA1,WA4(I),MODECH,ERR)
            IF (ERR .LT. CHKLIM) THEN
               WRITE (XERN1, '(I8)') I
               WRITE (XERN3, '(1PE15.6)') ERR
               CALL XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF FUNCTION '
     *            // XERN1 // ' MAY BE WRONG, ERR = ' // XERN3 //
     *            ' TOO CLOSE TO 0.', 7, 0)
            ENDIF
498         CONTINUE
C
            TEMP = FVEC(I)
            CALL DWUPDT(N,FJAC,LDFJAC,WA3,QTF,TEMP,WA1,WA2)
  500       CONTINUE
         NJEV = NJEV + 1
C
C        IF THE JACOBIAN IS RANK DEFICIENT, CALL DQRFAC TO
C        REORDER ITS COLUMNS AND UPDATE THE COMPONENTS OF QTF.
C
         SING = .FALSE.
         DO 510 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) SING = .TRUE.
            IPVT(J) = J
            WA2(J) = DENORM(J,FJAC(1,J))
  510       CONTINUE
         IF (.NOT.SING) GO TO 560
         CALL DQRFAC(N,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
         DO 550 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 540
            SUM = ZERO
            DO 520 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
  520         CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 530 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  530          CONTINUE
  540       CONTINUE
            FJAC(J,J) = WA1(J)
  550       CONTINUE
  560    CONTINUE
C
C        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
         IF (ITER .NE. 1) GO TO 80
         IF (MODE .EQ. 2) GO TO 60
         DO 50 J = 1, N
            DIAG(J) = WA2(J)
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   50       CONTINUE
   60    CONTINUE
C
C        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
C        AND INITIALIZE THE STEP BOUND DELTA.
C
         DO 70 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   70       CONTINUE
         XNORM = DENORM(N,WA3)
         DELTA = FACTOR*XNORM
         IF (DELTA .EQ. ZERO) DELTA = FACTOR
   80    CONTINUE
C
C        COMPUTE THE NORM OF THE SCALED GRADIENT.
C
         GNORM = ZERO
         IF (FNORM .EQ. ZERO) GO TO 170
         DO 160 J = 1, N
            L = IPVT(J)
            IF (WA2(L) .EQ. ZERO) GO TO 150
            SUM = ZERO
            DO 140 I = 1, J
               SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
  140          CONTINUE
            GNORM = MAX(GNORM,ABS(SUM/WA2(L)))
  150       CONTINUE
  160       CONTINUE
  170    CONTINUE
C
C        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
C
         IF (GNORM .LE. GTOL) INFO = 4
         IF (INFO .NE. 0) GO TO 300
C
C        RESCALE IF NECESSARY.
C
         IF (MODE .EQ. 2) GO TO 190
         DO 180 J = 1, N
            DIAG(J) = MAX(DIAG(J),WA2(J))
  180       CONTINUE
  190    CONTINUE
C
C        BEGINNING OF THE INNER LOOP.
C
  200    CONTINUE
C
C           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
C
            CALL DMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2,
     1                 WA3,WA4)
C
C           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
C
            DO 210 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  210          CONTINUE
            PNORM = DENORM(N,WA3)
C
C           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
            IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM)
C
C           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA2,WA4,FJAC,IJUNK)
            NFEV = NFEV + 1
            IF (IFLAG .LT. 0) GO TO 300
            FNORM1 = DENORM(M,WA4)
C
C           COMPUTE THE SCALED ACTUAL REDUCTION.
C
            ACTRED = -ONE
            IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C           COMPUTE THE SCALED PREDICTED REDUCTION AND
C           THE SCALED DIRECTIONAL DERIVATIVE.
C
            DO 230 J = 1, N
               WA3(J) = ZERO
               L = IPVT(J)
               TEMP = WA1(L)
               DO 220 I = 1, J
                  WA3(I) = WA3(I) + FJAC(I,J)*TEMP
  220             CONTINUE
  230          CONTINUE
            TEMP1 = DENORM(N,WA3)/FNORM
            TEMP2 = (SQRT(PAR)*PNORM)/FNORM
            PRERED = TEMP1**2 + TEMP2**2/P5
            DIRDER = -(TEMP1**2 + TEMP2**2)
C
C           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C           REDUCTION.
C
            RATIO = ZERO
            IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED
C
C           UPDATE THE STEP BOUND.
C
            IF (RATIO .GT. P25) GO TO 240
               IF (ACTRED .GE. ZERO) TEMP = P5
               IF (ACTRED .LT. ZERO)
     1            TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
               IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1
               DELTA = TEMP*MIN(DELTA,PNORM/P1)
               PAR = PAR/TEMP
               GO TO 260
  240       CONTINUE
               IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250
               DELTA = PNORM/P5
               PAR = P5*PAR
  250          CONTINUE
  260       CONTINUE
C
C           TEST FOR SUCCESSFUL ITERATION.
C
            IF (RATIO .LT. P0001) GO TO 290
C
C           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
            DO 270 J = 1, N
               X(J) = WA2(J)
               WA2(J) = DIAG(J)*X(J)
  270          CONTINUE
            DO 280 I = 1, M
               FVEC(I) = WA4(I)
  280          CONTINUE
            XNORM = DENORM(N,WA2)
            FNORM = FNORM1
            ITER = ITER + 1
  290       CONTINUE
C
C           TESTS FOR CONVERGENCE.
C
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     1          .AND. P5*RATIO .LE. ONE) INFO = 1
            IF (DELTA .LE. XTOL*XNORM) INFO = 2
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     1          .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3
            IF (INFO .NE. 0) GO TO 300
C
C           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C
            IF (NFEV .GE. MAXFEV) INFO = 5
            IF (ABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH
     1          .AND. P5*RATIO .LE. ONE) INFO = 6
            IF (DELTA .LE. EPSMCH*XNORM) INFO = 7
            IF (GNORM .LE. EPSMCH) INFO = 8
            IF (INFO .NE. 0) GO TO 300
C
C           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
C
            IF (RATIO .LT. P0001) GO TO 200
C
C        END OF THE OUTER LOOP.
C
         GO TO 30
  300 CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'DNLS1',
     +   'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DNLS1',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO .EQ. 4) CALL XERMSG ('SLATEC', 'DNLS1',
     +   'THIRD CONVERGENCE CONDITION, CHECK RESULTS BEFORE ACCEPTING.',
     +   1, 1)
      IF (INFO .EQ. 5) CALL XERMSG ('SLATEC', 'DNLS1',
     +   'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
      IF (INFO .GE. 6) CALL XERMSG ('SLATEC', 'DNLS1',
     +   'TOLERANCES TOO SMALL, NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
      RETURN
C
C     LAST CARD OF SUBROUTINE DNLS1.
C
      END
*/
