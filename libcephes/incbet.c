/*							incbet.c
 *
 *	Incomplete beta integral
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, md_incbet();
 *
 * y = md_incbet( a, b, x );
 *
 *
 * DESCRIPTION:
 *
 * Returns incomplete beta integral of the arguments, evaluated
 * from zero to x.  The function is defined as
 *
 *                  x
 *     -            -
 *    | (a+b)      | |  a-1     b-1
 *  -----------    |   t   (1-t)   dt.
 *   -     -     | |
 *  | (a) | (b)   -
 *                 0
 *
 * The domain of definition is 0 <= x <= 1.  In this
 * implementation a and b are restricted to positive values.
 * The integral from x to 1 may be obtained by the symmetry
 * relation
 *
 *    1 - md_incbet( a, b, x )  =  md_incbet( b, a, 1-x ).
 *
 * The integral is evaluated by a continued fraction expansion
 * or, when b*x is small, by a power series.
 *
 * ACCURACY:
 *
 * Tested at uniformly distributed random points (a,b,x) with a and b
 * in "domain" and x between 0 and 1.
 *                                        Relative error
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,5         10000       6.9e-15     4.5e-16
 *    IEEE      0,85       250000       2.2e-13     1.7e-14
 *    IEEE      0,1000      30000       5.3e-12     6.3e-13
 *    IEEE      0,10000    250000       9.3e-11     7.1e-12
 *    IEEE      0,100000    10000       8.7e-10     4.8e-11
 * Outputs smaller than the IEEE gradual underflow threshold
 * were excluded from these statistics.
 *
 * ERROR MESSAGES:
 *   message         condition      value returned
 * incbet domain      x<0, x>1          0.0
 * incbet underflow                     0.0
 */


/*
Cephes Math Library, Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*/

#include "mconf.h"

#ifdef DEC
#define MAXGAM 34.84425627277176174
#else
#define MAXGAM 171.624376956302725
#endif

extern double MACHEP, MINLOG, MAXLOG;
#ifdef ANSIPROT
extern double md_gamma ( double );
extern double md_lgam ( double );
extern double md_exp ( double );
extern double md_log ( double );
extern double md_pow ( double, double );
extern double md_fabs ( double );
static double md_incbcf(double, double, double);
static double md_incbd(double, double, double);
static double md_pseries(double, double, double);
#else
double md_gamma(), md_lgam(), md_exp(), md_log(), md_pow(), md_fabs();
static double md_incbcf(), md_incbd(), md_pseries();
#endif

static double big = 4.503599627370496e15;
static double biginv =  2.22044604925031308085e-16;


double md_incbet( aa, bb, xx )
double aa, bb, xx;
{
double a, b, t, x, xc, w, y;
int flag;

if( aa <= 0.0 || bb <= 0.0 )
	goto domerr;

if( (xx <= 0.0) || ( xx >= 1.0) )
	{
	if( xx == 0.0 )
		return(0.0);
	if( xx == 1.0 )
		return( 1.0 );
domerr:
	mtherr( "incbet", DOMAIN );
	return( 0.0 );
	}

flag = 0;
if( (bb * xx) <= 1.0 && xx <= 0.95)
	{
	t = md_pseries(aa, bb, xx);
		goto done;
	}

w = 1.0 - xx;

/* Reverse a and b if x is greater than the mean. */
if( xx > (aa/(aa+bb)) )
	{
	flag = 1;
	a = bb;
	b = aa;
	xc = xx;
	x = w;
	}
else
	{
	a = aa;
	b = bb;
	xc = w;
	x = xx;
	}

if( flag == 1 && (b * x) <= 1.0 && x <= 0.95)
	{
	t = md_pseries(a, b, x);
	goto done;
	}

/* Choose expansion for better convergence. */
y = x * (a+b-2.0) - (a-1.0);
if( y < 0.0 )
	w = md_incbcf( a, b, x );
else
	w = md_incbd( a, b, x ) / xc;

/* Multiply w by the factor
     a      b   _             _     _
    x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

y = a * md_log(x);
t = b * md_log(xc);
if( (a+b) < MAXGAM && md_fabs(y) < MAXLOG && md_fabs(t) < MAXLOG )
	{
	t = md_pow(xc,b);
	t *= md_pow(x,a);
	t /= a;
	t *= w;
	t *= md_gamma(a+b) / (md_gamma(a) * md_gamma(b));
	goto done;
	}
/* Resort to logarithms.  */
y += t + md_lgam(a+b) - md_lgam(a) - md_lgam(b);
y += md_log(w/a);
if( y < MINLOG )
	t = 0.0;
else
	t = md_exp(y);

done:

if( flag == 1 )
	{
	if( t <= MACHEP )
		t = 1.0 - MACHEP;
	else
		t = 1.0 - t;
	}
return( t );
}

/* Continued fraction expansion #1
 * for incomplete beta integral
 */

static double md_incbcf( a, b, x )
double a, b, x;
{
double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
double k1, k2, k3, k4, k5, k6, k7, k8;
double r, t, ans, thresh;
int n;

k1 = a;
k2 = a + b;
k3 = a;
k4 = a + 1.0;
k5 = 1.0;
k6 = b - 1.0;
k7 = k4;
k8 = a + 2.0;

pkm2 = 0.0;
qkm2 = 1.0;
pkm1 = 1.0;
qkm1 = 1.0;
ans = 1.0;
r = 1.0;
n = 0;
thresh = 3.0 * MACHEP;
do
	{
	
	xk = -( x * k1 * k2 )/( k3 * k4 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	xk = ( x * k5 * k6 )/( k7 * k8 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	if( qk != 0 )
		r = pk/qk;
	if( r != 0 )
		{
		t = md_fabs( (ans - r)/r );
		ans = r;
		}
	else
		t = 1.0;

	if( t < thresh )
		goto cdone;

	k1 += 1.0;
	k2 += 1.0;
	k3 += 2.0;
	k4 += 2.0;
	k5 += 1.0;
	k6 -= 1.0;
	k7 += 2.0;
	k8 += 2.0;

	if( (md_fabs(qk) + md_fabs(pk)) > big )
		{
		pkm2 *= biginv;
		pkm1 *= biginv;
		qkm2 *= biginv;
		qkm1 *= biginv;
		}
	if( (md_fabs(qk) < biginv) || (md_fabs(pk) < biginv) )
		{
		pkm2 *= big;
		pkm1 *= big;
		qkm2 *= big;
		qkm1 *= big;
		}
	}
while( ++n < 300 );

cdone:
return(ans);
}


/* Continued fraction expansion #2
 * for incomplete beta integral
 */

static double md_incbd( a, b, x )
double a, b, x;
{
double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
double k1, k2, k3, k4, k5, k6, k7, k8;
double r, t, ans, z, thresh;
int n;

k1 = a;
k2 = b - 1.0;
k3 = a;
k4 = a + 1.0;
k5 = 1.0;
k6 = a + b;
k7 = a + 1.0;;
k8 = a + 2.0;

pkm2 = 0.0;
qkm2 = 1.0;
pkm1 = 1.0;
qkm1 = 1.0;
z = x / (1.0-x);
ans = 1.0;
r = 1.0;
n = 0;
thresh = 3.0 * MACHEP;
do
	{
	
	xk = -( z * k1 * k2 )/( k3 * k4 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	xk = ( z * k5 * k6 )/( k7 * k8 );
	pk = pkm1 +  pkm2 * xk;
	qk = qkm1 +  qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	if( qk != 0 )
		r = pk/qk;
	if( r != 0 )
		{
		t = md_fabs( (ans - r)/r );
		ans = r;
		}
	else
		t = 1.0;

	if( t < thresh )
		goto cdone;

	k1 += 1.0;
	k2 -= 1.0;
	k3 += 2.0;
	k4 += 2.0;
	k5 += 1.0;
	k6 += 1.0;
	k7 += 2.0;
	k8 += 2.0;

	if( (md_fabs(qk) + md_fabs(pk)) > big )
		{
		pkm2 *= biginv;
		pkm1 *= biginv;
		qkm2 *= biginv;
		qkm1 *= biginv;
		}
	if( (md_fabs(qk) < biginv) || (md_fabs(pk) < biginv) )
		{
		pkm2 *= big;
		pkm1 *= big;
		qkm2 *= big;
		qkm1 *= big;
		}
	}
while( ++n < 300 );
cdone:
return(ans);
}

/* Power series for incomplete beta integral.
   Use when b*x is small and x not too close to 1.  */

static double md_pseries( a, b, x )
double a, b, x;
{
double s, t, u, v, n, t1, z, ai;

ai = 1.0 / a;
u = (1.0 - b) * x;
v = u / (a + 1.0);
t1 = v;
t = u;
n = 2.0;
s = 0.0;
z = MACHEP * ai;
while( md_fabs(v) > z )
	{
	u = (n - b) * x / n;
	t *= u;
	v = t / (a + n);
	s += v; 
	n += 1.0;
	}
s += t1;
s += ai;

u = a * md_log(x);
if( (a+b) < MAXGAM && md_fabs(u) < MAXLOG )
	{
	t = md_gamma(a+b)/(md_gamma(a)*md_gamma(b));
	s = s * t * md_pow(x,a);
	}
else
	{
	t = md_lgam(a+b) - md_lgam(a) - md_lgam(b) + u + md_log(s);
	if( t < MINLOG )
		s = 0.0;
	else
	s = md_exp(t);
	}
return(s);
}
