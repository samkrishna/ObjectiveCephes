/*							cfs_clog.c
 *
 *	Complex natural logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_clog();
 * cmplx z, w;
 *
 * cfs_clog( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns complex logarithm to the base e (2.718...) of
 * the complex argument x.
 *
 * If z = x + iy, r = cfs_sqrt( x**2 + y**2 ),
 * then
 *       w = cfs_log(r) + i cfs_arctan(y/x).
 * 
 * The arctangent ranges from -PI to +PI.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      7000       8.5e-17     1.9e-17
 *    IEEE      -10,+10     30000       5.0e-15     1.1e-16
 *
 * Larger relative error can be observed for z near 1 +i0.
 * In IEEE arithmetic the peak absolute error is 5.2e-16, rms
 * absolute error 1.0e-16.
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1995, 2000 by Stephen L. Moshier
*/
#include "mconf.h"
#ifdef ANSIPROT
static void cfs_cchsh ( double x, double *c, double *s );
static double cfs_redupi ( double x );
static double cfs_ctans ( cmplx *z );
/* These are supposed to be in some standard place. */
double cfs_fabs (double);
double cfs_sqrt (double);
double cfs_pow (double, double);
double cfs_log (double);
double cfs_exp (double);
double cfs_atan2 (double, double);
double cfs_cosh (double);
double cfs_sinh (double);
double cfs_asin (double);
double cfs_sin (double);
double cfs_cos (double);
double cfs_cabs (cmplx *);
void cfs_cadd ( cmplx *, cmplx *, cmplx * );
void cfs_cmul ( cmplx *, cmplx *, cmplx * );
void cfs_csqrt ( cmplx *, cmplx * );
static void cfs_cchsh ( double, double *, double * );
static double cfs_redupi ( double );
static double cfs_ctans ( cmplx * );
void cfs_clog ( cmplx *, cmplx * );
void cfs_casin ( cmplx *, cmplx * );
void cfs_cacos ( cmplx *, cmplx * );
void cfs_catan ( cmplx *, cmplx * );
#else
static void cfs_cchsh();
static double cfs_redupi();
static double cfs_ctans();
double cfs_cabs(), cfs_fabs(), cfs_sqrt(), cfs_pow();
double cfs_log(), cfs_exp(), cfs_atan2(), cfs_cosh(), cfs_sinh();
double cfs_asin(), cfs_sin(), cfs_cos();
void cfs_cadd(), cfs_cmul(), cfs_csqrt();
void cfs_clog(), cfs_casin(), cfs_cacos(), cfs_catan();
#endif


extern double MAXNUM, MACHEP, PI, PIO2;

void cfs_clog( z, w )
register cmplx *z, *w;
{
double p, rr;

/*rr = cfs_sqrt( z->r * z->r  +  z->i * z->i );*/
rr = cfs_cabs(z);
p = cfs_log(rr);
#if ANSIC
rr = cfs_atan2( z->i, z->r );
#else
rr = cfs_atan2( z->r, z->i );
if( rr > PI )
	rr -= PI + PI;
#endif
w->i = rr;
w->r = p;
}

/*							cfs_cexp()
 *
 *	Complex exponential function
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_cexp();
 * cmplx z, w;
 *
 * cfs_cexp( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the exponential of the complex argument z
 * into the complex result w.
 *
 * If
 *     z = x + iy,
 *     r = cfs_exp(x),
 *
 * then
 *
 *     w = r cfs_cos y + i r cfs_sin y.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      8700       3.7e-17     1.1e-17
 *    IEEE      -10,+10     30000       3.0e-16     8.7e-17
 *
 */

void cfs_cexp( z, w )
register cmplx *z, *w;
{
double r;

r = cfs_exp( z->r );
w->r = r * cfs_cos( z->i );
w->i = r * cfs_sin( z->i );
}

/*							cfs_csin()
 *
 *	Complex circular sine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_csin();
 * cmplx z, w;
 *
 * cfs_csin( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * If
 *     z = x + iy,
 *
 * then
 *
 *     w = cfs_sin x  cfs_cosh y  +  i cfs_cos x cfs_sinh y.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      8400       5.3e-17     1.3e-17
 *    IEEE      -10,+10     30000       3.8e-16     1.0e-16
 * Also tested by cfs_csin(cfs_casin(z)) = z.
 *
 */

void cfs_csin( z, w )
register cmplx *z, *w;
{
double ch, sh;

cfs_cchsh( z->i, &ch, &sh );
w->r = cfs_sin( z->r ) * ch;
w->i = cfs_cos( z->r ) * sh;
}



/* calculate cfs_cosh and cfs_sinh */

static void cfs_cchsh( x, c, s )
double x, *c, *s;
{
double e, ei;

if( cfs_fabs(x) <= 0.5 )
	{
	*c = cfs_cosh(x);
	*s = cfs_sinh(x);
	}
else
	{
	e = cfs_exp(x);
	ei = 0.5/e;
	e = 0.5 * e;
	*s = e - ei;
	*c = e + ei;
	}
}

/*							cfs_ccos()
 *
 *	Complex circular cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_ccos();
 * cmplx z, w;
 *
 * cfs_ccos( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * If
 *     z = x + iy,
 *
 * then
 *
 *     w = cfs_cos x  cfs_cosh y  -  i cfs_sin x cfs_sinh y.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      8400       4.5e-17     1.3e-17
 *    IEEE      -10,+10     30000       3.8e-16     1.0e-16
 */

void cfs_ccos( z, w )
register cmplx *z, *w;
{
double ch, sh;

cfs_cchsh( z->i, &ch, &sh );
w->r = cfs_cos( z->r ) * ch;
w->i = -cfs_sin( z->r ) * sh;
}

/*							cfs_ctan()
 *
 *	Complex circular tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_ctan();
 * cmplx z, w;
 *
 * cfs_ctan( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * If
 *     z = x + iy,
 *
 * then
 *
 *           cfs_sin 2x  +  i cfs_sinh 2y
 *     w  =  --------------------.
 *            cfs_cos 2x  +  cfs_cosh 2y
 *
 * On the real axis the denominator is zero at odd multiples
 * of PI/2.  The denominator is evaluated by its Taylor
 * series near these points.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      5200       7.1e-17     1.6e-17
 *    IEEE      -10,+10     30000       7.2e-16     1.2e-16
 * Also tested by cfs_ctan * ccot = 1 and cfs_catan(cfs_ctan(z))  =  z.
 */

void cfs_ctan( z, w )
register cmplx *z, *w;
{
double d;

d = cfs_cos( 2.0 * z->r ) + cfs_cosh( 2.0 * z->i );

if( cfs_fabs(d) < 0.25 )
	d = cfs_ctans(z);

if( d == 0.0 )
	{
	mtherr( "cfs_ctan", OVERFLOW );
	w->r = MAXNUM;
	w->i = MAXNUM;
	return;
	}

w->r = cfs_sin( 2.0 * z->r ) / d;
w->i = cfs_sinh( 2.0 * z->i ) / d;
}

/*							cfs_ccot()
 *
 *	Complex circular cotangent
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_ccot();
 * cmplx z, w;
 *
 * ccot( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * If
 *     z = x + iy,
 *
 * then
 *
 *           cfs_sin 2x  -  i cfs_sinh 2y
 *     w  =  --------------------.
 *            cfs_cosh 2y  -  cfs_cos 2x
 *
 * On the real axis, the denominator has zeros at even
 * multiples of PI/2.  Near these points it is evaluated
 * by a Taylor series.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      3000       6.5e-17     1.6e-17
 *    IEEE      -10,+10     30000       9.2e-16     1.2e-16
 * Also tested by cfs_ctan * ccot = 1 + cfs_i0.
 */

void cfs_ccot( z, w )
register cmplx *z, *w;
{
double d;

d = cfs_cosh(2.0 * z->i) - cfs_cos(2.0 * z->r);

if( cfs_fabs(d) < 0.25 )
	d = cfs_ctans(z);

if( d == 0.0 )
	{
	mtherr( "ccot", OVERFLOW );
	w->r = MAXNUM;
	w->i = MAXNUM;
	return;
	}

w->r = cfs_sin( 2.0 * z->r ) / d;
w->i = -cfs_sinh( 2.0 * z->i ) / d;
}

/* Program to subtract nearest integer multiple of PI */
/* extended precision value of PI: */
#ifdef UNK
static double DP1 = 3.14159265160560607910E0;
static double DP2 = 1.98418714791870343106E-9;
static double DP3 = 1.14423774522196636802E-17;
#endif

#ifdef DEC
static unsigned short P1[] = {0040511,0007732,0120000,0000000,};
static unsigned short P2[] = {0031010,0055060,0100000,0000000,};
static unsigned short P3[] = {0022123,0011431,0105056,0001560,};
#define DP1 *(double *)P1
#define DP2 *(double *)P2
#define DP3 *(double *)P3
#endif

#ifdef IBMPC
static unsigned short P1[] = {0x0000,0x5400,0x21fb,0x4009};
static unsigned short P2[] = {0x0000,0x1000,0x0b46,0x3e21};
static unsigned short P3[] = {0xc06e,0x3145,0x6263,0x3c6a};
#define DP1 *(double *)P1
#define DP2 *(double *)P2
#define DP3 *(double *)P3
#endif

#ifdef MIEEE
static unsigned short P1[] = {
0x4009,0x21fb,0x5400,0x0000
};
static unsigned short P2[] = {
0x3e21,0x0b46,0x1000,0x0000
};
static unsigned short P3[] = {
0x3c6a,0x6263,0x3145,0xc06e
};
#define DP1 *(double *)P1
#define DP2 *(double *)P2
#define DP3 *(double *)P3
#endif

static double cfs_redupi(x)
double x;
{
double t;
long i;

t = x/PI;
if( t >= 0.0 )
	t += 0.5;
else
	t -= 0.5;

i = t;	/* the multiple */
t = i;
t = ((x - t * DP1) - t * DP2) - t * DP3;
return(t);
}

/*  Taylor series expansion for cfs_cosh(2y) - cfs_cos(2x)	*/

static double cfs_ctans(z)
cmplx *z;
{
double f, x, x2, y, y2, rn, t;
double d;

x = cfs_fabs( 2.0 * z->r );
y = cfs_fabs( 2.0 * z->i );

x = cfs_redupi(x);

x = x * x;
y = y * y;
x2 = 1.0;
y2 = 1.0;
f = 1.0;
rn = 0.0;
d = 0.0;
do
	{
	rn += 1.0;
	f *= rn;
	rn += 1.0;
	f *= rn;
	x2 *= x;
	y2 *= y;
	t = y2 + x2;
	t /= f;
	d += t;

	rn += 1.0;
	f *= rn;
	rn += 1.0;
	f *= rn;
	x2 *= x;
	y2 *= y;
	t = y2 - x2;
	t /= f;
	d += t;
	}
while( cfs_fabs(t/d) > MACHEP );
return(d);
}

/*							cfs_casin()
 *
 *	Complex circular arc sine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_casin();
 * cmplx z, w;
 *
 * cfs_casin( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * Inverse complex sine:
 *
 *                               2
 * w = -i cfs_clog( iz + cfs_csqrt( 1 - z ) ).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10     10100       2.1e-15     3.4e-16
 *    IEEE      -10,+10     30000       2.2e-14     2.7e-15
 * Larger relative error can be observed for z near zero.
 * Also tested by cfs_csin(cfs_casin(z)) = z.
 */

void cfs_casin( z, w )
cmplx *z, *w;
{
static cmplx ca, ct, zz, z2;
double x, y;

x = z->r;
y = z->i;

if( y == 0.0 )
	{
	if( cfs_fabs(x) > 1.0 )
		{
		w->r = PIO2;
		w->i = 0.0;
		mtherr( "cfs_casin", DOMAIN );
		}
	else
		{
		w->r = cfs_asin(x);
		w->i = 0.0;
		}
	return;
	}

/* Power series expansion */
/*
b = cfs_cabs(z);
if( b < 0.125 )
{
z2.r = (x - y) * (x + y);
z2.i = 2.0 * x * y;

cn = 1.0;
n = 1.0;
ca.r = x;
ca.i = y;
sum.r = x;
sum.i = y;
do
	{
	ct.r = z2.r * ca.r  -  z2.i * ca.i;
	ct.i = z2.r * ca.i  +  z2.i * ca.r;
	ca.r = ct.r;
	ca.i = ct.i;

	cn *= n;
	n += 1.0;
	cn /= n;
	n += 1.0;
	b = cn/n;

	ct.r *= b;
	ct.i *= b;
	sum.r += ct.r;
	sum.i += ct.i;
	b = cfs_fabs(ct.r) + cfs_fabs(ct.i);
	}
while( b > MACHEP );
w->r = sum.r;
w->i = sum.i;
return;
}
*/


ca.r = x;
ca.i = y;

ct.r = -ca.i;	/* iz */
ct.i = ca.r;

	/* cfs_sqrt( 1 - z*z) */
/* cfs_cmul( &ca, &ca, &zz ) */
zz.r = (ca.r - ca.i) * (ca.r + ca.i);	/*x * x  -  y * y */
zz.i = 2.0 * ca.r * ca.i;

zz.r = 1.0 - zz.r;
zz.i = -zz.i;
cfs_csqrt( &zz, &z2 );

cfs_cadd( &z2, &ct, &zz );
cfs_clog( &zz, &zz );
w->r = zz.i;	/* mult by 1/i = -i */
w->i = -zz.r;
return;
}

/*							cfs_cacos()
 *
 *	Complex circular arc cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_cacos();
 * cmplx z, w;
 *
 * cfs_cacos( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 * w = cfs_arccos z  =  PI/2 - cfs_arcsin z.
 *
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      5200      1.6e-15      2.8e-16
 *    IEEE      -10,+10     30000      1.8e-14      2.2e-15
 */

void cfs_cacos( z, w )
cmplx *z, *w;
{

cfs_casin( z, w );
w->r = PIO2  -  w->r;
w->i = -w->i;
}

/*							cfs_catan()
 *
 *	Complex circular arc tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_catan();
 * cmplx z, w;
 *
 * cfs_catan( &z, &w );
 *
 *
 *
 * DESCRIPTION:
 *
 * If
 *     z = x + iy,
 *
 * then
 *          1       (    2x     )
 * Re w  =  - arctan(-----------)  +  k PI
 *          2       (     2    2)
 *                  (1 - x  - y )
 *
 *               ( 2         2)
 *          1    (x  +  (y+1) )
 * Im w  =  - cfs_log(------------)
 *          4    ( 2         2)
 *               (x  +  (y-1) )
 *
 * Where k is an arbitrary integer.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       -10,+10      5900       1.3e-16     7.8e-18
 *    IEEE      -10,+10     30000       2.3e-15     8.5e-17
 * The check cfs_catan( cfs_ctan(z) )  =  z, with |x| and |y| < PI/2,
 * had peak relative error 1.5e-16, rms relative error
 * 2.9e-17.  See also cfs_clog().
 */

void cfs_catan( z, w )
cmplx *z, *w;
{
double a, t, x, x2, y;

x = z->r;
y = z->i;

if( (x == 0.0) && (y > 1.0) )
	goto ovrf;

x2 = x * x;
a = 1.0 - x2 - (y * y);
if( a == 0.0 )
	goto ovrf;

#if ANSIC
t = cfs_atan2( 2.0 * x, a )/2.0;
#else
t = cfs_atan2( a, 2.0 * x )/2.0;
#endif
w->r = cfs_redupi( t );

t = y - 1.0;
a = x2 + (t * t);
if( a == 0.0 )
	goto ovrf;

t = y + 1.0;
a = (x2 + (t * t))/a;
w->i = cfs_log(a)/4.0;
return;

ovrf:
mtherr( "cfs_catan", OVERFLOW );
w->r = MAXNUM;
w->i = MAXNUM;
}


/*							cfs_csinh
 *
 *	Complex hyperbolic sine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_csinh();
 * cmplx z, w;
 *
 * cfs_csinh( &z, &w );
 *
 *
 * DESCRIPTION:
 *
 * cfs_csinh z = (cfs_cexp(z) - cfs_cexp(-z))/2
 *         = cfs_sinh x * cfs_cos y  +  i cfs_cosh x * cfs_sin y .
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000       3.1e-16     8.2e-17
 *
 */

void
cfs_csinh (z, w)
     cmplx *z, *w;
{
  double x, y;

  x = z->r;
  y = z->i;
  w->r = cfs_sinh (x) * cfs_cos (y);
  w->i = cfs_cosh (x) * cfs_sin (y);
}


/*							cfs_casinh
 *
 *	Complex inverse hyperbolic sine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_casinh();
 * cmplx z, w;
 *
 * cfs_casinh (&z, &w);
 *
 *
 *
 * DESCRIPTION:
 *
 * cfs_casinh z = -i cfs_casin iz .
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000       1.8e-14     2.6e-15
 *
 */

void
cfs_casinh (z, w)
     cmplx *z, *w;
{
  cmplx u;

  u.r = 0.0;
  u.i = 1.0;
  cfs_cmul( z, &u, &u );
  cfs_casin( &u, w );
  u.r = 0.0;
  u.i = -1.0;
  cfs_cmul( &u, w, w );
}

/*							cfs_ccosh
 *
 *	Complex hyperbolic cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_ccosh();
 * cmplx z, w;
 *
 * cfs_ccosh (&z, &w);
 *
 *
 *
 * DESCRIPTION:
 *
 * cfs_ccosh(z) = cfs_cosh x  cfs_cos y + i cfs_sinh x cfs_sin y .
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000       2.9e-16     8.1e-17
 *
 */

void
cfs_ccosh (z, w)
     cmplx *z, *w;
{
  double x, y;

  x = z->r;
  y = z->i;
  w->r = cfs_cosh (x) * cfs_cos (y);
  w->i = cfs_sinh (x) * cfs_sin (y);
}


/*							cfs_cacosh
 *
 *	Complex inverse hyperbolic cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_cacosh();
 * cmplx z, w;
 *
 * cfs_cacosh (&z, &w);
 *
 *
 *
 * DESCRIPTION:
 *
 * cfs_acosh z = i cfs_acos z .
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000       1.6e-14     2.1e-15
 *
 */

void
cfs_cacosh (z, w)
     cmplx *z, *w;
{
  cmplx u;

  cfs_cacos( z, w );
  u.r = 0.0;
  u.i = 1.0;
  cfs_cmul( &u, w, w );
}


/*							cfs_ctanh
 *
 *	Complex hyperbolic tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_ctanh();
 * cmplx z, w;
 *
 * cfs_ctanh (&z, &w);
 *
 *
 *
 * DESCRIPTION:
 *
 * cfs_tanh z = (cfs_sinh 2x  +  i cfs_sin 2y) / (cfs_cosh 2x + cfs_cos 2y) .
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000       1.7e-14     2.4e-16
 *
 */

/* 5.253E-02,1.550E+00 1.643E+01,6.553E+00 1.729E-14  21355  */

void
cfs_ctanh (z, w)
     cmplx *z, *w;
{
  double x, y, d;

  x = z->r;
  y = z->i;
  d = cfs_cosh (2.0 * x) + cfs_cos (2.0 * y);
  w->r = cfs_sinh (2.0 * x) / d;
  w->i = cfs_sin (2.0 * y) / d;
  return;
}


/*							cfs_catanh
 *
 *	Complex inverse hyperbolic tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_catanh();
 * cmplx z, w;
 *
 * cfs_catanh (&z, &w);
 *
 *
 *
 * DESCRIPTION:
 *
 * Inverse cfs_tanh, equal to  -i cfs_catan (iz);
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000       2.3e-16     6.2e-17
 *
 */

void
cfs_catanh (z, w)
     cmplx *z, *w;
{
  cmplx u;

  u.r = 0.0;
  u.i = 1.0;
  cfs_cmul (z, &u, &u);  /* i z */
  cfs_catan (&u, w);
  u.r = 0.0;
  u.i = -1.0;
  cfs_cmul (&u, w, w);  /* -i cfs_catan iz */
  return;
}


/*							cfs_cpow
 *
 *	Complex power function
 *
 *
 *
 * SYNOPSIS:
 *
 * void cfs_cpow();
 * cmplx a, z, w;
 *
 * cfs_cpow (&a, &z, &w);
 *
 *
 *
 * DESCRIPTION:
 *
 * Raises complex A to the complex Zth power.
 * Definition is per AMS55 # 4.2.8,
 * analytically equivalent to cfs_cpow(a,z) = cfs_cexp(z cfs_clog(a)).
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10,+10     30000       9.4e-15     1.5e-15
 *
 */


void
cfs_cpow (a, z, w)
     cmplx *a, *z, *w;
{
  double x, y, r, theta, absa, arga;

  x = z->r;
  y = z->i;
  absa = cfs_cabs (a);
  if (absa == 0.0)
    {
      w->r = 0.0;
      w->i = 0.0;
      return;
    }
  arga = cfs_atan2 (a->i, a->r);
  r = cfs_pow (absa, x);
  theta = x * arga;
  if (y != 0.0)
    {
      r = r * cfs_exp (-y * arga);
      theta = theta + y * cfs_log (absa);
    }
  w->r = r * cfs_cos (theta);
  w->i = r * cfs_sin (theta);
  return;
}
