/*							struve.c
 *
 *      Struve function
 *
 *
 *
 * SYNOPSIS:
 *
 * double v, x, y, cfs_struve();
 *
 * y = cfs_struve( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the Struve function Hv(x) of order v, argument x.
 * Negative x is rejected unless v is an integer.
 *
 * This module also contains the hypergeometric functions 1F2
 * and 3F0 and a routine for the Bessel function Yv(x) with
 * noninteger v.
 *
 *
 *
 * ACCURACY:
 *
 * Not accurately characterized, but spot checked against tables.
 *
 */


/*
Cephes Math Library Release 2.81:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/
#include "mconf.h"

#ifdef DEBUG
#undef DEBUG
#define DEBUG 0

#ifdef ANSIPROT
extern double cfs_gamma ( double );
extern double cfs_pow ( double, double );
extern double cfs_sqrt ( double );
extern double cfs_yn ( int, double );
extern double cfs_jv ( double, double );
extern double cfs_fabs ( double );
extern double cfs_floor ( double );
extern double cfs_sin ( double );
extern double cfs_cos ( double );
double cfs_yv ( double, double );
double cfs_onef2 (double, double, double, double, double * );
double cfs_threef0 (double, double, double, double, double * );
#else
double cfs_gamma(), cfs_pow(), cfs_sqrt(), cfs_yn(), cfs_yv(), cfs_jv(), cfs_fabs(), cfs_floor();
double cfs_sin(), cfs_cos();
double cfs_onef2(), cfs_threef0();
#endif
static double stop = 1.37e-17;
extern double MACHEP;

double cfs_onef2( a, b, c, x, err )
double a, b, c, x;
double *err;
{
double n, a0, sum, t;
double an, bn, cn, max, z;

an = a;
bn = b;
cn = c;
a0 = 1.0;
sum = 1.0;
n = 1.0;
t = 1.0;
max = 0.0;

do
	{
	if( an == 0 )
		goto done;
	if( bn == 0 )
		goto error;
	if( cn == 0 )
		goto error;
	if( (a0 > 1.0e34) || (n > 200) )
		goto error;
	a0 *= (an * x) / (bn * cn * n);
	sum += a0;
	an += 1.0;
	bn += 1.0;
	cn += 1.0;
	n += 1.0;
	z = cfs_fabs( a0 );
	if( z > max )
		max = z;
	if( sum != 0 )
		t = cfs_fabs( a0 / sum );
	else
		t = z;
	}
while( t > stop );

done:

*err = cfs_fabs( MACHEP*max /sum );

#if DEBUG
	printf(" onef2 cancellation error %.5E\n", *err );
#endif

goto xit;

error:
#if DEBUG
printf("onef2 does not converge\n");
#endif
*err = 1.0e38;

xit:

#if DEBUG
printf("onef2( %.2E %.2E %.2E %.5E ) =  %.3E  %.6E\n", a, b, c, x, n, sum);
#endif
return(sum);
}




double cfs_threef0( a, b, c, x, err )
double a, b, c, x;
double *err;
{
double n, a0, sum, t, conv, conv1;
double an, bn, cn, max, z;

an = a;
bn = b;
cn = c;
a0 = 1.0;
sum = 1.0;
n = 1.0;
t = 1.0;
max = 0.0;
conv = 1.0e38;
conv1 = conv;

do
	{
	if( an == 0.0 )
		goto done;
	if( bn == 0.0 )
		goto done;
	if( cn == 0.0 )
		goto done;
	if( (a0 > 1.0e34) || (n > 200) )
		goto error;
	a0 *= (an * bn * cn * x) / n;
	an += 1.0;
	bn += 1.0;
	cn += 1.0;
	n += 1.0;
	z = cfs_fabs( a0 );
	if( z > max )
		max = z;
	if( z >= conv )
		{
		if( (z < max) && (z > conv1) )
			goto done;
		}
	conv1 = conv;
	conv = z;
	sum += a0;
	if( sum != 0 )
		t = cfs_fabs( a0 / sum );
	else
		t = z;
	}
while( t > stop );

done:

t = cfs_fabs( MACHEP*max/sum );
#if DEBUG
	printf(" threef0 cancellation error %.5E\n", t );
#endif

max = cfs_fabs( conv/sum );
if( max > t )
	t = max;
#if DEBUG
	printf(" threef0 convergence %.5E\n", max );
#endif

goto xit;

error:
#if DEBUG
printf("threef0 does not converge\n");
#endif
t = 1.0e38;

xit:

#if DEBUG
printf("threef0( %.2E %.2E %.2E %.5E ) =  %.3E  %.6E\n", a, b, c, x, n, sum);
#endif

*err = t;
return(sum);
}




extern double PI;

double cfs_struve( v, x )
double v, x;
{
double y, ya, f, g, h, t;
double onef2err, threef0err;

f = cfs_floor(v);
if( (v < 0) && ( v-f == 0.5 ) )
	{
	y = cfs_jv( -v, x );
	f = 1.0 - f;
	g =  2.0 * cfs_floor(f/2.0);
	if( g != f )
		y = -y;
	return(y);
	}
t = 0.25*x*x;
f = cfs_fabs(x);
g = 1.5 * cfs_fabs(v);
if( (f > 30.0) && (f > g) )
	{
	onef2err = 1.0e38;
	y = 0.0;
	}
else
	{
	y = cfs_onef2( 1.0, 1.5, 1.5+v, -t, &onef2err );
	}

if( (f < 18.0) || (x < 0.0) )
	{
	threef0err = 1.0e38;
	ya = 0.0;
	}
else
	{
	ya = cfs_threef0( 1.0, 0.5, 0.5-v, -1.0/t, &threef0err );
	}

f = cfs_sqrt( PI );
h = cfs_pow( 0.5*x, v-1.0 );

if( onef2err <= threef0err )
	{
	g = cfs_gamma( v + 1.5 );
	y = y * h * t / ( 0.5 * f * g );
	return(y);
	}
else
	{
	g = cfs_gamma( v + 0.5 );
	ya = ya * h / ( f * g );
	ya = ya + cfs_yv( v, x );
	return(ya);
	}
}




/* Bessel function of noninteger order
 */

double cfs_yv( v, x )
double v, x;
{
double y, t;
int n;

y = cfs_floor( v );
if( y == v )
	{
	n = v;
	y = cfs_yn( n, x );
	return( y );
	}
t = PI * v;
y = (cfs_cos(t) * cfs_jv( v, x ) - cfs_jv( -v, x ))/cfs_sin(t);
return( y );
}

#endif // DEBUG

/* Crossover points between ascending series and asymptotic series
 * for Struve function
 *
 *	 v	 x
 * 
 *	 0	19.2
 *	 1	18.95
 *	 2	19.15
 *	 3	19.3
 *	 5	19.7
 *	10	21.35
 *	20	26.35
 *	30	32.31
 *	40	40.0
 */
