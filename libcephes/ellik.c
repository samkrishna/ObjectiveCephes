/*							ellik.c
 *
 *	Incomplete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, md_ellik();
 *
 * y = md_ellik( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi_\m)  =    |    ------------------
 *                |                   2
 *              | |    md_sqrt( 1 - m md_sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points with m in [0, 1] and phi as indicated.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -10,10       200000      7.4e-16     1.0e-16
 *
 *
 */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/*	Incomplete elliptic integral of first kind	*/

#include "mconf.h"
#ifdef ANSIPROT
extern double md_sqrt ( double );
extern double md_fabs ( double );
extern double md_log ( double );
extern double md_tan ( double );
extern double md_atan ( double );
extern double md_floor ( double );
extern double md_ellpk ( double );
double md_ellik ( double, double );
#else
double md_sqrt(), md_fabs(), md_log(), md_tan(), md_atan(), md_floor(), md_ellpk();
double md_ellik();
#endif
extern double PI, PIO2, MACHEP, MAXNUM;

double md_ellik( phi, m )
double phi, m;
{
double a, b, c, e, temp, t, K;
int d, mod, sign, npio2;

if( m == 0.0 )
	return( phi );
a = 1.0 - m;
if( a == 0.0 )
	{
	if( md_fabs(phi) >= PIO2 )
		{
		mtherr( "ellik", SING );
		return( MAXNUM );
		}
	return(  md_log(  md_tan( (PIO2 + phi)/2.0 )  )   );
	}
npio2 = md_floor( phi/PIO2 );
if( npio2 & 1 )
	npio2 += 1;
if( npio2 )
	{
	K = md_ellpk( a );
	phi = phi - npio2 * PIO2;
	}
else
	K = 0.0;
if( phi < 0.0 )
	{
	phi = -phi;
	sign = -1;
	}
else
	sign = 0;
b = md_sqrt(a);
t = md_tan( phi );
if( md_fabs(t) > 10.0 )
	{
	/* Transform the amplitude */
	e = 1.0/(b*t);
	/* ... but avoid multiple recursions.  */
	if( md_fabs(e) < 10.0 )
		{
		e = md_atan(e);
		if( npio2 == 0 )
			K = md_ellpk( a );
		temp = K - md_ellik( e, m );
		goto done;
		}
	}
a = 1.0;
c = md_sqrt(m);
d = 1;
mod = 0;

while( md_fabs(c/a) > MACHEP )
	{
	temp = b/a;
	phi = phi + md_atan(t*temp) + mod * PI;
	mod = (phi + PIO2)/PI;
	t = t * ( 1.0 + temp )/( 1.0 - temp * t * t );
	c = ( a - b )/2.0;
	temp = md_sqrt( a * b );
	a = ( a + b )/2.0;
	b = temp;
	d += d;
	}

temp = (md_atan(t) + mod * PI)/(d * a);

done:
if( sign < 0 )
	temp = -temp;
temp += npio2 * K;
return( temp );
}
