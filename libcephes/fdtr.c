/*							fdtr.c
 *
 *	F distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int df1, df2;
 * double x, y, cfs_fdtr();
 *
 * y = cfs_fdtr( df1, df2, x );
 *
 * DESCRIPTION:
 *
 * Returns the area from zero to x under the F density
 * function (also known as Snedcor's density or the
 * variance ratio density).  This is the density
 * of x = (u1/df1)/(u2/df2), where u1 and u2 are random
 * variables having Chi square distributions with df1
 * and df2 degrees of freedom, respectively.
 *
 * The incomplete beta integral is used, according to the
 * formula
 *
 *	P(x) = cfs_incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).
 *
 *
 * The arguments a and b are greater than zero, and x is
 * nonnegative.
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,x).
 *
 *                x     a,b                     Relative error:
 * arithmetic  domain  domain     # trials      peak         rms
 *    IEEE      0,1    0,100       100000      9.8e-15     1.7e-15
 *    IEEE      1,5    0,100       100000      6.5e-15     3.5e-16
 *    IEEE      0,1    1,10000     100000      2.2e-11     3.3e-12
 *    IEEE      1,5    1,10000     100000      1.1e-11     1.7e-13
 * See also incbet.c.
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * fdtr domain     a<0, b<0, x<0         0.0
 *
 */
/*							fdtrc()
 *
 *	Complemented F distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int df1, df2;
 * double x, y, cfs_fdtrc();
 *
 * y = cfs_fdtrc( df1, df2, x );
 *
 * DESCRIPTION:
 *
 * Returns the area from x to infinity under the F density
 * function (also known as Snedcor's density or the
 * variance ratio density).
 *
 *
 *                      inf.
 *                       -
 *              1       | |  a-1      b-1
 * 1-P(x)  =  ------    |   t    (1-t)    dt
 *            B(a,b)  | |
 *                     -
 *                      x
 *
 *
 * The incomplete beta integral is used, according to the
 * formula
 *
 *	P(x) = cfs_incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).
 *
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,x) in the indicated intervals.
 *                x     a,b                     Relative error:
 * arithmetic  domain  domain     # trials      peak         rms
 *    IEEE      0,1    1,100       100000      3.7e-14     5.9e-16
 *    IEEE      1,5    1,100       100000      8.0e-15     1.6e-15
 *    IEEE      0,1    1,10000     100000      1.8e-11     3.5e-13
 *    IEEE      1,5    1,10000     100000      2.0e-11     3.0e-12
 * See also incbet.c.
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * fdtrc domain    a<0, b<0, x<0         0.0
 *
 */
/*							fdtri()
 *
 *	Inverse of complemented F distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int df1, df2;
 * double x, p, cfs_fdtri();
 *
 * x = cfs_fdtri( df1, df2, p );
 *
 * DESCRIPTION:
 *
 * Finds the F density argument x such that the integral
 * from x to infinity of the F density is equal to the
 * given probability p.
 *
 * This is accomplished using the inverse beta integral
 * function and the relations
 *
 *      z = cfs_incbi( df2/2, df1/2, p )
 *      x = df2 (1-z) / (df1 z).
 *
 * Note: the following relations hold for the inverse of
 * the uncomplemented F distribution:
 *
 *      z = cfs_incbi( df1/2, df2/2, p )
 *      x = df2 z / (df1 (1-z)).
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,p).
 *
 *              a,b                     Relative error:
 * arithmetic  domain     # trials      peak         rms
 *  For p between .001 and 1:
 *    IEEE     1,100       100000      8.3e-15     4.7e-16
 *    IEEE     1,10000     100000      2.1e-11     1.4e-13
 *  For p between 10^-6 and 10^-3:
 *    IEEE     1,100        50000      1.3e-12     8.4e-15
 *    IEEE     1,10000      50000      3.0e-12     4.8e-14
 * See also fdtrc.c.
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * fdtri domain   p <= 0 or p > 1       0.0
 *                     v < 1
 *
 */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*/


#include "mconf.h"
#ifdef ANSIPROT
extern double cfs_incbet ( double, double, double );
extern double cfs_incbi ( double, double, double );
#else
double cfs_incbet(), cfs_incbi();
#endif

double cfs_fdtrc( ia, ib, x )
int ia, ib;
double x;
{
double a, b, w;

if( (ia < 1) || (ib < 1) || (x < 0.0) )
	{
	mtherr( "fdtrc", DOMAIN );
	return( 0.0 );
	}
a = ia;
b = ib;
w = b / (b + a * x);
return( cfs_incbet( 0.5*b, 0.5*a, w ) );
}



double cfs_fdtr( ia, ib, x )
int ia, ib;
double x;
{
double a, b, w;

if( (ia < 1) || (ib < 1) || (x < 0.0) )
	{
	mtherr( "fdtr", DOMAIN );
	return( 0.0 );
	}
a = ia;
b = ib;
w = a * x;
w = w / (b + w);
return( cfs_incbet(0.5*a, 0.5*b, w) );
}


double cfs_fdtri( ia, ib, y )
int ia, ib;
double y;
{
double a, b, w, x;

if( (ia < 1) || (ib < 1) || (y <= 0.0) || (y > 1.0) )
	{
	mtherr( "fdtri", DOMAIN );
	return( 0.0 );
	}
a = ia;
b = ib;
/* Compute probability for x = 0.5.  */
w = cfs_incbet( 0.5*b, 0.5*a, 0.5 );
/* If that is greater than y, then the solution w < .5.
   Otherwise, solve at 1-y to remove cancellation in (b - b*w).  */
if( w > y || y < 0.001)
	{
	w = cfs_incbi( 0.5*b, 0.5*a, y );
	x = (b - b*w)/(a*w);
	}
else
	{
	w = cfs_incbi( 0.5*a, 0.5*b, 1.0-y );
	x = b*w/(a*(1.0-w));
	}
return(x);
}
