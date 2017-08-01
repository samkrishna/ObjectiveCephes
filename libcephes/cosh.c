/*							cfs_cosh.c
 *
 *	Hyperbolic cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, cfs_cosh();
 *
 * y = cfs_cosh( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns hyperbolic cosine of argument in the range MINLOG to
 * MAXLOG.
 *
 * cfs_cosh(x)  =  ( cfs_exp(x) + cfs_exp(-x) )/2.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       +- 88       50000       4.0e-17     7.7e-18
 *    IEEE     +-MAXLOG     30000       2.6e-16     5.7e-17
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * cfs_cosh overflow    |x| > MAXLOG       MAXNUM
 *
 *
 */

/*							cfs_cosh.c */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1985, 1995, 2000 by Stephen L. Moshier
*/

#include "mconf.h"
#ifdef ANSIPROT
extern double cfs_exp ( double );
extern int cfs_isnan ( double );
extern int cfs_isfinite ( double );
#else
double cfs_exp();
int cfs_isnan(), cfs_isfinite();
#endif
extern double MAXLOG, INFINITY, LOGE2;

double cfs_cosh(x)
double x;
{
double y;

#ifdef NANS
if( cfs_isnan(x) )
	return(x);
#endif
if( x < 0 )
	x = -x;
if( x > (MAXLOG + LOGE2) )
	{
	mtherr( "cfs_cosh", OVERFLOW );
	return( INFINITY );
	}	
if( x >= (MAXLOG - LOGE2) )
	{
	y = cfs_exp(0.5 * x);
	y = (0.5 * y) * y;
	return(y);
	}
y = cfs_exp(x);
y = 0.5 * (y + 1.0 / y);
return( y );
}
