/*							pdtr.c
 *
 *	Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, cfs_pdtr();
 *
 * y = cfs_pdtr( k, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the first k terms of the Poisson
 * distribution:
 *
 *   k         j
 *   --   -m  m
 *   >   e    --
 *   --       j!
 *  j=0
 *
 * The terms are not summed directly; instead the incomplete
 * cfs_gamma integral is employed, according to the relation
 *
 * y = cfs_pdtr( k, m ) = cfs_igamc( k+1, m ).
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 * See cfs_igamc().
 *
 */

/*							cfs_pdtrc()
 *
 *	Complemented poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, cfs_pdtrc();
 *
 * y = cfs_pdtrc( k, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 to infinity of the Poisson
 * distribution:
 *
 *  inf.       j
 *   --   -m  m
 *   >   e    --
 *   --       j!
 *  j=k+1
 *
 * The terms are not summed directly; instead the incomplete
 * cfs_gamma integral is employed, according to the formula
 *
 * y = cfs_pdtrc( k, m ) = cfs_igam( k+1, m ).
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 * See igam.c.
 *
 */

/*							cfs_pdtri()
 *
 *	Inverse Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, cfs_pdtri();
 *
 * m = cfs_pdtri( k, y );
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the Poisson variable x such that the integral
 * from 0 to x of the Poisson density is equal to the
 * given probability y.
 *
 * This is accomplished using the inverse cfs_gamma integral
 * function and the relation
 *
 *    m = cfs_igami( k+1, y ).
 *
 *
 *
 *
 * ACCURACY:
 *
 * See igami.c.
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * pdtri domain    y < 0 or y >= 1       0.0
 *                     k < 0
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*/

#include "mconf.h"
#ifdef ANSIPROT
extern double cfs_igam ( double, double );
extern double cfs_igamc ( double, double );
extern double cfs_igami ( double, double );
#else
double cfs_igam(), cfs_igamc(), cfs_igami();
#endif

double cfs_pdtrc( k, m )
int k;
double m;
{
double v;

if( (k < 0) || (m <= 0.0) )
	{
	mtherr( "pdtrc", DOMAIN );
	return( 0.0 );
	}
v = k+1;
return( cfs_igam( v, m ) );
}



double cfs_pdtr( k, m )
int k;
double m;
{
double v;

if( (k < 0) || (m <= 0.0) )
	{
	mtherr( "pdtr", DOMAIN );
	return( 0.0 );
	}
v = k+1;
return( cfs_igamc( v, m ) );
}


double cfs_pdtri( k, y )
int k;
double y;
{
double v;

if( (k < 0) || (y < 0.0) || (y >= 1.0) )
	{
	mtherr( "pdtri", DOMAIN );
	return( 0.0 );
	}
v = k+1;
v = cfs_igami( v, y );
return( v );
}
