/*							pdtr.c
 *
 *	Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, md_pdtr();
 *
 * y = md_pdtr( k, m );
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
 * md_gamma integral is employed, according to the relation
 *
 * y = md_pdtr( k, m ) = md_igamc( k+1, m ).
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 * See md_igamc().
 *
 */

/*							md_pdtrc()
 *
 *	Complemented poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, md_pdtrc();
 *
 * y = md_pdtrc( k, m );
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
 * md_gamma integral is employed, according to the formula
 *
 * y = md_pdtrc( k, m ) = md_igam( k+1, m ).
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

/*							md_pdtri()
 *
 *	Inverse Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, md_pdtri();
 *
 * m = md_pdtri( k, y );
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
 * This is accomplished using the inverse md_gamma integral
 * function and the relation
 *
 *    m = md_igami( k+1, y ).
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
extern double md_igam ( double, double );
extern double md_igamc ( double, double );
extern double md_igami ( double, double );
#else
double md_igam(), md_igamc(), md_igami();
#endif

double md_pdtrc( k, m )
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
return( md_igam( v, m ) );
}



double md_pdtr( k, m )
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
return( md_igamc( v, m ) );
}


double md_pdtri( k, y )
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
v = md_igami( v, y );
return( v );
}
