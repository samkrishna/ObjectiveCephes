/*							gdtr.c
 *
 *	Gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, cfs_gdtr();
 *
 * y = cfs_gdtr( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from zero to x of the cfs_gamma probability
 * density function:
 *
 *
 *                x
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               0
 *
 *  The incomplete cfs_gamma integral is used, according to the
 * relation
 *
 * y = cfs_igam( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See cfs_igam().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtr domain         x < 0            0.0
 *
 */

/*							gdtrc.c
 *
 *	Complemented gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, cfs_gdtrc();
 *
 * y = cfs_gdtrc( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from x to infinity of the cfs_gamma
 * probability density function:
 *
 *
 *               inf.
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               x
 *
 *  The incomplete cfs_gamma integral is used, according to the
 * relation
 *
 * y = cfs_igamc( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See cfs_igamc().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtrc domain         x < 0            0.0
 *
 */

/*							gdtr()  */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*/

#include "mconf.h"
#ifdef ANSIPROT
extern double cfs_igam ( double, double );
extern double cfs_igamc ( double, double );
#else
double cfs_igam(), cfs_igamc();
#endif

double cfs_gdtr( a, b, x )
double a, b, x;
{

if( x < 0.0 )
	{
	mtherr( "gdtr", DOMAIN );
	return( 0.0 );
	}
return(  cfs_igam( b, a * x )  );
}



double cfs_gdtrc( a, b, x )
double a, b, x;
{

if( x < 0.0 )
	{
	mtherr( "gdtrc", DOMAIN );
	return( 0.0 );
	}
return(  cfs_igamc( b, a * x )  );
}
