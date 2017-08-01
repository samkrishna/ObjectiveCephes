/*							gdtr.c
 *
 *	Gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, md_gdtr();
 *
 * y = md_gdtr( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from zero to x of the md_gamma probability
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
 *  The incomplete md_gamma integral is used, according to the
 * relation
 *
 * y = md_igam( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See md_igam().
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
 * double a, b, x, y, md_gdtrc();
 *
 * y = md_gdtrc( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from x to infinity of the md_gamma
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
 *  The incomplete md_gamma integral is used, according to the
 * relation
 *
 * y = md_igamc( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See md_igamc().
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
extern double md_igam ( double, double );
extern double md_igamc ( double, double );
#else
double md_igam(), md_igamc();
#endif

double md_gdtr( a, b, x )
double a, b, x;
{

if( x < 0.0 )
	{
	mtherr( "gdtr", DOMAIN );
	return( 0.0 );
	}
return(  md_igam( b, a * x )  );
}



double md_gdtrc( a, b, x )
double a, b, x;
{

if( x < 0.0 )
	{
	mtherr( "gdtrc", DOMAIN );
	return( 0.0 );
	}
return(  md_igamc( b, a * x )  );
}
