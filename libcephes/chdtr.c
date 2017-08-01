/*							chdtr.c
 *
 *	Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double df, x, y, md_chdtr();
 *
 * y = md_chdtr( df, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the left hand tail (from 0 to x)
 * of the Chi square probability density function with
 * v degrees of freedom.
 *
 *
 *                                  inf.
 *                                    -
 *                        1          | |  v/2-1  -t/2
 *  P( x | v )   =   -----------     |   t      e     dt
 *                    v/2  -       | |
 *                   2    | (v/2)   -
 *                                   x
 *
 * where x is the Chi-square variable.
 *
 * The incomplete md_gamma integral is used, according to the
 * formula
 *
 *	y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).
 *
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 * See igam().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtr domain   x < 0 or v < 1        0.0
 */
/*							chdtrc()
 *
 *	Complemented Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double v, x, y, md_chdtrc();
 *
 * y = md_chdtrc( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the right hand tail (from x to
 * infinity) of the Chi square probability density function
 * with v degrees of freedom:
 *
 *
 *                                  inf.
 *                                    -
 *                        1          | |  v/2-1  -t/2
 *  P( x | v )   =   -----------     |   t      e     dt
 *                    v/2  -       | |
 *                   2    | (v/2)   -
 *                                   x
 *
 * where x is the Chi-square variable.
 *
 * The incomplete md_gamma integral is used, according to the
 * formula
 *
 *	y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).
 *
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 * See md_igamc().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtrc domain  x < 0 or v < 1        0.0
 */
/*							chdtri()
 *
 *	Inverse of complemented Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double df, x, y, md_chdtri();
 *
 * x = md_chdtri( df, y );
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given cumulative probability y.
 *
 * This is accomplished using the inverse md_gamma integral
 * function and the relation
 *
 *    x/2 = md_igami( df/2, y );
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
 *   message            condition      value returned
 * md_chdtri domain   y < 0 or y > 1        0.0
 *                        v < 1
 *
 */

/*								md_chdtr() */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

#include "mconf.h"
#ifdef ANSIPROT
extern double md_igamc ( double, double );
extern double md_igam ( double, double );
extern double md_igami ( double, double );
#else
double md_igamc(), md_igam(), md_igami();
#endif

double md_chdtrc(df,x)
double df, x;
{

if( (x < 0.0) || (df < 1.0) )
	{
	mtherr( "chdtrc", DOMAIN );
	return(0.0);
	}
return( md_igamc( df/2.0, x/2.0 ) );
}



double md_chdtr(df,x)
double df, x;
{

if( (x < 0.0) || (df < 1.0) )
	{
	mtherr( "chdtr", DOMAIN );
	return(0.0);
	}
return( md_igam( df/2.0, x/2.0 ) );
}



double md_chdtri( df, y )
double df, y;
{
double x;

if( (y < 0.0) || (y > 1.0) || (df < 1.0) )
	{
	mtherr( "chdtri", DOMAIN );
	return(0.0);
	}

x = md_igami( 0.5 * df, y );
return( 2.0 * x );
}
