/*							beta.c
 *
 *	Beta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, y, cfs_beta();
 *
 * y = cfs_beta( a, b );
 *
 *
 *
 * DESCRIPTION:
 *
 *                      -     -
 *                     | (a) | (b)
 * cfs_beta( a, b )  =  -----------.
 *                        -
 *                       | (a+b)
 *
 * For large arguments the logarithm of the function is
 * evaluated using cfs_lgam(), then exponentiated.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC        0,30        1700       7.7e-15     1.5e-15
 *    IEEE       0,30       30000       8.1e-14     1.1e-14
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * beta overflow    cfs_log(beta) > MAXLOG       0.0
 *                     a or b <0 integer        0.0
 *
 */

/*							beta.c	*/


/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"

#ifdef UNK
#define MAXGAM 34.84425627277176174
#endif
#ifdef DEC
#define MAXGAM 34.84425627277176174
#endif
#ifdef IBMPC
#define MAXGAM 171.624376956302725
#endif
#ifdef MIEEE
#define MAXGAM 171.624376956302725
#endif

#ifdef ANSIPROT
extern double cfs_fabs ( double );
extern double cfs_gamma ( double );
extern double cfs_lgam ( double );
extern double cfs_exp ( double );
extern double cfs_log ( double );
extern double cfs_floor ( double );
#else
double cfs_fabs(), cfs_gamma(), cfs_lgam(), cfs_exp(), cfs_log(), cfs_floor();
#endif
extern double MAXLOG, MAXNUM;
extern int sgngam;

double cfs_beta( a, b )
double a, b;
{
double y;
int sign;

sign = 1;

if( a <= 0.0 )
	{
	if( a == cfs_floor(a) )
		goto over;
	}
if( b <= 0.0 )
	{
	if( b == cfs_floor(b) )
		goto over;
	}


y = a + b;
if( cfs_fabs(y) > MAXGAM )
	{
	y = cfs_lgam(y);
	sign *= sgngam; /* keep track of the sign */
	y = cfs_lgam(b) - y;
	sign *= sgngam;
	y = cfs_lgam(a) + y;
	sign *= sgngam;
	if( y > MAXLOG )
		{
over:
		mtherr( "beta", OVERFLOW );
		return( sign * MAXNUM );
		}
	return( sign * cfs_exp(y) );
	}

y = cfs_gamma(y);
if( y == 0.0 )
	goto over;

if( a > b )
	{
	y = cfs_gamma(a)/y;
	y *= cfs_gamma(b);
	}
else
	{
	y = cfs_gamma(b)/y;
	y *= cfs_gamma(a);
	}

return(y);
}



/* Natural cfs_log of |beta|.  Return the sign of beta in sgngam.  */

double cfs_lbeta( a, b )
double a, b;
{
double y;
int sign;

sign = 1;

if( a <= 0.0 )
	{
	if( a == cfs_floor(a) )
		goto over;
	}
if( b <= 0.0 )
	{
	if( b == cfs_floor(b) )
		goto over;
	}


y = a + b;
if( cfs_fabs(y) > MAXGAM )
	{
	y = cfs_lgam(y);
	sign *= sgngam; /* keep track of the sign */
	y = cfs_lgam(b) - y;
	sign *= sgngam;
	y = cfs_lgam(a) + y;
	sign *= sgngam;
	sgngam = sign;
	return( y );
	}

y = cfs_gamma(y);
if( y == 0.0 )
	{
over:
	mtherr( "lbeta", OVERFLOW );
	return( sign * MAXNUM );
	}

if( a > b )
	{
	y = cfs_gamma(a)/y;
	y *= cfs_gamma(b);
	}
else
	{
	y = cfs_gamma(b)/y;
	y *= cfs_gamma(a);
	}

if( y < 0 )
  {
    sgngam = -1;
    y = -y;
  }
else
  sgngam = 1;

return( cfs_log(y) );
}
