/*							igami()
 *
 *      Inverse of complemented imcomplete cfs_gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, p, cfs_igami();
 *
 * x = cfs_igami( a, p );
 *
 * DESCRIPTION:
 *
 * Given p, the function finds x such that
 *
 * It is valid in the right-hand-tail of the distribution, p < 0.5.
 *  cfs_igamc( a, x ) = p.
 *
 * Starting with the approximate value
 *
 *         3
 *  x = a t
 *
 *  where
 *
 *  t = 1 - d - cfs_ndtri(p) cfs_sqrt(d)
 * 
 * and
 *
 *  d = 1/9a,
 *
 * the routine performs up to 10 Newton iterations to find the
 * root of cfs_igamc(a,x) - p = 0.
 *
 * ACCURACY:
 *
 * Tested at random a, p in the intervals indicated.
 *
 *                a        p                      Relative error:
 * arithmetic   domain   domain     # trials      peak         rms
 *    IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
 *    IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
 *    IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*/

#include "mconf.h"

extern double MACHEP, MAXNUM, MAXLOG, MINLOG;
#ifdef ANSIPROT
extern double cfs_igamc ( double, double );
extern double cfs_ndtri ( double );
extern double cfs_exp ( double );
extern double cfs_fabs ( double );
extern double cfs_log ( double );
extern double cfs_sqrt ( double );
extern double cfs_lgam ( double );
#else
double cfs_igamc(), cfs_ndtri(), cfs_exp(), cfs_fabs(), cfs_log(), cfs_sqrt(), cfs_lgam();
#endif

double cfs_igami( a, cfs_y0 )
double a, cfs_y0;
{
double x0, x1, x, yl, yh, y, d, lgm, dithresh;
int i, dir;

if( cfs_y0 > 0.5)
    mtherr( "igami", PLOSS);

/* bound the solution */
x0 = MAXNUM;
yl = 0;
x1 = 0;
yh = 1.0;
dithresh = 5.0 * MACHEP;

/* approximation to inverse function */
d = 1.0/(9.0*a);
y = ( 1.0 - d - cfs_ndtri(cfs_y0) * cfs_sqrt(d) );
x = a * y * y * y;

lgm = cfs_lgam(a);

for( i=0; i<10; i++ )
	{
	if( x > x0 || x < x1 )
		goto ihalve;
	y = cfs_igamc(a,x);
	if( y < yl || y > yh )
		goto ihalve;
	if( y < cfs_y0 )
		{
		x0 = x;
		yl = y;
		}
	else
		{
		x1 = x;
		yh = y;
		}
/* compute the derivative of the function at this point */
	d = (a - 1.0) * cfs_log(x) - x - lgm;
	if( d < -MAXLOG )
		goto ihalve;
	d = -cfs_exp(d);
/* compute the step to the next approximation of x */
	d = (y - cfs_y0)/d;
	if( cfs_fabs(d/x) < MACHEP )
		goto done;
	x = x - d;
	}

/* Resort to interval halving if Newton iteration did not converge. */
ihalve:

d = 0.0625;
if( x0 == MAXNUM )
	{
	if( x <= 0.0 )
		x = 1.0;
	while( x0 == MAXNUM )
		{
		x = (1.0 + d) * x;
		y = cfs_igamc( a, x );
		if( y < cfs_y0 )
			{
			x0 = x;
			yl = y;
			break;
			}
		d = d + d;
		}
	}
d = 0.5;
dir = 0;

for( i=0; i<400; i++ )
	{
	x = x1  +  d * (x0 - x1);
	y = cfs_igamc( a, x );
	lgm = (x0 - x1)/(x1 + x0);
	if( cfs_fabs(lgm) < dithresh )
		break;
	lgm = (y - cfs_y0)/cfs_y0;
	if( cfs_fabs(lgm) < dithresh )
		break;
	if( x <= 0.0 )
		break;
	if( y >= cfs_y0 )
		{
		x1 = x;
		yh = y;
		if( dir < 0 )
			{
			dir = 0;
			d = 0.5;
			}
		else if( dir > 1 )
			d = 0.5 * d + 0.5; 
		else
			d = (cfs_y0 - yl)/(yh - yl);
		dir += 1;
		}
	else
		{
		x0 = x;
		yl = y;
		if( dir > 0 )
			{
			dir = 0;
			d = 0.5;
			}
		else if( dir < -1 )
			d = 0.5 * d;
		else
			d = (cfs_y0 - yl)/(yh - yl);
		dir -= 1;
		}
	}
if( x == 0.0 )
	mtherr( "igami", UNDERFLOW );

done:
return( x );
}
