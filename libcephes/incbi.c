/*							incbi()
 *
 *      Inverse of imcomplete beta integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, md_incbi();
 *
 * x = md_incbi( a, b, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Given y, the function finds x such that
 *
 *  md_incbet( a, b, x ) = y .
 *
 * The routine performs interval halving or Newton iterations to find the
 * root of incbet(a,b,x) - y = 0.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 *                x     a,b
 * arithmetic   domain  domain  # trials    peak       rms
 *    IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13
 *    IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15
 *    IEEE      0,1     0,5       50000    1.1e-12   5.5e-15
 *    VAX       0,1    .5,100     25000    3.5e-14   1.1e-15
 * With a and b constrained to half-integer or integer values:
 *    IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13
 *    IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16
 * With a = .5, b constrained to half-integer or integer values:
 *    IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11
 */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1996, 2000 by Stephen L. Moshier
*/

#include "mconf.h"

extern double MACHEP, MAXNUM, MAXLOG, MINLOG;
#ifdef ANSIPROT
extern double md_ndtri ( double );
extern double md_exp ( double );
extern double md_fabs ( double );
extern double md_log ( double );
extern double md_sqrt ( double );
extern double md_lgam ( double );
extern double md_incbet ( double, double, double );
#else
double md_ndtri(), md_exp(), md_fabs(), md_log(), md_sqrt(), md_lgam(), md_incbet();
#endif

double md_incbi( aa, bb, yy0 )
double aa, bb, yy0;
{
double a, b, md_y0, d, y, x, x0, x1, lgm, yp, di, dithresh, yl, yh, xt;
int i, rflg, dir, nflg;


i = 0;
if( yy0 <= 0 )
	return(0.0);
if( yy0 >= 1.0 )
	return(1.0);
x0 = 0.0;
yl = 0.0;
x1 = 1.0;
yh = 1.0;
nflg = 0;

if( aa <= 1.0 || bb <= 1.0 )
	{
	dithresh = 1.0e-6;
	rflg = 0;
	a = aa;
	b = bb;
	md_y0 = yy0;
	x = a/(a+b);
	y = md_incbet( a, b, x );
	goto ihalve;
	}
else
	{
	dithresh = 1.0e-4;
	}
/* approximation to inverse function */

yp = -md_ndtri(yy0);

if( yy0 > 0.5 )
	{
	rflg = 1;
	a = bb;
	b = aa;
	md_y0 = 1.0 - yy0;
	yp = -yp;
	}
else
	{
	rflg = 0;
	a = aa;
	b = bb;
	md_y0 = yy0;
	}

lgm = (yp * yp - 3.0)/6.0;
x = 2.0/( 1.0/(2.0*a-1.0)  +  1.0/(2.0*b-1.0) );
d = yp * md_sqrt( x + lgm ) / x
	- ( 1.0/(2.0*b-1.0) - 1.0/(2.0*a-1.0) )
	* (lgm + 5.0/6.0 - 2.0/(3.0*x));
d = 2.0 * d;
if( d < MINLOG )
	{
	x = 1.0;
	goto under;
	}
x = a/( a + b * md_exp(d) );
y = md_incbet( a, b, x );
yp = (y - md_y0)/md_y0;
if( md_fabs(yp) < 0.2 )
	goto newt;

/* Resort to interval halving if not close enough. */
ihalve:

dir = 0;
di = 0.5;
for( i=0; i<100; i++ )
	{
	if( i != 0 )
		{
		x = x0  +  di * (x1 - x0);
		if( x == 1.0 )
			x = 1.0 - MACHEP;
		if( x == 0.0 )
			{
			di = 0.5;
			x = x0  +  di * (x1 - x0);
			if( x == 0.0 )
				goto under;
			}
		y = md_incbet( a, b, x );
		yp = (x1 - x0)/(x1 + x0);
		if( md_fabs(yp) < dithresh )
			goto newt;
		yp = (y-md_y0)/md_y0;
		if( md_fabs(yp) < dithresh )
			goto newt;
		}
	if( y < md_y0 )
		{
		x0 = x;
		yl = y;
		if( dir < 0 )
			{
			dir = 0;
			di = 0.5;
			}
		else if( dir > 3 )
			di = 1.0 - (1.0 - di) * (1.0 - di);
		else if( dir > 1 )
			di = 0.5 * di + 0.5; 
		else
			di = (md_y0 - y)/(yh - yl);
		dir += 1;
		if( x0 > 0.75 )
			{
			if( rflg == 1 )
				{
				rflg = 0;
				a = aa;
				b = bb;
				md_y0 = yy0;
				}
			else
				{
				rflg = 1;
				a = bb;
				b = aa;
				md_y0 = 1.0 - yy0;
				}
			x = 1.0 - x;
			y = md_incbet( a, b, x );
			x0 = 0.0;
			yl = 0.0;
			x1 = 1.0;
			yh = 1.0;
			goto ihalve;
			}
		}
	else
		{
		x1 = x;
		if( rflg == 1 && x1 < MACHEP )
			{
			x = 0.0;
			goto done;
			}
		yh = y;
		if( dir > 0 )
			{
			dir = 0;
			di = 0.5;
			}
		else if( dir < -3 )
			di = di * di;
		else if( dir < -1 )
			di = 0.5 * di;
		else
			di = (y - md_y0)/(yh - yl);
		dir -= 1;
		}
	}
mtherr( "incbi", PLOSS );
if( x0 >= 1.0 )
	{
	x = 1.0 - MACHEP;
	goto done;
	}
if( x <= 0.0 )
	{
under:
	mtherr( "incbi", UNDERFLOW );
	x = 0.0;
	goto done;
	}

newt:

if( nflg )
	goto done;
nflg = 1;
lgm = md_lgam(a+b) - md_lgam(a) - md_lgam(b);

for( i=0; i<8; i++ )
	{
	/* Compute the function at this point. */
	if( i != 0 )
		y = md_incbet(a,b,x);
	if( y < yl )
		{
		x = x0;
		y = yl;
		}
	else if( y > yh )
		{
		x = x1;
		y = yh;
		}
	else if( y < md_y0 )
		{
		x0 = x;
		yl = y;
		}
	else
		{
		x1 = x;
		yh = y;
		}
	if( x == 1.0 || x == 0.0 )
		break;
	/* Compute the derivative of the function at this point. */
	d = (a - 1.0) * md_log(x) + (b - 1.0) * md_log(1.0-x) + lgm;
	if( d < MINLOG )
		goto done;
	if( d > MAXLOG )
		break;
	d = md_exp(d);
	/* Compute the step to the next approximation of x. */
	d = (y - md_y0)/d;
	xt = x - d;
	if( xt <= x0 )
		{
		y = (x - x0) / (x1 - x0);
		xt = x0 + 0.5 * y * (x - x0);
		if( xt <= 0.0 )
			break;
		}
	if( xt >= x1 )
		{
		y = (x1 - x) / (x1 - x0);
		xt = x1 - 0.5 * y * (x1 - x);
		if( xt >= 1.0 )
			break;
		}
	x = xt;
	if( md_fabs(d/x) < 128.0 * MACHEP )
		goto done;
	}
/* Did not converge.  */
dithresh = 256.0 * MACHEP;
goto ihalve;

done:

if( rflg )
	{
	if( x <= MACHEP )
		x = 1.0 - MACHEP;
	else
		x = 1.0 - x;
	}
return( x );
}
