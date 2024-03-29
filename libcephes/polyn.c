/*							polyn.c
 *							polyr.c
 * Arithmetic operations on polynomials
 *
 * In the following descriptions a, b, c are polynomials of degree
 * na, nb, nc respectively.  The degree of a polynomial cannot
 * exceed a run-time value MAXPOL.  An operation that attempts
 * to use or generate a polynomial of higher degree may produce a
 * result that suffers truncation at degree MAXPOL.  The value of
 * MAXPOL is set by calling the function
 *
 *     polini( maxpol );
 *
 * where maxpol is the desired maximum degree.  This must be
 * done prior to calling any of the other functions in this module.
 * Memory for internal temporary polynomial storage is allocated
 * by polini().
 *
 * Each polynomial is represented by an array containing its
 * coefficients, together with a separately declared integer equal
 * to the degree of the polynomial.  The coefficients appear in
 * ascending order; that is,
 *
 *                                        2                      na
 * a(x)  =  a[0]  +  a[1] * x  +  a[2] * x   +  ...  +  a[na] * x  .
 *
 *
 *
 * sum = poleva( a, na, x );	Evaluate polynomial a(t) at t = x.
 * polprt( a, na, D );		Print the coefficients of a to D digits.
 * cfs_polclr( a, na );		Set a identically equal to zero, up to a[na].
 * cfs_polmov( a, na, b );		Set b = a.
 * cfs_poladd( a, na, b, nb, c );	c = b + a, nc = max(na,nb)
 * cfs_polsub( a, na, b, nb, c );	c = b - a, nc = max(na,nb)
 * cfs_polmul( a, na, b, nb, c );	c = b * a, nc = na+nb
 *
 *
 * Division:
 *
 * i = cfs_poldiv( a, na, b, nb, c );	c = b / a, nc = MAXPOL
 *
 * returns i = the degree of the first nonzero coefficient of a.
 * The computed quotient c must be divided by x^i.  An error message
 * is printed if a is identically zero.
 *
 *
 * Change of variables:
 * If a and b are polynomials, and t = a(x), then
 *     c(t) = b(a(x))
 * is a polynomial found by substituting a(x) for t.  The
 * subroutine call for this is
 *
 * cfs_polsbt( a, na, b, nb, c );
 *
 *
 * Notes:
 * cfs_poldiv() is an integer routine; poleva() is double.
 * Any of the arguments a, b, c may refer to the same array.
 *
 */

#include <stdio.h>
#include "mconf.h"
#ifdef ANSIPROT
void exit (int);
//extern void * malloc ( long );
extern void	*malloc(size_t __size) __result_use_check;
extern void free ( void * );
void cfs_polclr ( double *, int );
void cfs_polmov ( double *, int, double * );
void cfs_polmul ( double *, int, double *, int, double * );
int cfs_poldiv ( double *, int, double *, int, double * );
#else
void exit();
void * malloc();
void free ();
void cfs_polclr(), cfs_polmov(), cfs_polmul();
int cfs_poldiv();
#endif
#ifndef NULL
#define NULL 0
#endif

/* near pointer version of malloc() */
/*
#define malloc _nmalloc
#define free _nfree
*/

/* Pointers to internal arrays.  Note cfs_poldiv() allocates
 * and deallocates some temporary arrays every time it is called.
 */
static double *pt1 = 0;
static double *pt2 = 0;
static double *pt3 = 0;

/* Maximum degree of polynomial. */
int MAXPOL = 0;
extern int MAXPOL;

/* Number of bytes (chars) in maximum size polynomial. */
static int psize = 0;


/* Initialize max degree of polynomials
 * and allocate temporary storage.
 */
void cfs_polini( maxdeg )
int maxdeg;
{

MAXPOL = maxdeg;
psize = (maxdeg + 1) * sizeof(double);

/* Release previously allocated memory, if any. */
if( pt3 )
	free(pt3);
if( pt2 )
	free(pt2);
if( pt1 )
	free(pt1);

/* Allocate new arrays */
pt1 = (double * )malloc(psize); /* used by cfs_polsbt */
pt2 = (double * )malloc(psize); /* used by cfs_polsbt */
pt3 = (double * )malloc(psize); /* used by cfs_polmul */

/* Report if failure */
if( (pt1 == NULL) || (pt2 == NULL) || (pt3 == NULL) )
	{
	mtherr( "cfs_polini", ERANGE );
	exit(1);
	}
}



/* Print the coefficients of a, with d decimal precision.
 */
static char *form = "abcdefghijk";

void cfs_polprt( a, na, d )
double a[];
int na, d;
{
int i, j, d1;
char *p;

/* Create format descriptor string for the printout.
 * Do this partly by hand, since sprintf() may be too
 * bug-ridden to accomplish this feat by itself.
 */
p = form;
*p++ = '%';
d1 = d + 8;
sprintf( p, "%d ", d1 );
p += 1;
if( d1 >= 10 )
	p += 1;
*p++ = '.';
sprintf( p, "%d ", d );
p += 1;
if( d >= 10 )
	p += 1;
*p++ = 'e';
*p++ = ' ';
*p++ = '\0';


/* Now do the printing.
 */
d1 += 1;
j = 0;
for( i=0; i<=na; i++ )
	{
/* Detect end of available line */
	j += d1;
	if( j >= 78 )
		{
		printf( "\n" );
		j = d1;
		}
	printf( form, a[i] );
	}
printf( "\n" );
}



/* Set a = 0.
 */
void cfs_polclr( a, n )
register double *a;
int n;
{
int i;

if( n > MAXPOL )
	n = MAXPOL;

for( i=0; i<=n; i++ )
	*a++ = 0.0;
}



/* Set b = a.
 */
void cfs_polmov( a, na, b )
register double *a, *b;
int na;
{
int i;

if( na > MAXPOL )
	na = MAXPOL;

for( i=0; i<= na; i++ )
	{
	*b++ = *a++;
	}
}


/* c = b * a.
 */
void cfs_polmul( a, na, b, nb, c )
double a[], b[], c[];
int na, nb;
{
int i, j, k, nc;
double x;

nc = na + nb;
cfs_polclr( pt3, MAXPOL );

for( i=0; i<=na; i++ )
	{
	x = a[i];
	for( j=0; j<=nb; j++ )
		{
		k = i + j;
		if( k > MAXPOL )
			break;
		pt3[k] += x * b[j];
		}
	}

if( nc > MAXPOL )
	nc = MAXPOL;
for( i=0; i<=nc; i++ )
	c[i] = pt3[i];
}



 
/* c = b + a.
 */
void cfs_poladd( a, na, b, nb, c )
double a[], b[], c[];
int na, nb;
{
int i, n;


if( na > nb )
	n = na;
else
	n = nb;

if( n > MAXPOL )
	n = MAXPOL;

for( i=0; i<=n; i++ )
	{
	if( i > na )
		c[i] = b[i];
	else if( i > nb )
		c[i] = a[i];
	else
		c[i] = b[i] + a[i];
	}
}

/* c = b - a.
 */
void cfs_polsub( a, na, b, nb, c )
double a[], b[], c[];
int na, nb;
{
int i, n;


if( na > nb )
	n = na;
else
	n = nb;

if( n > MAXPOL )
	n = MAXPOL;

for( i=0; i<=n; i++ )
	{
	if( i > na )
		c[i] = b[i];
	else if( i > nb )
		c[i] = -a[i];
	else
		c[i] = b[i] - a[i];
	}
}



/* c = b/a
 */
int cfs_poldiv( a, na, b, nb, c )
double a[], b[], c[];
int na, nb;
{
double quot;
double *ta, *tb, *tq;
int i, j, k, sing;

sing = 0;

/* Allocate temporary arrays.  This would be quicker
 * if done automatically on the stack, but stack space
 * may be hard to obtain on a small computer.
 */
ta = (double * )malloc( psize );
cfs_polclr( ta, MAXPOL );
cfs_polmov( a, na, ta );

tb = (double * )malloc( psize );
cfs_polclr( tb, MAXPOL );
cfs_polmov( b, nb, tb );

tq = (double * )malloc( psize );
cfs_polclr( tq, MAXPOL );

/* What to do if leading (constant) coefficient
 * of denominator is zero.
 */
if( a[0] == 0.0 )
	{
	for( i=0; i<=na; i++ )
		{
		if( ta[i] != 0.0 )
			goto nzero;
		}
	mtherr( "cfs_poldiv", SING );
	goto done;

nzero:
/* Reduce the degree of the denominator. */
	for( i=0; i<na; i++ )
		ta[i] = ta[i+1];
	ta[na] = 0.0;

	if( b[0] != 0.0 )
		{
/* Optional message:
		printf( "cfs_poldiv singularity, divide quotient by x\n" );
*/
		sing += 1;
		}
	else
		{
/* Reduce degree of numerator. */
		for( i=0; i<nb; i++ )
			tb[i] = tb[i+1];
		tb[nb] = 0.0;
		}
/* Call self, using reduced polynomials. */
	sing += cfs_poldiv( ta, na, tb, nb, c );
	goto done;
	}

/* Long division algorithm.  ta[0] is nonzero.
 */
for( i=0; i<=MAXPOL; i++ )
	{
	quot = tb[i]/ta[0];
	for( j=0; j<=MAXPOL; j++ )
		{
		k = j + i;
		if( k > MAXPOL )
			break;
		tb[k] -= quot * ta[j];
		}
	tq[i] = quot;
	}
/* Send quotient to output array. */
cfs_polmov( tq, MAXPOL, c );

done:

/* Restore allocated memory. */
free(tq);
free(tb);
free(ta);
return( sing );
}




/* Change of variables
 * Substitute a(y) for the variable x in b(x).
 * x = a(y)
 * c(x) = b(x) = b(a(y)).
 */

void cfs_polsbt( a, na, b, nb, c )
double a[], b[], c[];
int na, nb;
{
int i, j, k, n2;
double x;

/* 0th degree term:
 */
cfs_polclr( pt1, MAXPOL );
pt1[0] = b[0];

cfs_polclr( pt2, MAXPOL );
pt2[0] = 1.0;
n2 = 0;

for( i=1; i<=nb; i++ )
	{
/* Form ith power of a. */
	cfs_polmul( a, na, pt2, n2, pt2 );
	n2 += na;
	x = b[i];
/* Add the ith coefficient of b times the ith power of a. */
	for( j=0; j<=n2; j++ )
		{
		if( j > MAXPOL )
			break;
		pt1[j] += x * pt2[j];
		}
	}

k = n2 + nb;
if( k > MAXPOL )
	k = MAXPOL;
for( i=0; i<=k; i++ )
	c[i] = pt1[i];
}




/* Evaluate polynomial a(t) at t = x.
 */
double poleva( a, na, x )
double a[];
int na;
double x;
{
double s;
int i;

s = a[na];
for( i=na-1; i>=0; i-- )
	{
	s = s * x + a[i];
	}
return(s);
}

