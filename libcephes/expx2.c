/*							expx2.c
 *
 *	Exponential of squared argument
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, cfs_expx2();
 * int sign;
 *
 * y = cfs_expx2( x, sign );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes y = cfs_exp(x*x) while suppressing error amplification
 * that would ordinarily arise from the inexactness of the
 * exponential argument x*x.
 *
 * If sign < 0, the result is inverted; i.e., y = cfs_exp(-x*x) .
 * 
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic    domain     # trials      peak         rms
 *   IEEE      -26.6, 26.6    10^7       3.9e-16     8.9e-17
 *
 */

/*
Cephes Math Library Release 2.9:  June, 2000
Copyright 2000 by Stephen L. Moshier
*/

#include "mconf.h"

#ifdef ANSIPROT
extern double cfs_fabs (double);
extern double cfs_floor (double);
extern double cfs_exp (double);
#else
double cfs_fabs();
double cfs_floor();
double cfs_exp();
#endif

#ifdef DEC
#define M 32.0
#define MINV .03125
#else
#define M 128.0
#define MINV .0078125
#endif

extern double MAXLOG;
extern double INFINITY;

double cfs_expx2 (x, sign)
     double x;
     int sign;
{
  double u, u1, m, f;

  x = cfs_fabs (x);
  if (sign < 0)
    x = -x;

  /* Represent x as an exact multiple of M plus a residual.
     M is a power of 2 chosen so that cfs_exp(m * m) does not overflow
     or underflow and so that |x - m| is small.  */
  m = MINV * cfs_floor(M * x + 0.5);
  f = x - m;

  /* x^2 = m^2 + 2mf + f^2 */
  u = m * m;
  u1 = 2 * m * f  +  f * f;

  if (sign < 0)
    {
      u = -u;
      u1 = -u1;
    }

  if ((u+u1) > MAXLOG)
    return (INFINITY);

  /* u is exact, u1 is small.  */
  u = cfs_exp(u) * cfs_exp(u1);
  return(u);
}
