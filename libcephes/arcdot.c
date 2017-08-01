/*							arcdot.c
 *
 *	Angle between two vectors
 *
 *
 *
 *
 * SYNOPSIS:
 *
 * double p[3], q[3], cfs_arcdot();
 *
 * y = cfs_arcdot( p, q );
 *
 *
 *
 * DESCRIPTION:
 *
 * For two vectors p, q, the angle A between them is given by
 *
 *      p.q / (|p| |q|)  = cfs_cos A  .
 *
 * where "." represents inner product, "|x|" the length of vector x.
 * If the angle is small, an expression in cfs_sin A is preferred.
 * Set r = q - p.  Then
 *
 *     p.q = p.p + p.r ,
 *
 *     |p|^2 = p.p ,
 *
 *     |q|^2 = p.p + 2 p.r + r.r ,
 *
 *                  p.p^2 + 2 p.p p.r + p.r^2
 *     cfs_cos^2 A  =  ----------------------------
 *                    p.p (p.p + 2 p.r + r.r)
 *
 *                  p.p + 2 p.r + p.r^2 / p.p
 *              =  --------------------------- ,
 *                     p.p + 2 p.r + r.r
 *
 *     cfs_sin^2 A  =  1 - cfs_cos^2 A
 *
 *                   r.r - p.r^2 / p.p
 *              =  --------------------
 *                  p.p + 2 p.r + r.r
 *
 *              =   (r.r - p.r^2 / p.p) / q.q  .
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -1, 1        10^6       1.7e-16     4.2e-17
 *
 */

/*
Cephes Math Library Release 2.3:  November, 1995
Copyright 1995 by Stephen L. Moshier
*/

#include "mconf.h"
#ifdef ANSIPROT
extern double cfs_sqrt ( double );
extern double cfs_acos ( double );
extern double cfs_asin ( double );
extern double cfs_atan ( double );
#else
double cfs_sqrt(), cfs_acos(), cfs_asin(), cfs_atan();
#endif
extern double PI;

double cfs_arcdot(p,q)
double p[], q[];
{
double pp, pr, qq, rr, rt, pt, qt, pq;
int i;

pq = 0.0;
qq = 0.0;
pp = 0.0;
pr = 0.0;
rr = 0.0;
for (i=0; i<3; i++)
  {
    pt = p[i];
    qt = q[i];
    pq += pt * qt;
    qq += qt * qt;
    pp += pt * pt;
    rt = qt - pt;
    pr += pt * rt;
    rr += rt * rt;
  }
if (rr == 0.0 || pp == 0.0 || qq == 0.0)
  return 0.0;
rt = (rr - (pr * pr) / pp) / qq;
if (rt <= 0.75)
  {
    rt = cfs_sqrt(rt);
    qt = cfs_asin(rt);
    if (pq < 0.0)
      qt = PI - qt;
  }
else
  {
    pt = pq / cfs_sqrt(pp*qq);
    qt = cfs_acos(pt);
  }
return qt;
}
