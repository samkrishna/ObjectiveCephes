#include "mconf.h"

typedef struct {
    double n;
    double d;
} fract;

extern double MACHEP;
extern double MAXLOG;
extern double MINLOG;
extern double MAXNUM;
extern double PI;
extern double PIO2;
extern double PIO4;
extern double SQRT2;
extern double SQRTH;
extern double LOG2E;
extern double SQ2OPI;
extern double LOGE2; 
extern double LOGSQ2;
extern double THPIO4;
extern double TWOOPI;

extern int MAXPOL;

extern double cfs_acosh ( double x );
extern int cfs_airy ( double x, double *y, double *z, double *u, double *v );
extern double cfs_asin ( double x );
extern double cfs_acos ( double x );
extern double cfs_arcdot ( double *p, double *q );
extern double cfs_asinh ( double x );
extern double cfs_atan ( double x );
extern double cfs_atan2 ( double y, double x );
extern double cfs_atanh ( double x );
extern double cfs_bdtrc ( int k, int n, double p );
extern double cfs_bdtr ( int k, int n, double p );
extern double cfs_bdtri ( int k, int n, double y );
extern double cfs_beta ( double a, double b );
extern double cfs_lbeta ( double a, double b );
extern double cfs_btdtr ( double a, double b, double x );
extern double cfs_cbrt ( double x );
extern double cfs_chbevl ( double x, void *P, int n );
extern double cfs_chdtrc ( double df, double x );
extern double cfs_chdtr ( double df, double x );
extern double cfs_chdtri ( double df, double y );
extern void cfs_clog ( cmplx *z, cmplx *w );
extern void cfs_cexp ( cmplx *z, cmplx *w );
extern void cfs_csin ( cmplx *z, cmplx *w );
extern void cfs_ccos ( cmplx *z, cmplx *w );
extern void cfs_ctan ( cmplx *z, cmplx *w );
extern void cfs_ccot ( cmplx *z, cmplx *w );
extern void cfs_casin ( cmplx *z, cmplx *w );
extern void cfs_cacos ( cmplx *z, cmplx *w );
extern void cfs_catan ( cmplx *z, cmplx *w );
extern void cfs_csinh ( cmplx *z, cmplx *w );
extern void cfs_casinh ( cmplx *z, cmplx *w );
extern void cfs_ccosh ( cmplx *z, cmplx *w );
extern void cfs_cacosh ( cmplx *z, cmplx *w );
extern void cfs_ctanh ( cmplx *z, cmplx *w );
extern void cfs_catanh ( cmplx *z, cmplx *w );
extern void cfs_cpow ( cmplx *a, cmplx *z, cmplx *w );
extern void cfs_radd ( fract *a, fract *b, fract *c );
extern void cfs_rsub ( fract *a, fract *b, fract *c );
extern void cfs_rmul ( fract *a, fract *b, fract *c );
extern void cfs_rdiv ( fract *a, fract *b, fract *c );
extern double cfs_euclid ( double *x, double *y);
extern void cfs_cadd ( cmplx *a, cmplx *b, cmplx *c );
extern void cfs_csub ( cmplx *a, cmplx *b, cmplx *c );
extern void cfs_cmul ( cmplx *a, cmplx *b, cmplx *c );
extern void cfs_cdiv ( cmplx *a, cmplx *b, cmplx *c );
extern void cfs_cmov ( void *a, void *b );
extern void cfs_cneg ( cmplx *a );
extern double cfs_cabs ( cmplx *z );
extern void cfs_csqrt ( cmplx *z, cmplx *w );
extern double cfs_hypot ( double x, double y );
extern double cfs_cosh ( double x );
extern double cfs_dawsn ( double xx );
extern double cfs_ellie ( double phi, double m );
extern double cfs_ellik ( double phi, double m );
extern double cfs_ellpe ( double x );
extern int cfs_ellpj ( double u, double m, double *x, double *y, double *z, double *a );
extern double cfs_ellpk ( double x );
extern double cfs_exp ( double x );
extern double cfs_exp10 ( double x );
/* extern double cfs_exp1m ( double x ); */
extern double cfs_exp2 ( double x );
extern double cfs_expn ( int n, double x );
extern double cfs_ei ( double x );
extern double cfs_fabs ( double x );
extern double cfs_fac ( int i );
extern double cfs_fdtrc ( int ia, int ib, double x );
extern double cfs_fdtr ( int ia, int ib, double x );
extern double cfs_fdtri ( int ia, int ib, double y );
extern double cfs_ceil ( double x );
extern double cfs_floor ( double x );
extern double cfs_frexp ( double x, int *n);
/* extern double cfs_frexp ( double x, int *pw2 ); */
extern double cfs_ldexp ( double x, int pw2 );
/* extern int cfs_signbit ( double x ); */
/* extern int cfs_isnan ( double x ); */
/* extern int cfs_isfinite ( double x ); */
extern int cfs_fresnl ( double xxa, double *x, double *y);
extern double cfs_gamma ( double x );
extern double cfs_lgam ( double x );
extern double cfs_gdtr ( double a, double b, double x );
extern double cfs_gdtrc ( double a, double b, double x );
extern double cfs_hyp2f1 ( double a, double b, double c, double x );
extern double cfs_hyperg ( double a, double b, double x );
extern double cfs_hyp2f0 ( double a, double b, double x, int type, double *y );
extern double cfs_i0 ( double x );
extern double cfs_i0e ( double x );
extern double cfs_i1 ( double x );
extern double cfs_i1e ( double x );
extern double cfs_igamc ( double a, double x );
extern double cfs_igam ( double a, double x );
extern double cfs_igami ( double a, double y0 );
extern double cfs_incbet ( double aa, double bb, double xx );
extern double cfs_incbi ( double aa, double bb, double yy0 );
extern double cfs_iv ( double v, double x );
extern double cfs_j0 ( double x );
extern double cfs_y0 ( double x );
extern double cfs_j1 ( double x );
extern double cfs_y1 ( double x );
extern double cfs_jn ( int n, double x );
extern double cfs_jv ( double n, double x );
extern double cfs_k0 ( double x );
extern double cfs_k0e ( double x );
extern double cfs_k1 ( double x );
extern double cfs_k1e ( double x );
extern double cfs_kn ( int nn, double x );
extern double cfs_log ( double x );
extern double cfs_log10 ( double x );
extern double cfs_log2 ( double x );
extern long cfs_lrand ( void );
extern long cfs_lsqrt ( long x );
extern int cfs_mtherr ( char *name, int code );
extern double cfs_polevl ( double x, void *P, int N );
extern double cfs_p1evl ( double x, void *P, int N );
extern double cfs_nbdtrc ( int k, int n, double p );
extern double cfs_nbdtr ( int k, int n, double p );
extern double cfs_nbdtri ( int k, int n, double p );
extern double cfs_ndtr ( double a );
extern double cfs_erfc ( double a );
extern double cfs_erf ( double x );
extern double cfs_ndtri ( double y0 );
extern double cfs_pdtrc ( int k, double m );
extern double cfs_pdtr ( int k, double m );
extern double cfs_pdtri ( int k, double y );
extern double cfs_pow ( double x, double y );
extern double cfs_powi ( double x, int nn );
extern double cfs_psi ( double x );
extern double cfs_rgamma ( double x );
extern double cfs_round ( double x );
extern int cfs_shichi ( double x, double *y, double *z );
extern int cfs_sici ( double x, double *y, double *z );
extern double cfs_sin ( double x );
extern double cfs_cos ( double x );
extern double cfs_radian ( double d, double m, double s );
/*
extern int cfs_sincos ( double x, double *y, double *z, int flg );
*/
extern double cfs_sindg ( double x );
extern double cfs_cosdg ( double x );
extern double cfs_sinh ( double x );
extern double cfs_spence ( double x );
extern double cfs_sqrt ( double x );
extern double cfs_stdtr ( int k, double t );
extern double cfs_stdtri ( int k, double p );
extern double cfs_onef2 ( double a, double b, double c, double x, double *y );
extern double cfs_threef0 ( double a, double b, double c, double x, double *y );
extern double cfs_struve ( double v, double x );
extern double cfs_tan ( double x );
extern double cfs_cot ( double x );
extern double cfs_tandg ( double x );
extern double cfs_cotdg ( double x );
extern double cfs_tanh ( double x );
extern double cfs_log1p ( double x );
extern double cfs_expm1 ( double x );
extern double cfs_cosm1 ( double x );
extern double cfs_yn ( int n, double x );
extern double cfs_yv ( double n, double x );
extern double cfs_zeta ( double x, double q );
extern double cfs_zetac ( double x );
extern int cfs_drand ( double *x );

extern double cfs_plancki(double w, double T);

extern void cfs_polini( int maxdeg );
extern void cfs_polclr ( double * A, int n);
extern void cfs_polmov ( double * A, int na, double * B );
extern void cfs_polmul ( double * A, int na, double * B, int nb, double * C );
extern int cfs_poldiv ( double * A, int na, double * B, int nb, double * C);
extern void cfs_poladd ( double * A, int na, double * B, int nb, double * C );
extern void cfs_polsub ( double * A, int na, double * B, int nb, double * C );
extern void cfs_polsbt ( double * A, int na, double * B, int nb, double * C );
extern void cfs_polprt ( double * A, int na, int d );
extern double cfs_poleva (double * A, int na, double x);
extern void cfs_polatn(double * A, double * B, double * C, int n);
extern void cfs_polsqt(double * A, double * B, int n);
extern void cfs_polsin(double * A, double * B, int n);
extern void cfs_polcos(double * A, double * B, int n);
extern int cfs_polrt_wrap(double * xcof, double * cof, int m, double * r, double * i);

extern void cfs_bernum_wrap(double * num, double * den);
extern double cfs_simpsn_wrap(double * f, int n, double h);
extern int cfs_minv(double * A, double * X, int n, double * B, int * IPS);
extern void cfs_mtransp(int n, double * A, double * X);
extern void cfs_eigens(double * A, double * EV, double * E, int n);
extern int cfs_simq(double * A, double * B, double * X, int n, int flag, int * IPS);
extern double cfs_polylog(int n, double x);
extern double cfs_arcdot(double * p, double * q);
extern double cfs_expx2(double x, int sign);
