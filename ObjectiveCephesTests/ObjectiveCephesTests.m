//
//  ObjectiveCephesTests.m
//  ObjectiveCephesTests
//
//  Created by Sam Krishna on 7/25/17.
//  Copyright © 2017 Sam Krishna. All rights reserved.
//

#include "protos.h"
#import <XCTest/XCTest.h>

@interface ObjectiveCephesTests : XCTestCase

@end

@implementation ObjectiveCephesTests

- (void)testBesselCases {
    // This is modified from the bessel.t code from the Math-Cephes Perl Module

    int x = 2;
    int y = 20;
    int n = 5;
    double v = 3.3;

    XCTAssertEqualWithAccuracy(cfs_j0(x), .2238907791, .00000000005, @"result = %.14f", cfs_j0(x));
    XCTAssertEqualWithAccuracy(cfs_j0(y), .1670246643, .00000000005, @"result = %.14f", cfs_j0(y));
    XCTAssertEqualWithAccuracy(cfs_j1(x), .5767248078, .00000000005, @"result = %.14f", cfs_j0(x));
    XCTAssertEqualWithAccuracy(cfs_j1(y), .06683312418, .00000000005, @"result = %.14f", cfs_j0(y));
    XCTAssertEqualWithAccuracy(cfs_jn(n, x), .007039629756, .00000000005, @"result = %.14f", cfs_jn(n, x));
    XCTAssertEqualWithAccuracy(cfs_jn(n, y), .1511697680, .00000000005, @"result = %.14f", cfs_jn(n, y));
    XCTAssertEqualWithAccuracy(cfs_jv(v, x), .08901510322, .00000000005, @"result = %.14f", cfs_jn(v, x));
    XCTAssertEqualWithAccuracy(cfs_jv(v, y), -.02862625778, .00000000005, @"result = %.14f", cfs_jn(v, y));
    XCTAssertEqualWithAccuracy(cfs_jv(v, y), -.02862625778, .00000000005, @"result = %.14f", cfs_jv(v, y));
    XCTAssertEqualWithAccuracy(cfs_y0(x), .5103756726, .00000000005, @"result = %.14f", cfs_y0(x));
    XCTAssertEqualWithAccuracy(cfs_y0(y), .06264059681, .00000000005, @"result = %.14f", cfs_y0(y));
    XCTAssertEqualWithAccuracy(cfs_y1(x), -.1070324315, .00000000005, @"result = %.14f", cfs_y1(x));
    XCTAssertEqualWithAccuracy(cfs_y1(y), -.1655116144 , .00000000005, @"result = %.14f", cfs_y1(y));
    XCTAssertEqualWithAccuracy(cfs_yn(n, x), -9.935989128 , .0000000005, @"result = %.14f", cfs_yn(n, x));
    XCTAssertEqualWithAccuracy(cfs_yn(n, y), -.1000357679, .00000000005, @"result = %.14f", cfs_yn(n, y));
    XCTAssertEqualWithAccuracy(cfs_yv(v, x), -1.412002815 , .00000000005, @"result = %.14f", cfs_yv(v, x));
    XCTAssertEqualWithAccuracy(cfs_yv(v, y), .1773183649, .00000000005, @"result = %.14f", cfs_yv(v, y));
    XCTAssertEqualWithAccuracy(cfs_i0(x), 2.279585302, .0000000005, @"result = %.14f", cfs_i0(x));
    XCTAssertEqualWithAccuracy(cfs_i0e(y), .08978031187, .00000000005, @"result = %.14f", cfs_i0e(y));
    XCTAssertEqualWithAccuracy(cfs_i1(x), 1.590636855 , .0000000005, @"result = %.14f", cfs_i1(x));
    XCTAssertEqualWithAccuracy(cfs_i1e(y), .08750622217, .0000000001, @"result = %.14f", cfs_i1e(y));
    XCTAssertEqualWithAccuracy(cfs_iv(v, x), .1418012924, .00000000005, @"result = %.14f", cfs_iv(v, x));
    XCTAssertEqualWithAccuracy(cfs_k0(x), .1138938727, .00000000005, @"result = %.14f", cfs_k0(x));
    XCTAssertEqualWithAccuracy(cfs_k0e(y), .2785448766 , .0000000001, @"result = %.14f", cfs_k0e(y));
    XCTAssertEqualWithAccuracy(cfs_k1(x), .1398658818, .00000000005, @"result = %.14f", cfs_k1(x));
    XCTAssertEqualWithAccuracy(cfs_k1e(y), .2854254970, .0000000001, @"result = %.14f", cfs_k1e(y));
    XCTAssertEqualWithAccuracy(cfs_kn(n, x), 9.431049101, .0000000005, @"result = %.14f", cfs_kn(n, x));
}

- (void)testBetaCases
{
    double x = 5.57;
    double y = 2.2;
    double u = 0.3;
    double z = cfs_beta(x, y);
    double gammaResult = cfs_gamma(x) * cfs_gamma(y) / cfs_gamma(7.77);
    XCTAssertEqualWithAccuracy(z, gammaResult, 0.0000000000000005, @"result is %.20f", gammaResult);
    XCTAssertEqual(cfs_lbeta(x, y), cfs_log(z));
    z = cfs_incbet(x, y, u);
    XCTAssertEqualWithAccuracy(z, 0.00761009624, 0.00000000001, @"result = %0.12f", z);
    XCTAssertEqual(cfs_incbi(x, y, z), u);

}

- (void)testDistributionCases
{
    double k = 2;
    double n = 10;
    double p = 0.5;
    double y = 0.6;
    XCTAssertEqualWithAccuracy(cfs_bdtr(k, n, p), cfs_incbet(n-k, k+1, 1-p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_bdtrc(k, n, p), cfs_incbet(k+1, n-k, p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_bdtri(k, n, y), 1-cfs_incbi(n-k, k+1, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_btdtr(k, n, y), cfs_incbet(k, n, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_chdtr(k, y), cfs_igam(k/2, y/2), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_chdtrc(k, y), cfs_igamc(k/2, y/2), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_chdtri(k, y), 2 * cfs_igami(k / 2, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_fdtr(k, n, y), cfs_incbet(k/2, n/2,k*y/(n + k*y)), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_fdtrc(k, n, y), cfs_incbet(n/2, k/2, n/(n + k*y)), 0.0000000000005);
    double z = cfs_incbi( n/2, k/2, p);
    XCTAssertEqualWithAccuracy(cfs_fdtri(k, n, p), n*(1-z)/(k*z), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_gdtr(k, n, y), cfs_igam(n, k*y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_gdtrc(k, n, y), cfs_igamc(n, k*y), 0.0000000000005);
    double w = cfs_nbdtr(k, n, p);
    XCTAssertEqualWithAccuracy(w, cfs_incbet(n, k+1, p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_nbdtrc(k, n, p), cfs_incbet(k+1, n, 1-p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_nbdtri(k, n, w), p, 0.0000000000005);
    w = cfs_ndtr(y);
    XCTAssertEqualWithAccuracy(w, (1+cfs_erf(y/cfs_sqrt(2)))/2, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_ndtri(w), y, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_pdtr(k, n), cfs_igamc(k+1, n), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_pdtrc(k, n), cfs_igam(k+1, n), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_pdtri(k, y), cfs_igami(k+1, y), 0.0000000000005);
    w = cfs_stdtr( k, y);
    z = k/(k + y*y);
    double zed = 1 - 0.5 * cfs_incbet(k / 2, (1.0 / 2.0), z);
    XCTAssertEqualWithAccuracy(w, zed, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_stdtri(k, w), y, 0.0000000000005);
}

- (void)testEllipticCases
{
    double x = 0.3;
    XCTAssertEqualWithAccuracy(cfs_ellpk(1-x*x), 1.608048620, 0.0000000005);
    XCTAssertEqualWithAccuracy(cfs_ellik(cfs_asin(0.2), x*x), .2014795901, 0.0000000005);
    XCTAssertEqualWithAccuracy(cfs_ellpe(1-x*x), 1.534833465, 0.0000000005);
    XCTAssertEqualWithAccuracy(cfs_ellie(cfs_asin(0.2), x*x), .2012363833, 0.0000000005);
    double phi = PIO4;
    double m = 0.3;
    double u = cfs_ellik(phi, m);

    // 2017-07-27 14:12:00 -0700
    // There's some weird Clang shit going on that causes the declaration of the double pointers
    // to inconsistently be initialized NULL/not NULL. As a result, this causes an EXC_BAD_ACCESS in the
    // guts of the cfs_ellpj() code. I had to explicitly malloc the pointers to get past this issue.
    double *sn = (double *)malloc(sizeof(double));
    double *cn = (double *)malloc(sizeof(double));
    double *dn = (double *)malloc(sizeof(double));
    double *phi_out = (double *)malloc(sizeof(double));
    int flag = cfs_ellpj(u, m, sn, cn, dn, phi_out);

    XCTAssertEqualWithAccuracy(flag, 0, 0.0000000000005);
    XCTAssertEqualWithAccuracy(phi, *phi_out, 0.0000000000005);
    XCTAssertEqualWithAccuracy(*sn, cfs_sin(*phi_out), 0.0000000000005);
    XCTAssertEqualWithAccuracy(*cn, cfs_cos(*phi_out), 0.0000000000005);
    XCTAssertEqualWithAccuracy(*dn, cfs_sqrt(1-m * cfs_sin(*phi_out) * cfs_sin(*phi_out)), 0.0000000000005);

    free(sn);
    free(cn);
    free(dn);
    free(phi_out);
}

- (void)testExpLogCases
{
    double e = cfs_exp(1);
    XCTAssertEqualWithAccuracy(cfs_log(cfs_pow(e, e)), e, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_log(e*e), 2, 0.0000000000005);
    XCTAssertEqualWithAccuracy(1 / cfs_log(2), LOG2E, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_exp(-1), 1/e, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_exp(LOGE2), 2, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_log10(10000), 4, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_log10(cfs_sqrt(10)), 0.5, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_exp2(8), 256, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_log2(SQRT2), 0.5, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_log2(256), 8, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_log1p(0.5), cfs_log(1.5), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_expm1(0.5), cfs_exp(0.5)-1, 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_expx2(2, -1), cfs_exp(-4), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_expx2(0.5, 1), cfs_exp(0.25), 0.0000000000005);
    XCTAssertEqualWithAccuracy(cfs_exp2(-0.5), SQRTH, 0.01);
}

- (void)testComplexCases
{
    cmplx x = {5, 6};
    cmplx y = {1, 3};
    cmplx z;
    double tolerance = 0.00000000005;

    cfs_cadd(&x, &y, &z);
    XCTAssertTrue(z.r == 6);
    XCTAssertTrue(z.i == 9);

    cfs_csub(&y, &x, &z);
    XCTAssertEqualWithAccuracy(z.r, 4.0, tolerance);
    XCTAssertEqualWithAccuracy(z.i, 3.0, tolerance);

    cfs_cmul(&y, &x, &z);
    XCTAssertEqualWithAccuracy(z.r, -13.0, tolerance);
    XCTAssertEqualWithAccuracy(z.i, 21.0, tolerance);

    cfs_cdiv(&y, &x, &z);
    XCTAssertEqualWithAccuracy(z.r, 2.3, tolerance);
    XCTAssertEqualWithAccuracy(z.i, -0.9, tolerance);

    cfs_cneg(&z);
    XCTAssertEqualWithAccuracy(z.r, -2.3, tolerance);
    XCTAssertEqualWithAccuracy(z.i, 0.9, tolerance);

    cfs_cmov(&x, &z);
    XCTAssertEqualWithAccuracy(z.r, 5.0, tolerance);
    XCTAssertEqualWithAccuracy(z.i, 6, tolerance);
    XCTAssertEqualWithAccuracy(cfs_cabs(&z), sqrt(61), tolerance);

    cfs_clog(&x, &z);
    XCTAssertEqualWithAccuracy(z.r, log(hypot(5, 6)), tolerance);
    XCTAssertEqualWithAccuracy(z.i, atan2(6, 5), tolerance);

    cfs_cexp(&x, &z);
    XCTAssertEqualWithAccuracy(z.r, exp(5)*cos(6), tolerance);
    XCTAssertEqualWithAccuracy(z.i, exp(5)*sin(6), tolerance);

    cfs_csin(&x, &z);
    cmplx d = { (sin(5) * cosh(6)), (cos(5) * sinh(6)) };
    XCTAssertEqualWithAccuracy(z.r, d.r, tolerance);
    XCTAssertEqualWithAccuracy(z.i, d.i, tolerance);

    cfs_casin(&d, &z);
    XCTAssertEqualWithAccuracy(z.r, 5-2*PI, tolerance);
    XCTAssertEqualWithAccuracy(z.i, 6.0, tolerance);

    d.r = cos(5)*cosh(6);
    d.i = -sin(5)*sinh(6);
    cfs_ccos(&x, &z);
    XCTAssertEqualWithAccuracy(z.r, d.r, tolerance);
    XCTAssertEqualWithAccuracy(z.i, d.i, tolerance);

    cfs_cacos(&d, &z);
    XCTAssertEqualWithAccuracy(z.r, 5-2*PI, tolerance);
    XCTAssertEqualWithAccuracy(z.i, 6.0, tolerance);

    double den = cos(10) + cosh(12);
    d.r = sin(10) / den;
    d.i = sinh(12) / den;
    cfs_ctan(&x, &z);
    XCTAssertEqualWithAccuracy(z.r, d.r, tolerance);
    XCTAssertEqualWithAccuracy(z.i, d.i, tolerance);

    cfs_catan(&d, &z);
    XCTAssertEqualWithAccuracy(z.r, 5-2*PI, tolerance);
    XCTAssertEqualWithAccuracy(z.i, 6.0, tolerance);

    cfs_ccot(&x, &z);
    den = cosh(12) - cos(10);
    XCTAssertEqualWithAccuracy(z.r, (sin(10) / den), tolerance);
    XCTAssertEqualWithAccuracy(z.i, (-sinh(12) / den), tolerance);

    cfs_csqrt(&x, &z);
    XCTAssertEqualWithAccuracy(z.r, (3/z.i), tolerance);
    XCTAssertEqualWithAccuracy(z.i, sqrt((sqrt(61) - 5) / 2), tolerance);

    d.r = 2;
    d.i = 3;
    cfs_csinh(&d, &z);
    XCTAssertEqualWithAccuracy(z.r, sinh(2)*cos(3), tolerance);
    XCTAssertEqualWithAccuracy(z.i, cosh(2)*sin(3), tolerance);

    cfs_casinh(&z, &y);
    XCTAssertEqualWithAccuracy(y.r, 2.0, tolerance);
    XCTAssertEqualWithAccuracy(y.i, 3.0, tolerance);

    cfs_ccosh(&d, &z);
    XCTAssertEqualWithAccuracy(z.r, cosh(2)*cos(3), tolerance);
    XCTAssertEqualWithAccuracy(z.i, sinh(2)*sin(3), tolerance);

    cfs_cacosh(&z, &y);
    XCTAssertEqualWithAccuracy(y.r, 2.0, tolerance);
    XCTAssertEqualWithAccuracy(y.i, 3.0, tolerance);

    den = cosh(4) + cos(6);
    cfs_ctanh(&d, &z);
    XCTAssertEqualWithAccuracy(z.r, sinh(4)/den, tolerance);
    XCTAssertEqualWithAccuracy(z.i, sin(6)/den, tolerance);

    cfs_catanh(&z, &y);
    XCTAssertEqualWithAccuracy(y.r, 2.0, tolerance);
    XCTAssertEqualWithAccuracy(y.i, 3-PI, tolerance);

    d.r = 4;
    d.i = 5;
    cfs_cpow(&d, &y, &z);
    cmplx c, f, g;
    cfs_clog(&d, &c);
    cfs_cmul(&y, &c, &f);
    cfs_cexp(&f, &g);
    XCTAssertEqualWithAccuracy(z.r, g.r, tolerance);
    XCTAssertEqualWithAccuracy(z.i, g.i, tolerance);

    x.r = 55;
    x.i = 66;
    XCTAssertEqualWithAccuracy(x.r, 55.0, tolerance);
    XCTAssertEqualWithAccuracy(x.i, 66.0, tolerance);
}

@end
