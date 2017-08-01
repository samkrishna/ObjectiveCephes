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

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (void)testBesselCases {
    // This is modified from the bessel.t code from the Math-Cephes Perl Module

    int x = 2;
    int y = 20;
    int n = 5;
    double v = 3.3;

    XCTAssertEqualWithAccuracy(md_j0(x), .2238907791, .00000000005, @"result = %.14f", md_j0(x));
    XCTAssertEqualWithAccuracy(md_j0(y), .1670246643, .00000000005, @"result = %.14f", md_j0(y));
    XCTAssertEqualWithAccuracy(md_j1(x), .5767248078, .00000000005, @"result = %.14f", md_j0(x));
    XCTAssertEqualWithAccuracy(md_j1(y), .06683312418, .00000000005, @"result = %.14f", md_j0(y));
    XCTAssertEqualWithAccuracy(md_jn(n, x), .007039629756, .00000000005, @"result = %.14f", md_jn(n, x));
    XCTAssertEqualWithAccuracy(md_jn(n, y), .1511697680, .00000000005, @"result = %.14f", md_jn(n, y));
    XCTAssertEqualWithAccuracy(md_jv(v, x), .08901510322, .00000000005, @"result = %.14f", md_jn(v, x));
    XCTAssertEqualWithAccuracy(md_jv(v, y), -.02862625778, .00000000005, @"result = %.14f", md_jn(v, y));
    XCTAssertEqualWithAccuracy(md_jv(v, y), -.02862625778, .00000000005, @"result = %.14f", md_jv(v, y));
    XCTAssertEqualWithAccuracy(md_y0(x), .5103756726, .00000000005, @"result = %.14f", md_y0(x));
    XCTAssertEqualWithAccuracy(md_y0(y), .06264059681, .00000000005, @"result = %.14f", md_y0(y));
    XCTAssertEqualWithAccuracy(md_y1(x), -.1070324315, .00000000005, @"result = %.14f", md_y1(x));
    XCTAssertEqualWithAccuracy(md_y1(y), -.1655116144 , .00000000005, @"result = %.14f", md_y1(y));
    XCTAssertEqualWithAccuracy(md_yn(n, x), -9.935989128 , .0000000005, @"result = %.14f", md_yn(n, x));
    XCTAssertEqualWithAccuracy(md_yn(n, y), -.1000357679, .00000000005, @"result = %.14f", md_yn(n, y));
    XCTAssertEqualWithAccuracy(md_yv(v, x), -1.412002815 , .00000000005, @"result = %.14f", md_yv(v, x));
    XCTAssertEqualWithAccuracy(md_yv(v, y), .1773183649, .00000000005, @"result = %.14f", md_yv(v, y));
    XCTAssertEqualWithAccuracy(md_i0(x), 2.279585302, .0000000005, @"result = %.14f", md_i0(x));
    XCTAssertEqualWithAccuracy(md_i0e(y), .08978031187, .00000000005, @"result = %.14f", md_i0e(y));
    XCTAssertEqualWithAccuracy(md_i1(x), 1.590636855 , .0000000005, @"result = %.14f", md_i1(x));
    XCTAssertEqualWithAccuracy(md_i1e(y), .08750622217, .0000000001, @"result = %.14f", md_i1e(y));
    XCTAssertEqualWithAccuracy(md_iv(v, x), .1418012924, .00000000005, @"result = %.14f", md_iv(v, x));
    XCTAssertEqualWithAccuracy(md_k0(x), .1138938727, .00000000005, @"result = %.14f", md_k0(x));
    XCTAssertEqualWithAccuracy(md_k0e(y), .2785448766 , .0000000001, @"result = %.14f", md_k0e(y));
    XCTAssertEqualWithAccuracy(md_k1(x), .1398658818, .00000000005, @"result = %.14f", md_k1(x));
    XCTAssertEqualWithAccuracy(md_k1e(y), .2854254970, .0000000001, @"result = %.14f", md_k1e(y));
    XCTAssertEqualWithAccuracy(md_kn(n, x), 9.431049101, .0000000005, @"result = %.14f", md_kn(n, x));
}

- (void)testBetaCases
{
    double x = 5.57;
    double y = 2.2;
    double u = 0.3;
    double z = md_beta(x, y);
    double gammaResult = md_gamma(x) * md_gamma(y) / md_gamma(7.77);
    XCTAssertEqualWithAccuracy(z, gammaResult, 0.0000000000000005, @"result is %.20f", gammaResult);
    XCTAssertEqual(md_lbeta(x, y), md_log(z));
    z = md_incbet(x, y, u);
    XCTAssertEqualWithAccuracy(z, 0.00761009624, 0.00000000001, @"result = %0.12f", z);
    XCTAssertEqual(md_incbi(x, y, z), u);

}

- (void)testDistributionCases
{
    double k = 2;
    double n = 10;
    double p = 0.5;
    double y = 0.6;
    XCTAssertEqualWithAccuracy(md_bdtr(k, n, p), md_incbet(n-k, k+1, 1-p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_bdtrc(k, n, p), md_incbet(k+1, n-k, p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_bdtri(k, n, y), 1-md_incbi(n-k, k+1, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_btdtr(k, n, y), md_incbet(k, n, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_chdtr(k, y), md_igam(k/2, y/2), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_chdtrc(k, y), md_igamc(k/2, y/2), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_chdtri(k, y), 2 * md_igami(k / 2, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_fdtr(k, n, y), md_incbet(k/2, n/2,k*y/(n + k*y)), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_fdtrc(k, n, y), md_incbet(n/2, k/2, n/(n + k*y)), 0.0000000000005);
    double z = md_incbi( n/2, k/2, p);
    XCTAssertEqualWithAccuracy(md_fdtri(k, n, p), n*(1-z)/(k*z), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_gdtr(k, n, y), md_igam(n, k*y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_gdtrc(k, n, y), md_igamc(n, k*y), 0.0000000000005);
    double w = md_nbdtr(k, n, p);
    XCTAssertEqualWithAccuracy(w, md_incbet(n, k+1, p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_nbdtrc(k, n, p), md_incbet(k+1, n, 1-p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_nbdtri(k, n, w), p, 0.0000000000005);
    w = md_ndtr(y);
    XCTAssertEqualWithAccuracy(w, (1+md_erf(y/md_sqrt(2)))/2, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_ndtri(w), y, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_pdtr(k, n), md_igamc(k+1, n), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_pdtrc(k, n), md_igam(k+1, n), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_pdtri(k, y), md_igami(k+1, y), 0.0000000000005);
    w = md_stdtr( k, y);
    z = k/(k + y*y);
    double zed = 1 - 0.5 * md_incbet(k / 2, (1.0 / 2.0), z);
    XCTAssertEqualWithAccuracy(w, zed, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_stdtri(k, w), y, 0.0000000000005);
}

- (void)testEllipticCases
{
    double x = 0.3;
    XCTAssertEqualWithAccuracy(md_ellpk(1-x*x), 1.608048620, 0.0000000005);
    XCTAssertEqualWithAccuracy(md_ellik(md_asin(0.2), x*x), .2014795901, 0.0000000005);
    XCTAssertEqualWithAccuracy(md_ellpe(1-x*x), 1.534833465, 0.0000000005);
    XCTAssertEqualWithAccuracy(md_ellie(md_asin(0.2), x*x), .2012363833, 0.0000000005);
    double phi = PIO4;
    double m = 0.3;
    double u = md_ellik(phi, m);

    // 2017-07-27 14:12:00 -0700
    // There's some weird Clang shit going on that causes the declaration of the double pointers
    // to inconsistently be initialized NULL/not NULL. As a result, this causes an EXC_BAD_ACCESS in the
    // guts of the md_ellpj() code. I had to explicitly malloc the pointers to get past this issue.
    double *sn = (double *)malloc(sizeof(double));
    double *cn = (double *)malloc(sizeof(double));
    double *dn = (double *)malloc(sizeof(double));
    double *phi_out = (double *)malloc(sizeof(double));
    int flag = md_ellpj(u, m, sn, cn, dn, phi_out);

    XCTAssertEqualWithAccuracy(flag, 0, 0.0000000000005);
    XCTAssertEqualWithAccuracy(phi, *phi_out, 0.0000000000005);
    XCTAssertEqualWithAccuracy(*sn, md_sin(*phi_out), 0.0000000000005);
    XCTAssertEqualWithAccuracy(*cn, md_cos(*phi_out), 0.0000000000005);
    XCTAssertEqualWithAccuracy(*dn, md_sqrt(1-m * md_sin(*phi_out) * md_sin(*phi_out)), 0.0000000000005);

    free(sn);
    free(cn);
    free(dn);
    free(phi_out);
}

- (void)testExpLogCases
{
    double e = md_exp(1);
    XCTAssertEqualWithAccuracy(md_log(md_pow(e, e)), e, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_log(e*e), 2, 0.0000000000005);
    XCTAssertEqualWithAccuracy(1 / md_log(2), LOG2E, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_exp(-1), 1/e, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_exp(LOGE2), 2, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_log10(10000), 4, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_log10(md_sqrt(10)), 0.5, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_exp2(8), 256, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_log2(SQRT2), 0.5, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_log2(256), 8, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_log1p(0.5), md_log(1.5), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_expm1(0.5), md_exp(0.5)-1, 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_expx2(2, -1), md_exp(-4), 0.0000000000005);
    XCTAssertEqualWithAccuracy(md_expx2(0.5, 1), md_exp(0.25), 0.0000000000005);

    // No idea what to do with this since SQRTH has the unusual unsigned short array
//    XCTAssertEqualWithAccuracy(md_exp2(-1/2), SQRTH, 0.01);
}

@end
