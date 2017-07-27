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
    XCTAssertEqualWithAccuracy(jv(v, x), .08901510322, .00000000005, @"result = %.14f", md_jn(v, x));
    XCTAssertEqualWithAccuracy(jv(v, y), -.02862625778, .00000000005, @"result = %.14f", md_jn(v, y));
    XCTAssertEqualWithAccuracy(jv(v, y), -.02862625778, .00000000005, @"result = %.14f", jv(v, y));
    XCTAssertEqualWithAccuracy(md_y0(x), .5103756726, .00000000005, @"result = %.14f", md_y0(x));
    XCTAssertEqualWithAccuracy(md_y0(y), .06264059681, .00000000005, @"result = %.14f", md_y0(y));
    XCTAssertEqualWithAccuracy(md_y1(x), -.1070324315, .00000000005, @"result = %.14f", md_y1(x));
    XCTAssertEqualWithAccuracy(md_y1(y), -.1655116144 , .00000000005, @"result = %.14f", md_y1(y));
    XCTAssertEqualWithAccuracy(md_yn(n, x), -9.935989128 , .0000000005, @"result = %.14f", md_yn(n, x));
    XCTAssertEqualWithAccuracy(md_yn(n, y), -.1000357679, .00000000005, @"result = %.14f", md_yn(n, y));
    XCTAssertEqualWithAccuracy(yv(v, x), -1.412002815 , .00000000005, @"result = %.14f", yv(v, x));
    XCTAssertEqualWithAccuracy(yv(v, y), .1773183649, .00000000005, @"result = %.14f", yv(v, y));
    XCTAssertEqualWithAccuracy(i0(x), 2.279585302, .0000000005, @"result = %.14f", i0(x));
    XCTAssertEqualWithAccuracy(i0e(y), .08978031187, .00000000005, @"result = %.14f", i0e(y));
    XCTAssertEqualWithAccuracy(i1(x), 1.590636855 , .0000000005, @"result = %.14f", i1(x));
    XCTAssertEqualWithAccuracy(i1e(y), .08750622217, .0000000001, @"result = %.14f", i1e(y));
    XCTAssertEqualWithAccuracy(iv(v, x), .1418012924, .00000000005, @"result = %.14f", iv(v, x));
    XCTAssertEqualWithAccuracy(k0(x), .1138938727, .00000000005, @"result = %.14f", k0(x));
    XCTAssertEqualWithAccuracy(k0e(y), .2785448766 , .0000000001, @"result = %.14f", k0e(y));
    XCTAssertEqualWithAccuracy(k1(x), .1398658818, .00000000005, @"result = %.14f", k1(x));
    XCTAssertEqualWithAccuracy(k1e(y), .2854254970, .0000000001, @"result = %.14f", k1e(y));
    XCTAssertEqualWithAccuracy(kn(n, x), 9.431049101, .0000000005, @"result = %.14f", kn(n, x));
}

- (void)testBetaCases
{
    double x = 5.57;
    double y = 2.2;
    double u = 0.3;
    double z = beta(x, y);
    double gammaResult = md_gamma(x) * md_gamma(y) / md_gamma(7.77);
    XCTAssertEqualWithAccuracy(z, gammaResult, 0.0000000000000005, @"result is %.20f", gammaResult);
    XCTAssertEqual( lbeta(x, y), md_log(z));
    z = incbet(x, y, u);
    XCTAssertEqualWithAccuracy(z, 0.00761009624, 0.00000000001, @"result = %0.12f", z);
    XCTAssertEqual( incbi(x, y, z), u);

}

- (void)testDistributionCases
{
    double k = 2;
    double n = 10;
    double p = 0.5;
    double y = 0.6;
    XCTAssertEqualWithAccuracy(bdtr(k, n, p), incbet(n-k, k+1, 1-p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(bdtrc(k, n, p), incbet(k+1, n-k, p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(bdtri(k, n, y), 1-incbi(n-k, k+1, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(btdtr(k, n, y), incbet(k, n, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(chdtr(k, y), igam(k/2, y/2), 0.0000000000005);
    XCTAssertEqualWithAccuracy(chdtrc(k, y), igamc(k/2, y/2), 0.0000000000005);
    XCTAssertEqualWithAccuracy(chdtri(k, y), 2 * igami(k / 2, y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(fdtr(k, n, y), incbet(k/2, n/2,k*y/(n + k*y)), 0.0000000000005);
    XCTAssertEqualWithAccuracy(fdtrc(k, n, y), incbet(n/2, k/2, n/(n + k*y)), 0.0000000000005);
    double z = incbi( n/2, k/2, p);
    XCTAssertEqualWithAccuracy(fdtri(k, n, p), n*(1-z)/(k*z), 0.0000000000005);
    XCTAssertEqualWithAccuracy(gdtr(k, n, y), igam(n, k*y), 0.0000000000005);
    XCTAssertEqualWithAccuracy(gdtrc(k, n, y), igamc(n, k*y), 0.0000000000005);
    double w = nbdtr(k, n, p);
    XCTAssertEqualWithAccuracy(w, incbet(n, k+1, p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(nbdtrc(k, n, p), incbet(k+1, n, 1-p), 0.0000000000005);
    XCTAssertEqualWithAccuracy(nbdtri(k, n, w), p, 0.0000000000005);
    w = ndtr(y);
    XCTAssertEqualWithAccuracy(w, (1+erf(y/sqrt(2)))/2, 0.0000000000005);
    XCTAssertEqualWithAccuracy(ndtri(w), y, 0.0000000000005);
    XCTAssertEqualWithAccuracy(pdtr(k, n), igamc(k+1, n), 0.0000000000005);
    XCTAssertEqualWithAccuracy(pdtrc(k, n), igam(k+1, n), 0.0000000000005);
    XCTAssertEqualWithAccuracy(pdtri(k, y), igami(k+1, y), 0.0000000000005);
    w = stdtr( k, y);
    z = k/(k + y*y);
    double zed = 1 - 0.5 * incbet(k / 2, (1.0 / 2.0), z);
    XCTAssertEqualWithAccuracy(w, zed, 0.0000000000005);
    XCTAssertEqualWithAccuracy(stdtri(k, w), y, 0.0000000000005);
}
@end
