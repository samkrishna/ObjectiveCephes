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
@end
