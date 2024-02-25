//
//  CFSPolynomial.h
//  ObjectiveCephes
//
//  Created by Sam Krishna on 2/25/24.
//  Copyright © 2024 Sam Krishna. All rights reserved.
//

@import Foundation;

NS_ASSUME_NONNULL_BEGIN

@interface NSArray (Conversion)

- (double *)convertNSArrayToCArray;
+ (NSArray<NSDecimalNumber *> *)convertCArrayToNSArray:(double *)cArray size:(NSUInteger)size;

@end

@interface CFSPolynomial : NSObject

+ (NSArray<NSDecimalNumber *> *)cos:(NSArray<NSNumber*> *)coefficients;

@end

NS_ASSUME_NONNULL_END
