//
//  CFSPolynomial.m
//  ObjectiveCephes
//
//  Created by Sam Krishna on 2/25/24.
//  Copyright © 2024 Sam Krishna. All rights reserved.
//

#import "protos.h"
#import "NSDecimalNumber+LJSAdditions.h"
#import "CFSPolynomial.h"

@implementation NSArray (Conversion)

- (double *)convertNSArrayToCArray {
    // Allocate memory for the C array
    double *cArray = (double *)malloc(sizeof(double) * self.count);

    // Check if memory allocation was successful
    if (cArray == NULL) {
        return NULL;
    }

    // Convert each NSNumber to double and store in the C array
    for (NSUInteger i =  0; i < self.count; i++) {
        cArray[i] = [self[i] doubleValue];
    }

    return cArray;
}

+ (NSArray<NSDecimalNumber *> *)convertCArrayToNSArray:(double *)cArray size:(NSUInteger)size {
    NSMutableArray *numbers = [NSMutableArray array];

    for (NSUInteger i = 0; i < size; i++) {
        NSDecimalNumber *num = [LjsDn dnWithDouble:cArray[i]];
        [numbers addObject:num];
    }

    return [numbers copy];
}

@end

@implementation CFSPolynomial

+ (NSArray<NSDecimalNumber *> *)cos:(NSArray<NSNumber *> *)coefficients {
    MAXPOL = 16;
    double *a = [coefficients convertNSArrayToCArray];
    double *b = (double *)malloc(sizeof(double) * MAXPOL + 1);

    for (int i = 0; i < MAXPOL + 1; i++) {
        b[i] = 0;
    }

    cfs_polcos(a, b, 16);
    NSArray *finalProduct = [NSArray convertCArrayToNSArray:b size:16];
    return finalProduct;
}

+ (NSArray<NSDecimalNumber *> *)sin:(NSArray<NSNumber *> *)coefficients {
    MAXPOL = 16;
    double *a = [coefficients convertNSArrayToCArray];
    double *b = (double *)malloc(sizeof(double) * MAXPOL + 1);

    for (int i = 0; i < MAXPOL + 1; i++) {
        b[i] = 0;
    }

    cfs_polsin(a, b, 16);
    NSArray *finalProduct = [NSArray convertCArrayToNSArray:b size:16];
    return finalProduct;
}

+ (NSArray<NSDecimalNumber *> *)sqrt:(NSArray<NSNumber *> *)coefficients {
    MAXPOL = 16;
    double *a = [coefficients convertNSArrayToCArray];
    double *b = (double *)malloc(sizeof(double) * MAXPOL + 1);

    for (int i = 0; i < MAXPOL + 1; i++) {
        b[i] = 0;
    }

    cfs_polsqt(a, b, 16);
    NSArray *finalProduct = [NSArray convertCArrayToNSArray:b size:16];
    return finalProduct;
}

@end
