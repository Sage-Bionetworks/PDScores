//
//  PDArray.h
//  PDScores
//
//  Created by Erin Mounts on 3/2/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <sys/types.h>
@import Accelerate;

@class PDArray;
@class PDIntArray;

typedef double (^applyRealBlock)(const double element);
typedef size_t (^applyRealIntBlock)(const double element);
typedef double (^applyRealArrayBlock)(const double element, const double otherArrayElement);
typedef size_t (^applyRealArrayIntBlock)(const double element, const double otherArrayElement);

@interface PDArray : NSObject

@property (nonatomic, assign) void *data;
@property (nonatomic, assign) size_t rows;
@property (nonatomic, assign) size_t cols;
@property (nonatomic, assign) size_t allocatedSize;

- (size_t)typeSize;

- (BOOL)setRows:(size_t)rows columns:(size_t)columns;

// all arrays must be of the same PDXxxArray subclass, and must have cols == 1
- (BOOL)concatenateColumnVectors:(NSArray *)arrays;

// all arrays must be of the same PDXxxArray subclass, must have cols == 1, and must have the same number of rows (or 0)
- (BOOL)addColumns:(NSArray *)arrays;

- (instancetype)subarrayWithRows:(NSRange)rows columns:(NSRange)columns;

- (PDIntArray *)find;
- (PDIntArray *)findFirst:(size_t)howMany;
- (PDIntArray *)findLast:(size_t)howMany;

- (instancetype)transpose;

// internal method used by find--overridden by subclasses
- (BOOL)isZero:(void *)valPtr;

@end

@interface PDComplexArray : PDArray

@property (nonatomic, assign) DSPDoubleComplex *data;

@end

@interface PDRealArray : PDArray

@property (nonatomic, assign) double *data;

+ (PDRealArray *)rowVectorWithStart:(double)start step:(double)step cap:(double)cap;
- (PDRealArray *)applyReal:(applyRealBlock)block;
- (PDRealArray *)applyReal:(applyRealArrayBlock)block withRealArray:(PDRealArray *)array;
- (PDIntArray *)applyInt:(applyRealIntBlock)block;
- (PDIntArray *)applyInt:(applyRealArrayIntBlock)block withRealArray:(PDRealArray *)array;
- (PDRealArray *)elementsWithIndices:(PDIntArray *)indexArray;
- (PDRealArray *)min;
- (PDRealArray *)minAndIndices:(PDIntArray **)indices;
- (PDRealArray *)max;
- (PDRealArray *)maxAndIndices:(PDIntArray **)indices;
- (PDRealArray *)mean;
- (PDRealArray *)median;
- (PDRealArray *)iqr;
- (PDRealArray *)var;
- (PDRealArray *)diff;
- (PDRealArray *)cumsum;
- (PDRealArray *)sum;
- (PDRealArray *)sum2; // sums along rows instead of down cols. MATLAB: sum(x, 2)
- (PDRealArray *)sqrt;
- (PDRealArray *)diag;
- (PDRealArray *)matmult:(PDRealArray *)matrix;

@end

@interface PDIntArray : PDArray

@property (nonatomic, assign) size_t *data;

@end
