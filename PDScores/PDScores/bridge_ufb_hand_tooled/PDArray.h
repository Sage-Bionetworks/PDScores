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
@class PDRealArray;
@class PDIntArray;

typedef double (^applyRealBlock)(const double element);
typedef size_t (^applyRealIntBlock)(const double element);
typedef double (^applyRealArrayBlock)(const double element, const double otherArrayElement);
typedef size_t (^applyRealArrayIntBlock)(const double element, const double otherArrayElement);
typedef double (^applyIntArrayRealBlock)(const double element, const size_t otherArrayElement);

@interface PDArray : NSObject

@property (nonatomic, assign) void *data;
@property (nonatomic, assign) size_t rows;
@property (nonatomic, assign) size_t cols;
@property (nonatomic, assign) size_t allocatedSize;

- (size_t)typeSize;

- (BOOL)setRows:(size_t)rows columns:(size_t)columns;

- (BOOL)isVector;

// all arrays must be of the same PDXxxArray subclass, and must have cols == 1
- (BOOL)concatenateColumnVectors:(NSArray *)arrays;

// all arrays must be of the same PDXxxArray subclass, must have cols == 1, and must have the same number of rows (or 0)
- (BOOL)addColumns:(NSArray *)arrays;

- (instancetype)subarrayWithRows:(NSRange)rows columns:(NSRange)columns;
- (instancetype)subarrayWithRowIndices:(PDIntArray *)rows columnIndices:(PDIntArray *)columns;

// self and array types must match, and array rows and columns must match lengths of rows and columns ranges (and fit within self)
- (void)setSubarrayRows:(NSRange)rows columns:(NSRange)columns fromArray:(PDArray *)array;

// number of elements in rows * columns must match size of array (need not be the same dimensions, but data will be treated
// as if it were)
- (void)setElementsWithRowIndices:(PDIntArray *)rows columnIndices:(PDIntArray *)columns fromArray:(PDArray *)array;

- (PDIntArray *)find;
- (PDIntArray *)findFirst:(size_t)howMany;
- (PDIntArray *)findLast:(size_t)howMany;
- (instancetype)elementsWithIndices:(PDIntArray *)indexArray;

// array type must match self, and its size must match that of indexArray
- (void)setElementsWithIndices:(PDIntArray *)indexArray fromArray:(PDArray *)array;

- (instancetype)transpose;
- (instancetype)flipud;

// internal method used by find--overridden by subclasses
- (BOOL)isZero:(void *)valPtr;

@end

@interface PDComplexArray : PDArray

@property (nonatomic, assign) DSPDoubleComplex *data;

- (PDRealArray *)abs;

@end

@interface PDRealArray : PDArray

@property (nonatomic, assign) double *data;

+ (PDRealArray *)rowVectorWithStart:(double)start step:(double)step cap:(double)cap;
- (PDRealArray *)applyReal:(applyRealBlock)block;
- (PDRealArray *)applyReal:(applyRealArrayBlock)block withRealArray:(PDRealArray *)array;
- (PDRealArray *)applyReal:(applyIntArrayRealBlock)block withIntArray:(PDIntArray *)array;
- (PDIntArray *)applyInt:(applyRealIntBlock)block;
- (PDIntArray *)applyInt:(applyRealArrayIntBlock)block withRealArray:(PDRealArray *)array;
//- (PDRealArray *)elementsWithIndices:(PDIntArray *)indexArray;
- (PDRealArray *)abs;
- (PDRealArray *)round;
- (PDRealArray *)min;
- (PDRealArray *)minAndIndices:(PDIntArray **)indices;
- (PDRealArray *)max;
- (PDRealArray *)maxAndIndices:(PDIntArray **)indices;
- (PDRealArray *)mean;
- (PDRealArray *)median;
- (double)norm;
- (PDRealArray *)iqr;
- (PDRealArray *)var;
- (PDRealArray *)diff;
- (PDRealArray *)cumsum;
- (PDRealArray *)sum;
- (PDRealArray *)sum2; // sums along rows instead of down cols. MATLAB: sum(x, 2)
- (PDRealArray *)square;
- (PDRealArray *)sqrt;
- (PDRealArray *)sin;
- (PDRealArray *)sinpi;
- (PDRealArray *)sin:(applyRealBlock)block;
- (PDRealArray *)cos;
- (PDRealArray *)cospi;
- (PDRealArray *)cos:(applyRealBlock)block;
- (PDRealArray *)atan2:(PDRealArray *)x;
- (PDRealArray *)log;
- (PDRealArray *)log2;
- (PDRealArray *)log10;
- (PDRealArray *)exp2;
- (PDRealArray *)pow:(double)exp;
- (PDRealArray *)oneOverX;
- (PDRealArray *)diag;
- (PDComplexArray *)fft;
- (PDRealArray *)matmult:(PDRealArray *)matrix;
- (PDRealArray *)multiply:(double)factor;
- (PDRealArray *)divide:(double)denominator;
- (PDRealArray *)divideElementByElement:(PDRealArray *)denominators;
- (PDRealArray *)divideRows:(NSRange)rows byRow:(size_t)row ofRealArray:(PDRealArray *)denominators;
- (PDRealArray *)under:(double)numerator;
- (PDRealArray *)add:(double)addend;
- (PDRealArray *)subtract:(double)subtrahend;
- (PDRealArray *)subtractFrom:(double)minuend;

@end

@interface PDIntArray : PDArray

@property (nonatomic, assign) size_t *data;

+ (PDIntArray *)rowVectorFrom:(size_t)start to:(size_t)end;

@end
