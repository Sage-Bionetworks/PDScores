//
//  PDArray.m
//  PDScores
//
//  Created by Erin Mounts on 3/2/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "PDArray.h"

#pragma mark - PDArray

@implementation PDArray

- (instancetype)init
{
    if (self = [super init]) {
        _size[0] = 0;
        _size[1] = 0;
        _numDimensions = 0;
        _allocatedSize = 0;
        _data = NULL;
    }
}

- (size_t)typeSize
{
    return 0; // subclasses *must* override
}

- (size_t)requiredCount
{
    size_t required = _size[0];
    if (_numDimensions > 1) {
        required *= _size[1];
    }
    
    return required;
}

- (BOOL)ensureCapacity
{
    size_t required = [self requiredCount];
    if (_allocatedSize < required) {
        size_t newSize = _allocatedSize * 2;
        while (newSize < required) {
            newSize *= 2;
        }
        _data = reallocf(_data, newSize * [self typeSize]);
    }
    
    return (_data != NULL);
}

@end

#pragma mark - PDComplexArray

@implementation PDComplexArray

- (size_t)typeSize
{
    return sizeof(DSPDoubleComplex);
}

- (DSPDoubleComplex *)data
{
    return (DSPDoubleComplex *)_data;
}

@end

#pragma mark - PDRealArray

@implementation PDRealArray

- (size_t)typeSize
{
    return sizeof(double);
}

- (double *)data
{
    return (double *)_data;
}

@end
