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

@interface PDArray : NSObject

@property (nonatomic, assign) void *data;
@property (nonatomic, assign) size_t size[2];
@property (nonatomic, assign) size_t allocatedSize;
@property (nonatomic, assign) size_t numDimensions;

- (size_t)typeSize;

- (BOOL)ensureCapacity;

@end

@interface PDComplexArray : PDArray

@property (nonatomic, strong) DSPDoubleComplex *data;

@end

@interface PDRealArray : PDArray

@property (nonatomic, strong) double *data;

@end
