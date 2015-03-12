//
//  PDScoresTestPerformanceTests.m
//  PDScoresTestPerformanceTests
//
//  Created by Erin Mounts on 2/20/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#import "signalprocessing.h"
@import Accelerate;

@interface PDScoresTestPerformanceTests : XCTestCase

@end

@implementation PDScoresTestPerformanceTests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

static const unsigned long kMinWindowSize = 256;
static const unsigned long kMaxWindowSize = 16384;

- (void)testMyHanningPerformance {
    [self measureBlock:^{
        double outBuf[kMaxWindowSize];
        for (int i = kMinWindowSize; i < kMaxWindowSize; ++i) {
            sp_hanning(outBuf, i);
        }
    }];
}

- (void)testApplesHanningPerformance {
    [self measureBlock:^{
        double outBuf[kMaxWindowSize];
        for (int i = kMinWindowSize; i < kMaxWindowSize; ++i) {
            vDSP_hann_windowD(outBuf, i, vDSP_HANN_DENORM);
        }
    }];
}

@end
