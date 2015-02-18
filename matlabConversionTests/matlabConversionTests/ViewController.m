//
//  ViewController.m
//  matlabConversionTests
//
//  Created by Erin Mounts on 2/13/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "ViewController.h"
#import "bridge_ufb_types.h"
#import "bridge_ufb_emxutil.h"
#import "interp1.h"

@implementation ViewController

- (NSURL *)inputURLForVar:(NSString *)varName
{
    NSFileManager *fm = [NSFileManager defaultManager];
    
    static NSURL *dirURL = nil;
    
    static dispatch_once_t once;
    dispatch_once(&once, ^{
        dirURL = [fm URLForDirectory:NSDesktopDirectory inDomain:NSUserDomainMask appropriateForURL:nil create:NO error:nil];
        dirURL = [dirURL URLByAppendingPathComponent:@"From Matlab" isDirectory:YES];
        [fm createDirectoryAtURL:dirURL withIntermediateDirectories:YES attributes:nil error:nil];
    });
    
    NSString *inName = [NSString stringWithFormat:@"%@.csv", varName];
    NSURL *inputURL = [dirURL URLByAppendingPathComponent:inName];
    
    return inputURL;
}

- (NSMutableArray *)numbersByColumnFromStringRowsArray:(NSArray *)rows
{
    NSMutableArray *byColumns = [NSMutableArray array];
    for (NSString *rowString in rows) {
        NSArray *rowOfStrings = [rowString componentsSeparatedByString:@","];
        for (NSUInteger i = 0; i < rowOfStrings.count; ++i) {
            NSMutableArray *column = nil;
            if (i == byColumns.count) {
                column = [NSMutableArray arrayWithCapacity:rows.count];
                byColumns[i] = column;
            } else {
                column = byColumns[i];
            }
            double value = [rowOfStrings[i] doubleValue];
            NSNumber *num = [NSNumber numberWithDouble:value];
            [column addObject:num];
        }
    }
    
    return byColumns;
}

- (emxArray_real_T *)createMatlabMatrixFromNumbersByColumn:(NSArray *)numbersByColumn
{
    NSUInteger columns = numbersByColumn.count;
    NSArray *firstCol = numbersByColumn[0];
    NSUInteger rows = firstCol.count;
    int dimensions = (columns > 1) ? 2 : 1;
    emxArray_real_T *matrix;
    emxInit_real_T(&matrix, dimensions);
    double startSize = matrix->size[0];
    if (dimensions > 1) {
        startSize *= matrix->size[1];
    }
    matrix->size[0] = (int)rows;
    if (dimensions > 1) {
        matrix->size[1] = (int)columns;
    }
    emxEnsureCapacity((emxArray__common *)matrix, startSize, (int)sizeof(double));
    
    double *p = matrix->data;
    for (NSArray *column in numbersByColumn) {
        for (NSNumber *number in column) {
            *p++ = [number doubleValue];
        }
    }
    
    return matrix;
}

- (emxArray_real_T *)createMatlabMatrixForVarName:(NSString *)varName
{
    NSURL *varURL = [self inputURLForVar:varName];
    NSString *varStr = [NSString stringWithContentsOfURL:varURL encoding:NSUTF8StringEncoding error:nil];
    NSArray *varRows = [varStr componentsSeparatedByString:@"\n"];
    NSPredicate *notEmptyPred = [NSPredicate predicateWithFormat:@"length > 0"];
    varRows = [varRows filteredArrayUsingPredicate:notEmptyPred];
    NSArray *varNumbersByColumn = [self numbersByColumnFromStringRowsArray:varRows];
    
    return [self createMatlabMatrixFromNumbersByColumn:varNumbersByColumn];
}

- (NSUInteger)dataSize:(emxArray_real_T *)matrix
{
    NSUInteger size = 1;
    for (NSUInteger dim = 0; dim < matrix->numDimensions; dim++) {
        size *= matrix->size[dim];
    }
    
    return size;
}

- (void)testInterp1
{
    emxArray_real_T *f = [self createMatlabMatrixForVarName:@"f"];
    emxArray_real_T *absX = [self createMatlabMatrixForVarName:@"abs(X)"];
    emxArray_real_T *fERBs = [self createMatlabMatrixForVarName:@"fERBs"];
    emxArray_real_T *matlabOutputForInterp1 = [self createMatlabMatrixForVarName:@"interp1Out"];
    
    emxArray_real_T *myOutputForInterp1;
    emxInit_real_T(&myOutputForInterp1, matlabOutputForInterp1->numDimensions);
    
    interp1(f, absX, fERBs, myOutputForInterp1);
    
    NSUInteger sizeMatlab = [self dataSize:matlabOutputForInterp1];
    NSUInteger sizeMine = [self dataSize:myOutputForInterp1];
    
    NSMutableArray *deltas = [NSMutableArray arrayWithCapacity:sizeMine];
    NSMutableArray *compare = [NSMutableArray arrayWithCapacity:sizeMine];
    NSNumber *numDelta;
    NSString *comparison;
    
    if (sizeMatlab != sizeMine) {
        NSLog(@"Failure: sizes differ!");
        goto cleanup;
    }
    
    double *pMine = myOutputForInterp1->data;
    double *pMineEnd = pMine + sizeMine;
    double *pMatlab = matlabOutputForInterp1->data;
    
    while (pMine < pMineEnd) {
        double mine = *pMine++;
        double matlabs = *pMatlab++;
        double delta = mine - matlabs;
        numDelta = [NSNumber numberWithDouble:delta];
        [deltas addObject:numDelta];
        
        comparison = [NSString stringWithFormat:@"%g\t%g", mine, matlabs];
        [compare addObject:comparison];
    }
    
    for (NSNumber *numDelta in deltas) {
        double delta = [numDelta doubleValue];
        if (delta != 0.0) {
            NSLog(@"Failure: values differ!");
            goto cleanup;
        }
    }
    
    NSLog(@"Success!");
    
cleanup:
    emxFree_real_T(&matlabOutputForInterp1);
    emxFree_real_T(&fERBs);
    emxFree_real_T(&absX);
    emxFree_real_T(&f);
}

- (void)viewDidLoad {
    [super viewDidLoad];

    // Do any additional setup after loading the view.
    dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
        [self testInterp1];
    });
}

- (void)setRepresentedObject:(id)representedObject {
    [super setRepresentedObject:representedObject];

    // Update the view, if already loaded.
}

@end
