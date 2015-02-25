//
//  PDScoresTestTests.m
//  PDScoresTestTests
//
//  Created by Erin Mounts on 2/11/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
#import "bridge_ufb_types.h"
#import "bridge_ufb_emxutil.h"
#import "interp1.h"
#import "swipep.h"
#import "signalprocessing.h"

/*static*/ void pitchStrengthAllCandidates(const emxArray_real_T *f, const
                                           emxArray_real_T *L, const emxArray_real_T *pc, emxArray_real_T *S);
/*static*/ void pitchStrengthOneCandidate(const emxArray_real_T *f, const
                                          emxArray_real_T *NL, double pc, emxArray_real_T *S);

@interface PDScoresTestTests : XCTestCase

@end

@implementation PDScoresTestTests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (NSURL *)inputURLForVar:(NSString *)varName
{
//    NSFileManager *fm = [NSFileManager defaultManager];
//    
//    static NSURL *dirURL = nil;
//    
//    static dispatch_once_t once;
//    dispatch_once(&once, ^{
//        dirURL = [fm URLForDirectory:NSDesktopDirectory inDomain:NSUserDomainMask appropriateForURL:nil create:NO error:nil];
//        dirURL = [dirURL URLByAppendingPathComponent:@"From Matlab" isDirectory:YES];
//        [fm createDirectoryAtURL:dirURL withIntermediateDirectories:YES attributes:nil error:nil];
//    });
//    
//    NSString *inName = [NSString stringWithFormat:@"%@.csv", varName];
//    NSURL *inputURL = [dirURL URLByAppendingPathComponent:inName];
    
    NSURL *inputURL = [[NSBundle bundleForClass:[self class]] URLForResource:varName withExtension:@".csv"];
    
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

- (emxArray__common *)createMatlabMatrixFromNumbersByColumn:(NSArray *)numbersByColumn imaginaryByColumn:(NSArray *)imByColumn
{
    NSUInteger columns = numbersByColumn.count;
    NSArray *firstCol = numbersByColumn[0];
    NSUInteger rows = firstCol.count;
    int dimensions = (columns > 1) ? 2 : 1;
    emxArray__common *matrix;
    size_t dataSize;
    BOOL complex = (imByColumn != nil);
    if (complex) {
        emxArray_creal_T *cmatrix;
        emxInit_creal_T(&cmatrix, dimensions);
        matrix = cmatrix;
        dataSize = sizeof(creal_T);
    } else {
        emxArray_real_T *rmatrix;
        emxInit_real_T(&rmatrix, dimensions);
        matrix = rmatrix;
        dataSize = sizeof(real_T);
    }
    double startSize = matrix->size[0];
    if (dimensions > 1) {
        startSize *= matrix->size[1];
    }
    matrix->size[0] = (int)rows;
    matrix->size[1] = (int)columns;

    emxEnsureCapacity(matrix, startSize, (int)dataSize);
    
    double *p = (double *)matrix->data;
    for (NSUInteger i = 0; i < numbersByColumn.count; ++i) {
        NSArray *columnReal = [numbersByColumn objectAtIndex:i];
        NSArray *columnIm;
        if (complex) {
            columnIm = [imByColumn objectAtIndex:i];
        }
        for (NSUInteger j = 0; j < columnReal.count; ++j) {
            NSNumber *real = [columnReal objectAtIndex:j];
            *p++ = [real doubleValue];
            if (complex) {
                NSNumber *im = [columnIm objectAtIndex:j];
                *p++ = [im doubleValue];
            }
        }
    }
    
    return matrix;
}

- (NSArray *)numbersByColumnForVarName:(NSString *)varName imaginaryPart:(BOOL)imaginaryPart
{
    if (imaginaryPart) {
        varName = [varName stringByAppendingString:@"im"];
    }
    NSURL *varURL = [self inputURLForVar:varName];
    NSString *varStr = [NSString stringWithContentsOfURL:varURL encoding:NSUTF8StringEncoding error:nil];
    NSArray *varRows = [varStr componentsSeparatedByString:@"\n"];
    NSPredicate *notEmptyPred = [NSPredicate predicateWithFormat:@"length > 0"];
    varRows = [varRows filteredArrayUsingPredicate:notEmptyPred];
    NSArray *varNumbersByColumn = [self numbersByColumnFromStringRowsArray:varRows];
    
    return varNumbersByColumn;
}

- (emxArray_real_T *)createMatlabMatrixForVarName:(NSString *)varName
{
    NSArray *varNumbersByColumn = [self numbersByColumnForVarName:varName imaginaryPart:NO];
    
    return [self createMatlabMatrixFromNumbersByColumn:varNumbersByColumn imaginaryByColumn:nil];
}

- (emxArray_creal_T *)createMatlabComplexMatrixForVarName:(NSString *)varName
{
    NSArray *varNumbersByColumn = [self numbersByColumnForVarName:varName imaginaryPart:NO];
    NSArray *varImNumbersByColumn = [self numbersByColumnForVarName:varName imaginaryPart:YES];
    
    return [self createMatlabMatrixFromNumbersByColumn:varNumbersByColumn imaginaryByColumn:varImNumbersByColumn];
}

- (NSUInteger)dataSize:(emxArray__common *)matrix
{
    NSUInteger size = 1;
    for (NSUInteger dim = 0; dim < matrix->numDimensions; dim++) {
        size *= matrix->size[dim];
    }
    
    return size;
}

- (double)compareDouble:(double)mine toDouble:(double)matlabs trackMax:(double *)max saveResultInArray:(NSMutableArray *)compare
{
    double denom = fabs(matlabs);
    double diff = fabs(mine - matlabs);
    double delta = (denom == 0) ? diff : diff / denom;
    if (max) {
        if (delta > *max) {
            *max = delta;
        }
    }

    NSString *comparison = [NSString stringWithFormat:@"%.8g\t%.8g\t%.5g", mine, matlabs, delta];
    [compare addObject:comparison];
    
    if (delta > 1e-10) {
        NSLog(@"%@", comparison);
    }
    return delta;
}

- (BOOL)compare:(emxArray__common *)myOutput to:(emxArray__common *)matlabsOutput isComplex:(BOOL)complex tolerance:(double)epsilon varName:(NSString *)varName
{
    NSUInteger sizeMatlab = [self dataSize:matlabsOutput];
    NSUInteger sizeMine = [self dataSize:myOutput];
    
    NSMutableArray *compare = [NSMutableArray arrayWithCapacity:sizeMine];
    NSMutableArray *imCompare = [NSMutableArray arrayWithCapacity:sizeMine];
    
    XCTAssert(sizeMine == sizeMatlab, @"Failure: sizes of %@ differ!", varName);
    if (sizeMatlab != sizeMine) {
        return NO;
    }
    
    size_t typeSize = complex ? 2 : 1;
    
    double *pMine = myOutput->data;
    double *pMineEnd = pMine + sizeMine * typeSize;
    double *pMatlab = matlabsOutput->data;
    
    double maxDelta = 0.0;
    double maxImDelta = 0.0;
    while (pMine < pMineEnd) {
        double mine = *pMine++;
        double matlabs = *pMatlab++;
        double delta = [self compareDouble:mine toDouble:matlabs trackMax:&maxDelta saveResultInArray:compare];
        
        if (complex) {
            mine = *pMine++;
            matlabs = *pMatlab++;
            delta = [self compareDouble:mine toDouble:matlabs trackMax:&maxImDelta saveResultInArray:imCompare];
        }
    }
    
    NSLog(@"MaxDelta for %@: %g", varName, maxDelta);
    if (complex) {
        NSLog(@"MaxImDelta for %@: %g", varName, maxImDelta);
        XCTAssert(fabs(maxImDelta) < epsilon, @"Failure: imaginary values of %@ differ!", varName);
    }
    XCTAssert(fabs(maxDelta) < epsilon, @"Failure: values of %@ differ!", varName);
    
    if (fabs(maxDelta) >= epsilon || fabs(maxImDelta) >= epsilon) {
        return NO;
    }
    
    return YES;
}

- (void)testInterp1
{
    for (int ex = 9; ex < 14; ++ex) {
        int pow2 = 1 << ex;
        
        NSString *varfName = [NSString stringWithFormat:@"f%d", pow2];
        NSString *varabsXName = [NSString stringWithFormat:@"abs(X)%d", pow2];
        NSString *varfERBsName = [NSString stringWithFormat:@"fERBs%d", pow2];
        NSString *varinterp1Name = [NSString stringWithFormat:@"interp1%d", pow2];
        NSString *varLName = [NSString stringWithFormat:@"L%d", pow2];
        emxArray_real_T *f = [self createMatlabMatrixForVarName:varfName];
        emxArray_real_T *absX = [self createMatlabMatrixForVarName:varabsXName];
        emxArray_real_T *fERBs = [self createMatlabMatrixForVarName:varfERBsName];
        emxArray_real_T *matlabinterp1 = [self createMatlabMatrixForVarName:varinterp1Name];
        emxArray_real_T *matlabL = [self createMatlabMatrixForVarName:varLName];
        
        emxArray_real_T *myOutputForInterp1;
        emxInit_real_T(&myOutputForInterp1, matlabL->numDimensions);
        
        interp1(f, absX, fERBs, myOutputForInterp1);
        
        [self compare:myOutputForInterp1 to:matlabinterp1 isComplex:NO tolerance:5e-13 varName:varinterp1Name];
        
        emxArray_real_T *myL;
        emxInit_real_T(&myL, myOutputForInterp1->numDimensions);
        [self setSizes:myOutputForInterp1->size ofMatrix:myL elementSize:sizeof(double)];
        
        NSUInteger size = [self dataSize:myL];
        double *p1 = myL->data;
        double *p1end = p1 + size;
        double *p2 = myOutputForInterp1->data;
        while (p1 < p1end) {
            *p1++ = sqrt(MAX(0.0, *p2));
            p2++; // don't do ++ in a macro in case it gets evaluated multiple times
        }
        
        [self compare:myL to:matlabL isComplex:NO tolerance:5e-13 varName:varLName];
        
        emxFree_real_T(&myL);
        emxFree_real_T(&myOutputForInterp1);
        emxFree_real_T(&matlabL);
        emxFree_real_T(&fERBs);
        emxFree_real_T(&absX);
        emxFree_real_T(&f);
    }
}

- (void)transposeMatrix:(emxArray_real_T *)matrix intoMatrix:(emxArray_real_T *)transpose
{
    int oldSize = transpose->size[0] * transpose->size[1];
    transpose->size[0] = matrix->size[1];
    transpose->size[1] = matrix->size[0];
    emxEnsureCapacity((emxArray__common *)transpose, oldSize, (int)sizeof(double));
    int matrixRows = matrix->size[0];
    for (int matrixRow = 0; matrixRow < matrixRows; matrixRow++) {
        int matrixCols = matrix->size[1];
        for (int matrixCol = 0; matrixCol < matrixCols; matrixCol++) {
            transpose->data[matrixCol + transpose->size[0] * matrixRow] = matrix->data[matrixRow + matrix->size[0]
                                                                       * matrixCol];
        }
    }
}

- (void)testPitchStrengthAllCandidates
{
    for (int ex = 9; ex < 14; ++ex) {
        int pow2 = 1 << ex;
        
        NSString *varfERBsName = [NSString stringWithFormat:@"fERBs%d", pow2];
        NSString *varLName = [NSString stringWithFormat:@"L%d", pow2];
        NSString *varpcName = [NSString stringWithFormat:@"pc(j)%d", pow2];
        NSString *varSiName = [NSString stringWithFormat:@"Si%d", pow2];
        NSString *varInterp1LinearName = [NSString stringWithFormat:@"interp1linear%d", pow2];
        NSString *vartName = @"tForPitchStrength";
        NSString *vartiName = [NSString stringWithFormat:@"specgramti(%d)", pow2];
        emxArray_real_T *fERBs = [self createMatlabMatrixForVarName:varfERBsName];
        emxArray_real_T *L = [self createMatlabMatrixForVarName:varLName];
        emxArray_real_T *pc = [self createMatlabMatrixForVarName:varpcName];
        emxArray_real_T *Si = [self createMatlabMatrixForVarName:varSiName];
        emxArray_real_T *interp1linear = [self createMatlabMatrixForVarName:varInterp1LinearName];
        emxArray_real_T *t = [self createMatlabMatrixForVarName:vartName];
        emxArray_real_T *ti = [self createMatlabMatrixForVarName:vartiName];
        
        emxArray_real_T *mySi;
        emxInit_real_T(&mySi, Si->numDimensions);
        
        pitchStrengthAllCandidates(fERBs, L, pc, mySi);
        
        [self compare:mySi to:Si isComplex:NO tolerance:5e-10 varName:varSiName];
        
        emxArray_real_T *b_interp1out;
        emxInit_real_T(&b_interp1out, interp1linear->numDimensions);
        
        // Take the matrix transpose of Si
        emxArray_real_T *SiPrime;
        emxInit_real_T(&SiPrime, Si->numDimensions);
        [self transposeMatrix:Si intoMatrix:SiPrime];
        
        b_interp1(ti, SiPrime, t, b_interp1out);
        
        emxArray_real_T *myInterp1Linear;
        emxInit_real_T(&myInterp1Linear, b_interp1out->numDimensions);
        [self transposeMatrix:b_interp1out intoMatrix:myInterp1Linear];
        
        [self compare:myInterp1Linear to:interp1linear isComplex:NO tolerance:1e-16 varName:varInterp1LinearName];
        
        emxFree_real_T(&myInterp1Linear);
        emxFree_real_T(&SiPrime);
        emxFree_real_T(&b_interp1out);
        emxFree_real_T(&mySi);
        emxFree_real_T(&ti);
        emxFree_real_T(&t);
        emxFree_real_T(&interp1linear);
        emxFree_real_T(&Si);
        emxFree_real_T(&pc);
        emxFree_real_T(&L);
        emxFree_real_T(&fERBs);
    }
}

- (void)testHanning
{
    emxArray_real_T *myHanning;
    emxInit_real_T(&myHanning, 1);
    
    for (int ex = 9; ex < 14; ++ex) {
        int pow2 = 1 << ex;
        NSString *varName = [NSString stringWithFormat:@"hanning(%d)", pow2];
        
        emxArray_real_T *matlabHanning = [self createMatlabMatrixForVarName:varName];
        
        double oldSize = myHanning->size[0];
        myHanning->size[0] = pow2;
        emxEnsureCapacity((emxArray__common *)myHanning, oldSize, (int)sizeof(double));
        
        hanning(myHanning->data, pow2);
        
        BOOL ok = [self compare:myHanning to:matlabHanning isComplex:NO tolerance:5e-10 varName:varName];
        
        emxFree_real_T(&matlabHanning);
        
        if (!ok) {
            break;
        }
    }

    
cleanup:
    emxFree_real_T(&myHanning);
}

- (void)zeroPadSignal:(emxArray_real_T *)signal intoMatrix:(emxArray_real_T *)paddedSignal windowSize:(double)windowSize
{
    double hopSize = windowSize / 2.0;

    double halfWin = windowSize / 2.0;
    int i15 = paddedSignal->size[0];
    paddedSignal->size[0] = ((int)halfWin + signal->size[0]) + (int)(hopSize + halfWin);
    emxEnsureCapacity((emxArray__common *)paddedSignal, i15, (int)sizeof(double));
    int b_ndbl = (int)halfWin;
    for (i15 = 0; i15 < b_ndbl; i15++) {
        paddedSignal->data[i15] = 0.0;
    }
    
    b_ndbl = signal->size[0];
    for (i15 = 0; i15 < b_ndbl; i15++) {
        paddedSignal->data[i15 + (int)halfWin] = signal->data[i15];
    }
    
    b_ndbl = (int)(hopSize + halfWin);
    for (i15 = 0; i15 < b_ndbl; i15++) {
        paddedSignal->data[(i15 + (int)halfWin) + signal->size[0]] = 0.0;
    }
}

- (void)setSizes:(int *)sizes ofMatrix:(emxArray__common *)m elementSize:(size_t)elementSize
{
    int oldSize = 1;
    int dims = MAX(2, m->numDimensions);
    for (int i = 0; i < dims; ++i) {
        oldSize *= m->size[i];
        
        // matrices always have at least 2 sizes, even if only one dimension
        if (i >= m->numDimensions) {
            m->size[i] = 1;
        } else {
            m->size[i] = sizes[i];
        }
    }
    emxEnsureCapacity((emxArray__common *)m, oldSize, (int)elementSize);
}

- (void)testSpectrogram
{
    emxArray_real_T *audiotrim = [self createMatlabMatrixForVarName:@"audiotrim"];
    emxArray_real_T *xzp;
    emxInit_real_T(&xzp, 1);
    emxArray_creal_T *X;
    emxInit_creal_T(&X, 2);
    emxArray_real_T *f;
    emxInit_real_T(&f, 1);
    emxArray_real_T *ti;
    emxInit_real_T(&ti, 1);
    emxArray_real_T *myHanning;
    emxInit_real_T(&myHanning, 1);
    
    for (int ex = 9; ex < 14; ++ex) {
        int windowSize = 1 << ex;
        int halfWindow = windowSize / 2;
        int hopSize = windowSize / 2;
        
        int myHanningSizes[2] = {windowSize, 1};
        [self setSizes:myHanningSizes ofMatrix:myHanning elementSize:sizeof(double)];
        
        hanning(myHanning->data, windowSize);
        
        [self zeroPadSignal:audiotrim intoMatrix:xzp windowSize:windowSize];
        
        int overlap = windowSize - hopSize;
        int framestep = windowSize - overlap;
        int bins = halfWindow + 1;
        int frames = (xzp->size[0] - overlap) / framestep;
        
        int Xsizes[2] = {bins, frames};
        int fsizes[2] = {bins, 1};
        int tisizes[2] = {frames, 1};
        [self setSizes:Xsizes ofMatrix:X elementSize:sizeof(creal_T)];
        [self setSizes:fsizes ofMatrix:f elementSize:sizeof(double)];
        [self setSizes:tisizes ofMatrix:ti elementSize:sizeof(double)];

        spectrogram(X->data, f->data, ti->data, xzp->data, xzp->size[0], myHanning->data, overlap, windowSize, 44100.0);
        
        NSString *varX = [NSString stringWithFormat:@"specgramX(%d)", windowSize];
        NSString *varf = [NSString stringWithFormat:@"specgramf(%d)", windowSize];
        NSString *varti = [NSString stringWithFormat:@"specgramti(%d)", windowSize];
        emxArray_real_T *matlabX = [self createMatlabComplexMatrixForVarName:varX];
        emxArray_real_T *matlabf = [self createMatlabMatrixForVarName:varf];
        emxArray_real_T *matlabti = [self createMatlabMatrixForVarName:varti];
        
        [self compare:X to:matlabX isComplex:YES tolerance:5e-5 varName:varX];
        [self compare:f to:matlabf isComplex:NO tolerance:1e-16 varName:varf];
        [self compare:ti to:matlabti isComplex:NO tolerance:3.2e-14 varName:varti];
        
        emxFree_real_T(&matlabti);
        emxFree_real_T(&matlabf);
        emxFree_real_T(&matlabX);
    }
    
    
cleanup:
    emxFree_real_T(&myHanning);
    emxFree_real_T(&ti);
    emxFree_real_T(&f);
    emxFree_creal_T(&X);
    emxFree_real_T(&xzp);
    emxFree_real_T(&audiotrim);
}

- (void)testB_swipep
{
    emxArray_real_T *audiotrim = [self createMatlabMatrixForVarName:@"audiotrim"];
    emxArray_real_T *matlabSwipepOutF0 = [self createMatlabMatrixForVarName:@"swipepOutF0"];
    emxArray_real_T *matlabSwipepOutF0t = [self createMatlabMatrixForVarName:@"swipepOutF0t"];
    emxArray_real_T *matlabSwipepOutF0p = [self createMatlabMatrixForVarName:@"swipepOutF0p"];
    
    emxArray_real_T *mySwipepOutF0;
    emxArray_real_T *mySwipepOutF0t;
    emxArray_real_T *mySwipepOutF0p;
    emxInit_real_T(&mySwipepOutF0, matlabSwipepOutF0->numDimensions);
    emxInit_real_T(&mySwipepOutF0t, matlabSwipepOutF0t->numDimensions);
    emxInit_real_T(&mySwipepOutF0p, matlabSwipepOutF0p->numDimensions);
    
    b_swipep(audiotrim, 44100.0, mySwipepOutF0, mySwipepOutF0t, mySwipepOutF0p);
    
    [self compare:mySwipepOutF0 to:matlabSwipepOutF0 isComplex:NO tolerance:5e-8 varName:@"swipeOutF0"];
    [self compare:mySwipepOutF0t to:matlabSwipepOutF0t isComplex:NO tolerance:1e-16 varName:@"swipeOutF0t"];
    [self compare:mySwipepOutF0p to:matlabSwipepOutF0p isComplex:NO tolerance:5e-8 varName:@"swipeOutF0p"];
    
cleanup:
    emxFree_real_T(&mySwipepOutF0p);
    emxFree_real_T(&mySwipepOutF0t);
    emxFree_real_T(&mySwipepOutF0t);
    emxFree_real_T(&matlabSwipepOutF0p);
    emxFree_real_T(&matlabSwipepOutF0t);
    emxFree_real_T(&matlabSwipepOutF0);
    emxFree_real_T(&audiotrim);
}


@end
