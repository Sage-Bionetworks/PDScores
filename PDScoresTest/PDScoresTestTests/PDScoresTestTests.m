//
//  PDScoresTestTests.m
//  PDScoresTestTests
//
//  Created by Erin Mounts on 2/11/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>
//#import "bridge_ufb_types.h"
//#import "bridge_ufb_emxutil.h"
//#import "interp1.h"
//#import "swipep.h"
#import "bridge_ufb_hand_tooled.h"
#import "signalprocessing.h"

@interface PDArray(with_data)

@property (nonatomic, readonly) void *data;

@end

@implementation PDArray(with_data)

- (void *)data
{
    return (void *)[(id)self bytes];
}

@end

///*static*/ void pitchStrengthAllCandidatesX(const emxArray_real_T *f, const
//                                           emxArray_real_T *L, const emxArray_real_T *pc, emxArray_real_T *S);
///*static*/ void pitchStrengthOneCandidateX(const emxArray_real_T *f, const
//                                          emxArray_real_T *NL, double pc, emxArray_real_T *S);

extern PDRealArray *pitchStrengthAllCandidates(PDRealArray *f, PDRealArray *L, PDRealArray *pc);

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

- (PDArray *)createMatlabMatrixFromNumbersByColumn:(NSArray *)numbersByColumn imaginaryByColumn:(NSArray *)imByColumn
{
    NSUInteger columns = numbersByColumn.count;
    NSArray *firstCol = numbersByColumn[0];
    NSUInteger rows = firstCol.count;
//    int dimensions = (columns > 1) ? 2 : 1;
//    emxArray__common *matrix;
    PDArray *matrix;
//    size_t dataSize;
    BOOL complex = (imByColumn != nil);
    if (complex) {
//        emxArray_creal_T *cmatrix;
//        emxInit_creal_T(&cmatrix, dimensions);
//        matrix = cmatrix;
//        dataSize = sizeof(creal_T);
        matrix = [PDComplexArray new];
    } else {
//        emxArray_real_T *rmatrix;
//        emxInit_real_T(&rmatrix, dimensions);
//        matrix = rmatrix;
//        dataSize = sizeof(real_T);
        matrix = [PDRealArray new];
    }
//    double startSize = matrix->size[0];
//    if (dimensions > 1) {
//        startSize *= matrix->size[1];
//    }
//    matrix->size[0] = (int)rows;
//    matrix->size[1] = (int)columns;
//
//    emxEnsureCapacity(matrix, startSize, (int)dataSize);
    [matrix setRows:rows columns:columns];
    
    double *p = (double *)matrix.data;
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

- (PDArray *)createMatlabMatrixForVarName:(NSString *)varName
{
    NSArray *varNumbersByColumn = [self numbersByColumnForVarName:varName imaginaryPart:NO];
    
    return [self createMatlabMatrixFromNumbersByColumn:varNumbersByColumn imaginaryByColumn:nil];
}

- (PDComplexArray *)createMatlabComplexMatrixForVarName:(NSString *)varName
{
    NSArray *varNumbersByColumn = [self numbersByColumnForVarName:varName imaginaryPart:NO];
    NSArray *varImNumbersByColumn = [self numbersByColumnForVarName:varName imaginaryPart:YES];
    
    return (PDComplexArray *)[self createMatlabMatrixFromNumbersByColumn:varNumbersByColumn imaginaryByColumn:varImNumbersByColumn];
}

- (NSUInteger)dataSize:(PDArray *)matrix
{
//    NSUInteger size = 1;
//    for (NSUInteger dim = 0; dim < matrix->numDimensions; dim++) {
//        size *= matrix->size[dim];
//    }
//    
//    return size;
    return matrix.rows * matrix.cols;
}

- (double)compareDouble:(double)mine toDouble:(double)matlabs trackMax:(double *)max saveResultInArray:(NSMutableArray *)compare
{
    static const double epsilon = 1e-10;
    double denom = fabs(matlabs);
    double diff = fabs(mine - matlabs);
    double delta = (denom < epsilon) ? diff : diff / denom;
    if (max) {
        if (delta > *max) {
            *max = delta;
        }
    }

    NSString *comparison = [NSString stringWithFormat:@"%.8g\t%.8g\t%.5g", mine, matlabs, delta];
    [compare addObject:comparison];
    
//    if (delta > 1e-7) {
//        NSLog(@"%@", comparison);
//    }
    return delta;
}

- (BOOL)compare:(PDArray *)myOutput to:(PDArray *)matlabsOutput isComplex:(BOOL)complex tolerance:(double)epsilon varName:(NSString *)varName
{
    NSUInteger sizeMatlab = [self dataSize:matlabsOutput];
    NSUInteger sizeMine = [self dataSize:myOutput];
    
    NSMutableArray *compare = [NSMutableArray arrayWithCapacity:sizeMine];
    NSMutableArray *imCompare = [NSMutableArray arrayWithCapacity:sizeMine];
    
    XCTAssert(sizeMine == sizeMatlab, @"Failure: sizes of %@ differ: mine is %zux%zu, matlabs is %zux%zu", varName, myOutput.rows, myOutput.cols, matlabsOutput.rows, matlabsOutput.cols);
    if (sizeMatlab != sizeMine) {
        return NO;
    }
    
    size_t typeSize = complex ? 2 : 1;
    
    double *pMine = myOutput.data;
    double *pMineEnd = pMine + sizeMine * typeSize;
    double *pMatlab = matlabsOutput.data;
    
    double maxDelta = 0.0;
    double maxImDelta = 0.0;
    size_t outOfTolerance = 0;
    while (pMine < pMineEnd) {
        double mine = *pMine++;
        double matlabs = *pMatlab++;
        double delta = [self compareDouble:mine toDouble:matlabs trackMax:&maxDelta saveResultInArray:compare];
        
        if (complex) {
            mine = *pMine++;
            matlabs = *pMatlab++;
            delta = [self compareDouble:mine toDouble:matlabs trackMax:&maxImDelta saveResultInArray:imCompare];
        }
        
        if (fabs(delta) > epsilon) {
            ++outOfTolerance;
        }
    }
    
    NSLog(@"MaxDelta for %@: %g", varName, maxDelta);
    if (complex) {
        NSLog(@"MaxImDelta for %@: %g", varName, maxImDelta);
        XCTAssert(fabs(maxImDelta) <= epsilon, @"Failure: imaginary values of %@ differ!", varName);
    }
    XCTAssert(fabs(maxDelta) <= epsilon, @"Failure: %lu values of %@ (%d%%)differ by up to %g with tolerance given as %g", outOfTolerance, varName, (int)round(100.0 * outOfTolerance / (myOutput.rows * myOutput.cols)), maxDelta, epsilon);
    
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
        PDRealArray *f = (PDRealArray *)[self createMatlabMatrixForVarName:varfName];
        PDRealArray *absX = (PDRealArray *)[self createMatlabMatrixForVarName:varabsXName];
        PDRealArray *fERBs = (PDRealArray *)[self createMatlabMatrixForVarName:varfERBsName];
        PDRealArray *matlabinterp1 = (PDRealArray *)[self createMatlabMatrixForVarName:varinterp1Name];
        PDRealArray *matlabL = (PDRealArray *)[self createMatlabMatrixForVarName:varLName];
        
        // first test it against itself--does spline(x) == y for all x, y pairs in original data?
        PDRealArray *myOutputForInterp1 = interp1(f, absX, f, PDInterp1MethodSpline, 0.0);
        [self compare:myOutputForInterp1 to:absX isComplex:NO tolerance:4e-13 varName:[NSString stringWithFormat:@"%@ self-check", varinterp1Name]];

        // now compare it to Matlab's data
        myOutputForInterp1 = interp1(f, absX, fERBs, PDInterp1MethodSpline, 0.0);
        
        [self compare:myOutputForInterp1 to:matlabinterp1 isComplex:NO tolerance:2e-12 varName:varinterp1Name];
        
        PDRealArray *myL = [myOutputForInterp1 applyReal:^double(const double element) {
            return sqrt(MAX(0.0, element));
        }];
//        [self setSizes:myOutputForInterp1->size ofMatrix:myL elementSize:sizeof(double)];
        
//        NSUInteger size = [self dataSize:myL];
//        double *p1 = myL.data;
//        double *p1end = p1 + size;
//        double *p2 = myOutputForInterp1.data;
//        while (p1 < p1end) {
//            *p1++ = sqrt(MAX(0.0, *p2));
//            p2++; // don't do ++ in a macro in case it gets evaluated multiple times
//        }
        
        [self compare:myL to:matlabL isComplex:NO tolerance:1e-12 varName:varLName];
        
//        emxFree_real_T(&myL);
//        emxFree_real_T(&myOutputForInterp1);
//        emxFree_real_T(&matlabL);
//        emxFree_real_T(&fERBs);
//        emxFree_real_T(&absX);
//        emxFree_real_T(&f);
    }
}

//- (void)transposeMatrix:(emxArray_real_T *)matrix intoMatrix:(emxArray_real_T *)transpose
//{
//    int oldSize = transpose->size[0] * transpose->size[1];
//    transpose->size[0] = matrix->size[1];
//    transpose->size[1] = matrix->size[0];
//    emxEnsureCapacity((emxArray__common *)transpose, oldSize, (int)sizeof(double));
//    int matrixRows = matrix->size[0];
//    for (int matrixRow = 0; matrixRow < matrixRows; matrixRow++) {
//        int matrixCols = matrix->size[1];
//        for (int matrixCol = 0; matrixCol < matrixCols; matrixCol++) {
//            transpose->data[matrixCol + transpose->size[0] * matrixRow] = matrix->data[matrixRow + matrix->size[0]
//                                                                       * matrixCol];
//        }
//    }
//}

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
        PDRealArray *fERBs = (PDRealArray *)[self createMatlabMatrixForVarName:varfERBsName];
        PDRealArray *L = (PDRealArray *)[self createMatlabMatrixForVarName:varLName];
        PDRealArray *pc = (PDRealArray *)[self createMatlabMatrixForVarName:varpcName];
        PDRealArray *Si = (PDRealArray *)[self createMatlabMatrixForVarName:varSiName];
        PDRealArray *interp1linear = (PDRealArray *)[self createMatlabMatrixForVarName:varInterp1LinearName];
        PDRealArray *t = (PDRealArray *)[self createMatlabMatrixForVarName:vartName];
        PDRealArray *ti = (PDRealArray *)[self createMatlabMatrixForVarName:vartiName];
        
        PDRealArray *mySi = pitchStrengthAllCandidates(fERBs, L, pc);
        
        [self compare:mySi to:Si isComplex:NO tolerance:6e-10 varName:varSiName];
        
        PDRealArray *interp1out;
        
        // Take the matrix transpose of Si
        PDRealArray *SiPrime = [Si transpose];
        
        interp1out = interp1(ti, SiPrime, t, PDInterp1MethodLinear, NAN);
        
        PDRealArray *myInterp1Linear = [interp1out transpose];
        
        [self compare:myInterp1Linear to:interp1linear isComplex:NO tolerance:0.0 varName:varInterp1LinearName];
        
    }
}

- (void)testHanning
{
    PDRealArray *myHanning;
    
    for (int ex = 9; ex < 14; ++ex) {
        int pow2 = 1 << ex;
        NSString *varName = [NSString stringWithFormat:@"hanning(%d)", pow2];
        
        PDRealArray *matlabHanning = (PDRealArray *)[self createMatlabMatrixForVarName:varName];
        
        myHanning = hanning(pow2);
        
        BOOL ok = [self compare:myHanning to:matlabHanning isComplex:NO tolerance:6e-11 varName:varName];
        
        if (!ok) {
            break;
        }
    }
}

- (void)zeroPadSignal:(PDRealArray *)signal intoMatrix:(PDRealArray *)paddedSignal windowSize:(double)windowSize
{
    double hopSize = windowSize / 2.0;

    double halfWin = windowSize / 2.0;
//    int i15 = paddedSignal.rows;
//    paddedSignal->size[0] = ((int)halfWin + signal.rows) + (int)(hopSize + halfWin);
//    emxEnsureCapacity((emxArray__common *)paddedSignal, i15, (int)sizeof(double));
    [paddedSignal setRows:((int)halfWin + signal.rows) + (int)(hopSize + halfWin) columns:1];
    int b_ndbl = (int)halfWin;
    for (int i = 0; i < b_ndbl; i++) {
        paddedSignal.data[i] = 0.0;
    }
    
    b_ndbl = (int)signal.rows;
    for (int i = 0; i < b_ndbl; i++) {
        paddedSignal.data[i + (int)halfWin] = signal.data[i];
    }
    
    b_ndbl = (int)(hopSize + halfWin);
    for (int i = 0; i < b_ndbl; i++) {
        paddedSignal.data[(i + (int)halfWin) + signal.rows] = 0.0;
    }
}

//- (void)setSizes:(int *)sizes ofMatrix:(emxArray__common *)m elementSize:(size_t)elementSize
//{
//    int oldSize = 1;
//    int dims = MAX(2, m->numDimensions);
//    for (int i = 0; i < dims; ++i) {
//        oldSize *= m->size[i];
//        
//        // matrices always have at least 2 sizes, even if only one dimension
//        if (i >= m->numDimensions) {
//            m->size[i] = 1;
//        } else {
//            m->size[i] = sizes[i];
//        }
//    }
//    emxEnsureCapacity((emxArray__common *)m, oldSize, (int)elementSize);
//}

- (void)testSpectrogram
{
    PDRealArray *audiotrim = (PDRealArray *)[self createMatlabMatrixForVarName:@"audiotrim"];
    PDRealArray *xzp = [PDRealArray new];
    PDComplexArray *X;
    PDRealArray *f;
    PDRealArray *ti;
    PDRealArray *myHanning;
    
    for (int ex = 9; ex < 14; ++ex) {
        int windowSize = 1 << ex;
//        int halfWindow = windowSize / 2;
        int hopSize = windowSize / 2;
        
//        int myHanningSizes[2] = {windowSize, 1};
//        [self setSizes:myHanningSizes ofMatrix:myHanning elementSize:sizeof(double)];
        
        myHanning = hanning(windowSize);
        
        [self zeroPadSignal:audiotrim intoMatrix:xzp windowSize:windowSize];
        
        int overlap = windowSize - hopSize;
        
        X = specgram(xzp, windowSize, 44100.0, myHanning, overlap, &f, &ti);
//        spectrogram(X->data, f->data, ti->data, xzp->data, xzp->size[0], myHanning->data, overlap, windowSize, 44100.0);
        
        NSString *varX = [NSString stringWithFormat:@"specgramX(%d)", windowSize];
        NSString *varf = [NSString stringWithFormat:@"specgramf(%d)", windowSize];
        NSString *varti = [NSString stringWithFormat:@"specgramti(%d)", windowSize];
        PDRealArray *matlabX = (PDRealArray *)[self createMatlabComplexMatrixForVarName:varX];
        PDRealArray *matlabf = (PDRealArray *)[self createMatlabMatrixForVarName:varf];
        PDRealArray *matlabti = (PDRealArray *)[self createMatlabMatrixForVarName:varti];
        
        [self compare:X to:matlabX isComplex:YES tolerance:7e-8 varName:varX];
        [self compare:f to:matlabf isComplex:NO tolerance:0.0 varName:varf];
        [self compare:ti to:matlabti isComplex:NO tolerance:3e-14 varName:varti];
    }
}

- (void)testSwipep
{
    PDRealArray *audiotrim = (PDRealArray *)[self createMatlabMatrixForVarName:@"audiotrim"];
    PDRealArray *matlabSwipepOutF0 = (PDRealArray *)[self createMatlabMatrixForVarName:@"swipepOutF0"];
    PDRealArray *matlabSwipepOutF0t = (PDRealArray *)[self createMatlabMatrixForVarName:@"swipepOutF0t"];
    PDRealArray *matlabSwipepOutF0p = (PDRealArray *)[self createMatlabMatrixForVarName:@"swipepOutF0p"];
    
    PDRealArray *mySwipepOutF0;
    PDRealArray *mySwipepOutF0t;
    PDRealArray *mySwipepOutF0p;
    
    double plim[2] = { 50, 500 };
    swipep(audiotrim, 44100.0, plim, 0.02, 1.0 / 48.0, 0.1, 0.5, -INFINITY, &mySwipepOutF0, &mySwipepOutF0t, &mySwipepOutF0p);
    
    [self compare:mySwipepOutF0 to:matlabSwipepOutF0 isComplex:NO tolerance:2e-15 varName:@"swipeOutF0"];
    [self compare:mySwipepOutF0t to:matlabSwipepOutF0t isComplex:NO tolerance:3e-16 varName:@"swipeOutF0t"];
    [self compare:mySwipepOutF0p to:matlabSwipepOutF0p isComplex:NO tolerance:6e-12 varName:@"swipeOutF0p"];
}


@end
