//
//  PDArray.m
//  PDScores
//
//  Created by Erin Mounts on 3/2/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "PDArray.h"
#import "PDMatlab.h"
@import Accelerate;

#pragma mark - PDArray

@implementation PDArray

- (instancetype)init
{
    if (self = [super init]) {
        _rows = 0;
        _cols = 0;
        _allocatedSize = 0;
        _data = NULL;
    }

    return self;
}

- (size_t)typeSize
{
    return 0; // subclasses *must* override
}

- (BOOL)setRows:(size_t)rows columns:(size_t)columns
{
    size_t required = rows * columns;
    if (_allocatedSize < required) {
        size_t newSize = _allocatedSize * 2;
        if (!newSize) {
            newSize = 1;
        }
        while (newSize < required) {
            newSize <<= 1;
        }
        size_t typeSize = self.typeSize;
        size_t oldBytes = _allocatedSize * typeSize;
        size_t newBytes = newSize * typeSize;
        _data = reallocf(_data, newBytes);
        if  (_data) {
            int addedSize = (int)(newBytes - oldBytes);
            if (addedSize > 0) {
                // zero out the newly allocated part
                memset(_data + oldBytes, 0, addedSize);
            }
            _allocatedSize = newSize;
        }
    }
    
    if (_data) {
        _rows = rows;
        _cols = columns;
    }
    
    return (_data != NULL);
}

- (BOOL)concatenateColumnVectors:(NSArray *)arrays
{
    if (self.cols > 1) {
        return NO;
    }
    
    size_t totalSize = _rows;
    for (PDArray *array in arrays) {
        if (array.cols > 1) {
            return NO;
        }
        totalSize +=  array.rows * array.cols;
    }
    
    size_t offsetToNew = _rows * _cols * self.typeSize;
    BOOL canDo = [self setRows:totalSize columns:1];
    char *p = (char *)_data + offsetToNew;
    if (canDo) {
        for (PDArray *array in arrays) {
            size_t addedBytes = array.rows * array.cols * array.typeSize;
            memcpy(p, (char *)array.data, addedBytes);
            p += addedBytes;
        }
        _rows = totalSize;
    }
    
    return canDo;
}

- (BOOL)addColumns:(NSArray *)arrays
{
    // make sure they're all the right size
    size_t rows = _rows;
    for (PDArray *array in arrays) {
        if (!rows) {
            rows = array.rows;
        }
        if (array.rows != rows) {
            return NO;
        }
        if (array.cols != 1) {
            return NO;
        }
    }
    
    size_t totalCols = _cols + arrays.count;
    size_t totalRows = rows;
    
    // Matlab uses column-major order; we're going to follow their convention to avoid confusion in the conversion
    BOOL canDo = [self concatenateColumnVectors:arrays];
    _rows = totalRows;
    _cols = totalCols;
    
    return canDo;
}

- (instancetype)copy
{
    PDArray *copy = [[self class] new];
    [copy setRows:_rows columns:_cols];
    memcpy(copy.data, _data, _rows * _cols * self.typeSize);
    
    return copy;
}

- (instancetype)subarrayWithRows:(NSRange)rows columns:(NSRange)columns
{
    PDArray *subarray = [[self class] new];
    [subarray setRows:rows.length columns:columns.length];
    char *pDest = (char *)subarray.data;
    size_t endCol = columns.length + columns.location;
    for (size_t col = columns.location; col < endCol; ++col) {
        char *p = (char *)_data + self.typeSize * (_rows * col + rows.location);
        char *pEnd = p + self.typeSize * rows.length;
        while (p < pEnd) {
            *pDest++ = *p++;
        }
    }
    
    return subarray;
}

- (PDIntArray *)find
{
    PDIntArray *found = [PDIntArray new];
    
    // make it plenty big enough to start with
    size_t totals = _rows * _cols;
    size_t typeSize = self.typeSize;
    [found setRows:totals columns:1];
    char *p = _data;
    char *pEnd = p + totals * typeSize;
    size_t *pDest = found.data;
    size_t numFound = 0;
    while (p < pEnd) {
        if (![self isZero:p]) {
            *pDest++ = (p - (char *)_data) / typeSize + 1;
            ++numFound;
        }
        p += typeSize;
    }
    
    [found setRows:numFound columns:1];
    return found;
}

- (PDIntArray *)findFirst:(size_t)howMany
{
    PDIntArray *found = [PDIntArray new];
    [found setRows:howMany columns:1];
    
    size_t totals = _rows * _cols;
    size_t typeSize = self.typeSize;
    char *p = _data;
    char *pEnd = p + totals * typeSize;
    size_t *pDest = found.data;
    size_t numFound = 0;
    while (p < pEnd && numFound < howMany) {
        if (![self isZero:p]) {
            *pDest++ = (p - (char *)_data) / typeSize + 1;
            ++numFound;
        }
        p += typeSize;
    }
    
    [found setRows:numFound columns:1];
    return found;
}

- (PDIntArray *)findLast:(size_t)howMany
{
    PDIntArray *found = [PDIntArray new];
    [found setRows:howMany columns:1];
    
    size_t totals = _rows * _cols;
    size_t typeSize = self.typeSize;
    char *pEnd = _data ;
    char *p = pEnd + (totals - 1) * typeSize;
    size_t *pDest = found.data;
    size_t numFound = 0;
    while (p >= pEnd && numFound < howMany) {
        if (![self isZero:p]) {
            *pDest++ = (p - (char *)_data) / typeSize + 1;
            ++numFound;
        }
        p -= typeSize;
    }
    
    [found setRows:numFound columns:1];
    return found;
}

static inline void blit(char *dest, const char *src, size_t count)
{
    char *end = dest + count;
    while (dest < end) {
        *dest++ = *src++;
    }
}

- (instancetype)transpose
{
    PDArray *transpose = [[self class] new];
    [transpose setRows:_cols columns:_rows];
    if (_rows == 1 || _cols == 1) {
        // it's just a vector--copy it and we're done
        memcpy(transpose.data, _data, _rows * _cols * self.typeSize);
    } else {
        size_t typeSize = self.typeSize;
        size_t tColStride = transpose.rows;
        const char *pSrc = _data;
        const char *pEnd = pSrc + _rows * _cols * typeSize;
        char *rowStart = transpose.data;
        char *pDest = rowStart;
        while (pSrc < pEnd) {
            for (size_t srcRow = 0; srcRow < _rows; ++srcRow) {
                blit(pDest, pSrc, typeSize);
                pSrc += typeSize;
                pDest += tColStride;
            }
            rowStart += typeSize;
            pDest = rowStart;
        }
    }
    return transpose;
}

- (BOOL)isZero:(void *)valPtr
{
    return YES; // must override in all subclasses
}

- (void)dealloc
{
    free(_data);
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
    return (DSPDoubleComplex *)super.data;
}

- (BOOL)isZero:(void *)valPtr
{
    DSPDoubleComplex cval = *(DSPDoubleComplex *)valPtr;
    return (cval.real == 0.0 &&  cval.imag == 0.0);
}

@end

#pragma mark - PDRealArray

@implementation PDRealArray

+ (PDRealArray *)rowVectorWithStart:(double)start step:(double)step cap:(double)cap
{
    PDRealArray *rv = [PDRealArray new];
    if (step == 0.0 || (cap - start) / step < 0.0) {
        return rv;
    }
    size_t cols = round((cap - start) / step) + 1;
    [rv setRows:1 columns:cols];
    double *p = rv.data;
    double *pEnd = p + cols;
    double current = start;
    while (p < pEnd) {
        *p++ = current;
        current += step;
    }
    
    return rv;
}

- (size_t)typeSize
{
    return sizeof(double);
}

- (double *)data
{
    return (double *)super.data;
}

- (PDRealArray *)applyReal:(applyRealBlock)block
{
    PDRealArray *applied = [PDRealArray new];
    [applied setRows:self.rows columns:self.cols];
    double *p = self.data;
    double *pEnd = p + self.rows * self.cols;
    double *pDest = applied.data;
    while (p < pEnd) {
        *pDest++ = block(*p++);
    }
    return applied;
}

- (PDRealArray *)applyReal:(applyRealArrayBlock)block withRealArray:(PDRealArray *)array
{
    // sizes must match if there's another array
    if (array.rows != self.rows || array.cols != self.cols) {
        return nil;
    }
    
    PDRealArray *applied = [PDRealArray new];
    [applied setRows:self.rows columns:self.cols];
    const double *p = self.data;
    const double *pOther = array.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = applied.data;
    while (p < pEnd) {
        *pDest++ = block(*p++, *pOther++);
    }
    return applied;
}

- (PDIntArray *)applyInt:(applyRealIntBlock)block
{
    PDIntArray *applied = [PDIntArray new];
    [applied setRows:self.rows columns:self.cols];
    double *p = self.data;
    double *pEnd = p + self.rows * self.cols;
    size_t *pDest = applied.data;
    while (p < pEnd) {
        *pDest++ = block(*p++);
    }
    return applied;
}

- (PDIntArray *)applyInt:(applyRealArrayIntBlock)block withRealArray:(PDRealArray *)array
{
    // sizes must match if there's another array
    if (array.rows != self.rows || array.cols != self.cols) {
        return nil;
    }
    
    PDIntArray *applied = [PDIntArray new];
    [applied setRows:self.rows columns:self.cols];
    const double *p = self.data;
    const double *pOther = array.data;
    const double *pEnd = p + self.rows * self.cols;
    size_t *pDest = applied.data;
    while (p < pEnd) {
        *pDest++ = block(*p++, *pOther++);
    }
    return applied;
}

- (PDRealArray *)elementsWithIndices:(PDIntArray *)indexArray
{
    PDRealArray *elements = [PDRealArray new];
    size_t rows = indexArray.rows;
    size_t cols = 1;
    if (rows == 1) {
        rows = indexArray.cols;
    }
    [elements setRows:rows columns:cols];
    
    const size_t *pIdx = indexArray.data;
    const size_t *pEnd = pIdx + rows;
    const double *source = self.data;
    double *pDest = elements.data;
    while (pIdx < pEnd) {
        size_t index = *pIdx++;
        *pDest++ = source[index];
    }
    
    return elements;
}

double minForColumn(const double *colStart, size_t rows, size_t *index)
{
    const double *p = colStart;
    const double *pEnd = p + rows;
    double min = DBL_MAX;
    while (p < pEnd) {
        if (*p < min) {
            if (index) {
                *index = p - colStart;
            }
            min = *p;
        }
        ++p;
    }
    
    return min;
}

- (PDRealArray *)min
{
    PDRealArray *min = [PDRealArray new];
    [min setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = min.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [min setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = minForColumn(p, self.cols, NULL);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = minForColumn(p, self.rows, NULL);
            p += self.rows;
        }
    }
    
    return min;
}

- (PDRealArray *)minAndIndices:(PDIntArray *__autoreleasing *)indices
{
    PDRealArray *min = [PDRealArray new];
    [min setRows:1 columns:self.cols];
    *indices = [PDIntArray new];
    [*indices setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = min.data;
    size_t *pIDest = (*indices).data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [min setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = minForColumn(p, self.cols, pIDest++);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = minForColumn(p, self.rows, pIDest++);
            p += self.rows;
        }
    }
    
    return min;
}

double maxForColumn(const double *colStart, size_t rows, size_t *index)
{
    const double *p = colStart;
    const double *pEnd = p + rows;
    double max = -DBL_MAX;
    while (p < pEnd) {
        if (*p > max) {
            if (index) {
                *index = p - colStart;
            }
            max = *p;
        }
        ++p;
    }
    
    return max;
}

- (PDRealArray *)max
{
    PDRealArray *max = [PDRealArray new];
    [max setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = max.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [max setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = maxForColumn(p, self.cols, NULL);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = maxForColumn(p, self.rows, NULL);
            p += self.rows;
        }
    }
    
    return max;
}

- (PDRealArray *)maxAndIndices:(PDIntArray *__autoreleasing *)indices
{
    PDRealArray *max = [PDRealArray new];
    [max setRows:1 columns:self.cols];
    *indices = [PDIntArray new];
    [*indices setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = max.data;
    size_t *pIDest = (*indices).data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [max setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = maxForColumn(p, self.cols, pIDest++);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = maxForColumn(p, self.rows, pIDest++);
            p += self.rows;
        }
    }
    
    return max;
}

double meanForColumn(const double *colStart, size_t rows)
{
    const double *p = colStart;
    const double *pEnd = p + rows;
    double sum = 0.0;
    while (p < pEnd) {
        sum += *p++;
    }
    
    return sum / (double)rows;
}

- (PDRealArray *)mean
{
    PDRealArray *mean = [PDRealArray new];
    [mean setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = mean.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [mean setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = meanForColumn(p, self.cols);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = meanForColumn(p, self.rows);
            p += self.rows;
        }
    }
    
    return mean;
}

double medianForColumn(const double *colStart, size_t rows)
{
    double sortedCol[rows];
    memcpy(sortedCol, colStart, rows * sizeof(double));
    qsort_b(sortedCol, rows, sizeof(double), ^int(const void *ptr1, const void *ptr2) {
        double v1 = *(double *)ptr1;
        double v2 = *(double *)ptr2;
        return v1 < v2 ? -1 : v1 > v2 ? 1 : 0;
    });
    
    double median = NAN;
    size_t mid = rows >> 1;
    if (rows % 2) {
        // odd--take the middle one
        median = sortedCol[mid];
    } else {
        // even--split the difference
        double low = sortedCol[mid - 1];
        double high = sortedCol[mid];
        median = (low + high) / 2.0;
    }
    
    return median;
}

- (PDRealArray *)median
{
    PDRealArray *median = [PDRealArray new];
    [median setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = median.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [median setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = medianForColumn(p, self.cols);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = medianForColumn(p, self.rows);
            p += self.rows;
        }
    }
    
    return median;
}

- (PDRealArray *)iqr
{
    PDRealArray *quartile = quantile(self, 0.25);
    PDRealArray *threequartile = quantile(self, 0.75);
    PDRealArray *iqr = [threequartile applyReal:^double(const double element, const double otherArrayElement) {
        return element - otherArrayElement;
    } withRealArray:quartile];
    
    return iqr;
}

double varForColumn(const double *colStart, size_t rows)
{
    double mean = meanForColumn(colStart, rows);
    
    const double *p = colStart;
    const double *pEnd = p + rows;
    double sumsq = 0.0;
    while (p < pEnd) {
        double dev = *p++ - mean;
        sumsq += dev * dev;
    }
    
    // see docs for Matlab var() function
    double N = (rows > 1) ? (double)(rows - 1) : 1.0;
    
    return sumsq / N;
}


- (PDRealArray *)var
{
    PDRealArray *var = [PDRealArray new];
    [var setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = var.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [var setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = varForColumn(p, self.cols);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = varForColumn(p, self.rows);
            p += self.rows;
        }
    }
    
    return var;
}

void diffsForColumn(double *destStart, const double *colStart, size_t rows)
{
    const double *p = colStart;
    const double *pNext = p + 1;
    const double *pEnd = p + rows;
    double *pDest = destStart;
    while (pNext < pEnd) {
        *pDest++ = *pNext++ - *p++;
    }
}

- (PDRealArray *)diff
{
    PDRealArray *diff = [PDRealArray new];
    if (self.rows == 0 || self.cols == 0) {
        return diff;
    }
    [diff setRows:self.rows - 1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = diff.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [diff setRows:1 columns:self.cols - 1];
        while (p < pEnd) {
            diffsForColumn(pDest, p, self.rows);
            p += self.cols;
            pDest += diff.cols;
        }
    } else {
        while (p < pEnd) {
            diffsForColumn(pDest, p, self.rows);
            p += self.rows;
            pDest += diff.rows;
        }
    }
    
    return diff;
}

void cumsumsForColumn(double *destStart, const double *colStart, size_t rows)
{
    const double *pPrev = destStart;
    const double *p = colStart + 1;
    const double *pEnd = colStart + rows;
    double *pDest = destStart;
    *pDest++ = *colStart;
    while (p < pEnd) {
        *pDest++ = *p++ + *pPrev++;
    }
}

- (PDRealArray *)cumsum
{
    PDRealArray *cumsum = [PDRealArray new];
    [cumsum setRows:self.rows columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = cumsum.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        while (p < pEnd) {
            cumsumsForColumn(pDest, p, self.cols);
            p += self.cols;
            pDest += cumsum.cols;
        }
    } else {
        while (p < pEnd) {
            cumsumsForColumn(pDest, p, self.rows);
            p += self.rows;
            pDest += cumsum.rows;
        }
    }
    
    return cumsum;
}

double sumForColumn(const double *colStart, size_t rows)
{
    const double *p = colStart;
    const double *pEnd = p + rows;
    double sum = 0.0;
    while (p < pEnd) {
        sum += *p++;
    }
    
    return sum;
}

- (PDRealArray *)sum
{
    PDRealArray *sum = [PDRealArray new];
    [sum setRows:1 columns:self.cols];
    const double *p = self.data;
    const double *pEnd = p + self.rows * self.cols;
    double *pDest = sum.data;
    if (self.rows == 1) {
        // pretend it's a column vector by swapping indices
        [sum setRows:1 columns:1];
        while (p < pEnd) {
            *pDest++ = sumForColumn(p, self.cols);
            p += self.cols;
        }
    } else {
        while (p < pEnd) {
            *pDest++ = sumForColumn(p, self.rows);
            p += self.rows;
        }
    }
    
    return sum;
}

void addRowsFromColumn(double *destStart, const double *colStart, size_t rows)
{
    const double *p = colStart;
    const double *pEnd = p + rows;
    double *pDest = destStart;
    while (p < pEnd) {
        *pDest++ += *p++;
    }
}

- (PDRealArray *)sum2
{
    PDRealArray *sum2 = zeros(self.rows, 1);
    for (size_t column = 0; column < self.cols; ++column) {
        addRowsFromColumn(sum2.data, self.data + column * self.rows, self.rows);
    }
    
    return sum2;
}

- (PDRealArray *)sqrt
{
    // assumes all values in self are non-negative; seems to be true for Max Little's code
    return [self applyReal:^double(const double element) {
        return sqrt(element);
    }];
}

- (PDRealArray *)diag
{
    PDRealArray *diag = [PDRealArray new];
    size_t length = self.rows * self.cols;
    size_t srcStride;
    size_t destStride;

    if (self.rows == 1 || self.cols == 1) {
        // create a square diagonal matrix from this vector
        [diag setRows:length columns:length];
        srcStride = 1;
        destStride = length + 1;
    } else {
        // create a column vector from the main diagonal of this matrix
        [diag setRows:length columns:1];
        srcStride = length + 1;
        destStride = 1;
    }
    
    const double *pSrc = self.data;
    const double *pEnd = pSrc + length;
    double *pDest = diag.data;
    while (pSrc < pEnd) {
        *pDest = *pSrc;
        pSrc += srcStride;
        pDest += destStride;
    }
    
    return diag;
}

- (PDRealArray *)matmult:(PDRealArray *)matrix
{
    // inner dimensions must match
    if (self.cols != matrix.rows) {
        return nil;
    }
    PDRealArray *result = [PDRealArray new];
    [result setRows:self.rows columns:matrix.cols];
    
    // reverse the order of the inputs because of Matlab's column-major vs everyone else's row-major thingy
    // also swap rows for cols and vice-versa for the same stupid reason
    vDSP_mmulD(matrix.data, 1, self.data, 1, result.data, 1, matrix.cols, self.rows, self.cols);
    return result;
}

- (BOOL)isZero:(void *)valPtr
{
    return *(double *)valPtr == 0.0;
}

@end

#pragma mark - PDIntArray

@implementation PDIntArray

- (size_t)typeSize
{
    return sizeof(size_t);
}

- (BOOL)isZero:(void *)valPtr
{
    return *(size_t *)valPtr == 0;
}

@end
