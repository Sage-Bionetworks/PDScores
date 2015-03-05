//
//  PDMatlab.m
//  PDScores
//
//  Created by Erin Mounts on 3/3/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "PDMatlab.h"
#include <stdlib.h>
#include <math.h>

PDRealArray *zeros(size_t rows, size_t columns)
{
    PDRealArray *ra = [PDRealArray new];
    [ra setRows:rows columns:columns];
    
    return ra;
}

PDRealArray *NaN(size_t rows, size_t columns)
{
    PDRealArray *ra = zeros(rows, columns);
    double *p = ra.data;
    double *pEnd = ra.data + rows * columns;
    while (p < pEnd) {
        *p++ = NAN;
    }
    
    return ra;
}

PDRealArray *sortrows(PDRealArray *table, size_t column)
{
    PDRealArray *sorted = [table copy];
    
    if (sorted.cols == 1) {
        // just sort it in place and be done with it
        qsort_b(sorted.data, sorted.rows, sorted.typeSize, ^int(const void *ptr1, const void *ptr2) {
            double v1 = *(double *)ptr1;
            double v2 = *(double *)ptr2;
            return v1 < v2 ? -1 : v1 > v2 ? 1 : 0;
        });
    } else {
        // create a permutation array of row indices--we'll sort this first, then permute the rows
        PDIntArray *permute = [PDIntArray new];
        [permute setRows:table.rows columns:1];
        size_t *p = permute.data;
        size_t *pEnd = p + permute.rows;
        size_t index = 0;
        while (p < pEnd) {
            *p++ = index++;
        }
        
        const double *refColumn = sorted.data + sorted.rows * column;
        qsort_b(permute.data, permute.rows, permute.typeSize, ^int(const void *ptr1, const void *ptr2) {
            size_t index1 = *(size_t *)ptr1;
            size_t index2 = *(size_t *)ptr2;
            double v1 = refColumn[index1];
            double v2 = refColumn[index2];
            return v1 < v2 ? -1 : v1 > v2 ? 1 : 0;
        });
        
        // now permute the rows of each column of sorted according to the permutation array
        size_t rows = permute.rows;
        for (size_t column = 0; column < sorted.cols; ++column) {
            // make a copy of the permutation array to work with--we are going to mess with it
            size_t indices[rows];
            memcpy(indices, permute.data, rows * permute.typeSize);
            
            double *coldata = sorted.data + column * rows;
            for (size_t i = 0; i < rows; ++i) {
                // Check if this element needs to be permuted
                size_t iSrc = indices[i];
                if (iSrc == i) {
                    continue; // already where it needs to be--skip
                }
                
                size_t iDst = i;
                double temp = coldata[iDst];
                
                // Follow the permutation cycle
                do {
                    coldata[iDst] = coldata[iSrc];
                    indices[iSrc] = iDst;
                    
                    iDst = iSrc;
                    iSrc = indices[iSrc];
                } while (iSrc != i);
                
                coldata[iDst] = temp;
                indices[iDst] = iDst;
            }
        }
    }
    
    return sorted;
}

inline static double sqr(double x) {
    return x*x;
}

BOOL quadraticFit(int n, const double x[], const double y[], double coeffs[])
{
    double   sumx = 0.0;                        /* sum of x                      */
    double   sumx2 = 0.0;                       /* sum of x**2                   */
    double   sumx3 = 0.0;                       /* sum of x**3                   */
    double   sumx4 = 0.0;                       /* sum of x**4                   */
    double   sumxy = 0.0;                       /* sum of x * y                  */
    double   sumx2y = 0.0;                      /* sum of x**2 * y               */
    double   sumy = 0.0;                        /* sum of y                      */
    
    for (int i=0;i<n;i++)
    {
        sumx   += x[i];
        sumx2  += sqr(x[i]);
        sumx3  += sqr(x[i]) * x[i];
        sumx4  += sqr(sqr(x[i]));
        sumxy  += x[i] * y[i];
        sumx2y += sqr(x[i]) * y[i];
        sumy   += y[i];
    }
    
    double sxx = sumx2 - (sqr(sumx) / n);
    double sxy = sumxy - (sumx * sumy / n);
    double sxx2 = sumx3 - (sumx * sumx2 / n);
    double sx2y = sumx2y - (sumx2 * sumy / n);
    double sx2x2 = sumx4 - (sqr(sumx2) / n);
    
    double denom = sxx * sx2x2 - sqr(sxx2);
    if (denom == 0) {
        // singular matrix. can't solve the problem.
        coeffs[0] = 0;
        coeffs[1] = 0;
        coeffs[3] = 0;
        return NO;
    }
    
    coeffs[0] = (sx2y * sxx - sxy * sxx2) / denom;
    coeffs[1] = (sxy * sx2x2 - sx2y * sxx2) / denom;
    coeffs[3] = sumy / n - coeffs[1] * sumx / n - coeffs[0] * sumx2 / n;
    
    return YES;
}

BOOL linearFit(int n, const double x[], const double y[], double coeffs[])
{
    double   sumx = 0.0;                        /* sum of x                      */
    double   sumx2 = 0.0;                       /* sum of x**2                   */
    double   sumxy = 0.0;                       /* sum of x * y                  */
    double   sumy = 0.0;                        /* sum of y                      */
    double   sumy2 = 0.0;                       /* sum of y**2                   */
    
    for (int i=0;i<n;i++)
    {
        sumx  += x[i];
        sumx2 += sqr(x[i]);
        sumxy += x[i] * y[i];
        sumy  += y[i];
        sumy2 += sqr(y[i]);
    }
    
    double denom = (n * sumx2 - sqr(sumx));
    if (denom == 0) {
        // singular matrix. can't solve the problem.
        coeffs[0] = 0;
        coeffs[1] = 0;
        return NO;
    }
    
    coeffs[0] = (n * sumxy  -  sumx * sumy) / denom;
    coeffs[1] = (sumy * sumx2  -  sumx * sumxy) / denom;
    
    return YES;
}

PDRealArray *polyfit(PDRealArray *x, PDRealArray *y, int order)
{
    // both have to be single dimension vectors (row or col vector doesn't matter)
    if (!(x.rows == 1 || x.cols == 1) || !(y.rows == 1 || y.cols == 1)) {
        return nil;
    }
    
    // both have to be the same size
    size_t numPoints = x.rows * x.cols;
    if (y.rows * y.cols != numPoints) {
        return nil;
    }
    
    PDRealArray *coeffs = [PDRealArray new];
    [coeffs setRows:order + 1 columns:1];
    if (order == 1) {
        linearFit((int)numPoints, x.data, y.data, coeffs.data);
    } else if (order == 2) {
        quadraticFit((int)numPoints, x.data, y.data, coeffs.data);
    }
    
    return coeffs;
}

PDRealArray *repmat(PDRealArray *x, size_t rowsreps, size_t colsreps)
{
    PDRealArray *reppedmat = [PDRealArray new];
    [reppedmat setRows:x.rows * rowsreps columns:x.cols * colsreps];
    
    size_t originalRows = x.rows;

    // extend each column rowsreps times
    const double *pSrc = x.data;
    const double *pEnd = pSrc + originalRows;
    double *pDest = reppedmat.data;
    for (size_t col = 0; col < x.cols; ++col) {
        for (size_t rowrep = 0; rowrep < rowsreps; ++rowrep) {
            // copy the original column rowsreps times sequentially.
            // ptr blit is way faster than memcpy function call for smaller numbers of rows,
            // probably reasonably close for larger numbers.
            const double *p = pSrc;
            while (p < pEnd) {
                *pDest++ = *p++;
            }
        }
        // move pSrc, pEnd for the next column
        pSrc = pEnd;
        pEnd += originalRows;
    }
    
    // now replicate the whole shebang colsreps times
    double *pSrcAll = reppedmat.data;
    const size_t reppedRowsSize = reppedmat.rows * x.cols;
    const size_t wholeShebang = reppedRowsSize * reppedmat.typeSize;
    double *pDestAll = pSrcAll + reppedRowsSize;
    for (size_t colrep = 1; colrep < colsreps; ++colrep) {
        // this is more likely to involve copying a larger amount of data a smaller number of times,
        // which is where memcpy really shines.
        memcpy(pDestAll, pSrcAll, wholeShebang);
        pDestAll += reppedRowsSize;
    }
    
    return reppedmat;
}

double quantileForColumn(PDRealArray *column, double p)
{
    PDRealArray *sorted = sortrows(column, 1);
    size_t n = column.rows;
    double indexish = p * n - 0.5;
    if (indexish < 0.0) {
        return sorted.data[0];
    }
    if (indexish >= n - 1.0) {
        return sorted.data[n - 1];
    }
    if (indexish == floor(indexish)) {
        return sorted.data[(size_t)indexish];
    }
    
    // it's between two values, so do a linear interpolation
    size_t lower = floor(indexish);
    size_t higher = indexish + 1;
    double betwixt = indexish - lower;
    double xLower = sorted.data[lower];
    double xHigher = sorted.data[higher];
    double value = xLower + betwixt * (xHigher - xLower);
    return value;
}

PDRealArray *quantile(PDRealArray *x, double p)
{
    PDRealArray *quantile = [PDRealArray new];
    PDRealArray *array = x;
    
    if (x.rows == 1) {
        // just transpose it to a column vector by switching the dimensions
        [quantile setRows:1 columns:1];
        PDRealArray *xPrime = [x copy];
        xPrime.rows = x.cols;
        xPrime.cols = x.rows;
        array = xPrime;
    }
    [quantile setRows:1 columns:array.cols];
    
    double *pDest = quantile.data;
    for (size_t col = 0; col < array.cols; ++col) {
        PDRealArray *column = [array subarrayWithRows:NSMakeRange(0, array.rows) columns:NSMakeRange(col, 1)];
        *pDest++ = quantileForColumn(column, p);
    }
    
    return quantile;
}

