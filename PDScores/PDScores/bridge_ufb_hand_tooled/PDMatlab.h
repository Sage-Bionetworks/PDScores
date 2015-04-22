//
//  PDMatlab.h
//  PDScores
//
//  Created by Erin Mounts on 3/3/15.
//  (CC BY-SA 3.0) Sage Bionetworks, 2015.
//

#import "PDArray.h"

typedef NS_ENUM(int, PDInterp1Method) {
    PDInterp1MethodLinear = 0,
    PDInterp1MethodSpline
};

extern PDRealArray *zeros(size_t rows, size_t columns);
extern PDRealArray *ones(size_t rows, size_t columns);
extern PDRealArray *NaN(size_t rows, size_t columns);
extern PDRealArray *sortrows(PDRealArray *table, size_t column);
extern PDRealArray *polyfit(PDRealArray *x, PDRealArray *y, int order);
extern PDRealArray *polyval(PDRealArray *c, PDRealArray *x);
extern PDRealArray *repmat(PDRealArray *x, size_t rowsreps, size_t colsreps);
extern PDRealArray *quantile(PDRealArray *x, double p);
extern PDRealArray *buffer(PDRealArray *x, size_t n, size_t p);
extern PDRealArray *linspace(double start, double end, size_t n);
extern PDRealArray *hamming(size_t windowSize);
extern PDRealArray *hanning(size_t windowSize);
extern PDComplexArray *specgram(PDRealArray *x, size_t windowSize, double samplingRate, PDRealArray *window, size_t overlap, PDRealArray **freqs, PDRealArray **times);
extern PDRealArray *interp1(PDRealArray *x, PDRealArray *v, PDRealArray *xq, PDInterp1Method method, double extrapolation);
