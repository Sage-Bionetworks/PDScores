//% Performs fast detrended fluctuation analysis on a nonstationary input signal to
//% obtain an estimate for the scaling exponent.
//%
//% Useage:
//% [alpha, intervals, flucts] = fastdfa(x)
//% [alpha, intervals, flucts] = fastdfa(x, intervals)
//% Inputs
//%    x          - input signal: must be a row vector
//% Optional inputs
//%    intervals  - List of sample interval widths at each scale
//%                 (If not specified, then a binary subdivision is constructed)
//%
//% Outputs:
//%    alpha      - Estimated scaling exponent
//%    intervals  - List of sample interval widths at each scale
//%    flucts     - List of fluctuation amplitudes at each scale
//%
//% (CC BY-SA 3.0) Max Little, 2006-2014.
//% If you use this code for academic publication, please cite:
//% M. Little, P. McSharry, I. Moroz, S. Roberts (2006),
//% Nonlinear, biophysically-informed speech pathology detection
//% in Proceedings of ICASSP 2006, IEEE Publishers: Toulouse, France.
//
// Objective-C port (CC BY-SA 3.0) Sage Bionetworks, 2015.

#import "bridge_ufb_hand_tooled.h"
#import "fastdfa_core_nomex.h"

//function [alpha, intervals, flucts] = fastdfa(x, varargin)
void fastdfa(const PDRealArray *x, double *alpha, PDRealArray **intervals, PDRealArray **flucts)
{
    //[xpts, ypts] = fastdfa_core(x, varargin{:});
    double elements = x.rows * x.cols;
    size_t nscales = floor(log2(elements));
    if (1 << (nscales - 1) > elements / 2.5) {
        nscales = nscales - 1;
    }
    
    PDRealArray *xpts = zeros(nscales,1);
    PDRealArray *ypts = zeros(nscales,1);
    fastdfa_core_nomex(xpts.data, ypts.data, x.data, x.rows, x.cols);

    //% Sort the intervals, and produce a log-log straight line fit
    //datapts   = sortrows([xpts ypts],1);
    PDRealArray *both = [PDRealArray new];
    [both addColumns:@[xpts, ypts]];
    PDRealArray *datapts = sortrows(both, 0);
    
    //intervals = datapts(:,1);
    *intervals = [datapts subarrayWithRows:NSMakeRange(0, datapts.rows) columns:NSMakeRange(0, 1)];
    //flucts    = datapts(:,2);
    *flucts = [datapts subarrayWithRows:NSMakeRange(0, datapts.rows) columns:NSMakeRange(1, 1)];

    //coeffs    = polyfit(log10(xpts), log10(ypts), 1);
    PDRealArray *log10xpts = [xpts applyReal:^double(const double element) {
        return log10(element);
    }];
    PDRealArray *log10ypts = [ypts applyReal:^double(const double element) {
        return log10(element);
    }];
    PDRealArray *coeffs = polyfit(log10xpts, log10ypts, 1);
    
    //alpha     = coeffs(1);
    *alpha = coeffs.data[0];
}