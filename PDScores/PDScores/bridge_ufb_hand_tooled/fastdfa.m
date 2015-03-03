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
//function [alpha, intervals, flucts] = fastdfa(x, varargin)
//
//[xpts, ypts] = fastdfa_core(x, varargin{:});
//
//% Sort the intervals, and produce a log-log straight line fit
//datapts   = sortrows([xpts ypts],1);
//intervals = datapts(:,1);
//flucts    = datapts(:,2);
//
//coeffs    = polyfit(log10(xpts), log10(ypts), 1);
//alpha     = coeffs(1);

#import "PDArray.h"
#import "fastdfa_core_nomex.h"

void fastdfa(const PDRealArray *x, double *alpha, PDRealArray *intervals, PDRealArray *flucts)
{
    PDRealArray *outIntervals = [[PDRealArray alloc] initWithDimensions:<#(size_t)#>]
    fastdfa_core_nomex(<#double *outIntervals#>, <#double *outFlucts#>, <#double *inSignal#>, <#unsigned long rowsIn#>, <#unsigned long colsIn#>)
}