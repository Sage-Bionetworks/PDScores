#import "bridge_ufb_hand_tooled.h"

//function [f,P,prob] = lomb(t,h,ofac,hifac)
void lomb(PDRealArray *t, PDRealArray *h, double ofac, double hifac, PDRealArray **f, PDRealArray **P, PDRealArray **prob)
{
    //% LOMB(T,H,OFAC,HIFAC) computes the Lomb normalized periodogram (spectral
    //% power as a function of frequency) of a sequence of N data points H,
    //% sampled at times T, which are not necessarily evenly spaced. T and H must
    //% be vectors of equal size. The routine will calculate the spectral power
    //% for an increasing sequence of frequencies (in reciprocal units of the
    //% time array T) up to HIFAC times the average Nyquist frequency, with an
    //% oversampling factor of OFAC (typically >= 4).
    //%
    //% The returned values are arrays of frequencies considered (f), the
    //% associated spectral power (P) and estimated significance of the power
    //% values (prob).  Note: the significance returned is the false alarm
    //% probability of the null hypothesis, i.e. that the data is composed of
    //% independent gaussian random variables.  Low probability values indicate a
    //% high degree of significance in the associated periodic signal.
    //%
    //% Although this implementation is based on that described in Press,
    //% Teukolsky, et al. Numerical Recipes  In C, section 13.8, rather than using
    //% trigonometric rercurrences, this takes advantage of MATALB's array
    //% operators to calculate the exact spectral power as defined in equation
    //% 13.8.4 on page 577.  This may cause memory issues for large data sets and
    //% frequency ranges.
    //%
    //% Example
    //%    [f,P,prob] = lomb(t,h,4,1);
    //%    plot(f,P)
    //%    [Pmax,jmax] = max(P)
    //%    disp(['Most significant period is ',num2str(1/f(jmax)),...
    //%         ' with FAP of ',num2str(prob(jmax))])
    //%
    //% Written by Dmitry Savransky 21 May 2008

    //%sample length and time span
    //N = length(h);
    size_t N = MAX(h.rows, h.cols);
    //T = max(t) - min(t);
    double T = [t max].data[0] - [t min].data[0];

    //%mean and variance
    //mu = mean(h);
    double mu = [h mean].data[0];
    //s2 = var(h);
    double s2 = [h var].data[0];

    //%calculate sampling frequencies
    double start = 1.0 / (T * ofac);
    double step = start;
    double cap = hifac * (double)N / (2.0 * T);
    //f = (1/(T*ofac):1/(T*ofac):hifac*N/(2*T)).';
    *f = [[PDRealArray rowVectorWithStart:start step:step cap:cap] transpose];

    //%angular frequencies and constant offsets
    //w = 2*pi*f;
    PDRealArray *w = [*f multiply:2.0 * M_PI];
    
    //tau = atan2(sum(sin(2*w*t.'),2),sum(cos(2*w*t.'),2))./(2*w);
    // naïve approach was excessively slow and memory intensive:
//    PDRealArray *slowTau;
//    PDRealArray *w_t_prime;
//    @autoreleasepool {
//        PDRealArray *t_prime = [t transpose];
//        w_t_prime = [w matmult:t_prime];
//        PDRealArray *two_w_t_prime = [w_t_prime multiply:2.0];
//        PDRealArray *sumsinSlow = [[two_w_t_prime sin] sum2];
//        PDRealArray *sumcosSlow = [[two_w_t_prime cos] sum2];
//        PDRealArray *atan2_of_sums_slow = [sumsinSlow atan2:sumcosSlow];
//        slowTau = [atan2_of_sums_slow applyReal:^double(const double element, const double otherArrayElement) {
//            return element / (2.0 * otherArrayElement);
//        } withRealArray:w];
//    }
    // drastically simplified matrix calculations save hundreds of megabytes and several seconds
    PDRealArray *sumsin = [w applyReal:^double(const double wElement) {
        PDRealArray *sumSinThisRow = [[t applyReal:^double(const double tElement) {
            return sin(2.0 * wElement * tElement);
        }] sum];
        return sumSinThisRow.data[0];
    }];
    PDRealArray *sumcos = [w applyReal:^double(const double wElement) {
        PDRealArray *sumCosThisRow = [[t applyReal:^double(const double tElement) {
            return cos(2.0 * wElement * tElement);
        }] sum];
        return sumCosThisRow.data[0];
    }];
    PDRealArray *atan2_of_sums = [sumsin atan2:sumcos];
    PDRealArray *tau = [atan2_of_sums applyReal:^double(const double element, const double otherArrayElement) {
        return element / (2.0 * otherArrayElement);
    } withRealArray:w];
    
    //%spectral power
    
    //P = (sum(cterm*diag(h-mu),2).^2./sum(cterm.^2,2) + ...
    //     sum(sterm*diag(h-mu),2).^2./sum(sterm.^2,2))/(2*s2);

    // see above about bloated and slow
//    PDRealArray *ctermSlow, *stermSlow;
//    @autoreleasepool {
//        double lengthT = MAX(t.rows, t.cols);
//        PDRealArray *reppedw_tau = repmat([w applyReal:^double(const double element, const double otherArrayElement) {
//            return element * otherArrayElement;
//        } withRealArray:tau], 1, lengthT);
//        slowTau = nil;
//        PDRealArray *w_t_prime_minus_repped_w_tau = [w_t_prime applyReal:^double(const double element, const double otherArrayElement) {
//            return element - otherArrayElement;
//        } withRealArray:reppedw_tau];
//        w_t_prime = nil;
//        
//        //cterm = cos(w*t.' - repmat(w.*tau,1,length(t)));
//        ctermSlow = [w_t_prime_minus_repped_w_tau cos];
//        
//        //sterm = sin(w*t.' - repmat(w.*tau,1,length(t)));
//        stermSlow = [w_t_prime_minus_repped_w_tau sin];
//    }
//    PDRealArray *numeratorSlow;
//    @autoreleasepool {
//        PDRealArray *diag_h_mu = [[h subtract:mu] diag];
//        PDRealArray *cterm_diag_sum_sqSlow = [[[ctermSlow matmult:diag_h_mu] sum2] square];
//        PDRealArray *sterm_diag_sum_sqSlow = [[[stermSlow matmult:diag_h_mu] sum2] square];
//        PDRealArray *cterm_sq_sumSlow = [[ctermSlow square] sum2];
//        PDRealArray *sterm_sq_sumSlow = [[stermSlow square] sum2];
//        ctermSlow = nil;
//        stermSlow = nil;
//        numeratorSlow = [[cterm_diag_sum_sqSlow applyReal:^double(const double element, const double otherArrayElement) {
//            return element / otherArrayElement;
//        } withRealArray:cterm_sq_sumSlow] applyReal:^double(const double element, const double otherArrayElement) {
//            return element + otherArrayElement;
//        } withRealArray:[sterm_diag_sum_sqSlow applyReal:^double(const double element, const double otherArrayElement) {
//            return element / otherArrayElement;
//        } withRealArray:sterm_sq_sumSlow]];
//    }
    
    // see above about lean and fast
    PDRealArray *cterm_diag_sum_sq = [[w applyReal:^double(const double wElement, const double tauElement) {
        PDRealArray *cterm_diag_sum = [[t applyReal:^double(const double tElement, const double hElement) {
            return cos(wElement *(tElement - tauElement)) * (hElement - mu);
        } withRealArray:h] sum];
        return cterm_diag_sum.data[0];
    } withRealArray:tau] square];
    PDRealArray *sterm_diag_sum_sq = [[w applyReal:^double(const double wElement, const double tauElement) {
        PDRealArray *sterm_diag_sum = [[t applyReal:^double(const double tElement, const double hElement) {
            return sin(wElement *(tElement - tauElement)) * (hElement - mu);
        } withRealArray:h] sum];
        return sterm_diag_sum.data[0];
    } withRealArray:tau] square];
    PDRealArray *cterm_sq_sum = [w applyReal:^double(const double wElement, const double tauElement) {
        PDRealArray *csqsum = [[t applyReal:^double(const double tElement) {
            double costerm = cos(wElement *(tElement - tauElement));
            return costerm * costerm;
        }] sum];
        return csqsum.data[0];
    } withRealArray:tau];
    PDRealArray *sterm_sq_sum = [w applyReal:^double(const double wElement, const double tauElement) {
        PDRealArray *ssqsum = [[t applyReal:^double(const double tElement) {
            double sinterm = sin(wElement *(tElement - tauElement));
            return sinterm * sinterm;
        }] sum];
        return ssqsum.data[0];
    } withRealArray:tau];
    PDRealArray *numerator = [[cterm_diag_sum_sq divideElementByElement:cterm_sq_sum] applyReal:^double(const double element, const double otherArrayElement) {
        return element + otherArrayElement;
    } withRealArray:[sterm_diag_sum_sq divideElementByElement:sterm_sq_sum]];
    *P = [numerator divide:(2.0 * s2)];

    //%estimate of the number of independent frequencies
    //M=2*length(f)/ofac;
    double lengthf = MAX((*f).rows, (*f).cols);
    double M = 2.0 * lengthf / ofac;

    //%statistical significane of power
    //prob = M*exp(-P);
    *prob = [*P applyReal:^double(const double element) {
        return M * exp(-element);
    }];
    //inds = prob > 0.01;
    PDRealArray *inds = [*prob applyReal:^double(const double element) {
        return element > 0.01 ? 1.0 : 0.0;
    }];
    //prob(inds) = 1-(1-exp(-P(inds))).^M;
    PDRealArray *Pinds = [*P applyReal:^double(const double element, const double otherArrayElement) {
        double value = NAN;
        if (otherArrayElement) {
            value = 1.0 - pow((1.0 - exp(-element)), M);
        }
        return value;
    } withRealArray:inds];
    
    *prob = [*prob applyReal:^double(const double element, const double otherArrayElement) {
        double value = element;
        if (!isnan(otherArrayElement)) {
            value = otherArrayElement;
        }
        return value;
    } withRealArray:Pinds];
    
}
