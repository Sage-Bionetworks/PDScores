//% Computes basic posture test features.
//% Inputs:
//%  post - posture accelerometry vector: post(:,1) - time points,
//%         post(:,2:4) - X,Y,Z acceleration data
//%
//% (CC BY-SA 3.0) Max Little, 2014

#import "bridge_ufb_hand_tooled.h"

//function ft = features_bpa(post)
void features_bpa(PDRealArray *post, PDRealArray **ft)
{
    //% Output feature vector
    //ft = NaN(1,3);
    *ft = NaN(1, 3);

    //% Ignore zero-length inputs
    //N = size(post,1);
    size_t N = post.rows;
    //if (N == 0)
    //    return;
    //end
    if (N == 0) {
        return;
    }

    //% Calculate relative time
    //t = post(:,1)-post(1,1);
    PDRealArray *times = [post subarrayWithRows:NSMakeRange(0, post.rows) columns:NSMakeRange(0, 1)];
    double startTime = times.data[0];
    PDRealArray *t = [times applyReal:^double(const double element) {
        return element - startTime;
    }];
    
    //Tstart = 3.0;
    //Tend   = 19.0;
    double Tstart = 3.0;
    double Tend = 19.0;
    
    //% Ignore posture tests which do not contain enough data
    //if (max(t) < Tend)
    //    return;
    //end
    if (t.max.data[0] < Tend) {
        return;
    }

    //% Trim away start/end of test
    //istart = find(t >= Tstart,1,'first');
    PDRealArray *fromStartTime = [t applyReal:^double(const double element) {
        return element >= Tstart ? 1.0 : 0.0;
    }];
    size_t istart = [fromStartTime findFirst:1].data[0];
    //iend   = find(t <= Tend,1,'last');
    PDRealArray *untilEndTime = [t applyReal:^double(const double element) {
        return element <= Tend ? 1.0 : 0.0;
    }];
    size_t iend = [untilEndTime findLast:1].data[0];

    //if (iend < istart)
    //    return;
    //end
    if (iend < istart) {
        return;
    }

    //posttrim = post(istart:iend,2:4);
    NSRange startToEnd = NSMakeRange(istart - 1, iend - istart + 1);
    PDRealArray *posttrim = [post subarrayWithRows:startToEnd columns:NSMakeRange(1, 3)];
    //ttrim = t(istart:iend);
    PDRealArray *ttrim = [t subarrayWithRows:startToEnd columns:NSMakeRange(0, 1)];

    //N = length(posttrim);
    //T = ttrim(end);
    //dT = T-ttrim(1);
    N = MAX(posttrim.rows, posttrim.cols);
    double T = ttrim.data[ttrim.rows - 1];
    double dT = T - ttrim.data[0];

    //% Orientation
    //mg = mean(posttrim);
    PDRealArray *mg = [posttrim mean];

    //% Orientation-corrected force signals
    //postforce = posttrim-repmat(mg,N,1);
    PDRealArray *reppedmg = repmat(mg, N, 1);
    PDRealArray *postforce = [posttrim applyReal:^double(const double element, const double otherArrayElement) {
        return element - otherArrayElement;
    } withRealArray:reppedmg];

    //% Scaled velocity signals
    //dt = diff(ttrim);
    PDRealArray *dt = [ttrim diff];
    //dt(end+1) = dt(end);
    [dt setRows:dt.rows + 1 columns:dt.cols];
    dt.data[dt.rows * dt.cols - 1] = dt.data[dt.rows * dt.cols - 2];
    //postvel = cumsum(postforce.*repmat(dt,1,3));
    PDRealArray *reppeddt = repmat(dt, 1, 3);
    PDRealArray *postvel = [[postforce applyReal:^double(const double element, const double otherArrayElement) {
        return element * otherArrayElement;
    } withRealArray:reppeddt] cumsum];

    //% Average scaled power X,Y,Z
    //postpower = mean(sum(0.5*70*postvel.^2)/dT)/1e4;
    PDRealArray *postpower = [[[[[postvel applyReal:^double(const double element) {
        return pow(element, 2.0) * 70.0 * 0.5;
    }] sum] applyReal:^double(const double element) {
        return element / dT;
    }] mean] applyReal:^double(const double element) {
        return element / 1e4;
    }];

    //% Force vector magnitude signal
    //postmag = sqrt(sum(postforce.^2,2));
    PDRealArray *pfsquared = [postforce applyReal:^double(const double element) {
        return pow(element, 2.0);
    }];
    PDRealArray *postmagsquared = [pfsquared sum2];
    PDRealArray *postmag = [postmagsquared sqrt];

    //% Maximum force
    //postpeak = quantile(postmag,0.95)/10;
    PDRealArray *postpeak = [quantile(postmag, 0.95) applyReal:^double(const double element) {
        return element / 10.0;
    }];

    //% Detrended fluctuation analysis scaling exponent
    //[alpha,iv,fl] = fastdfa(postmag);
    double alpha;
    PDRealArray *iv;
    PDRealArray *fl;
    fastdfa(postmag, &alpha, &iv, &fl);

    //% Output posture test feature vector
    //ft = [postpeak postpower alpha];
    (*ft).data[0] = postpeak.data[0];
    (*ft).data[1] = postpower.data[0];
    (*ft).data[2] = alpha;
}
