//% Computes basic gait test features.
//% Inputs:
//%  gait - gait accelerometry vector: gait(:,1) - time points,
//%         gait(:,2:4) - X,Y,Z acceleration data
//%
//% (CC BY-SA 3.0) Max Little, 2014

#import "bridge_ufb_hand_tooled.h"

//function ft = features_bga(gait)
void features_bga(PDRealArray *gait, PDRealArray **ft)
{
    //% Output feature vector
    //ft = NaN(1,7);
    *ft = NaN(1, 7);

    //% Ignore zero-length inputs
    //N = size(gait,1);
    size_t N = gait.rows;
    //if (N == 0)
    //    return;
    //end
    if (N == 0) {
        return;
    }

    //% Calculate relative time
    //t = gait(:,1)-gait(1,1);
    PDRealArray *times = [gait subarrayWithRows:NSMakeRange(0, gait.rows) columns:NSMakeRange(0, 1)];
    double startTime = times.data[0];
    PDRealArray *t = [times subtract:startTime];
    
    //Tstart = 3.0;
    //Tend   = 19.0;
    double Tstart = 3.0;
    double Tend = 19.0;
    
    //% Ignore gait tests which do not contain enough data
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
    
    //gaittrim = gait(istart:iend,2:4);
    NSRange startToEnd = NSMakeRange(istart - 1, iend - istart + 1);
    PDRealArray *gaittrim = [gait subarrayWithRows:startToEnd columns:NSMakeRange(1, 3)];
    //ttrim = t(istart:iend);
    PDRealArray *ttrim = [t subarrayWithRows:startToEnd columns:NSMakeRange(0, 1)];
    
    //N = length(gaittrim);
    //T = ttrim(end);
    //dT = T-ttrim(1);
    N = MAX(gaittrim.rows, gaittrim.cols);
    double T = ttrim.data[ttrim.rows - 1];
    double dT = T - ttrim.data[0];
    
    //% Orientation
    //mg = mean(gaittrim);
    PDRealArray *mg = [gaittrim mean];
    
    //% Orientation-corrected force signals
    //gaitforce = gaittrim-repmat(mg,N,1);
    PDRealArray *reppedmg = repmat(mg, N, 1);
    PDRealArray *gaitforce = [gaittrim applyReal:^double(const double element, const double otherArrayElement) {
        return element - otherArrayElement;
    } withRealArray:reppedmg];
    
    //% Scaled velocity signals
    //dt = diff(ttrim);
    PDRealArray *dt = [ttrim diff];
    //dt(end+1) = dt(end);
    [dt setRows:dt.rows + 1 columns:dt.cols];
    dt.data[dt.rows * dt.cols - 1] = dt.data[dt.rows * dt.cols - 2];
    //gaitvel = cumsum(gaitforce.*repmat(dt,1,3));
    PDRealArray *reppeddt = repmat(dt, 1, 3);
    PDRealArray *gaitvel = [[gaitforce applyReal:^double(const double element, const double otherArrayElement) {
        return element * otherArrayElement;
    } withRealArray:reppeddt] cumsum];
    
    //% Average scaled power X,Y,Z
    //gaitpower = mean(sum(0.5*70*gaitvel.^2)/dT)/1e4;
    PDRealArray *gaitpower = [[[[[gaitvel applyReal:^double(const double element) {
        return pow(element, 2.0) * 70.0 * 0.5;
    }] sum] divide:dT] mean] divide:1e4];
    
    //% Force vector magnitude signal
    //gaitmag = sqrt(sum(gaitforce.^2,2));
    PDRealArray *gfsquared = [gaitforce square];
    PDRealArray *gaitmagsquared = [gfsquared sum2];
    PDRealArray *gaitmag = [gaitmagsquared sqrt];
    
    //% Maximum force
    //gaitpeak = quantile(gaitmag,0.95)/10;
    PDRealArray *gaitpeak = [quantile(gaitmag, 0.95) divide:10.0];
    
    //% Zero crossing events
    //mg = mean(gaitmag);
    double gmm = [gaitmag mean].data[0];
    //izc = (gaitmag(1:end-1) < mg) & (gaitmag(2:end) >= mg);
    PDRealArray *gmBefore = [gaitmag subarrayWithRows:NSMakeRange(0, gaitmag.rows - 1) columns:NSMakeRange(0, 1)];
    PDRealArray *gmAfter = [gaitmag subarrayWithRows:NSMakeRange(1, gaitmag.rows - 1) columns:NSMakeRange(0, 1)];
    PDIntArray *izc = [gmBefore applyInt:^size_t(const double element, const double otherArrayElement) {
        return element < gmm && otherArrayElement >= gmm;
    } withRealArray:gmAfter];
    PDIntArray *indices = [izc find];
    //tzc = ttrim(izc);
    PDRealArray *tzc = [ttrim elementsWithIndices:indices];
    //dtzc = diff(tzc);
    PDRealArray *dtzc = [tzc diff];
    //zcr = median(dtzc);
    PDRealArray *zcr = [dtzc median];
    //zcv = iqr(dtzc);
    PDRealArray *zcv = [dtzc iqr];
    //
    //% Non-uniform frequency features
    //[F,P,prob] = lomb(ttrim,gaitmag-mg,4,1);
    PDRealArray *F;
    PDRealArray *P;
    PDRealArray *prob;
    PDRealArray *normGaitmag = [gaitmag subtract:gmm];
    lomb(ttrim, normGaitmag, 4.0, 1.0, &F, &P, &prob);
    
    //% Peak frequency
    //[Pmax,imax] = max(P);
    PDIntArray *imax;
    PDRealArray *Pmax = [P maxAndIndices:&imax];
    //F0 = F(imax);
    double F0 = F.data[imax.data[0] - 1];
    //probF0 = -log10(prob(imax));
    double probF0 = -log10(prob.data[imax.data[0] - 1]);

    //% SNR of peak frequency
    //i = (F <= F0);
    PDIntArray *i = [F applyInt:^size_t(const double element) {
        return element <= F0;
    }];
    //mP = median(log10(P(i)));
    PDRealArray *Psubi = [P elementsWithIndices:[i find]];
    PDRealArray *mP = [[Psubi log10] median];
    //snrF0 = log10(Pmax)-mP;
    PDRealArray *snrF0 = [Pmax applyReal:^double(const double element, const double otherArrayElement) {
        return log10(element) - otherArrayElement;
    } withRealArray:mP];

    //% Output gait test feature vector
    //ft = [gaitpeak gaitpower zcr zcv F0 snrF0 probF0];
    (*ft).data[0] = gaitpeak.data[0];
    (*ft).data[1] = gaitpower.data[0];
    (*ft).data[2] = zcr.data[0];
    (*ft).data[3] = zcv.data[0];
    (*ft).data[4] = F0;
    (*ft).data[5] = snrF0.data[0];
    (*ft).data[6] = probF0;
}
