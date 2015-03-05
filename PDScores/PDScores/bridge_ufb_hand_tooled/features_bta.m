//% Computes basic tapping test features.
//% Inputs:
//%  tap - tapping data vector: tap(:,1) - time points,
//%        tap(:,2:3) - X,Y touch screen coordinates
//%
//% (CC BY-SA 3.0) Max Little, 2014

#import "bridge_ufb_hand_tooled.h"

//function ft = features_bta(tap)
void features_bta(PDRealArray *tap, PDRealArray **ft)
{
    //% Output feature vector
    //ft = NaN(1,2);
    *ft = NaN(1, 2);
    
    //% Ignore zero-length inputs
    //N = size(tap,1);
    size_t N = tap.rows;
    //if (N == 0)
    //    return;
    //end
    if (N == 0) {
        return;
    }
    
    //% Calculate relative time
    //t = tap(:,1)-tap(1,1);
    //T = tap(end,1);
    NSRange tapRows = NSMakeRange(0, tap.rows);
    PDRealArray *times = [tap subarrayWithRows:tapRows columns:NSMakeRange(0, 1)];
    double startTime = times.data[0];
    PDRealArray *t = [times applyReal:^double(const double element) {
        return element - startTime;
    }];
    
    //% X,Y offset
    //tapx = tap(:,2);
    PDRealArray *tapx = [tap subarrayWithRows:tapRows columns:NSMakeRange(1, 1)];
    //tapy = tap(:,3);
    PDRealArray *tapy = [tap subarrayWithRows:tapRows columns:NSMakeRange(2, 1)];
    //tapx = tapx-mean(tapx);
    double tapxmean = [tapx mean].data[0];
    tapx = [tapx applyReal:^double(const double element) {
        return element - tapxmean;
    }];
    //tapy = tapy-mean(tapy);
    double tapymean = [tapy mean].data[0];
    tapy = [tapy applyReal:^double(const double element) {
        return element - tapymean;
    }];

    //% Find left/right finger 'depress' events
    //dx = diff(tapx);
    PDRealArray *dx = [tapx diff];
    //i = (abs(dx) > 20);
    PDIntArray *i = [dx applyInt:^size_t(const double element) {
        return abs(element) > 20;
    }];
    //ttx = t(i);
    PDRealArray *ttx = [t elementsWithIndices:[i find]];

    //% Find depress event intervals
    //dtx = diff(ttx);
    PDRealArray *dtx = [ttx diff];

    //% Median and spread of tapping rate
    //tapdt = median(dtx);
    PDRealArray *tapdt = [dtx median];
    //tapiqr = iqr(dtx);
    PDRealArray *tapiqr = [dtx iqr];

    //% Output tapping test feature vector
    //ft = [tapdt tapiqr];
    (*ft).data[0] = tapdt.data[0];
    (*ft).data[1] = tapiqr.data[0];
}
