//% Computes basic phonation test features.
//% Inputs:
//%  audio - mono floating-point, [-1,+1] normalized phonation
//%          signal
//%  srate - sample rate in Hz
//%
//% (CC BY-SA 3.0) Max Little, 2014
// Objective-C port (CC BY-SA 3.0) Sage Bionetworks, 2015.

#import "bridge_ufb_hand_tooled.h"

//function ft = features_bvav2(audio,srate)
void features_bvav2(PDRealArray *audio, double srate, PDRealArray **ft)
{
    //
    //% Output feature vector
    //ft = NaN(1,13);
    *ft = NaN(1, 13);
    
    //% Ignore zero-length inputs
    //N = length(audio);
    size_t N = MAX(audio.rows, audio.cols);
    //if (N == 0)
    //    return;
    //end
    if (N == 0) {
        return;
    }
    
    //% Trim away start/end of phonation
    //Tstart = 1.0;
    //Tend   = 10.0;
    double Tstart = 1.0;
    double Tend = 10.0;
    //istart = floor(Tstart*srate);
    size_t istart = floor(Tstart * srate);
    //iend   = ceil(Tend*srate);
    size_t iend = ceil(Tend * srate);
    //if (istart > N)
    //    istart = 1;
    //end
    if (istart > N) {
        istart = 1;
    }
    //if (iend > N)
    //    iend = N;
    //end
    if (iend > N) {
        iend = N;
    }
    //audiotrim = audio(istart:iend);
    PDRealArray *audiotrim = [audio elementsWithIndices:[PDIntArray rowVectorFrom:istart to:iend]];

    //N = length(audiotrim);
    N = MAX(audiotrim.rows, audiotrim.cols);
    //T = N/srate;
    double T = (double)N / srate;

    //% F0 features
    //f0dt = 0.02;
    //f0lo = 50;
    //f0hi = 500;
    double f0dt = 0.02;
    double f0lo = 50;
    double f0hi = 500;

    //% Get F0 estimates
    //[f0,f0t,f0p] = swipep(audiotrim,srate,[f0lo f0hi],f0dt);
    PDRealArray *f0, *f0t, *f0p;
    double plim[2] = { f0lo, f0hi };
    swipep(audiotrim, srate, plim, f0dt, 1.0 / 48.0, 0.1, 0.5, -INFINITY, &f0, &f0t, &f0p);

    //% Strip away all low-probability poor F0 estimates
    //i = (f0p >= 0.2);
    PDIntArray *i = [[f0p applyInt:^size_t(const double element) {
        return element >= 0.2;
    }] find];
    //f0 = f0(i);
    f0 = [f0 elementsWithIndices:i];
    //f0t = f0t(i);
    f0t = [f0t elementsWithIndices:i];
    //f0p = f0p(i);
    f0p = [f0p elementsWithIndices:i];

    //% Median F0
    //f0m = median(f0);
    PDRealArray *f0m = [f0 median];

    //% Two jitter estimates
    PDRealArray *abs_diff_f0_over_f0dt = [[[f0 diff] abs] divide:f0dt];
    //f0j = mean(abs(diff(f0))/f0dt);
    PDRealArray *f0j = [abs_diff_f0_over_f0dt mean];
    //f0jr = median(abs(diff(f0))/f0dt);
    PDRealArray *f0jr = [abs_diff_f0_over_f0dt median];

    //% Amplitude features
    //ampdt = 0.05;
    //ampwin = round(ampdt*srate);
    //abuf = buffer(audiotrim,ampwin);
    //l1amp = mean(abs(abuf));
    //l2amp = sqrt(mean(abuf.^2));
    double ampdt = 0.05;
    double ampwin = round(ampdt * srate);
    PDRealArray *abuf = buffer(audiotrim, ampwin, 0);
//    PDRealArray *l1amp = [[abuf abs] mean];
    PDRealArray *l2amp = [[[abuf square] mean] sqrt];

    //% Two shimmer estimates
    //ash = mean(abs(diff(l2amp))/ampdt)
    PDRealArray *l2amp_diff_abs_over_ampdt = [[[l2amp diff] abs] divide:ampdt];
    PDRealArray *ash = [l2amp_diff_abs_over_ampdt mean];
    //ashr = median(abs(diff(l2amp))/ampdt);
    PDRealArray *ashr = [l2amp_diff_abs_over_ampdt median];

    //% Compute MFCCs
    //cep = mfcc(audiotrim)';
    PDRealArray *cep;
    mfcc(audiotrim, &cep);
    cep = [cep transpose];

    //% MFCC summaries across time
    //cepm = median(cep);
    PDRealArray *cepm = [cep median];
    //cept = linspace(0,T,size(cep,1))';
    PDRealArray *cept = [linspace(0.0, T, cep.rows) transpose];
    //cepdt = cept(2)-cept(1);
    double cepdt = cept.data[1] - cept.data[0];

    //% MFCC jitter measure
    //cepj = mean(abs(diff(cep))/cepdt);
    PDRealArray *cepj = [[[[cep diff] abs] divide:cepdt] mean];

    //% Output phonation feature vector
    //ft = [f0m f0j f0jr ash ashr cepm(1:4) cepj(1:4)];
    (*ft).data[0] = f0m.data[0];
    (*ft).data[1] = f0j.data[0];
    (*ft).data[2] = f0jr.data[0];
    (*ft).data[3] = ash.data[0];
    (*ft).data[4] = ashr.data[0];
    memcpy((*ft).data + 5, cepm.data, sizeof(double) * 4);
    memcpy((*ft).data + 9, cepj.data, sizeof(double) * 4);
}
