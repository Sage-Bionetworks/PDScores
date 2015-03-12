//% Computes mel-frequency cepstral coefficients (MFCCs) for an audio signal.
//% Inputs:
//%  samples - mono audio signal
//% Outputs:
//%  cepstra - MFCC matrix
//%
//% (CC BY-SA 3.0) Max Little, 2014

#import "bridge_ufb_hand_tooled.h"

//fmel2hz = @(m) 700*(10.^(m/2595)-1);
static inline double fmel2hz(double m)
{
    return 700.0 * (__exp10(m / 2595.0) - 1.0);
}

//fhz2mel = @(f) 2595*log10(1+f/700);
static inline double fhz2mel(double f)
{
    return 2595.0 * log10(1.0 + f / 700.0);
}

//function cepstra = mfcc(samples)
void mfcc(PDRealArray *samples, PDRealArray **cepstra)
{
    //% Sample rate
    //sr = 44100;
    double sr = 44100.0;

    //% MFCC parameters
    //wintime = 0.020;        % Frame duration (sec)
    //numcep = 13;            % Number of cepstral coeffs to return
    //lifterexp = 0.6;        % Exponent for liftering
    //minfreq = 0;            % Lowest band edge of mel filters (Hz)
    //maxfreq = 8000;         % Highest band edge of mel filters (Hz)
    //nbands = 20;            % Number of warped spectral bands to use
    double wintime = 0.020;     // Frame duration (sec)
    double numcep = 13.0;       // Number of cepstral coeffs to return
    double lifterexp = 0.6;     // Exponent for liftering
    double minfreq = 0.0;       // Lowest band edge of mel filters (Hz)
    double maxfreq = 8000.0;    // Highest band edge of mel filters (Hz)
    size_t nbands = 20;         // Number of warped spectral bands to use

    //% Pre-compute data (assuming fixed sample rate and window length)
    //winpts = round(wintime*sr);
    double winpts = round(wintime * sr);
    //nfft = 2^ceil(log2(winpts));
    size_t nfft = 1 << (int)ceil(log2(winpts));
    //nfreqs = nfft/2+1;
    size_t nfreqs = nfft / 2 + 1;
    //wndw = hamming(nfft);
    PDRealArray *wndw = hamming(nfft);
    //nframes = ceil(length(samples)/nfft);
    size_t lengthsamples = MAX(samples.rows, samples.cols);
    size_t nframes = ceil((double)lengthsamples / (double)nfft);

    //% Conversion between mel/linear scales using O'Shaugnessy's formula
    //fmel2hz = @(m) 700*(10.^(m/2595)-1);
    //fhz2mel = @(f) 2595*log10(1+f/700);
    // --see C static inline functions above
    
    //% Pre-compute mel-scale auditory perceptual spectrum
    //melwts = zeros(nbands,nfft);
    PDRealArray *melwts = zeros(nbands, nfft);
    //fftfrqs = (0:nfft-1)/nfft*sr;           % Center freqs of each FFT bin
    PDRealArray *fftfrqs = linspace(0, (double)(nfft - 1) / (double)nfft * sr, nfft);

    //% Get center of mel bands uniformly spaced over given frequency range
    //minmel = fhz2mel(minfreq);
    double minmel = fhz2mel(minfreq);
    //maxmel = fhz2mel(maxfreq);
    double maxmel = fhz2mel(maxfreq);
    //binfrqs = fmel2hz(minmel+(0:(nbands+1))/(nbands+1)*(maxmel-minmel));
    PDRealArray *binfrqs = [linspace(minmel, maxmel, nbands + 2) applyReal:^double(const double element) {
        return fmel2hz(element);
    }];

    //for i = 1:nbands
    //  fs = binfrqs(i+[0 1 2]);
    //
    //  % Compute lower and upper slopes for all bins
    //  loslope = (fftfrqs-fs(1))/(fs(2)-fs(1));
    //  hislope = (fs(3)-fftfrqs)/(fs(3)-fs(2));
    //
    //  % Form triangular mel-scale band kernel function
    //  melwts(i,:) = max(0,min(loslope,hislope));
    //end

    for (size_t i = 0; i < nbands; ++i) {
        NSRange range = NSMakeRange(i, 3);
        NSRange otherRange = NSMakeRange(0, 1);
        PDRealArray *fs = [binfrqs subarrayWithRows:otherRange columns:range];
        
        PDRealArray *loslope = [[fftfrqs subtract:fs.data[0]] divide:fs.data[1] - fs.data[0]];
        PDRealArray *hislope = [[fftfrqs subtractFrom:fs.data[2]] divide:fs.data[2] - fs.data[1]];
        
        PDRealArray *maxminslopes = [loslope applyReal:^double(const double element, const double otherArrayElement) {
            return MAX(0, MIN(element, otherArrayElement));
        } withRealArray:hislope];
        [melwts setSubarrayRows:NSMakeRange(i, 1) columns:NSMakeRange(0, maxminslopes.cols) fromArray:maxminslopes];
    }

    //% Only use positive frequency part of Fourier transform
    //melwts = melwts(:,1:nfreqs);
    melwts = [melwts subarrayWithRows:NSMakeRange(0, melwts.rows) columns:NSMakeRange(0, nfreqs)];

    //% Make the unitary DCT matrix for log mel-power spectra to cepstra
    //dctm = zeros(numcep, nbands);
    PDRealArray *dctm = zeros(numcep, nbands);
    //for i = 1:numcep
    //    dctm(i,:) = cos((i-1)*(1:2:(2*nbands-1))/(2*nbands)*pi) * sqrt(2/nbands);
    //end
    NSRange dctmcols = NSMakeRange(0, dctm.cols);
    PDRealArray *dctmbands = [PDRealArray rowVectorWithStart:1.0 step:2.0 cap:2.0 * (double)nbands - 1.0];
    for (size_t i = 0; i < numcep; ++i) {
        NSRange rows = NSMakeRange(i, 1);
        [dctm setSubarrayRows:rows columns:dctmcols fromArray:[[[dctmbands multiply:(double)i / (2.0 * (double)nbands) * M_PI] cos] multiply:sqrt(2.0/(double)nbands)]];
    }
    //dctm(1,:) = dctm(1,:)/sqrt(2);
    NSRange firstRow = NSMakeRange(0, 1);
    [dctm setSubarrayRows:firstRow columns:dctmcols fromArray:[[dctm subarrayWithRows:firstRow columns:dctmcols] divide:M_SQRT2]];

    //% Construct liftering weight matrix (N.B. lifterexp > 0)
    //liftwts = [1 (1:(numcep-1)).^lifterexp];
    PDRealArray *numceps_to_the_lifterexp = [linspace(1, numcep - 1, numcep - 1) pow:lifterexp];
    PDRealArray *liftwts = zeros(1, numcep);
    liftwts.data[0] = 1.0;
    [liftwts setSubarrayRows:firstRow columns:NSMakeRange(1, numcep - 1) fromArray:numceps_to_the_lifterexp];
    //liftermat = diag(liftwts);
    PDRealArray *liftermat = [liftwts diag];

    //% Compute windowed power spectrum for each frame
    //sampbuff = buffer(samples,nfft);
    PDRealArray *sampbuff = buffer(samples, nfft, 0);
    //sampbuff = sampbuff.*repmat(wndw,1,nframes);
    sampbuff = [sampbuff applyReal:^double(const double element, const double otherArrayElement) {
        return element * otherArrayElement;
    } withRealArray:repmat(wndw, 1, nframes)];
    //pspectrum = abs(fft(sampbuff)).^2;
    PDRealArray *pspectrum = [[[sampbuff fft] abs] square];
    //pspectrum = pspectrum(1:(nfft/2+1),:);
    NSRange pspectrum_rows_to_use = NSMakeRange(0, nfft / 2 + 1);
    NSRange pspectrum_cols = NSMakeRange(0, pspectrum.cols);
    pspectrum = [pspectrum subarrayWithRows:pspectrum_rows_to_use columns:pspectrum_cols];

    //% Sum over FFT bins to form mel-scale bins
    //aspectrum = melwts*pspectrum;
    PDRealArray *aspectrum = [melwts matmult:pspectrum];

    //% Convert to cepstra via DCT
    //cepstra = dctm*log(aspectrum+1e-5);
    *cepstra = [dctm matmult:[[aspectrum add:1e-5] log]];

    //% Apply liftering weights
    //cepstra = liftermat*cepstra;
    *cepstra = [liftermat matmult:*cepstra];
}
