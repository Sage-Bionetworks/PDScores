% Computes mel-frequency cepstral coefficients (MFCCs) for an audio signal.
% Inputs:
%  samples - mono audio signal
% Outputs:
%  cepstra - MFCC matrix
%
% (CC BY-SA 3.0) Max Little, 2014

function cepstra = mfcc(samples)

% Sample rate
sr = 44100;

% MFCC parameters
wintime = 0.020;        % Frame duration (sec)
numcep = 13;            % Number of cepstral coeffs to return
lifterexp = 0.6;        % Exponent for liftering
minfreq = 0;            % Lowest band edge of mel filters (Hz)
maxfreq = 8000;         % Highest band edge of mel filters (Hz)
nbands = 20;            % Number of warped spectral bands to use

% Pre-compute data (assuming fixed sample rate and window length)
winpts = round(wintime*sr);
nfft = 2^ceil(log2(winpts));
nfreqs = nfft/2+1;
wndw = hamming(nfft);
nframes = ceil(length(samples)/nfft);

% Conversion between mel/linear scales using O'Shaugnessy's formula
% Anonymous functions are not supported for code generation, so these have
% to be defined as their own functions (below)
%fmel2hz = @(m) 700*(10.^(m/2595)-1);
%fhz2mel = @(f) 2595*log10(1+f/700);

% Pre-compute mel-scale auditory perceptual spectrum
melwts = zeros(nbands,nfft);
fftfrqs = (0:nfft-1)/nfft*sr;           % Center freqs of each FFT bin

% Get center of mel bands uniformly spaced over given frequency range
minmel = fhz2mel(minfreq);
maxmel = fhz2mel(maxfreq);
binfrqs = fmel2hz(minmel+(0:(nbands+1))/(nbands+1)*(maxmel-minmel));

for i = 1:nbands
  fs = binfrqs(i+[0 1 2]);

  % Compute lower and upper slopes for all bins
  loslope = (fftfrqs-fs(1))/(fs(2)-fs(1));
  hislope = (fs(3)-fftfrqs)/(fs(3)-fs(2));

  % Form triangular mel-scale band kernel function
  melwts(i,:) = max(0,min(loslope,hislope));
end

% Only use positive frequency part of Fourier transform
melwts = melwts(:,1:nfreqs);

% Make the unitary DCT matrix for log mel-power spectra to cepstra
dctm = zeros(numcep, nbands);
for i = 1:numcep
    dctm(i,:) = cos((i-1)*(1:2:(2*nbands-1))/(2*nbands)*pi) * sqrt(2/nbands);
end
dctm(1,:) = dctm(1,:)/sqrt(2);

% Construct liftering weight matrix (N.B. lifterexp > 0)
liftwts = [1 (1:(numcep-1)).^lifterexp];
liftermat = diag(liftwts);

% Compute windowed power spectrum for each frame
if coder.target('MATLAB')
    sampbuff = buffer(samples,nfft);
else
    buflen = (length(samples) + nfft - 1) / nfft;
    sampbuff = zeros(nfft,buflen);
    coder.ceval('buffer_nooverlap',coder.ref(sampbuff),samples,length(samples),nfft);
end
%sampbuff = buffer(samples,nfft);
sampbuff = sampbuff.*repmat(wndw,1,nframes);
pspectrum = abs(fft(sampbuff)).^2;
pspectrum = pspectrum(1:(nfft/2+1),:);

% Sum over FFT bins to form mel-scale bins
aspectrum = melwts*pspectrum;

% Convert to cepstra via DCT
cepstra = dctm*log(aspectrum+1e-5);

% Apply liftering weights
cepstra = liftermat*cepstra;


function hz = fmel2hz(m)
hz = 700*(10.^(m/2595)-1);

function mel = fhz2mel(f)
mel = 2595*log10(1+f/700);
