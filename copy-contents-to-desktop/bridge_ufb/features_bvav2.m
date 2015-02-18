% Computes basic phonation test features.
% Inputs:
%  audio - mono floating-point, [-1,+1] normalized phonation
%          signal
%  srate - sample rate in Hz
%
% (CC BY-SA 3.0) Max Little, 2014
function ft = features_bvav2(audio,srate)

% Output feature vector
ft = NaN(1,13);

% Ignore zero-length inputs
N = length(audio);
if (N == 0)
    return;
end

% Trim away start/end of phonation
Tstart = 1.0;
Tend   = 10.0;
istart = floor(Tstart*srate);
iend   = ceil(Tend*srate);
if (istart > N)
    istart = 1;
end
if (iend > N)
    iend = N;
end
audiotrim = audio(istart:iend);

N = length(audiotrim);
T = N/srate;

% F0 features
f0dt = 0.02;
f0lo = 50;
f0hi = 500;

% Get F0 estimates
[f0,f0t,f0p] = swipep(audiotrim,srate,[f0lo f0hi],f0dt,NaN,NaN,NaN,NaN);

% Strip away all low-probability poor F0 estimates
i = (f0p >= 0.2);
f0 = f0(i);
f0t = f0t(i);
f0p = f0p(i);

% Median F0
f0m = median(f0);

% Two jitter estimates
f0j = mean(abs(diff(f0))/f0dt);
f0jr = median(abs(diff(f0))/f0dt);

% Amplitude features
ampdt = 0.05;
ampwin = round(ampdt*srate);
%abuf = buffer(audiotrim,ampwin);
if coder.target('MATLAB')
    abuf = buffer(audiotrim,ampwin);
else
    buflen = (N + ampwin - 1) / ampwin;
    abuf = zeros(ampwin,buflen);
    coder.ceval('buffer_nooverlap',coder.ref(abuf),audiotrim,N,ampwin);
end

l1amp = mean(abs(abuf));
l2amp = sqrt(mean(abuf.^2));

% Two shimmer estimates
ash = mean(abs(diff(l2amp))/ampdt);
ashr = median(abs(diff(l2amp))/ampdt);

% Compute MFCCs
cep = mfcc(audiotrim)';

% MFCC summaries across time
cepm = median(cep);
cept = linspace(0,T,size(cep,1))';
cepdt = cept(2)-cept(1);

% MFCC jitter measure
cepj = mean(abs(diff(cep))/cepdt);

% Output phonation feature vector
ft = [f0m f0j f0jr ash ashr cepm(1:4) cepj(1:4)];
