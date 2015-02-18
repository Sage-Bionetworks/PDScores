% Computes basic posture test features.
% Inputs:
%  post - posture accelerometry vector: post(:,1) - time points,
%         post(:,2:4) - X,Y,Z acceleration data
%
% (CC BY-SA 3.0) Max Little, 2014
function ft = features_bpa(post)

% Output feature vector
ft = NaN(1,3);

% Ignore zero-length inputs
N = size(post,1);
if (N == 0)
    return;
end

% Calculate relative time
t = post(:,1)-post(1,1);

Tstart = 3.0;
Tend   = 19.0;

% Ignore posture tests which do not contain enough data
if (max(t) < Tend)
    return;
end

% Trim away start/end of test
istart = find(t >= Tstart,1,'first');
iend   = find(t <= Tend,1,'last');

if (iend < istart)
    return;
end

posttrim = post(istart:iend,2:4);
ttrim = t(istart:iend);

N = length(posttrim);
T = ttrim(end);
dT = T-ttrim(1);

% Orientation
mg = mean(posttrim);

% Orientation-corrected force signals
postforce = posttrim-repmat(mg,N,1);

% Scaled velocity signals
dt = diff(ttrim);
dt(end+1) = dt(end);
postvel = cumsum(postforce.*repmat(dt,1,3));

% Average scaled power X,Y,Z
postpower = mean(sum(0.5*70*postvel.^2)/dT)/1e4;

% Force vector magnitude signal
postmag = sqrt(sum(postforce.^2,2));

% Maximum force
postpeak = quantile(postmag,0.95)/10;

% Detrended fluctuation analysis scaling exponent
[alpha,iv,fl] = fastdfa(postmag);

% Output posture test feature vector
ft = [postpeak postpower alpha];
