% Computes basic gait test features.
% Inputs:
%  gait - gait accelerometry vector: gait(:,1) - time points,
%         gait(:,2:4) - X,Y,Z acceleration data
%
% (CC BY-SA 3.0) Max Little, 2014
function ft = features_bga(gait)

% Output feature vector
ft = NaN(1,7);

% Ignore zero-length inputs
N = size(gait,1);
if (N == 0)
    return;
end

% Calculate relative time
t = gait(:,1)-gait(1,1);

Tstart = 3.0;
Tend   = 19.0;

% Ignore gait tests which do not contain enough data
if (max(t) < Tend)
    return;
end

% Trim away start/end of test
istart = find(t >= Tstart,1,'first');
iend   = find(t <= Tend,1,'last');

gaittrim = gait(istart:iend,2:4);
ttrim = t(istart:iend);

N = length(gaittrim);
T = ttrim(end);
dT = T-ttrim(1);

% Orientation
mg = mean(gaittrim);

% Orientation-corrected force signals
gaitforce = gaittrim-repmat(mg,N,1);

% Scaled velocity signals
dt = diff(ttrim);
dt(end+1) = dt(end);
gaitvel = cumsum(gaitforce.*repmat(dt,1,3));

% Average scaled power X,Y,Z
gaitpower = mean(sum(0.5*70*gaitvel.^2)/dT)/1e4;

% Force vector magnitude signal
gaitmag = sqrt(sum(gaitforce.^2,2));

% Maximum force
gaitpeak = quantile(gaitmag,0.95)/10;

% Zero crossing events
mg = mean(gaitmag);
izc = (gaitmag(1:end-1) < mg) & (gaitmag(2:end) >= mg);
tzc = ttrim(izc);
dtzc = diff(tzc);
zcr = median(dtzc);
zcv = iqr(dtzc);

% Non-uniform frequency features
[F,P,prob] = lomb(ttrim,gaitmag-mg,4,1);

% Peak frequency
[Pmax,imax] = max(P);
F0 = F(imax);
probF0 = -log10(prob(imax));

% SNR of peak frequency
i = (F <= F0);
mP = median(log10(P(i)));
snrF0 = log10(Pmax)-mP;

% Output gait test feature vector
ft = [gaitpeak gaitpower zcr zcv F0 snrF0 probF0];
