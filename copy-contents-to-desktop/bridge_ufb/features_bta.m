% Computes basic tapping test features.
% Inputs:
%  tap - tapping data vector: tap(:,1) - time points,
%        tap(:,2:3) - X,Y touch screen coordinates
%
% (CC BY-SA 3.0) Max Little, 2014
function ft = features_bta(tap)

% Output feature vector
ft = NaN(1,2);

% Ignore zero-length inputs
N = size(tap,1);
if (N == 0)
    return;
end

% Calculate relative time
t = tap(:,1)-tap(1,1);
T = tap(end,1);

% X,Y offset
tapx = tap(:,2);
tapy = tap(:,3);
tapx = tapx-mean(tapx);
tapy = tapy-mean(tapy);

% Find left/right finger 'depress' events
dx = diff(tapx);
i = (abs(dx) > 20);
ttx = t(i);

% Find depress event intervals
dtx = diff(ttx);

% Median and spread of tapping rate
tapdt = median(dtx);
tapiqr = iqr(dtx);

% Output tapping test feature vector
ft = [tapdt tapiqr];
