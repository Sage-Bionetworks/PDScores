% Segment a phonation recording, based on an attack-sustain-release state
% machine model.
% Inputs:
%  x - mono floating-point, [-1,+1] normalized phonation recording
% Outputs:
%  y - phonation segments
%  i - sample index of start of each segment
%
% (CC BY-SA 3.0) Max Little, 2014
function [y,i] = vadsplitphon(x)

% Normalize audio input
x = x/max(abs(x));

srate = 44100;              % Sample rate
framedur = 0.04;            % Frame duration (sec)
framechange = 0.01;         % Frame change duration (sec)
maxsilencedur = 0.5;        % Post-voice activity silence dur (sec)
ampthresh = 0.05;           % Voice threshold

framesize  = ceil(framedur*srate);
frameinc   = ceil(framechange*srate);
maxsilenceframes = ceil(maxsilencedur/framechange);

% Mean absolute amplitude
fskip = ceil(framesize/frameinc);   % Ignore first few incomplete frames
X = buffer(x,framesize,framesize-frameinc);
amp = mean(abs(X(:,fskip:end)))';
N = length(amp);

% Loop pre-conditions
n = 1;
k = 0;
y = {};
i = [];
state = 0;  % 0-silent wait for voice, 1-voice active, 2-silent sustain

while (n <= N)
    
    switch (state)

        % Entry state: waiting for voice activity.
        case 0,
            if (amp(n) > ampthresh)
                activeframes = 1;
                tstart = n;
                state = 1;
            end

        % Tracking voice activity
        case 1,
            if (amp(n) > ampthresh)
                % Continue counting frames exceeding threshold
                activeframes = activeframes + 1;
            else
                % Amplitude below threshold: begin silence sustain
                activeframes = activeframes + 1;
                silenceframes = 1;
                state = 2;
            end

        % Silence sustain phase
        case 2,
            if (amp(n) <= ampthresh)
                % Amplitude below threshold: continue sustain
                silenceframes = silenceframes + 1;
                activeframes = activeframes + 1;
                
                if (silenceframes >= maxsilenceframes)
                    % Exceeded max silence sustain duration, capture return
                    % to wait for next activity.
                    tend = tstart + activeframes - 1;
                    i0 = (tstart-1)*frameinc+1;
                    i1 = tend*frameinc;
                    k = k + 1;
                    y{k} = x(i0:i1);
                    i(k,1) = i0;
                    i(k,2) = i1;
                    state = 0;
                end
            else
                % Ampltude goes above threshold again, cancel sustain
                activeframes = activeframes + 1;
                state = 1;
            end

    end
    n = n + 1;
    
end

% Deal with ending in non-entry state
if ((state == 1) || (state == 2))
    tend = tstart + activeframes - 1;
    i0 = (tstart-1)*frameinc+1;
    i1 = tend*frameinc;
    k = k + 1;
    y{k} = x(i0:i1);
    i(k,1) = i0;
    i(k,2) = i1;
end
