% Finds the best segment of a phonation recording.
% Inputs:
%  rawaudio - mono floating-point, [-1,+1] normalized phonation recording
% Outputs:
%  audiout - best phonation signal
%  indout  - segment index of phonation signal
%
% (CC BY-SA 3.0) Max Little, 2014
function [audioout, indout] = extractaudiophon(rawaudio)

srate = 44100;
mindur = floor(2.0*srate);

[aud,ind] = vadsplitphon(rawaudio);

audioout = [];
indout = [];
if (length(aud) >= 1)
    L = [];
    for i = 1:length(aud)
        L(i) = length(aud{i});
    end
    [N,i] = max(L);

    if (N >= mindur)
        audioout = aud{i};
        indout = ind(i,:);
    end
end
