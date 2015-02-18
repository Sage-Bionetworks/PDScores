% Receives raw features, normalization parameters and feature weights, and
% computes user feedback value.
% Inputs:
%  ftvec - feature vector (1xN)
%  wvec  - weight vector (Nx1)
%  ilog  - list of features to apply logarithm (variance stabilization)
%  ftmin - feature normalization vector minimum (1xN)
%  ftmax - feature normalization vector maximum (1xN)
%  fbmin - user feedback minimum
%  fbmax - user feedback maximum
% Outputs:
%  fbnorm - normalized user feedback value
%
% (CC BY-SA 3.0) Max Little, 2014
function fbnorm = features_ufb(ftvec,wvec,ilog,ftmin,ftmax,fbmin,fbmax)

ftvecnorm = ftvec;
ftvecnorm(:,ilog) = log10(ftvecnorm(:,ilog));
ftvecnorm = 2*((ftvecnorm-ftmin)./(ftmax-ftmin))-1;

% Projection of features onto 1st PC to create simple user feedback
fb = ftvecnorm*wvec;

% Normalize and re-scale 0-100, clamp
fbnorm = max(min((fb-fbmin)/(fbmax-fbmin)*100,100),0);
