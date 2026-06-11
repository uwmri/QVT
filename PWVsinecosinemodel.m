function modelWaveMat = PWVsinecosinemodel(inParams,distance,T_heart,ntimes)

% Computes model velocity waveforms
%
% n = number of cross-sections
% m = ntimes = number of timepoints in one velocity waveform
%
% Input:
% inparams: guesses for the wave coefficients and PWV, PWV last. 
%           Size [1 x m+1]
%
% distance: Vector containing the distanses from the seed-point
%           cross-sections of the arterial tree to the subsequent
%           cross-sections. Usually sorted after distance. In meters. 
%           Size [n x 1]
%
% T_heart: 60/bpm, i.e. heartbeat repetition time.
%
% ntimes: called m in a previous version
%

n = length(distance);

% time is a repeated row vector, m columns repeated n times
% position is a repeated column vector, n rows repeated m times
% result should be tM and xM have the same dimensions

tV = (T_heart/ntimes)*(0:(ntimes-1));
tM = repmat(tV,n,1);
xM = repmat(distance,1,ntimes);

npar = length(inParams);
ncoef = npar-1;
pwv = inParams(npar);

if mod(ncoef,2) ~= 0
    modelWaveMat = [];
    return
end

phaseArg = (2*pi/T_heart)*(tM - xM/pwv);
modelWaveMat = zeros(n,ntimes);

for k = 1:(ncoef/2)
    theta = k*phaseArg;
    modelWaveMat = modelWaveMat + inParams(2*k-1)*sin(theta) + inParams(2*k)*cos(theta);
end

end
