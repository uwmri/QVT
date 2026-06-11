function modelWaveMat = PWVceciliamodel(inParams,dv,T_heart,ntimes)

% Computes model velocity waveforms using an algorithm which should be
% identical to Cecilia's model.
%
% n = number of cross-sections
% m = number of timepoints in one velocity waveform
% ntimes = number of timepoints in the output waveform
%
% Input:
% inparams: guesses for the wave coefficients and PWV, PWV last. 
%           Size [1 x m+1]
%
% distance: Vector containing the distanses from the seed-point
%           cross-sections of the arterial tree to the subsequent
%           cross-sections. Usually sorted after distance. In meters. 
%           Size [n x 1]
%      be careful, this has to be a column vector!
%
% T_heart: 60/bpm, i.e. heartbeat repetition time.
%
% ntimes: called m in a previous version
%

m = size(inParams,2)-1;
n = length(dv);
tres = T_heart/ntimes;
pwv = inParams(m+1);

% these times are in units of tres
% we will interpolate query points given values inParams(1:m) at sample
% points
tpts = (0:(ntimes-1))*(m/ntimes);
tptsquery = repmat(tpts,n,1) - (1.0/(pwv*tres))*repmat(dv,1,ntimes);
tqmin = min(tptsquery,[],'all');
tqmax = max(tptsquery,[],'all');
tmin = m*floor(tqmin/m);
tmax = m*ceil((tqmax+1)/m);
ncycle = (tmax-tmin)/m;
tsample = tmin:(tmax-1);
vsample = repmat(inParams(1:m),1,ncycle);

modelWaveMat = zeros(size(tptsquery));
for kk = 1:n
    modelWaveMat(kk,:) = interp1(tsample,vsample,tptsquery(kk,:));
end

end
