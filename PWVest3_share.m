function qDiffSqr = PWVest3_share(inParams,distance,waveMat,tRes,weight)

% Estimates an "underlying" velocity waveform and a pwv
%
% n = number of cross-sections
% m = number of timepoints in one velocity waveform
%
% Input:
% inparams: Initial guesses for the velocity waveform and PWV, waveform
%           datapoints first and PWV last. 
%           Size [1 x m+1]
%
% distance: Vector containing the distanses from the seed-point
%           cross-sections of the arterial tree to the subsequent
%           cross-sections. Usually sorted after distance. In meters. 
%           Size [n x 1]
%
% waveMat:  Matrix containing the velocity waveform for each
%           cross-section. Same sorting as in "distance".
%           Usually normalized to have zero mean and unit std.
%           Size [n x m]
%
% tRes:     Time between consecutive frames in seconds.
%
% scaling:  Weight factor for each cross-section.
%           
%           "cross-sections-area / scalingFactor^2", where scalingFactor is
%           what was used to variance-normalize waveforms.
%
%           scaling can probably be set to ones(n,1) for a single vessel of
%           approximatly equal size along its path.
%                  
%           Size [n x 1]
%
% Example Usage:
% fun1=@(inParams)PWVest3_share(inParams,d,waveMat,tRes,w); 
% pwv0 = 10; %initial guess of pwv
% mean_flow = mean(waveMat); %initial guess of waveform
% initialGuess=[mean_flow, pwv0]; 
% options = optimset('Display','iter', 'TolCon', 1e-7, 'TolX', 1e-7, 'TolFun', 1e-7,'DiffMinChange', 1e-3);
% [params,exitflag,output] = fminunc(fun1,initialGuess, options);%,options1);
% pwv = params(end) % have a look at the PWV
%
% FOR REFERENCE SEE BJÖRNFOT et al JCBFM 2021

m = size(inParams,2)-1;

tV = 0:tRes:(tRes*3*m-tRes);

velocity = inParams(1:m);

pwv = inParams(m+1);

region = [m+1:m*2];

deltaT = repmat(tV,length(distance),1)-repmat(distance,1,length(tV))/pwv;

vShift=interp1(tV,repmat(velocity,1,3),deltaT,'linear');

weightm = repmat(weight,1,m);

vShift = vShift(:,region);

qDiffSqr=sum(weightm(:).*(vShift(:)-waveMat(:)).^2);

end
