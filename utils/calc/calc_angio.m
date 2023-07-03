function angio = calc_angio(MAG,vMean,Venc)
%CALC_ANGIO: Calculates angiogram from pcvipr header and velocity data
%   
%   Used by: loadpcvipr.m
%   Dependencies: NONE

Vmag = sqrt(sum(vMean.^2,4)); %get speed image
idx = find(Vmag > Venc); %find where flow velocity > VENC.
Vmag(idx) = Venc; %cap Vmag at VENC

%Create complex-difference angiogram
angio = MAG.*sin( (pi/2*Vmag) / Venc);
return
