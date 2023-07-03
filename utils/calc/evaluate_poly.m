function val = evaluate_poly(x,y,z,fit)
%EVALUATE_POLY: Evaluates polynomial data from background phase correction.
% ** See background_phase_correction.m for more information.
%
%   Kevin Johnson, UW-Madison 2014
%   Used by: load_pcvipr.m, background_phase_correction.m
%   Dependencies: NONE

val = 0;
for i = 1:numel(fit.px) %numel(fit.px) = number of polynomial combinations
    val = val + fit.vals(i)*( x.^fit.px(i).*y.^fit.py(i).*z.^fit.pz(i));
end
