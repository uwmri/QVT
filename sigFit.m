function [sigmoid,t,t1] = sigFit(meanROI,times)
%%% See the following article by Anas Dogui in JMRI:
% Measurement of Aortic Arch Pulse Wave Velocity in Cardiovascular MR:
% Comparison of Transit Time Estimators and Description of a New Approach

    [~,peak] = max(meanROI); %find max
    upslope = meanROI(1:peak); %find upslope region of flow curve
    t0 = times(1:peak);
    t = linspace(1,times(peak),1000); %interpolate even more
    upslope = interp1(t0,upslope,t); %interpolate upslope
    upslope = rescale(upslope); %normalize from 0 to 1
    [~,MIN] = min(upslope);
    upslope = upslope(MIN:end);
    t = t(MIN:end);
    dt = t(2)-t(1); %new temporal resolution (=0.1)
    midpoint = round(length(upslope)/2);
    
    % c1 = b, c2 = a, c3 = x0, c4 = dx
    % Note that we could assume the equation e^t/(1+e^(t-t0)) since c1=1 and
    % c2=0. However, will keep the same as the Dogui paper.
    sigmoidModel = @(c) c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ) - upslope;
    c0 = [0,1,t(midpoint),dt/2]; %initial params for upslope region
    opts = optimset('Display', 'off'); %turn off display output
    c = lsqnonlin(sigmoidModel,c0,[],[],opts); %get nonlinear LSQ solution
    sigmoid = c(1) + ( (c(2)-c(1)) ) ./ ( 1+exp((t-c(3))./c(4)) ); %calculate our sigmoid fit with our new params

    dy = diff(sigmoid,1);
    dy(end) = [];
    dx = diff(t,1);
    dx(end) = [];
    ddy = diff(sigmoid,2);
    curvature = ddy.*dx./(dx.^2 + dy.^2).^(3/2);
    [~,tIdx] = max(curvature);
    t1 = t(tIdx);
    