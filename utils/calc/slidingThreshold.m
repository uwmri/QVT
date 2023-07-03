function [optThresh,segment] = slidingThreshold(image2segment,step,UPthresh,SMf,shiftHM_flag,medFilt_flag)
%SLIDINGTHRESHOLD: Performs sliding threshold segmentation on input image 
%   Step 1: Remove noise by applying median filter
%   Step 2: Sum all active pixels in binary mask for every threshold in x
%   Step 3: Smooth L0 vs threshold (Sval) plot by convolving w/ mean filter
%   Step 4: Calculate curvature plot as a function of threshold
%   Step 5: Calculate FWHM of curvature
%   Step 6: Set optimal threshold as point of maximum curvature + 1 FWHM
%       Used by: loadpcvipr.m, paramMap_params_new.m
%       Dependencies: NONE

if medFilt_flag
    CDcrop = medfilt3(image2segment); %removes noise on surface of vessels
    maxVal = max(CDcrop(:));
else
    CDcrop = image2segment;
    maxVal = max(CDcrop(:));
end 

x = 0:step:UPthresh; %array of all possible thresholds
Sval = zeros(size(x),'single');

iter = 1;
for n = 0:step:UPthresh
    temp = CDcrop>(maxVal.*n); %mask all pixels > current threshold
    Sval(iter) = sum(temp(:)); %take sum of mask (L0 norm)
    iter = iter+1; 
end

% Perform smoothing of sliding-threshold L0-norm curve
y = conv(Sval,ones(SMf,1),'same')./SMf; %convolve w/ mean filter (size SMf)
y = y./max(y); 

% Calculate max curvature of sliding-threshold curve
dx  = gradient(x);  %1st derivative, note that 2nd derivative=0
dy  = gradient(y);  
ddy = gradient(dy); 
num   = dx.*ddy;
denom = dx.*dx + dy.*dy; 
curvatureSM = num ./ (sqrt(denom)).^3; %calculate curvature
curvatureSM(curvatureSM < 0) = 0; %throw negative curvatures out
curvatureSM = conv(curvatureSM,ones(SMf,1),'same')./SMf; %smooth out curvature plot
[~,Idx] = max(curvatureSM); %get index of max curvature

if shiftHM_flag
    halfMaxLeftIndex = find(curvatureSM>=max(curvatureSM)/2, 1, 'first'); 
    halfMaxRightIndex = find(curvatureSM>=max(curvatureSM)/2, 1, 'last');
    fullWidth = halfMaxRightIndex - halfMaxLeftIndex;
    optThresh = maxVal.*x(Idx+fullWidth); %shift max value by FWHM
else
    optThresh = maxVal.*x(Idx); %dont shift by FWHM, stick to max value
end 
segment = CDcrop>optThresh; %segment with optimized threshold

end

