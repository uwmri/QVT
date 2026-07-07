function analyzePWVendpoints(vessel, flow, timeres)

% vessel = branchList(vesselInds,:) where vesselInds is an array of the
% indices for points to be analyzed
% flow = flowPulsatile_val(vesselInds,:)
   
% Distance
positions = vessel(:,1:3);
distances = cumsum(vecnorm(diff(positions),2,2));
DX = distances(end);

nFrames = size(flow,2);
scale = 12;
smoothLevel = 15;
nFrames_interp = scale*nFrames;
times = (timeres/scale)*(1:nFrames_interp);
timeresInt = timeres/scale;

% interp params
x = 1:nFrames;
xq = linspace(1,nFrames,nFrames_interp);

wave1_smooth = smoothdata(flow(1,:),'gaussian',smoothLevel);
wave1_interp = rescale(interp1(x,wave1_smooth,xq,'linear'));
[maxFl, indMax1] = max(wave1_interp);
midPt = round(length(wave1_interp)/2); %midpoint of flow curve
wave1 = circshift(wave1_interp, midPt-indMax1);

wave2_smooth = smoothdata(flow(end,:),'gaussian',smoothLevel);
wave2_interp = rescale(interp1(x,wave2_smooth,xq,'linear'));
[maxFl, indMax2] = max(wave2_interp);
wave2 = circshift(wave2_interp, midPt-indMax1);
%figure; plot(times,wave1); hold on; plot(times,wave2);


% FIRST CURVE
indStart = max(find(wave1(1:midPt) < 0.2)) + 1; % 20% of max, first to the left
indEnd = max(find(wave1(indStart:midPt) < 0.8)) + indStart -1;
pts1 = indStart:indEnd;
upstroke1 = wave1(pts1);

[~,Idx50] = min(abs(upstroke1-0.5)); % TTP
ttp1 = times(Idx50+indStart-1);

[p,~] = polyfit(times(pts1),upstroke1,1); % TTF 
ttf1 = -p(2)/p(1); %y=0 intercept

[~,~,ttu1] = sigFit(wave1,times); %% TTU: see sigFit function below


% SECOND CURVE
indStart = max(find(wave2(1:indMax2) < 0.2*maxFl)) + 1;
indEnd   = max(find(wave2(indStart:indMax2) < 0.8*maxFl)) + indStart -1;
pts2 = indStart:indEnd;
upstroke2 = wave2(pts2);
    
[~,Idx50] = min(abs(upstroke2-0.5));
ttp2 = times(Idx50+indStart-1);
TTP = ttp2-ttp1;
PWV_ttp = DX/TTP; 
disp(['TTP: ' num2str(PWV_ttp) ' m/s']);

[p,~] = polyfit(times(pts2),upstroke2,1);
ttf2 = -p(2)/p(1); %y=0 intercept
TTF = ttf2-ttf1;
PWV_ttf = DX/TTF; 
disp(['TTF: ' num2str(PWV_ttf) ' m/s']);

[~,~,ttu2] = sigFit(wave2,times); %see sigFit function below
TTU = ttu2-ttu1;
PWV_ttu = DX/TTU; 
disp(['TTU: ' num2str(PWV_ttu) ' m/s']);

[Xcorrs,lags] = xcorr(wave2,wave1,'normalized'); %perform cross correlation between flow curves
[~,maxXcorrIdx] = max(Xcorrs); %get index of max Xcorr value
XCORR = lags(maxXcorrIdx)*timeresInt; %find time lag of Xcorr peak
PWV_xcorr = DX/XCORR; 
disp(['XCORR: ' num2str(PWV_xcorr) ' m/s']);

