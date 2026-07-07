function analyzePWVfull(vessel, flow, timeres)

% vessel = branchList(vesselInds,:) where vesselInds is an array of the
% indices for points to be analyzed
% flow = flowPulsatile_val(vesselInds,:)

% Distance
positions = vessel(:,1:3);
distances = cumsum(vecnorm(diff(positions),2,2));

nWaves = size(flow,1);
nFrames = size(flow,2);
scale = 12;
smoothLevel = 15;
nFrames_interp = scale*nFrames;
times = (timeres/scale)*(1:nFrames_interp);
timeresInt = timeres/scale;

% interp params
x = 1:nFrames;
xq = linspace(1,nFrames,nFrames_interp);

% find max of first wave curve
temp1 = interp1(x,smoothdata(flow(1,:),'gaussian',smoothLevel),xq,'cubic');
[maxFl, indMax1] = max(temp1);
midPt = round(length(temp1)/2); %midpoint of flow curve
for i=1:nWaves
    flow_smooth(i,:) = smoothdata(flow(i,:),'gaussian',smoothLevel);
    flow_interp(i,:) = rescale(interp1(x,flow_smooth(i,:),xq,'cubic'));
    waveforms(i,:) = circshift(flow_interp(i,:), midPt-indMax1);
end 
waveforms = rescale(waveforms);
wave1 = waveforms(1,:);
% figure; hold on
% f = 1:nFrames_interp;
% for i = 1:nWaves
%    plot3(f,i*ones(size(f)),waveforms(i,:))
% end

% 20% of max, first to the left
indStart = max(find(wave1(1:midPt) < 0.2)) + 1;
indEnd   = max(find(wave1(indStart:midPt) < 0.8)) + indStart -1;
pts1 = indStart:indEnd;
upstroke1 = wave1(pts1);

% TTP
[~,Idx50] = min(abs(upstroke1-0.5));
ttp1 = times(Idx50+indStart-1);

% TTF 
[p,~] = polyfit(times(pts1),upstroke1,1);
ttf1 = -p(2)/p(1); %y=0 intercept

% TTU
[~,~,ttu1] = sigFit(wave1,times); %see sigFit function below

%Cross correlation
% nothing needed here

% Wavelet 
% dT = 0.001;
% PAD = 1;
% DERIV = 4;
% [y1,PERIOD,~,~,~,~,~] = contwt(wave1,dT,PAD,[],[],[],'dog',DERIV);
% f1 = 1./PERIOD;
% commenting out the wavelet analysis here since it's commented out in the
% loop below, check the other branch to see if contwt is defined there
% (it doesn't seem to be a MATLAB function) -- RVC 6/15/26

%% Get Second Curves
for w = 2:size(waveforms,1) % loop over all flow curves
    wave2 = circshift(waveforms(w,:), midPt-indMax1); % next flow waveform,
    [maxFl, indMax2] = max(wave2);
    % 20% of max, first to the left
    indStart = max(find(wave2(1:indMax2) < 0.2*maxFl)) + 1;
    indEnd   = max(find(wave2(indStart:indMax2) < 0.8*maxFl)) + indStart -1;
    pts2 = indStart:indEnd;
    upstroke2 = wave2(pts2);
    
    %% TTP
    [~,Idx50] = min(abs(upstroke2-0.5));
    ttp2 = times(Idx50+indStart-1);
    TTP(w-1) = ttp2-ttp1;
    
    %% TTF
    [p,~] = polyfit(times(pts2),upstroke2,1);
    ttf2 = -p(2)/p(1); %y=0 intercept
    TTF(w-1) = ttf2-ttf1;
    
    %% TTU
    [~,~,ttu2] = sigFit(wave2,times); %see sigFit function below
    TTU(w-1) = ttu2-ttu1;
    
    %% Cross Correlation
    [Xcorrs,lags] = xcorr(wave2,wave1,'normalized'); %perform cross correlation between flow curves
    [~,maxXcorrIdx] = max(Xcorrs); %get index of max Xcorr value
    XCORR(w-1) = lags(maxXcorrIdx)*timeresInt; %find time lag of Xcorr peak
    
    %% WAVELET
%     [y2,~,~,~,~,~,~] = contwt(wave2,dT,PAD,[],[],[],'dog',DERIV);
% 
%     % cross-spectrum is only performed on 
%     % 1) points during systole
%     % 2) over frequencies between fc and 10 Hz
%     % the points between foot of flow_sl1 and peak of flow_sl2
%     pts = unique([pts1, pts2]);
%     fc = 1/(numel(pts)*1000);     % in Hz, the fundamental freq
%     ind_f = find(f1 >= fc & f1 <= 10);
% 
%     % the complex cross-spectrum
%     y = y1(ind_f,pts).*conj(y2(ind_f,pts));
% 
%     % calculate y_hat
%     y_hat = y/(sum(sum(abs(y))));
% 
%     psi =   atan2(imag(y),real(y))./repmat(f1(ind_f)'*2*pi,[1 length(pts)]);
%     aa =    dot(y_hat,psi);  % eqn 7 in Bargiotas paper
%     tempDelay = abs(sum(aa(:)))*1000 * scale;
% 
%     % re-scale delay from zscore
%     [~, m1, s1] = zscore(wave1(pts));
%     [~, m2, s2] = zscore(wave2(pts));
%     tempDelay = tempDelay/(s1/s2);        
%     WAVELET(w-1) = tempDelay;
end

[TTP, TF] = rmoutliers(TTP,'movmedian', 30, 'ThresholdFactor',2);
dist = distances;
dist(TF) = [];
figure; scatter(dist,TTP','filled'); title('TTP');
[p, ~] = polyfit(dist,TTP,1);
PWV_ttp = 1/p(1); 
disp(['TTP: ' num2str(PWV_ttp) ' m/s']);

[TTF, TF] = rmoutliers(TTF,'movmedian', 30, 'ThresholdFactor',2);
dist = distances;
dist(TF) = [];
figure; scatter(dist,TTF','filled'); title('TTF');
[p, ~] = polyfit(dist,TTF,1);
PWV_ttf = 1/p(1); 
disp(['TTF: ' num2str(PWV_ttf) ' m/s']);

[TTU, TF] = rmoutliers(TTU,'movmedian', 30, 'ThresholdFactor',2);
dist = distances;
dist(TF) = [];
figure; scatter(dist,TTU','filled'); title('TTU');
[p, ~] = polyfit(dist,TTU,1);
PWV_ttu = 1/p(1); 
disp(['TTU: ' num2str(PWV_ttu) ' m/s']);

[XCORR, TF] = rmoutliers(XCORR,'movmedian', 30, 'ThresholdFactor',2);
dist = distances;
dist(TF) = [];
figure; scatter(dist,XCORR','filled'); title('XCORR');
[p, ~] = polyfit(dist,XCORR,1);
PWV_xcorr = 1/p(1); 
disp(['XCORR: ' num2str(PWV_xcorr) ' m/s']);

% [WAVELET, TF] = rmoutliers(WAVELET,'movmedian', 30, 'ThresholdFactor',2);
% dist = distances;
% dist(TF) = [];
% figure; scatter(dist,WAVELET','filled'); title('WAVELET');
% [p, ~] = polyfit(dist,WAVELET,1);
% PWV_wavelet = 1/p(1); 
% disp(['WAVELET: ' num2str(PWV_wavelet) ' m/s']);