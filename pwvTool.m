function varargout = pwvTool(varargin)
% PWVTOOL MATLAB code for pwvTool.fig
%      PWVTOOL, by itself, creates a new PWVTOOL or raises the existing
%      singleton*.
%
%      H = PWVTOOL returns the handle to a new PWVTOOL or the handle to
%      the existing singleton*.
%
%      PWVTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PWVTOOL.M with the given input arguments.
%
%      PWVTOOL('Property','Value',...) creates a new PWVTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pwvTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pwvTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pwvTool

% Last Modified by GUIDE v2.5 03-Nov-2021 20:49:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pwvTool_OpeningFcn, ...
                   'gui_OutputFcn',  @pwvTool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pwvTool is made visible.
function pwvTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pwvTool (see VARARGIN)
global branchList res timeres flowPulsatile_val SavePath 

branchList = varargin{1}; %corners for all planes
res = varargin{2}; %locations/labels for all vessel points
timeres = varargin{3}; %segmentation mask
flowPulsatile_val = varargin{4}; %directory of pcviprData file (imageData)
SavePath = varargin{5}; %pixel resolution (mm)

% Choose default command line output for pwvTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pwvTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pwvTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function branchNumber_Callback(hObject, eventdata, handles)
% hObject    handle to branchNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of branchNumber as text
%        str2double(get(hObject,'String')) returns contents of branchNumber as a double


% --- Executes during object creation, after setting all properties.
function branchNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to branchNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function firstCLpoint_Callback(hObject, eventdata, handles)
% hObject    handle to firstCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of firstCLpoint as text
%        str2double(get(hObject,'String')) returns contents of firstCLpoint as a double


% --- Executes during object creation, after setting all properties.
function firstCLpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to firstCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lastCLpoint_Callback(hObject, eventdata, handles)
% hObject    handle to lastCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lastCLpoint as text
%        str2double(get(hObject,'String')) returns contents of lastCLpoint as a double


% --- Executes during object creation, after setting all properties.
function lastCLpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lastCLpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in analyzeEndpointsButton.
function analyzeEndpointsButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeEndpointsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global branchList timeres flowPulsatile_val SavePath 

text = get(handles.branchNumber,'String');
branchNumbers = sscanf(text, '%g,');
firstCLpoint = str2double(get(handles.firstCLpoint,'String'))+1;
lastCLpoint = str2double(get(handles.lastCLpoint,'String'))+1;

% Define orientation
start = find((branchList(:,4)==branchNumbers(1)&(branchList(:,5)==firstCLpoint)));
stop = find((branchList(:,4)==branchNumbers(end)&(branchList(:,5)==lastCLpoint)));
if branchList(start,3)>branchList(stop,3)
    orientSI = 1;
else
    orientSI = 0;
end 

% Get flow along segments
for b=1:length(branchNumbers)
    vesselInds = find(branchList(:,4)==branchNumbers(b));
    vessel_seg = branchList(vesselInds,:);
    flow_seg = flowPulsatile_val(vesselInds,:);
    
    if orientSI
        if vessel_seg(1,3)<vessel_seg(end,3)
            vessel_seg = flipud(vessel_seg);
            flow_seg = flipud(flow_seg);
        end 
    else
        if vessel_seg(1,3)>vessel_seg(end,3)
            vessel_seg = flipud(vessel_seg);
            flow_seg = flipud(flow_seg);
        end 
    end 
    
    if b==1
        firstPointIdx = find(vessel_seg(:,5)==firstCLpoint);
        vessel_seg = vessel_seg(firstPointIdx:end,:);  
        flow_seg = flow_seg(firstPointIdx:end,:);  
        vessel = vessel_seg;
        flow = flow_seg;
    elseif b==length(branchNumbers)
        lastPointIdx = find(vessel_seg(:,5)==lastCLpoint);
        vessel_seg = vessel_seg(1:lastPointIdx,:);  
        flow_seg = flow_seg(1:lastPointIdx,:); 
        vessel = [vessel; vessel_seg];
        flow = [flow; flow_seg];
    else
        vessel = [vessel; vessel_seg];
        flow = [flow; flow_seg];
    end 
end 
    
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



% --- Executes on button press in analyzeFullButton.
function analyzeFullButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeFullButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global branchList timeres flowPulsatile_val SavePath 

text = get(handles.branchNumber,'String');
branchNumbers = sscanf(text, '%g,');
firstCLpoint = str2double(get(handles.firstCLpoint,'String'))+1;
lastCLpoint = str2double(get(handles.lastCLpoint,'String'))+1;

% Define orientation
start = find((branchList(:,4)==branchNumbers(1)&(branchList(:,5)==firstCLpoint)));
stop = find((branchList(:,4)==branchNumbers(end)&(branchList(:,5)==lastCLpoint)));
if branchList(start,3)>branchList(stop,3)
    orientSI = 1;
else
    orientSI = 0;
end 

% Get flow along segments
for b=1:length(branchNumbers)
    vesselInds = find(branchList(:,4)==branchNumbers(b));
    vessel_seg = branchList(vesselInds,:);
    flow_seg = flowPulsatile_val(vesselInds,:);
    
    if orientSI
        if vessel_seg(1,3)<vessel_seg(end,3)
            vessel_seg = flipud(vessel_seg);
            flow_seg = flipud(flow_seg);
        end 
    else
        if vessel_seg(1,3)>vessel_seg(end,3)
            vessel_seg = flipud(vessel_seg);
            flow_seg = flipud(flow_seg);
        end 
    end 
    
    if b==1
        firstPointIdx = find(vessel_seg(:,5)==firstCLpoint);
        vessel_seg = vessel_seg(firstPointIdx:end,:);  
        flow_seg = flow_seg(firstPointIdx:end,:);  
        vessel = vessel_seg;
        flow = flow_seg;
    elseif b==length(branchNumbers)
        lastPointIdx = find(vessel_seg(:,5)==lastCLpoint);
        vessel_seg = vessel_seg(1:lastPointIdx,:);  
        flow_seg = flow_seg(1:lastPointIdx,:); 
        vessel = [vessel; vessel_seg];
        flow = [flow; flow_seg];
    else
        vessel = [vessel; vessel_seg];
        flow = [flow; flow_seg];
    end 
end 
    
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
dT = 0.001;
PAD = 1;
DERIV = 4;
[y1,PERIOD,~,~,~,~,~] = contwt(wave1,dT,PAD,[],[],[],'dog',DERIV);
f1 = 1./PERIOD;

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
    
    