function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
    maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
    VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
    timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean] = loadHDF5_py(directory,handles)
%LOADHDF5: loadhdf5 reads in python-reconstructed PCVIPR data.
%   Used by: paramMap.m
%   Dependencies: background_phase_correction.m, evaluate_poly.m, calc_angio.m,
%   feature_extraction.m, paramMap_params_new.m, makeITPlane.m, slidingThreshold.m

%% Read HDF5
filetype = 'hdf5_py';
set(handles.TextUpdate,'String','Loading .HDF5 Data'); drawnow;
%cd = h5read(fullfile(directory,'Flow.h5'),'/ANGIO');
mag = h5read(fullfile(directory,'Flow.h5'),'/MAG');
vx = h5read(fullfile(directory,'Flow.h5'),'/VX')*10; % units cm/s to mm/s
vy = h5read(fullfile(directory,'Flow.h5'),'/VY')*10;
vz = h5read(fullfile(directory,'Flow.h5'),'/VZ')*10;

matrix(1) = size(mag,1);                 
matrix(2) = size(mag,2);
matrix(3) = size(mag,3);
nframes = size(mag,4);

disp('Computing time-averaged data')
%CD = mean(cd,4)*32000;
MAG = mean(mag,4)*32000;
V(:,:,:,1) = mean(vx,4); 
V(:,:,:,2) = mean(vy,4);
V(:,:,:,3) = mean(vz,4);

clear mag vx vy vz

%% Reads PCVIPR Header (need to generate accurate header for python code)
%fid = fopen([directory filesep 'pcvipr_header.txt'], 'r');
%delimiter = ' ';
%formatSpec = '%s%s%[^\n\r]'; %read 2 strings(%s%s),end line(^\n),new row(r)
% Info from headers are placed in dataArray, 1x2 cell array.
%dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, ...
%    'MultipleDelimsAsOne', true, 'ReturnOnError', false);
%fclose(fid);

% Converts value column from strings to structure with nums.
%dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
%pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
pcviprHeader.a = 1;
pcviprHeader.b = {'This array is empty?'};
%%%% SPATIAL RESOLUTION ASSUMED TO BE ISOTROPIC (PCVIPR)
%hard-coded for now

%timeres = pcviprHeader.timeres; %temporal resolution (ms)
%res = nonzeros(abs([pcviprHeader.ix,pcviprHeader.iy,pcviprHeader.iz])); %spatial res (mm)
%VENC = pcviprHeader.VENC;
timeres = 50; %50*20 = 1000ms
res = 0.6875; %mm
VENC= 800; %mm/s

%% Reads Data Header
% Checks if automatic background phase correction was performed in recon
fid = fopen([directory filesep 'data_header.txt'], 'r');
fid = 0 % dont read for now
if fid>0
    dataHeader = textscan(fid, formatSpec, 'Delimiter', delimiter, ...
        'MultipleDelimsAsOne', true, 'ReturnOnError', false);
    fclose(fid);
    bgpcIdx = find(contains(dataHeader{1,1},'automatic_BGPC_flag'));
    if isempty(bgpcIdx)
        BGPCdone = 0;
    else
        BGPCdone = dataHeader{1,2}{bgpcIdx};
        BGPCdone = str2double(BGPCdone);
    end 
else
    BGPCdone = 0; %assume automatic backg. phase corr. wasnt done in recon
end 

%% Read Gating Track
find_gating_file = dir('*pcvipr_track');
if ~isempty(find_gating_file)
    name = find_gating_file.name;
    fid = fopen(name);
    gate = fread(fid,'int32','b');
    gate = reshape(gate,[numel(gate)/4 4]);
    gate = sortrows(gate,3);
    fclose(fid);

    % Put data into gatingTrack structure
    gatingTrack.ecg = gate(:,1);

    gatingTrack.unclean = 2*median(gatingTrack.ecg);% this is the RR interval
    within_rr = gatingTrack.ecg < gatingTrack.unclean;
    ecg_filtered = gatingTrack.ecg(within_rr);
    sum_within = size(ecg_filtered,1);
    sum_total = size(gatingTrack.ecg,1);
    gatingTrack.pct_within_rr = 100.0*sum_within / sum_total;  
    gatingTrack.uncleanbpm = 60000./(gatingTrack.unclean);
    gating_rr = nanmedian(gatingTrack.unclean);
    gating_hr = round(nanmedian(gatingTrack.uncleanbpm));
    gating_var = round(nanmedian(gatingTrack.pct_within_rr));
else
    gating_rr = "missing gating file";
    gating_hr = "missing gating file";
    gating_var = "missing gating file";
end 

%% Auto crop images (from MAG data)
% Done to save memory when loading in TR velocity data below.
SUMnumA = squeeze(sum(sum(MAG,1),2)); %1D axial projection
SUMnumS = squeeze(sum(sum(MAG,1),3))'; %1D sagittal projection
SUMnumC = squeeze(sum(sum(MAG,2),3)); %1D coronal projection

% Chop off edges of projections (usually noisy data)
SUMnumA(1:3) = 0; SUMnumA(end-2:end) = 0;
SUMnumS(1:3) = 0; SUMnumS(end-2:end) = 0;
SUMnumC(1:3) = 0; SUMnumC(end-2:end) = 0;

% Normalize values (from 0-1)
SUMnumC = rescale(SUMnumC,'InputMin',min(SUMnumC),'InputMax',max(SUMnumC)); 
BIN = SUMnumC>0.25; % Find where projection crosses thresh value of 0.25 
[~,IDXstart(1)] = max(BIN,[],1); %get first thresh crossing
[~,IDXend(1)] = max(flipud(BIN),[],1); 
IDXend(1) = matrix(1) - IDXend(1) + 1; %get last thresh crossing

SUMnumS = rescale(SUMnumS,'InputMin',min(SUMnumS),'InputMax',max(SUMnumS)); 
BIN = SUMnumS>0.25; % Find where projection crosses thresh value of 0.25 
[~,IDXstart(2)] = max(BIN,[],1); %get first thresh crossing
[~,IDXend(2)] = max(flipud(BIN),[],1); 
IDXend(2) = matrix(2) - IDXend(2) + 1; %get last thresh crossing

SUMnumA = rescale(SUMnumA,'InputMin',min(SUMnumA),'InputMax',max(SUMnumA)); 
BIN = SUMnumA>0.25; % Find where projection crosses thresh value of 0.25 
[~,IDXstart(3)] = max(BIN,[],1); %get first thresh crossing
[~,IDXend(3)] = max(flipud(BIN),[],1); 
IDXend(3) = matrix(3) - IDXend(3) + 1; %get last thresh crossing

% Crop data with new dimensions
MAG = MAG(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
newDIM = size(MAG);
    
%% Read Average Velocity
vMean = zeros(newDIM(1),newDIM(2),newDIM(3),3,'single'); 
% Looped reading of average velocity data.
for n = 1:3
    vMean(:,:,:,n) = V(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3),n);
end

%% Manual Background Phase Correction (if necessary)
back = zeros(size(vMean),'single');
if ~BGPCdone
    set(handles.TextUpdate,'String','Phase Correction with Polynomial'); drawnow;
    
    [poly_fitx,poly_fity,poly_fitz] = background_phase_correction(MAG,vMean(:,:,:,1),vMean(:,:,:,2),vMean(:,:,:,3));
    disp('Correcting data with polynomial');
    xrange = single(linspace(-1,1,size(MAG,1)));
    yrange = single(linspace(-1,1,size(MAG,2)));
    zrange = single(linspace(-1,1,size(MAG,3)));
    [Y,X,Z] = meshgrid(yrange,xrange,zrange);
    
    % Get poly data and correct average velocity for x,y,z dimensions
    back(:,:,:,1) = single(evaluate_poly(X,Y,Z,poly_fitx));
    back(:,:,:,2) = single(evaluate_poly(X,Y,Z,poly_fity));
    back(:,:,:,3) = single(evaluate_poly(X,Y,Z,poly_fitz));
    vMean = vMean - back;
    BGPCdone = 1;
    clear X Y Z poly_fitx poly_fity poly_fitz range1 range2 range3 temp
end

%% Create Angio
% Calculate complex difference angiogram for visualization.
set(handles.TextUpdate,'String','Creating Angiogram'); drawnow;
timeMIP = calc_angio(MAG, vMean, VENC);
% NOTE: timeMIP is an approximated complex difference image.
% The result is nearly equivalent to loading 'CD.dat'.

%%%% OPTION FOR ANGIOGRAM %%%%%
% CD = CD(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
% angio = max(MAG(:))./max(CD(:))*CD - MAG;
% v = sqrt(vx.^2 + vy.^2 + vz.^2);
% v = v(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
% SDv = std(v,0,4);
% threshed = adaptthresh(SDv);
% SDvBW = imcomplement(imbinarize(threshed));
% SE = strel('sphere',3);
% SDvBW = imdilate(SDvBW,SE);
% timeMIP = angio.*SDvBW;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timeMIP = timeMIP(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
% clear vx vy vz v angio v SDv threshed SDvBW SE

%% Find optimum global threshold 
step = 0.001; %step size for sliding threshold
UPthresh = 0.8; %max upper threshold when creating Sval curvature plot
SMf = 10; %smoothing factor
shiftHM_flag = 1; %flag to shift max curvature by FWHM
medFilt_flag = 1; %flag for median filtering of CD image
[~,segment] = slidingThreshold(timeMIP,step,UPthresh,SMf,shiftHM_flag,medFilt_flag);
areaThresh = round(sum(segment(:)).*0.005); %minimum area to keep
conn = 6; %connectivity (i.e. 6-pt)
segment = bwareaopen(segment,areaThresh,conn); %inverse fill holes

clear CDcrop x y SMf temp n halfMaxRightIndex halfMaxLeftIndex Idx BIN V
clear curvatureSM denom num ddy dy ddx dx areaThresh fullWidth conn
clear Sval iter maxThresh newDIM dataArray fid formatSpec delimiter
clear SUMnumA SUMnumC SUMnumS SUMnum step ans dataHeader UPthresh shiftHM_flag

imageData.MAG = MAG;
imageData.CD = timeMIP; 
imageData.V = vMean;
imageData.Segmented = segment;
imageData.pcviprHeader = pcviprHeader;
imageData.gating_rr = gating_rr;
imageData.gating_hr = gating_hr;
imageData.gating_var = gating_var;
clear step UPthresh SMf shiftHM_flag medFilt_flag areaThresh conn ans

%% Feature Extraction
% Get trim and create the centerline data
sortingCriteria = 3; %sorts branches by junctions/intersects 
spurLength = 15; %minimum branch length (removes short spurs)
[~,~,branchList,~] = feature_extraction(sortingCriteria,spurLength,vMean,segment,handles);

% Flow parameter calculation, bulk of code is in paramMap_parameters.m
SEG_TYPE = 'kmeans'; %kmeans or thresh
if strcmp(SEG_TYPE,'kmeans')
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val, ...
        velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
        vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
        = paramMap_params_kmeans(filetype,branchList,matrix,timeMIP,vMean,back, ...
        BGPCdone,directory,nframes,res,MAG,IDXstart,IDXend,handles);
elseif strcmp(SEG_TYPE,'thresh')
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
       velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
       vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
       = paramMap_params_new(filetype,branchList,matrix,timeMIP,vMean,back,...
      BGPCdone,directory,nframes,res,MAG,IDXstart,IDXend,handles);
else
    disp("Incorrect segmentation type selected, please select 'kmeans' or 'thresh'");
end 

set(handles.TextUpdate,'String','All Data Loaded'); drawnow;

return









