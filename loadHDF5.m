function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
    maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
    VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
    timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean] = loadHDF5(directory,handles)
%LOADHDF5: loadhdf5 reads in PCVIPR data saved in h5
%   Used by: paramMap.m
%   Dependencies: background_phase_correction.m, evaluate_poly.m, calc_angio.m,
%   feature_extraction.m, paramMap_params_new.m, slidingThreshold.m

%% Read HDF5
filetype = 'hdf5';
set(handles.TextUpdate,'String','Loading .HDF5 Data'); drawnow;
%cd = h5read(fullfile(directory,'Flow.h5'),'/ANGIO');
mag = h5read(fullfile(directory,'Flow.h5'),'/Data/MAG');
vx = h5read(fullfile(directory,'Flow.h5'),'/Data/comp_vd_1');
vy = h5read(fullfile(directory,'Flow.h5'),'/Data/comp_vd_2');
vz = h5read(fullfile(directory,'Flow.h5'),'/Data/comp_vd_3');

disp('Computing time-averaged data')
%CD = mean(cd,4)*32000;
MAG = single(mag);
V(:,:,:,1) = single(vx);
V(:,:,:,2) = single(vy);
V(:,:,:,3) = single(vz);

clear mag vx vy vz


%% Reads PCVIPR Header from H5 (only to get gating data)
headerh5 = h5info(fullfile(directory,'Flow.h5'),'/Header');
for i = 1:size(headerh5.Attributes,1)
    new_field = headerh5.Attributes(i).Name;
    if(strcmp(new_field,'2d_flag')) 
     %skip    
    else
        attr_val = headerh5.Attributes(i).Value;    
        if (isnumeric(attr_val))
            attr_val = double(attr_val);
        end
        pcviprHeader.(new_field) = attr_val; 
    end   
end

%%%% SPATIAL RESOLUTION ASSUMED TO BE ISOTROPIC (PCVIPR)
%timeres = pcviprHeader.timeres; %temporal resolution (ms)
%res = nonzeros(abs([pcviprHeader.ix,pcviprHeader.iy,pcviprHeader.iz])); %spatial res (mm)
%VENC = pcviprHeader.VENC;
%matrix(1) = pcviprHeader.matrixx; %number of pixels in rows (ASSUMED ISOTROPIC)
%matrix(2) = pcviprHeader.matrixy;
%matrix(3) = pcviprHeader.matrixz;
%nframes = pcviprHeader.frames;
% Checks if automatic background phase correction was performed in recon
BGPCdone = pcviprHeader.automatic_BGPC_flag;
BGPCdone = 0; %assume automatic backg. phase corr. wasnt done in recon

imageData.gating_rr = pcviprHeader.median_rr_interval_ms;
imageData.gating_hr = pcviprHeader.expected_hr_bpm;
imageData.gating_var = round(pcviprHeader.vals_within_expected_rr_pct);

% for some reason now the h5 header has the wrong matrix size, something
% change in PSD
 
%% Reads PCVIPR Header
set(handles.TextUpdate,'String','Loading .DAT Data'); drawnow;
fid = fopen([directory filesep 'pcvipr_header.txt'], 'r');
delimiter = ' ';
formatSpec = '%s%s%[^\n\r]'; %read 2 strings(%s%s),end line(^\n),new row(r)
% Info from headers are placed in dataArray, 1x2 cell array.
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, ...
    'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fid);

% Converts value column from strings to structure with nums.
dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);

%%%%% SPATIAL RESOLUTION ASSUMED TO BE ISOTROPIC (PCVIPR)
nframes = pcviprHeader.frames; %number of reconstructed frames
timeres = pcviprHeader.timeres; %temporal resolution (ms)
res = nonzeros(abs([pcviprHeader.ix,pcviprHeader.iy,pcviprHeader.iz])); %spatial res (mm)
matrix(1) = pcviprHeader.matrixx; %number of pixels in rows (ASSUMED ISOTROPIC)
matrix(2) = pcviprHeader.matrixy;
matrix(3) = pcviprHeader.matrixz;
VENC = pcviprHeader.VENC;
clear dataArray fid

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
    clear X Y Z poly_fitx poly_fity poly_fitz range1 range2 range3 temp
end

%% Create Angio
% Calculate complex difference angiogram for visualization.
set(handles.TextUpdate,'String','Creating Angiogram'); drawnow;
timeMIP = calc_angio(MAG, vMean, VENC);
% NOTE: timeMIP is an approximated complex difference image.
% The result is nearly equivalent to loading 'CD.dat'.

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

% save raw (cropped) images to imageData structure (for Visual Tool)
imageData.MAG = MAG;
imageData.CD = timeMIP; 
imageData.V = vMean;
imageData.Segmented = segment;
imageData.pcviprHeader = pcviprHeader;



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
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val, ...
        flowPulsatile_val,velMean_val,VplanesAllx,VplanesAlly,VplanesAllz, ...
        r,timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, ...
        bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
        = paramMap_params_new(filetype,branchList,matrix,timeMIP,vMean,...
        back,BGPCdone,directory,nframes,res,MAG,IDXstart,IDXend,handles);
else
    disp("Incorrect segmentation type selected, please select 'kmeans' or 'thresh'");
end 

set(handles.TextUpdate,'String','All Data Loaded'); drawnow;

return




