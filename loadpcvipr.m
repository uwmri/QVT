function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
    maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
    VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
    timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean] = loadpcvipr(directory,handles)
%LOADPCVIPR: loadpcvipr reads in header information and reconstructed data 
%(velocity, vmean, etc.) and transforms data into usable matlab variables.
%   Used by: paramMap.m
%   Dependencies: load_dat.m, background_phase_correction.m, evaluate_poly.m
%     calc_angio.m, feature_extraction.m, paramMap_params_new.m, slidingThreshold.m

%% Reads PCVIPR Header
filetype = 'dcm';
set(handles.TextUpdate,'String','Loading .DCM Data'); drawnow;

cd 'MR_4D_FLOW_PD_P'
dcms = dir('*.dcm');
P_info = dicominfo(dcms(1).name);

nframes = P_info.CardiacNumberOfImages; %number of reconstructed frames
timeres = P_info.NominalInterval/nframes; %temporal resolution (ms)
res = P_info.PixelSpacing(1); %spatial res (mm) (ASSUMED ISOTROPIC)
matrix(1) = P_info.Rows; %number of pixels in rows
matrix(2) = P_info.Columns; %number of pixels in columns
matrix(3) = 40; %number of slices
vdims = 3; %number of velocity dimensions
VENC = str2double(cell2mat(regexp(P_info.Private_0051_1014,'\d*','Match')));
BGPCdone = 0; %assumes background phase correction was done...

%% Import Phase
set(handles.TextUpdate,'String','Loading Phase/Velocity Data'); drawnow;
for d=1:length(dcms)
    temp(:,:,d) = dicomread(dcms(d).name);
end
temp = reshape(temp,[matrix(1),matrix(2),vdims,nframes,matrix(3)]);
phase = single(permute(temp,[1 2 5 3 4]));

% Conver to velocity
maxx = max(phase,[],'all');
minn = min(phase,[],'all');
v = (2 * (phase-minn)/(maxx-minn) - 1) * VENC; %range values to VENCs
vMean = mean(v,5);
clear temp maxx minn

%% Import Magnitude
set(handles.TextUpdate,'String','Loading Magnitude Data'); drawnow;
cd(basedir)
cd 'MR_4D_FLOW_PD_M'
dcms = dir('*.dcm');
M_info = dicominfo(dcms(1).name);
for d=1:length(dcms) %60 ims per slice at 20 frames
    temp(:,:,d) = dicomread(dcms(d).name); 
end
temp = reshape(temp,[matrix(1),matrix(2),vdims,nframes,matrix(3)]);
mag = single(permute(temp,[1 2 5 3 4]));
MAG = mean(mag,[4 5]); %mean across vdims and nframes
clear temp

%% Import Complex Difference
set(handles.TextUpdate,'String','Loading Complex Difference Data'); drawnow;
cd(basedir)
cd 'MR_4D_FLOW_PD_MSUM'
dcms = dir('*.dcm');
MSUM_info = dicominfo(dcms(1).name);
for d=1:length(dcms)
    temp(:,:,d) = dicomread(dcms(d).name);
end
temp = reshape(temp,[matrix(1),matrix(2),nframes,matrix(3)]);
compd = single(permute(temp,[1 2 4 3]));
timeMIP = mean(compd,4);
clear temp

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
    clear X Y Z poly_fitx poly_fity poly_fitz xrange yrange zrange
end

%% Find optimum global threshold 
step = 0.001; %step size for sliding threshold
UPthresh = 0.8; %max upper threshold when creating Sval curvature plot
SMf = 10;
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
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
        velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
        vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
        = paramMap_params_kmeans(filetype,branchList,matrix,timeMIP,vMean,back,...
        BGPCdone,directory,nframes,res,MAG,handles);
elseif strcmp(SEG_TYPE,'thresh')
    [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
        velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
        vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
        = paramMap_params_new(filetype,branchList,matrix,timeMIP,vMean,back,...
        BGPCdone,directory,nframes,res,MAG,handles);
else
    disp("Incorrect segmentation type selected, please select 'kmeans' or 'thresh'");
end 

set(handles.TextUpdate,'String','All Data Loaded'); drawnow;

return