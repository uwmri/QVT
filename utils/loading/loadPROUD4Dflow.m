function [nframes,matrix,res,timeres,VENC,area_val,diam_val,flowPerHeartCycle_val, ...
    maxVel_val,PI_val,RI_val,flowPulsatile_val,velMean_val, ...
    VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
    timeMIPcrossection,segmentFull,vTimeFrameave,MAGcrossection, imageData, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean] = loadPROUD4Dflow(directory,handles)
%loadPROUD4Dflow: loadPROUD4Dflow reads in header information and reconstructed data 
%(velocity, vmean, etc.) and transforms data into usable matlab variables.
%   Used by: paramMap.m
%   Dependencies: readrec_V4_2.m calc_angio.m, feature_extraction.m, 
%                   paramMap_params_new.m, slidingThreshold.m

fBase = dir([directory filesep '*.rec']);
% loop through and remove and PCMRA-based recons
name2Delete = [];
for ii = 1:length(fBase)
    if contains(fBase(ii).name, '_PCMRA_')
        name2Delete = [name2Delete,ii];
    end
end
fBase(name2Delete) = [];
fBase = fBase(1).name(1:end-5);

%% grab each parrec and save corresponding data
disp('Loading data')
PARRECFILE = fullfile(directory,[fBase, '1.rec']);
[IMG1,~] = readrec_V4_2(PARRECFILE);
IMG1 = double(IMG1);
% this is the readout direction
vx = squeeze(IMG1(:,:,:,:,:,2,:));  
mag1 = squeeze(IMG1(:,:,:,:,:,1,:));

PARRECFILE = fullfile(directory,[fBase, '2.rec']);
[IMG2,~] = readrec_V4_2(PARRECFILE);
IMG2 = double(IMG2);
% this is the phase direction
vy = squeeze(IMG2(:,:,:,:,:,2,:));
mag2 = squeeze(IMG2(:,:,:,:,:,1,:));

PARRECFILE = fullfile(directory,[fBase, '3.rec']);
[IMG3,header] = readrec_V4_2(PARRECFILE);
IMG3 = double(IMG3);
% this is the slice direction
vz = squeeze(IMG3(:,:,:,:,:,2,:));
mag3 = squeeze(IMG3(:,:,:,:,:,1,:));
warning('on','all');

MAG = mean(cat(5,mag1,mag2,mag3),5);
clear mag1 mag2 mag3 IMG1 IMG2 IMG3

nframes = header.nphases;                                       % number of reconstructed frames
timeres = max(header.tbl(:,header.tblcols.ttime))/nframes;      % temporal resolution, in ms
fov = header.fov;                                               % Field of view in cm
VENC = max(header.pevelocity)*10;                               % venc, in mm/s
res = header.pixdim;                                            % the reconstructed resolution, in mm
ori = header.tbl(1,26);                                         % orientation number (1 - axial, 2 - sagittal, 3 - coronal)
phasedir = header.prepdir;

% %% manually change velocity directions depending on scan orientations
% % sagittal aortic scan with phase direction = AP
% if (ori == 2 && strcmp('AP',phasedir))
%     vx = -vx;
%     vz = -vz;
% end
% 
% % coronal carotid scan with phase direction = RL
% if (ori == 3 && strcmp('RL',phasedir))
%     vx = -vx;
%     vz = -vz;
% end
% 
% % axial brain scan with phase direction = RL
% if (ori == 1 && strcmp('RL',phasedir))
%     vx = -vx;
%     vy = -vy;
%     vz = -vz;
% end
% 
% % axial whole-heart scan
% if (ori == 1 && strcmp('AP',phasedir))
%     % do nothing!
% end

%%
% scale velocity
vx = (vx./3142)*VENC;
vy = (vy./3142)*VENC;
vz = (vz./3142)*VENC;
v = cat(5,vx,vy,vz); v = permute(v, [1 2 3 5 4]);

matrix = size(MAG,1:3);
BGPCdone = 1;

% take the means
vMean = mean(v,5);
MAG = MAG./max(MAG(:));
MAG = mean(MAG,4);
%%%%%%%%%%%%%%%%%

%% Create Angio
% Calculate complex difference angiogram for visualization.
filetype = 'rec';
set(handles.TextUpdate,'String','Creating Angiogram'); drawnow;
timeMIP = calc_angio(mean(MAG,4), vMean, VENC);
% NOTE: timeMIP is an approximated complex difference image.
% The result is nearly equivalent to loading 'CD.dat'.

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
% imageData.Header = pcviprHeader;


%% Feature Extraction
% Get trim and create the centerline data
sortingCriteria = 3; %sorts branches by junctions/intersects 
spurLength = 15; %minimum branch length (removes short spurs)
[~,~,branchList,~] = feature_extraction(sortingCriteria,spurLength,vMean,segment,handles);

back = zeros(size(vMean),'single'); % skipping the background phase correction
IDXstart = [1 1 1]; IDXend = [matrix(2) matrix(1) matrix(3)]; % no cropping for PROUD
% Flow parameter calculation, bulk of code is in paramMap_parameters.m
% [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
%     velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
%     vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
%     = paramMap_params_kmeans(filetype,branchList,matrix,timeMIP,vMean,back,...
%     BGPCdone,directory,nframes,res,MAG,IDXstart,IDXend,handles,vx,vy,vz);

[area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
    velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
    vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
    = paramMap_params_new(filetype,branchList,matrix,timeMIP,vMean,back,...
    BGPCdone,directory,nframes,res,MAG,IDXstart,IDXend,handles,vx,vy,vz);

set(handles.TextUpdate,'String','All Data Loaded'); drawnow;

return