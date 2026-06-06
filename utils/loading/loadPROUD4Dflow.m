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
BGPCdone = 0;
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

warning('off','all');
% find all reconstructions in the same folder
files = dir([directory, fBase, '*.par']);
% quickly load in headers to parse data type
keepPar = zeros(size(files));
delayedReconFlag = 0;
for ii = 1:length(files)
    [~,h] = readrec_V4_2(fullfile(files(ii).folder,files(ii).name),'quiet','par');
    if unique(h.tbl(:,5)) == 0  % we only have type = MAG
        if sum(unique(h.tbl(:,6)) == [2;4]) == 2    % we have magnitude and pcmra, this is a delayed recon
            delayedReconFlag = 1; keepPar(ii) = 2;
        end
    else
        keepPar(ii) = 1;
    end
end

%% grab each parrec and save corresponding data
disp('Loading data')
if delayedReconFlag
    for ii = 1:length(files)
        if keepPar(ii) > 0
            if keepPar(ii) == 2
                PARRECFILE = fullfile(directory,files(ii).name);
                [IMG,h1] = readrec_V4_2(PARRECFILE);
                MAG = double(squeeze(IMG(:,:,:,:,:,1,:))); clear IMG;
            else
                PARRECFILE = fullfile(directory,files(ii).name);
                [IMG,h3] = readrec_V4_2(PARRECFILE);
                [~,I] = max(h3.pevelocity);
                switch I
                    case 1  % measurement direction velocity, in mm/s
                        vx = squeeze(double(IMG))*10;
                    case 2  % phase direction velocity, in mm/s
                        vy = squeeze(double(IMG))*10;
                    case 3  % slice direction velocity, in mm/s
                        vz = squeeze(double(IMG))*10;
                end
            end
        end
    end
else
    PARRECFILE = fullfile(directory,[fBase, '1.rec']);
    [IMG,h1] = readrec_V4_2(PARRECFILE);
    IMG = single(IMG);
    % this is the readout direction
    vx = squeeze(IMG(:,:,:,:,:,2,:));
    mag1 = squeeze(IMG(:,:,:,:,:,1,:)); clear IMG;

    PARRECFILE = fullfile(directory,[fBase, '2.rec']);
    [IMG,h2] = readrec_V4_2(PARRECFILE);
    IMG = single(IMG);
    % this is the phase direction
    vy = squeeze(IMG(:,:,:,:,:,2,:));
    mag2 = squeeze(IMG(:,:,:,:,:,1,:)); clear IMG;

    PARRECFILE = fullfile(directory,[fBase, '3.rec']);
    [IMG,h3] = readrec_V4_2(PARRECFILE);
    IMG = single(IMG);
    % this is the slice direction
    vz = squeeze(IMG(:,:,:,:,:,2,:));
    mag3 = squeeze(IMG(:,:,:,:,:,1,:)); clear IMG;
    warning('on','all');

    MAG = double(mean(cat(5,mag1,mag2,mag3),5));
    clear mag1 mag2 mag3
    vx = double(vx); vy = double(vy); vz = double(vz);
end

nframes = h3.nphases;                                       % number of reconstructed frames
timeres = max(h3.tbl(:,h3.tblcols.ttime))/nframes;          % temporal resolution, in ms
fov = h3.fov;                                               % Field of view in cm
matrix = round([h3.nrows h3.ncols h3.nslices]);                % number of pixels in row,col,slices
VENC = max(h3.pevelocity)*10;                               % venc, in mm/s
res = h3.pixdim;                                         % the reconstructed resolution
phasedir = h3.prepdir;
if strcmp(phasedir,'Anterior-Posterior')
    phasedir = 'AP';
end
if strcmp(phasedir,'Right-Left')
    phasedir = 'RL';
end
if strcmp(phasedir,'Foot-Head')
    phasedir = 'FH';
end
angulations = h3.tbl(1,17:19);                              % the rotated axis, ap, fh, rl

%% manually change velocity directions depending on scan orientations

% velocity directions correspond to the following:
% vx: in-plane up-down
% vy: in-plane right-left
% vz: through-plane
if delayedReconFlag
    switch h3.tbl(1,26) % orientation number (1 - axial, 2 - sagittal, 3 - coronal)
        case 1
            ori.label = 'axial';
            angulations = [angulations(3) angulations(1) angulations(2)]; % rl, ap, fh
            if strcmp('AP',phasedir)
                ori.vxlabel = 'R-L';
                ori.vylabel = 'A-P';
                tmp = vx; vx = vy; vy = tmp; clear tmp;
            else    % phasedir == 'RL'
                ori.vxlabel = 'A-P';
                ori.vylabel = 'R-L';
                angulations = [angulations(2) angulations(1) angulations(3)]; % ap, rl, fh
            end
            ori.vzlabel = 'F-H';
        case 2
            ori.label = 'sagittal';
            angulations = [angulations(2) angulations(1) angulations(3)]; % fh, ap, rl
            if strcmp('AP',phasedir)
                ori.vxlabel = 'H-F';   % rows: Head(top)→Foot(bottom) in HFS
                ori.vylabel = 'A-P';
            else    % phasedir == 'FH'
                ori.vxlabel = 'A-P';
                ori.vylabel = 'H-F';   % cols become rows after swap: Head→Foot
                angulations = [angulations(2) angulations(1) angulations(3)]; % ap, fh, rl
            end
            tmp = vx; vx = vz; vz = tmp; clear tmp;
            vx = -vx;
            % Slice direction: determine from per-slice RL offcentre
            rl_step = h3.tbl(2,22) - h3.tbl(1,22);
            if rl_step < 0
                ori.vzlabel = 'L-R';
            else
            ori.vzlabel = 'R-L';
            end
        case 3
            ori.label = 'coronal';
            angulations = [angulations(2) angulations(3) angulations(1)]; % fh, rl, ap
            if strcmp('RL',phasedir)
                ori.vxlabel = 'F-H';
                ori.vylabel = 'R-L';
            else    % phasedir == 'FH'
                ori.vxlabel = 'R-L';
                ori.vylabel = 'F-H';
                angulations = [angulations(2) angulations(1) angulations(3)]; % rl, fh, ap
            end
            tmp = vx; vx = vz; tmp2 = vy; vy = tmp; vz = tmp2; clear tmp tmp2;
            ori.vzlabel = 'A-P';
            vx = -vx;
%             vz = -vz;
    end
else
    switch h3.tbl(1,26) % orientation number (1 - axial, 2 - sagittal, 3 - coronal)
        case 1
            ori.label = 'axial';
            angulations = [angulations(3) angulations(1) angulations(2)]; % rl, ap, fh
            if strcmp('AP',phasedir)
                ori.vxlabel = 'R-L';
                ori.vylabel = 'A-P';
                tmp = vx; vx = vy; vy = tmp; clear tmp;
            else    % phasedir == 'RL'
                ori.vxlabel = 'A-P';
                ori.vylabel = 'R-L';
                angulations = [angulations(2) angulations(1) angulations(3)]; % ap, rl, fh
            end
            ori.vzlabel = 'F-H';
        case 2
            ori.label = 'sagittal';
            angulations = [angulations(2) angulations(1) angulations(3)]; % fh, ap, rl
            if strcmp('AP',phasedir)
                ori.vxlabel = 'H-F';   % rows: Head(top)→Foot(bottom) in HFS
                ori.vylabel = 'A-P';
            else    % phasedir == 'FH'
                ori.vxlabel = 'A-P';
                ori.vylabel = 'H-F';   % cols become rows after swap: Head→Foot
                tmp = vx; vx = vy; vy = tmp; clear tmp;
                angulations = [angulations(2) angulations(1) angulations(3)]; % ap, fh, rl
            end
            % Slice direction: determine from per-slice RL offcentre
            % RL+ = patient LEFT; if rl decreases (slice 1 left, slice N right) → L-R
            rl_step = h3.tbl(2,22) - h3.tbl(1,22); % offcentre_rl step slice1→slice2
            if rl_step < 0
                ori.vzlabel = 'L-R'; % slices increase toward Right
            else
                ori.vzlabel = 'R-L'; % slices increase toward Left
            end
            vx = -vx;
            vz = -vz;
        case 3
            ori.label = 'coronal';
            angulations = [angulations(2) angulations(3) angulations(1)]; % fh, rl, ap
            vx = -vx;
            if strcmp('RL',phasedir)
                ori.vxlabel = 'H-F';
                ori.vylabel = 'R-L';
            else    % phasedir == 'FH'
                ori.vxlabel = 'R-L';
                ori.vylabel = 'F-H';
                angulations = [angulations(2) angulations(1) angulations(3)]; % rl, fh, ap
            end
            ori.vzlabel = 'A-P';
    end
end
ori.angulations = angulations;

%% Refine vzlabel for axial/coronal using actual slice stepping direction
% and patient position (h3.position: HFS/FFS/HFP/FFP)
if strcmp(ori.label,'axial')
    % Axial: slices step in FH direction
    % HFS/HFP: feet first → fh increases (feet=low, head=high) → 'F-H'
    % FFS/FFP: head first → fh decreases → 'H-F'
    if any(strncmpi(h3.position, {'HFS','HFP'}, 3))
        ori.vzlabel = 'F-H';
    else
        ori.vzlabel = 'H-F';
    end
elseif strcmp(ori.label,'coronal')
    % Coronal: slices step in AP direction - derive from per-slice AP offcentre
    ap_step = h3.tbl(2,20) - h3.tbl(1,20); % offcentre_ap step
    if ap_step > 0
        ori.vzlabel = 'A-P'; % slices increase toward Posterior
    else
        ori.vzlabel = 'P-A'; % slices increase toward Anterior
    end
end
%%
v = cat(5,vx,vy,vz); v = permute(v, [1 2 3 5 4]);
if ~delayedReconFlag    % delayed recon is already scaled
    % scale velocity
    v = (v./3142)*VENC;
    vx = (vx./3142)*VENC;
    vy = (vy./3142)*VENC;
    vz = (vz./3142)*VENC;
end

% take the means
vMean = mean(v,5);
MAG = MAG./max(MAG(:));
MAG = mean(MAG,4);

% Calculate complex difference angiogram for visualization.
filetype = 'rec';
set(handles.TextUpdate,'String','Creating Angiogram'); drawnow;
timeMIP = calc_angio(MAG, vMean, VENC);

% NOTE: timeMIP is an approximated complex difference image.
% The result is nearly equivalent to loading 'CD.dat'.

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