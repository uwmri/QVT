function [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
    velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
    vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes,pixelSpace,segmentFullJS] ...
    = paramMap_params_threshS(varargin)
%filetype,branchList,matrix,timeMIP,vMean,back,...
  %  BGPCdone,directory,nframes,res,MAG,IDXstart,IDXend,handles)
%PARAMMAP_PARAMS_NEW: Create tangent planes and calculate hemodynamics
%   Based on a sliding threshold segmentation algorithm developed by Carson Hoffman
%   Used by: loadpcvipr.m
%   Dependencies: slidingThreshold.m
%% Argument Inputs
filetype = varargin{1};
branchList = varargin{2};
matrix = varargin{3};
timeMIP = varargin{4};
vMean = varargin{5};
back = varargin{6};
BGPCdone = varargin{7};
directory = varargin{8};
nframes = varargin{9};
res = varargin{10};
MAG = varargin{11};
handles = varargin{12};
v = varargin{13};
sliceSpace=varargin{14};
%% Tangent Plane Creation
r=10; %radius of plan generation (in pixels)
InterpVals = 4; %choose the interpolation between pixel values
width = r*InterpVals.*2+1; %width of plane in pixels
[x_full,y_full,z_full,x,y,z,Tangent_V,Planes] = create_planes(branchList,r,single(matrix),InterpVals,width);
segments=length(branchList); %number of segments

%% Calculate out of plane length distortion
xyNorm=[0 0 1];
pixelSpace=zeros([segments,1]);
for i=1:length(Tangent_V(:,1))
    T=Tangent_V(i,:)';
    DOT=xyNorm*T;
    RAD=acos(DOT);
    pixelSpace(i)=res+(sliceSpace-res)*sin(RAD); %New pixel space
end
clear DOT RAD T xyNorm
%% Interpolate Matrices Into the Planes
set(handles.TextUpdate,'String','Interpolating Data');drawnow;
% Get interpolated velocity from 3 directions, multipley w/ tangent vector
[v1] = interp_vol_to_planes(vMean(:,:,:,1),x,y,z,x_full,y_full,z_full,width,segments);
[v2] = interp_vol_to_planes(vMean(:,:,:,2),x,y,z,x_full,y_full,z_full,width,segments);
[v3] = interp_vol_to_planes(vMean(:,:,:,3),x,y,z,x_full,y_full,z_full,width,segments);
temp = zeros([size(v1),3]); % used to hold velocity data information
temp(:,:,1) = bsxfun(@times,v1,Tangent_V(:,1)); %dot product here
temp(:,:,2) = bsxfun(@times,v2,Tangent_V(:,2)); %make veloc. through-plane
temp(:,:,3) = bsxfun(@times,v3,Tangent_V(:,3)); %(mm/s)
% Through-plane mean SPEED for all points (tangent vector dotted with 3D vel)
vTimeFrameave = sqrt(temp(:,:,1).^2 + temp(:,:,2).^2 + temp(:,:,3).^2); %(mm/s)
% Interpolation for complex difference data
[timeMIPcrossection] = interp_vol_to_planes(timeMIP,x,y,z,x_full,y_full,z_full,width,segments);
% Interpolation for magnitude data
[MAGcrossection] = interp_vol_to_planes(MAG,x,y,z,x_full,y_full,z_full,width,segments);

clear v1 v2 v3 MAG timeMIP temp vtimeave
%% In-Plane Segmentation
set(handles.TextUpdate,'String','Performing In-Plane Segmentation');drawnow;
[area_val,diam_val,segmentFull] = segment_cross_section_thresh(segments,width,timeMIPcrossection,vTimeFrameave,MAGcrossection,res,r,InterpVals,pixelSpace);
if length(varargin) == 16
    Exseg=varargin{15};
    [Excross] = interp_vol_to_planes(single(Exseg),x,y,z,x_full,y_full,z_full,width,segments);
    [~,~,segmentFullJS] = segment_cross_section_ExternalSeg(segments,width,Excross,res,r,InterpVals);
else
    segmentFullJS=segmentFull;
end
%% Extract Time-Resolved Velocities
% Initialize time-resolved hemodynamic parameters 
% Sliding Threshold
flowPulsatile_val = zeros(size(area_val,1),nframes);
maxVelFrame = zeros(size(area_val,1),nframes);
velPulsatile_val = zeros(size(area_val,1),nframes);
bnumMeanFlow = zeros(max(branchList(:,4)),1);
bnumStdvFlow = zeros(max(branchList(:,4)),1);   
% Initialize time-resolved velocity matrix (not interpolated yet)
VplanesAllx = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
VplanesAlly = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
VplanesAllz = zeros([length(branchList),(r.*2+1).^2 nframes],'single');
% Extract single interp location Idx
ROW = repmat((1:InterpVals:width)',[1 r*2+1]); %replicate up-down
COL = repmat(1:InterpVals*(width):(width)^2,[r*2+1 1])-1; %rep. lf-rt
idCOL = reshape(ROW+COL,[1 numel(ROW)]); %interp query points
for j = 1:nframes
    if strcmp(filetype,'dat')
        set(handles.TextUpdate,'String',['Calculating Quantitative Params Frame: ' num2str(j) '/' num2str(nframes)]);drawnow;
        % Load x,y,z components of velocity - single frame
        vx = load_dat(fullfile(directory, ['ph_' num2str(j-1,'%03i') '_vd_1.dat']),[matrix(1) matrix(2) matrix(3)]);
        vy = load_dat(fullfile(directory, ['ph_' num2str(j-1,'%03i') '_vd_2.dat']),[matrix(1) matrix(2) matrix(3)]);
        vz = load_dat(fullfile(directory, ['ph_' num2str(j-1,'%03i') '_vd_3.dat']),[matrix(1) matrix(2) matrix(3)]);
        % Crop velocity using crop indices from load_pcvipr.m
        vz = vz(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3)); 
        vx = vx(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
        vy = vy(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
        % Correct background phase (if autoBGPC_flag is off)
        if ~BGPCdone
            vx = vx - back(:,:,:,1); %subtract off background phase in x-dir
            vy = vy - back(:,:,:,2);
            vz = vz - back(:,:,:,3);
        end 
    elseif strcmp(filetype,'hdf5')
        set(handles.TextUpdate,'String',['Calculating Quantitative - Parameters Time Frame: ' num2str(j) '/' num2str(nframes)]);drawnow;
        xvel_label = append('/Data/',['ph_' num2str(j-1,'%03i') '_vd_1']);
        yvel_label = append('/Data/',['ph_' num2str(j-1,'%03i') '_vd_2']);
        zvel_label = append('/Data/',['ph_' num2str(j-1,'%03i') '_vd_3']);
        % Load x,y,z components of velocity (cropped) - single frame
        vx = single(h5read(fullfile(directory,'Flow.h5'),xvel_label, ... 
            [IDXstart(1),IDXstart(2),IDXstart(3)], ...
            [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1]));
        vy = single(h5read(fullfile(directory,'Flow.h5'),yvel_label, ... 
            [IDXstart(1),IDXstart(2),IDXstart(3)], ...
            [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1]));
        vz = single(h5read(fullfile(directory,'Flow.h5'),zvel_label, ... 
            [IDXstart(1),IDXstart(2),IDXstart(3)], ...
            [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1]));
        % Correct background phase (if autoBGPC_flag is off)
        if ~BGPCdone
            vx = vx - back(:,:,:,1); %subtract off background phase in x-dir
            vy = vy - back(:,:,:,2);
            vz = vz - back(:,:,:,3);
        end
    elseif strcmp(filetype,'dcm')
        set(handles.TextUpdate,'String',['Calculating Quantitative Params Frame: ' num2str(j) '/' num2str(nframes)]);drawnow;
        vx = squeeze(v(:,:,:,1,j));
        vy = squeeze(v(:,:,:,2,j));
        vz = squeeze(v(:,:,:,3,j));
    else
        set(handles.TextUpdate,'String',['Calculating Quantitative (python H5 export) - Parameters Time Frame: ' num2str(j) '/' num2str(nframes)]);drawnow;
        % use for flow python
        % Load x,y,z components of velocity (cropped) - single frame
         vx = h5read(fullfile(directory,'Flow.h5'),'/VX', ... 
             [IDXstart(1),IDXstart(2),IDXstart(3),j], ...
             [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1,1])*10; % Flow units are cm/s, convert to mm/s (*10)
         vy = h5read(fullfile(directory,'Flow.h5'),'/VY', ... 
             [IDXstart(1),IDXstart(2),IDXstart(3),j], ...
             [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1,1])*10;
         vz = h5read(fullfile(directory,'Flow.h5'),'/VZ', ... 
             [IDXstart(1),IDXstart(2),IDXstart(3),j], ...
             [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1,1])*10;

        % Correct background phase (if autoBGPC_flag is off)
        if ~BGPCdone
            vx = vx - back(:,:,:,1); %subtract off background phase in x-dir
            vy = vy - back(:,:,:,2);
            vz = vz - back(:,:,:,3);
        end
    end 
        %use for flow python
        % Load x,y,z components of velocity (cropped) - single frame
%         vx = h5read(fullfile(directory,'Flow.h5'),'/VX', ... 
%             [IDXstart(1),IDXstart(2),IDXstart(3),j], ...
%             [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1,1]);
%         vy = h5read(fullfile(directory,'Flow.h5'),'/VY', ... 
%             [IDXstart(1),IDXstart(2),IDXstart(3),j], ...
%             [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1,1]);
%         vz = h5read(fullfile(directory,'Flow.h5'),'/VZ', ... 
%             [IDXstart(1),IDXstart(2),IDXstart(3),j], ...
%             [IDXend(1)-IDXstart(1)+1,IDXend(2)-IDXstart(2)+1,IDXend(3)-IDXstart(3)+1,1]);
    % Interpolation of time-resolved velocities
    [v1] = interp_vol_to_planes(vx,x,y,z,x_full,y_full,z_full,width,segments);
    [v2] = interp_vol_to_planes(vy,x,y,z,x_full,y_full,z_full,width,segments);
    [v3] = interp_vol_to_planes(vz,x,y,z,x_full,y_full,z_full,width,segments);
    v1 = bsxfun(@times,v1,Tangent_V(:,1)); %dot product here
    v2 = bsxfun(@times,v2,Tangent_V(:,2)); %make velocity through-plane
    v3 = bsxfun(@times,v3,Tangent_V(:,3)); %mm/s)
    VplanesAllx(:,:,j) = v1(:,idCOL); %uninterpolated TR vel. (mm/s)
    VplanesAlly(:,:,j) = v2(:,idCOL);
    VplanesAllz(:,:,j) = v3(:,idCOL);
    % Sliding Threshold
    vTimeFrame = segmentFull.*(0.1*(v1 + v2 + v3)); %masked velocity (cm/s)
    %size(vTimeFrame)
    %SegPlanes(:,:,j)=vTimeFrame;
    vTimeFramerowMean = sum(vTimeFrame,2) ./ sum(vTimeFrame~=0,2); %mean vel
    flowPulsatile_val(:,j) = vTimeFramerowMean.*area_val; %TR flow (ml/s)
    maxVelFrame(:,j) = max(vTimeFrame,[],2); %max vel. each frame (cm/s)
    velPulsatile_val(:,j) = vTimeFramerowMean;%mean vel. each frame (cm/s)  
end 
%save('SegPlanes.mat','SegPlanes','-v7.3' );
clear COL ROW idCOL Tangent_V v1 v2 v3 vx vy vz x_full y_full z_full x y z SegPlanes
%% Compute Hemodynamic Parameters
maxVel_val = max(maxVelFrame,[],2); %max in-plane veloc. for all frames
flowPerHeartCycle_val = sum(flowPulsatile_val,2)./(nframes); %TA flow (ml/s)
velMean_val = sum(velPulsatile_val,2)./(nframes); %TA in-plane velocities
% Pulsatility Index (PI) = (systolic vel - diastolic vel)/(mean vel)
PI_val = abs(max(flowPulsatile_val,[],2) - min(flowPulsatile_val,[],2))./mean(flowPulsatile_val,2);
% Resistivity Index (RI) = (systolic vel - diastolic vel)/(systolic vel)
RI_val = abs(max(flowPulsatile_val,[],2) - min(flowPulsatile_val,[],2))./max(flowPulsatile_val,[],2);
% Mean and standard deviation of flow along all branches
for i=1:max(branchList(:,4))
    idx1 = branchList(:,4)==i; %find all points along branch
    bnumMeanFlow(i) = mean(flowPerHeartCycle_val(idx1)); %mean TA flow
    bnumStdvFlow(i) = std(flowPerHeartCycle_val(idx1)); %stdv TA flow
end 
% Get coefficient of variation (stdv from mean) for all points along branch
% Looks at local stdv and mean (window width of 5).
StdvFromMean = flowPerHeartCycle_val;
for n = 1:max(branchList(:,4))
    IDbranch = find(branchList(:,4)== n); %extract points for branch n
    len=length(IDbranch);
    L_ID=ones([1 len]);
    L_ID(4:end)=L_ID(4:end)+1:(len-2);
    R_ID=3:len+2;
    R_ID((end-2):end)=len;
    for m = 1:len
        QV_meanflow_var = 1-std(flowPerHeartCycle_val(IDbranch(L_ID(m):R_ID(m))))./abs(mean(flowPerHeartCycle_val(IDbranch(L_ID(m):R_ID(m)))));
        QV_area_var = 1-(std(area_val(IDbranch(L_ID(m):R_ID(m))))./abs(mean(area_val(IDbranch(L_ID(m):R_ID(m)))))); %ideal 1, range -inf,1
        QV_circularity = mean(diam_val(IDbranch(L_ID(m):R_ID(m)))); %ideal =1, range 0,1
        for kk=1:length(flowPulsatile_val(1,:))
            flows_phase=flowPulsatile_val(IDbranch(L_ID(m):R_ID(m)),kk);
            %std_phase(kk)=std(flows_phase);
            minmax_phase(kk)=max(flows_phase)-min(flows_phase);
        end
        QV_tightness=1-mean(minmax_phase)./abs(mean(flowPerHeartCycle_val(IDbranch(L_ID(m):R_ID(m)))));
        StdvFromMean(IDbranch(m)) = QV_meanflow_var + QV_area_var + QV_circularity+QV_tightness;        
    end
end