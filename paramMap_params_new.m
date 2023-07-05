function [area_val,diam_val,flowPerHeartCycle_val,maxVel_val,PI_val,RI_val,flowPulsatile_val,...
    velMean_val,VplanesAllx,VplanesAlly,VplanesAllz,r,timeMIPcrossection,segmentFull,...
    vTimeFrameave,MAGcrossection,bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
    = paramMap_params_new(varargin)
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

%% Tangent Plane Creation
%set(handles.TextUpdate,'String','Creating Tangent Planes');drawnow;
d = 2; %dist. behind/ahead of current pt for tangent plane calc (d=2->5pts)
Tangent_V = zeros(0,3);
for n = 1:max(branchList(:,4))
    branchActual = branchList(branchList(:,4)==n,:);
    dir_temp = zeros(size(branchActual,1),3);
    for i = 1:size(branchActual,1)
        % Extract normal to cross-section
        if i < d+1 %if near 1st endpoint
            dir = (branchActual(i+d,1:3) - branchActual(i,1:3));
        elseif i >= size(branchActual,1)-d %if near 2nd endpoint 
            dir = (branchActual(i,1:3) - branchActual(i-d,1:3));
        else %calculate tangent from d points ahead/behind curr point
            dir = (branchActual(i+d,1:3) - branchActual(i-d,1:3));
        end
        dir_temp(i,:) = dir/norm(dir); %tangent vector with magnitude of 1
    end
    Tangent_V = [Tangent_V;dir_temp]; %add all tangents to large list
end

% This will find a normalized vector perpendicular to the tangent vector
[~,idx_max] = max(abs(Tangent_V),[],2); %get max unit along rows
idx_max(idx_max==2) = 1; %flatten to 2D
max_pts = sub2ind(size(Tangent_V),(1:size(Tangent_V,1))',idx_max);
temp = zeros(size(Tangent_V));
temp(max_pts) = 1; %binary matrix of location of max unit vectors
[~,idx_shift] = max(abs(circshift(temp,1,2)),[],2); %rotate (ie x->y,z->x)
shift_pts = sub2ind(size(Tangent_V),(1:size(Tangent_V,1))',idx_shift);
V2 = zeros(size(Tangent_V));
V2(max_pts) = Tangent_V(shift_pts);
V2(shift_pts) = -Tangent_V(max_pts);
N = repmat(sqrt(sum(abs(V2).^2,2)),[1 3]); %repeat vel. magnitude as Nx3
V2 = V2./N;
V3 = cross(Tangent_V,V2); %Third vector that is normalized
% V3,V2,Tangent_V are all orthogonal (i.e. dot( V3(1,:),Tangent_V(1,:) )=0)

%% Interpolate
% Get the full tangent plane for all the points
r = 10; %size of plane to select from non interpolated data is r*2+1
InterpVals = 4; %choose the interpolation between points
Side = r*InterpVals; %creates correct number of points for interpolation
width = Side.*2+1; %width of plane in pixels
Mid = zeros(length(branchList),1);

% Find x values on line
temp = repmat(V2(:,1)./InterpVals,[1 Side]);
temp = cumsum(temp,2); %runs from 0 to +(r*interpVals) by unit dist/interp
temp2 = -fliplr(temp); %runs from -(r*interpVals) to 0 by unit dist/interp
x_val = [temp2 Mid temp]; %combine temps--size = N x (r*interpVals*2)+1
x_val = bsxfun(@plus,x_val,branchList(:,1)); %pointwise addition
x_val = reshape(x_val,[numel(x_val) 1]); %stretch into vector

% Find y values on line
temp = repmat(V2(:,2)./InterpVals,[1 Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
y_val = [temp2 Mid temp];
y_val = bsxfun(@plus,y_val,branchList(:,2));
y_val = reshape(y_val,[numel(y_val) 1]);

% Find z values on the line
temp = repmat(V2(:,3)./InterpVals,[1 Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
z_val = [temp2 Mid temp];
z_val = bsxfun(@plus,z_val,branchList(:,3));
z_val = reshape(z_val,[numel(z_val) 1]);

% At this point x,y,z values have created a tangent line perpendicular to
% the normal vector for all centerline points.
% Now, we begin filling out the other perpendicular line to create a plane.

% Find x values on plane
Mid = zeros(length(branchList)*(width),1);
temp = repmat(V3(:,1)./InterpVals,[width Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
x_full = [temp2 Mid temp];
x_full = bsxfun(@plus,x_full,x_val);
x_full = reshape(x_full,[length(branchList)*(width).^2,1]);

% Find y values on plane
temp = repmat(V3(:,2)./InterpVals,[(width) Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
y_full = [temp2 Mid temp];
y_full = bsxfun(@plus,y_full,y_val);
y_full = reshape(y_full,[length(branchList)*(width).^2,1]);

% Find z values on plane
temp = repmat(V3(:,3)./InterpVals,[(width) Side]);
temp = cumsum(temp,2);
temp2 = -fliplr(temp);
z_full = [temp2 Mid temp];
z_full = bsxfun(@plus,z_full,z_val);
z_full = reshape(z_full,[length(branchList)*(width).^2,1]);

% Typecast to single and reshape
x_full = reshape(single(x_full),[length(branchList),(width).^2]);
y_full = reshape(single(y_full),[length(branchList),(width).^2]);
z_full = reshape(single(z_full),[length(branchList),(width).^2]);

% Get corners of UNINTERPOLATED planes
Planes = zeros(size(branchList,1),4,3);
Planes(:,:,1) = [x_full(:,1),x_full(:,width-InterpVals),x_full(:,end),x_full(:,end-width+1)];
Planes(:,:,2) = [y_full(:,1),y_full(:,width-InterpVals),y_full(:,end),y_full(:,end-width+1)];
Planes(:,:,3) = [z_full(:,1),z_full(:,width-InterpVals),z_full(:,end),z_full(:,end-width+1)];

dimIM = size(timeMIP);
x = 1:dimIM(1);
y = 1:dimIM(2);
z = 1:dimIM(3);

clear V2 V3 shift_pts temp temp2 idx_max idx_shift Mid x_val y_val z_val 
clear N max_pts d dimIM

% Might use to speed up interpolation time if needed for large matrix sizes
%SE = ones(10,10,10);
%CD_bin_new = imdilate(CD_bin,SE);

%% Interpolation
%set(handles.TextUpdate,'String','Interpolating Data');drawnow;
% Get interpolated velocity from 3 directions, multipley w/ tangent vector
v1 = interp3(y,x,z,vMean(:,:,:,1),y_full(:),x_full(:),z_full(:),'linear',0);
v2 = interp3(y,x,z,vMean(:,:,:,2),y_full(:),x_full(:),z_full(:),'linear',0);
v3 = interp3(y,x,z,vMean(:,:,:,3),y_full(:),x_full(:),z_full(:),'linear',0);
v1 = reshape(v1,[length(branchList),(width).^2]);
v2 = reshape(v2,[length(branchList),(width).^2]);
v3 = reshape(v3,[length(branchList),(width).^2]);
temp = zeros([size(v1),3]); % used to hold velocity data information
temp(:,:,1) = bsxfun(@times,v1,Tangent_V(:,1)); %dot product here
temp(:,:,2) = bsxfun(@times,v2,Tangent_V(:,2)); %make veloc. through-plane
temp(:,:,3) = bsxfun(@times,v3,Tangent_V(:,3)); %(mm/s)

% Through-plane SPEED for all points (tangent vector dotted with 3D vel)
vTimeFrameave = sqrt(temp(:,:,1).^2 + temp(:,:,2).^2 + temp(:,:,3).^2); %(mm/s)

% Interpolation for complex difference data
CD_int = interp3(y,x,z,timeMIP,y_full(:),x_full(:),z_full(:),'linear',0);
timeMIPcrossection = reshape(CD_int,[length(branchList),(width).^2]);

% Interpolation for magnitude data
Mag_int = interp3(y,x,z,MAG,y_full(:),x_full(:),z_full(:),'linear',0);
MAGcrossection = reshape(Mag_int,[length(branchList),(width).^2]);

clear v1 v2 v3 MAG timeMIP temp CD_int Mag_int vtimeave 

%% In-Plane Segmentation
%set(handles.TextUpdate,'String','Performing In-Plane Segmentation');drawnow;
area_val = zeros(size(Tangent_V,1),1);
diam_val = zeros(size(Tangent_V,1),1);
segmentFull = zeros([length(branchList),(width).^2]);

for n = 1:size(Tangent_V,1)
    %%%%%% SLIDING THRESHOLD %%%%%%
    % Get Planes and normalize
    cdSLICE = reshape(timeMIPcrossection(n,:),[(width),(width)]);
    temp = cdSLICE - min(cdSLICE); %shift the minimum to 0
    cdSLICE = temp./max(temp(:)); %now normalize from 0 to 1
    
    velSLICE = reshape(vTimeFrameave(n,:),[(width),(width)]);
    temp = velSLICE - min(velSLICE);
    velSLICE = temp./max(temp(:));
    
    magSLICE = reshape(MAGcrossection(n,:),[(width),(width)]);
    temp = magSLICE - min(magSLICE);
    magSLICE = temp./max(temp(:));
    
    weightIMS = [.2 .8 .2]; % Weights = [Mag CD Vel]
    weightIMAGE = (weightIMS(1).*magSLICE) + (weightIMS(2).*cdSLICE) + (weightIMS(3).*velSLICE);
    
    step = 0.001;
    UPthresh = 0.8;
    SMf = 90; %smoothing factor
    shiftHM_flag = 0; %do not shift by FWHM
    medFilt_flag = 1; %flag for median filtering of CD image
    [~,segment] = slidingThreshold(weightIMAGE,step,UPthresh,SMf,shiftHM_flag,medFilt_flag);
    areaThresh = round(sum(segment(:)).*0.05); %minimum area to keep
    conn = 6; %connectivity (i.e. 6-pt)
    segment = bwareaopen(segment,areaThresh,conn); %inverse fill holes
    % Can compare in-plane segmentation to initial global segmentation. 
    % To do this, the 'segment' variable from 'loadpcvipr' needs to be 
    % passed as an arg. I did this by adding 'segment_old' as 2nd input
    %segment_old = interp3(y,x,z,single(segment_old),y_full(:),x_full(:),z_full(:),'linear',0);
    %segment_old = reshape(segment_old,[length(branchList),(width).^2]);
    %segSlice = reshape(segment_old(n,:),[(width),(width)]);
    %figure; imshow(imbinarize(segSlice));
    %figure; imshow(segment);
    
    % Remove all segments not closest to the center
    s = regionprops(segment,'centroid'); %centroids of unique lbls  
    CenterIm = [size(segment,1)/2,size(segment,2)/2]; %loc image center
    Centroids = reshape([s(:).Centroid],[2,length([s(:).Centroid])/2])';
    DisCen = sqrt(sum((Centroids - repmat(CenterIm,[size(Centroids,1),1])).^2,2));
    [~,CenIdx]  = min(DisCen); %find centroid closest to center

    % Fill in the holes and clean up
    [L,Num] = bwlabel(segment); %find centroid index
    LabUse = 1:Num;
    segment = L==LabUse(CenIdx); %cut out other centroids
    %segment = imopen(segment,ones(3,3)); %morphological opening
    
    % Vessel area measurements
    dArea = (res/10).^2; %pixel size (cm^2)
    area_val(n,1) = sum(segment(:))*dArea*((2*r+1)/(2*r*InterpVals+1))^2;
   
    
    segmentFull(n,:) = segment(:);
    
    % New with ratios of areas. Ratio of smallest inner circle over
    % largest encompassing outer circle (assume circular area). Measure of
    % circularity of vessel (ratio =1 is circle,ratio<1 is irregular shape)
    D = bwdist(~segment); %euclidean distance transform
    Rin = max(D(:)); %distance from center to closest non-zero entry
    [xLoc,yLoc] = find(bwperim(segment)); %get perimeter
    D = pdist2([xLoc,yLoc],[xLoc,yLoc]); %distance b/w perimeter points
    Rout = max(D(:))/2; %radius of largest outer circle
    diam_val(n,1) = Rin^2/Rout^2; %ratio of areas
    diam_val(diam_val==inf) = 0;
    %diam_val(n) = 2*sqrt(area_val(n)/pi); %equivalent diameter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
clear cdSLICE magSLICE temp segment weightIMAGE dArea L LabUse CenIdx Num

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
        %set(handles.TextUpdate,'String',['Calculating Quantitative Params Frame: ' num2str(j) '/' num2str(nframes)]);drawnow;

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
        %set(handles.TextUpdate,'String',['Calculating Quantitative - Parameters Time Frame: ' num2str(j) '/' num2str(nframes)]);drawnow;
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
        vx = squeeze(v(:,:,:,1,j));
        vy = squeeze(v(:,:,:,2,j));
        vz = squeeze(v(:,:,:,3,j));
    else
        %set(handles.TextUpdate,'String',['Calculating Quantitative (python H5 export) - Parameters Time Frame: ' num2str(j) '/' num2str(nframes)]);drawnow;
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
    v1 = interp3(y,x,z,vx,y_full(:),x_full(:),z_full(:),'linear',0);
    v2 = interp3(y,x,z,vy,y_full(:),x_full(:),z_full(:),'linear',0);
    v3 = interp3(y,x,z,vz,y_full(:),x_full(:),z_full(:),'linear',0);
    v1 = reshape(v1,[length(branchList),(width).^2]);
    v2 = reshape(v2,[length(branchList),(width).^2]);
    v3 = reshape(v3,[length(branchList),(width).^2]);
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
    % Calculate near branch start
    StdvFromMean(IDbranch(1)) = std(flowPerHeartCycle_val(IDbranch(1:3))) ./ abs(mean(flowPerHeartCycle_val(IDbranch(1:3))));
    StdvFromMean(IDbranch(2)) = std(flowPerHeartCycle_val(IDbranch(1:4))) ./ abs(mean(flowPerHeartCycle_val(IDbranch(1:4))));
    % Calculate for middle of branch (window width of 5)
    for m = 1:numel(IDbranch)-4
        StdvFromMean(IDbranch(m+2)) = std(flowPerHeartCycle_val(IDbranch(m:m+4)))./abs(mean(flowPerHeartCycle_val(IDbranch(m:m+4))));
    end
    % Calculate near branch end
    StdvFromMean(IDbranch(end-1)) = std(flowPerHeartCycle_val(IDbranch(end-3:end)))./abs(mean(flowPerHeartCycle_val(IDbranch(end-3:end))));
    StdvFromMean(IDbranch(end)) = std(flowPerHeartCycle_val(IDbranch(end-2:end)))./abs(mean(flowPerHeartCycle_val(IDbranch(end-2:end))));
end
StdvFromMean = StdvFromMean - min(StdvFromMean(:)); %shift the minimum to 0
StdvFromMean = StdvFromMean./max(StdvFromMean(:)); %normalize range 0-1
end
