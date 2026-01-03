function [cl,branchMat,branchList,branchTextList,jListStruct] = feature_extraction( ...
    sortingCriteria,spurLength,vMean,segment,handles)
%FEATURE_EXTRACTION: Create vessel centerlines and label branches
%   Used by: loadpcvipr.m
%   Dependencies: centerlineX.m, centerline_new.m

set(handles.TextUpdate,'String','Completing Centerline Extraction and Labeling'); drawnow;

%% Skeletonization - Vascular Tree Construction
% New Skeleton built in functions (completes skeleton and trimming)
SkelBin = bwskel(logical(segment),'MinBranchLength',spurLength);
zeroEdger = padarray(ones(size(SkelBin)-2,'logical'),[1 1 1] ,0);
SkelBin = SkelBin.*zeroEdger; %make edges 0 (also turns SkelBin to double)

% specify sortingCriteria as either
% = 2 to get all branches connected to each other (few branches,no junctions)
% = 3 to get branch by branch sorting (many branches)
[ cl,branchMat,junctionMat,branchTextList,junctionList,~ ] = centerlineX(SkelBin, 1, sortingCriteria);
Cbin4CL = imbinarize(zeroEdger.*segment);

% Prepare clData structure and settings for centerline_new
clData.branchMat = branchMat;
clData.branchTextList = branchTextList;
settings.cl = struct();
settings.cl.branchMinLength = 5; %this is trimming of junctions
CLsettings = settings.cl;
branch = centerline_new(Cbin4CL,clData,CLsettings);

%% Branch List Sorting
% sort the branchList so that
% 1. all the same labels are connected along the rows
% 2. low -> high row index is in the same direction as the flow
branchListSorted = zeros(0,5);
segmentCutoff = 8;

for nbr = 1:length(branch)
    % Find branch
    branchActual = [branch(nbr).y',branch(nbr).x',branch(nbr).z',ones(numel(branch(nbr).x),1).*nbr];

    % Check if a -> b is in the direction of the flow...
    if size(branchActual,1) < segmentCutoff
        v0x = 0; v0y = 0; v0z = 0;
        for j = 1:size(branchActual, 1)
            v0x = v0x + vMean(branchActual(j,1), branchActual(j,2), branchActual(j,3), 1);
            v0y = v0y + vMean(branchActual(j,1), branchActual(j,2), branchActual(j,3), 2);
            v0z = v0z + vMean(branchActual(j,1), branchActual(j,2), branchActual(j,3), 3);
        end
        isReverse = dot(double(branchActual(end, 1:3) - branchActual(1, 1:3)), double([v0x v0y v0z]));
    else %cutoff segment after some amount of points (just to analyze flow)
        v0x = 0; v0y = 0; v0z = 0;
        for j = 1:segmentCutoff %iteratively add velocities along segment
            v0x = v0x + vMean(branchActual(j,1), branchActual(j,2), branchActual(j,3), 1);
            v0y = v0y + vMean(branchActual(j,1), branchActual(j,2), branchActual(j,3), 2);
            v0z = v0z + vMean(branchActual(j,1), branchActual(j,2), branchActual(j,3), 3);
        end
        % If velocity and xyz coords run in same direction, isReverse>0
        isReverse = dot(double(branchActual(segmentCutoff, 1:3) - branchActual(1, 1:3)), double([v0x v0y v0z]));
    end

    if isReverse < 0 % ...if not, reverse segment indices of xyz locations
        branchActual = flipud(branchActual);
    end
    
    branchActual = [branchActual,(1:numel(branch(nbr).x))'];
    branchListSorted = [branchListSorted; branchActual];
end
branchList = branchListSorted;

%% Centerline smoothing
% Smooths labeled centerline w/ splenic spline fit
branchListSmooth = ones([size(branchList,1),size(branchList,2)]);   
smoothParameter = 0.3750; %user-defined degree of smoothing
for n = 1:max(branchList(:,4))
    branchActual = branchList(branchList(:,4)==n,:); %branch locations(xyz)
    xyz = [branchActual(:,1)';branchActual(:,2)';branchActual(:,3)'];
    [ndim,npts] = size(xyz);
    xyzp = zeros(size(xyz)); %initialize spline xyz matrix
    
    % Cubic spline smoothing (see function details)
    % Default smoothParameter is = 1/(1 + spacing^3/6) = 0.8571
    for k=1:ndim %for each dimension (xyz)
        pp = csaps(1:npts,xyz(k,:),smoothParameter);
        xyzp(k,:)=ppval(pp,1:npts); %apply spline fit params to new matrix
    end
    branchListSmooth(branchList(:,4)==n,1:3) = xyzp'; %reassign xyz locs
    branchListSmooth(branchList(:,4)==n,4:5) = branchList(branchList(:,4)==n,4:5);
end
branchList = branchListSmooth;

%% match junctions to branches
% RVC December 2025
imgsize = size(cl);
% branchJunctions = - ones(2, max(branchList(:,4))); % are we using this?
njunc = max(junctionList(:,4)); % is this used more than once?

% preallocate array using dummy structure
junctemp.pos = zeros(1,3);
junctemp.branches = [];
jListStruct = repmat(junctemp, 1, njunc);

% for k=1:max(branchList(:,4))
%     br = branchList(branchList(:,4)==k,:);
% end
% 
% for k=1:njunc
%     thisjunction = junctionList(junctionList(:,4)==k,1:3);
%     if isempty(thisjunction)
%         continue
%     end
% 
%     % define a box around the junction
%     d = size(thisjunction);
%     if d(1) > 1
%         boxmin = min(thisjunction) - 1; % lower edge is one less than minimum value in each column
%         boxmax = max(thisjunction) + 1; % upper edge
%     else
%         boxmin = thisjunction - 1;
%         boxmax = thisjunction + 1;
%     end
%     boxmin = max(boxmin,5);         % at least five voxels from edge of image (may not be necessary)
%     boxmax = min(boxmax,imgsize - 5);
%     naybrhood = cl(boxmin(1):boxmax(1),boxmin(2):boxmax(2),boxmin(3):boxmax(3));
%     [ rowvec, colvec, pgvec ] = ind2sub(size(naybrhood),find(naybrhood==2)); % branch "coordinates"
% 
%     % look at all of the branch points around this junction
%     for idxnbr=1:length(rowvec)
%         lbptmatch(:,3) = abs(boxmin(3) + pgvec(idxnbr) - 1 - branchList(:,3)) < 1;
%         lbptmatch(:,2) = abs(boxmin(2) + colvec(idxnbr) - 1 - branchList(:,2)) < 1;
%         lbptmatch(:,1) = abs(boxmin(1) + rowvec(idxnbr) - 1 - branchList(:,1)) < 1;
% 
%         idxbrpts = find(all(lbptmatch,2));
%         if isempty(idxbrpts) % maybe the branch was trimmed away?
% 	        continue
%         end
% 
%         fprintf('branch %d attaches to junction %d:\n', branchList(idxbrpts(1),4), k)
%         disp(branchList(idxbrpts,:))
%     end
% end

end

