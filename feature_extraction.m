function [branchList,branchJunctions,jListStruct] = feature_extraction( ...
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
imgsize = size(cl); % needed ?
nbrnch = max(branchList(:,4));
branchJunctions = - ones(nbrnch, 2); % are we using this?
% preallocate array using dummy structure
njunc = max(junctionList(:,4)); % is this used more than once?
junctemp.pos = zeros(1,3);
junctemp.idbrs = [];
junctemp.brposmat = [];
jListStruct = repmat(junctemp, 1, njunc);

for k=1:nbrnch
    thisbrlist = branchList(branchList(:,4)==k,:);
    ri = thisbrlist(1,1:3);
    d = size(thisbrlist);
    rf = thisbrlist(d(1),1:3);
    idbri = findJunctions(ri, junctionMat, imgsize);
    idbrf = findJunctions(rf, junctionMat, imgsize);
    if idbri < 0 && idbrf < 0
        continue
    end
    if idbri == idbrf
        fprintf('WARNING: branch %d appears to be a loop, discarding\n',k);
        continue
    end
    if idbri > 0
        branchJunctions(k,1) = idbri;
        n = 1 + length(jListStruct(idbri).idbrs);
        jListStruct(idbri).idbrs(n) = k;
        jListStruct(idbri).brposmat(n,:) = ri;
    end
    if idbrf > 0
        branchJunctions(k,2) = idbrf;
        n = 1 + length(jListStruct(idbrf).idbrs); % ugh, repitition
        jListStruct(idbrf).idbrs(n) = k;
        jListStruct(idbrf).brposmat(n,:) = rf;
    end
end

for k=1:njunc
    n = length(jListStruct(k).idbrs);
    if n > 1
        jListStruct(k).pos = mean(jListStruct(k).brposmat);
    elseif n == 1
        jListStruct(k).pos = jListStruct(k).brposmat;
    end
end

end

function idbr = findJunctions(r, jMat, imgsz)
% return the id of the junction nearest to r

idbr = -1;
for ii=1:3
    if r(ii) < 2 || r(ii) + 2 > imgsz(ii)
        return % don't look at edges of image
    end
end
submat = jMat(juncsearch(r(1)), juncsearch(r(2)), juncsearch(r(3)));
[ ~, ~, brlist ] = find(submat);
if isempty(brlist)
    return
end
ulist = unique(brlist);
if length(ulist) == 1
    idbr = ulist(1);
    return
end
disp('found multiple junction candidates, end coordinate is')
disp(r)
disp('candidate junctions are')
disp(ulist)

end

function idxarray = juncsearch(x)
    idxarray = ceil(x-2):floor(x+2);
end