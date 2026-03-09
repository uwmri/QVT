function [branchList,branchJunctions,junctionList] = feature_extraction( ...
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
% RVC starting December 2025
imgsize = size(cl); % needed ?
nbrnch = max(branchList(:,4));
branchJunctions = cell(nbrnch, 1);
% preallocate array using dummy structure
njunc = max(junctionList(:,4)); % yes this is used more than once
% junctionList: first column is arrays of branch IDs for the junction nx2,
% second column is 1x3 arrays representing junction position
jtmp.idbranches = [];
jtmp.brpos = zeros(0,3);
jtmp.rjunc = zeros(1,3);
jtmp.distbr = zeros(0,1);
jtmp.dseedsigned = 0;
jtmp.dseedabsolute = 0;
junctionList = repmat(jtmp,njunc,1); 

for k=1:nbrnch
    thisbrlist = branchList(branchList(:,4)==k,:);
    ri = thisbrlist(1,1:3);
    d = size(thisbrlist);
    rf = thisbrlist(d(1),1:3);
    idbri = findJunctions(ri, junctionMat, imgsize); % ids at start/source of branch
    idbrf = findJunctions(rf, junctionMat, imgsize); % ids at branch sink
    if isempty(idbri) && isempty(idbrf)
        fprintf('debugging: no junctions for branch %d\n',k)
        continue
    end
    if any(ismember(idbri, idbrf))
        fprintf('WARNING: branch %d appears to be a loop, discarding\n',k);
        continue
    end

    nji = length(idbri);
    njf = length(idbrf);
    idxlisti = 1:nji;
    idxlistf = nji + (1:njf);
    branchJunctions{k} = zeros(nji+njf,2);
    branchJunctions{k}(idxlisti,1) = idbri;
    branchJunctions{k}(idxlistf,1) = idbrf;
    branchJunctions{k}(idxlisti,2) = -1;
    branchJunctions{k}(idxlistf,2) = 1;
    for kk = 1:nji
        junctionList(idbri(kk)).idbranches(end+1) = k;
        junctionList(idbri(kk)).brpos(end+1,:) = ri;
    end
    for kk = 1:njf
        junctionList(idbrf(kk)).idbranches(end+1) = k;
        junctionList(idbrf(kk)).brpos(end+1,:) = rf;
    end
end

% compute position of each junction
for k=1:njunc
    nbr = length(junctionList(k).idbranches);
    if nbr < 1
        continue
    elseif nbr == 1
        junctionList(k).rjunc = junctionList(k).brpos(1,:);
    else
        junctionList(k).rjunc = mean(junctionList(k).brpos);
    end
    for kk=1:nbr
        brjuncmat = branchJunctions{junctionList(k).idbranches(kk)};
        indx = find(brjuncmat(:,1) == k);
        if length(indx) ~= 1
            disp('PROGRAMMING ERROR: junction assigned to branch multiple times')
            return
        end
        dvec = junctionList(k).rjunc - junctionList(k).brpos(kk,:);
        % distbr should be positive for flow into the junction
        junctionList(k).distbr(kk) = sqrt(dvec*dvec')*brjuncmat(indx,2);
    end
end

end

function idjuncs = findJunctions(r, jMat, imgsz)
% return the ids of the junctions near r

xrng = limitnbrs(r(1), imgsz(1));
yrng = limitnbrs(r(2), imgsz(2));
zrng = limitnbrs(r(3), imgsz(3));
[ rvec, idxvec, jlist ] = find(jMat(xrng,yrng,zrng));
idjuncs = unique(jlist);
if isempty(idjuncs)
    return
end

roffset = [ r(1)-xrng(1)+1, r(2)-yrng(1)+1, r(3)-zrng(1)+1 ];
% idxvec is linear indices so we need to convert to c(olumn)vec and
% p(age)vec
[ cvec, pvec ] = ind2sub([ length(yrng), length(zrng) ], idxvec);
mindists = 20*ones(size(idjuncs));
for k = 1:length(jlist)
    indxthisj = find(idjuncs == jlist(k));
    rdiff = roffset - [ rvec(k), cvec(k), pvec(k) ];
    d = sqrt(rdiff*rdiff');
    if d < mindists(indxthisj)
        mindists(indxthisj) = d;
    end
end

idfinal = [];
for k = 1:length(mindists)
    if mindists(k) < 2.0001 % this, or <= 2 ?? to go higher, limitnbrs() should return longer arrays
        idfinal(end+1) = idjuncs(k); %#ok<AGROW>
    end
end
idjuncs = idfinal;

end

function idxarray = limitnbrs(x, lim)
ii = ceil(x-2);
if ii < 1
    ii = 1;
end
jj = floor(x+2);
if jj > lim
    jj = lim;
end
idxarray = ii:jj;
end
