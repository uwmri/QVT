function branch = centerline_new(Cbin4CL,clData,CLsettings)
%CENTERLINE_NEW: Given 3D volume, compute centerline branches, parameterize 
%centerline branches and give branch indices by depth from main branch
%   Centerline data (clData) computed via centerlineX (Erik Spaak)
%
%   INPUT
%   1. Cbin - 3D static volume binarized
%   2. C - 3D static volume
%   3. clData - struct of centerlineX results
%   4. settings - struct containing:
%       window - 2D window size for calculating Area
%       branchMinLength - minimum branch length
%       branchMinDist - min. dist. b/w branches for adjacency consideration
%       brStartOri - starting point for branch adjacency orientation (1,-1)
%       brStartAdj - start point of adjacency analysis
%
%   OUTPUT
%   1. branch - Vector of struct for each branch
%       each branch struct contains, by parameterized index, the following:
%           branch(i).x - coordinates of centerline in x
%           branch(i).y - coordinates of centerline in y
%           branch(i).z - coordinates of centerline in z
%           branch(i).dS - incremental displacement along centerline
%           branch(i).S - total displacement along centerline from start
%
%       Gabe Shaughnessy, UW-Madison 9/23/2015
%       Carson Hoffman, UW-Madison 02/05/2019
%       Used by: feature_extraction.m
%       Dependencies: NONE


%% Setup and Initialize
volDim = size(Cbin4CL);
branchMinLength = CLsettings.branchMinLength;

if ~isfield(CLsettings,'brStartOri')
    CLsettings.brStartOri = 1;
end

branchMat = clData.branchMat;
branchTextList = clData.branchTextList;

%% Setup Branch Volume Indices
% Based on branch id, jbranch takes some time depending on tree.
Nbranch = size(branchTextList,1);
brindx = cell(Nbranch,1); %list all branches in cell array
indBrMat = find(branchMat~=0); %get all active pixels

% Cut branches shorter than branchMinLength (note 'centerlineX' cuts spurs)
nbrdel = 0;
for jbranch = 1:Nbranch
    brindx{jbranch} = indBrMat(find(branchMat(indBrMat)==jbranch));
    if length(brindx{jbranch}) < branchMinLength %if branch is too small
        nbrdel = nbrdel + 1;
        brdelList(nbrdel) = jbranch; %add to small branch list
        branchMat(brindx{jbranch})=0;
    end
end

for i = 1:nbrdel
    brindx{brdelList(i)} = []; %remove short branches from index list
end
brindx = brindx(~cellfun('isempty',brindx)); %chop entries from array
Nbranch = length(brindx);

%% Endpoint, Midpoint, and Branch Parameterization
nskip = 0;
endpt = zeros(3,2,Nbranch); %initialize endpoints/midpoints of each branch
for jbranch = 1:Nbranch 
    clear pos
    ns = length(brindx{jbranch}); %number of points in branch
    if ns <= 1  %neglect trivial branches
        nskip = nskip + 1;
        iskip(nskip) = jbranch;
        continue
    end
    [iy,ix,iz] = ind2sub(volDim,brindx{jbranch}); %index to image subscript

    ixTmp = ix; %copy indices
    iyTmp = iy;
    izTmp = iz;
    mid = brindx{jbranch}(round(ns/2)); %find middle of branch (index)
    [iymid,ixmid,izmid] = ind2sub(volDim,mid);
    iymidTmp = iymid; %copy midpoints
    ixmidTmp = ixmid;
    izmidTmp = izmid;
    midpt = [iymid,ixmid,izmid]; %middle of branch in image coordinates
    minds = 0;
    
    % Find first endpoint
    while minds < 2
        % ds = euclidean distance between successive points along branch
        ds = sqrt( (iyTmp-iymidTmp).^2 + (ixTmp-ixmidTmp).^2 + (izTmp-izmidTmp).^2 );
        if isempty(ds)
            break
        end
        ds(ds==0) = 1e9; %make true middle point gigantic
        [minds,I] = min(ds); %get minimum of distance array
        if minds >= 2
            continue
        end
        ixmidTmp = ixTmp(I); %fill next point in parameterized centerline
        iymidTmp = iyTmp(I);
        izmidTmp = izTmp(I);
        endpt(:,1,jbranch) = [iymidTmp ixmidTmp izmidTmp];
        ixTmp(I)=[]; %remove point for future consideration
        iyTmp(I)=[];
        izTmp(I)=[];
    end
    
    % Remove true middle point
    ds = sqrt( (iyTmp-iymid).^2 + (ixTmp-ixmid).^2 + (izTmp-izmid).^2 );
    [~,I] = min(ds);
    iymidTmp = iyTmp(I);
    ixmidTmp = ixTmp(I);
    izmidTmp = izTmp(I);
    iyTmp(I) = [];
    ixTmp(I) = [];
    izTmp(I) = [];
    minds = 0;
    
    % Find other endpoint
    while minds < 2
        ds = sqrt( (iyTmp-iymidTmp).^2 + (ixTmp-ixmidTmp).^2 + (izTmp-izmidTmp).^2 );
        if isempty(ds)
            break
        end
        ds(ds==0) = 1e9; %make true middle point gigantic
        [minds,I] = min(ds);
        if minds >= 2
            continue
        end
        ixmidTmp = ixTmp(I); %fill next point in parameterized centerline
        iymidTmp = iyTmp(I);
        izmidTmp = izTmp(I);
        endpt(:,2,jbranch) = [iymidTmp ixmidTmp izmidTmp];
        ixTmp(I) = []; %remove point for future consideration
        iyTmp(I) = [];
        izTmp(I) = [];
    end
    
    % Set starting point as first endpoint encountered
    ns = 1;
    if sum(endpt(:,1,jbranch),1) == 0 %if we didnt set an endpoint
        endpt(:,1,jbranch) = midpt;
    end
    if sum(endpt(:,2,jbranch),1) == 0
        endpt(:,2,jbranch) = midpt;
    end
    pos.y(ns) = endpt(1,1,jbranch);
    pos.x(ns) = endpt(2,1,jbranch);
    pos.z(ns) = endpt(3,1,jbranch);
    
    % Find distance from previous point in parameterized cl to all
    % others in cl volume with same cl index
    ds = sqrt( (iy-pos.y(ns)).^2 + (ix-pos.x(ns)).^2 + (iz-pos.z(ns)).^2 );
    [~,I] = min(ds);
    % Remove point for future consideration
    ix(I) = [];
    iy(I) = [];
    iz(I) = [];
    
    % Parameterize centerline between junctions
    while length(ix)>0
        ns = ns + 1;
        ds = sqrt( (iy-pos.y(ns-1)).^2 + (ix-pos.x(ns-1)).^2 + (iz-pos.z(ns-1)).^2 );
        [~,I] = min(ds);
        
        pos.x(ns) = ix(I); %fill next point in parameterized cl
        pos.y(ns) = iy(I); 
        pos.z(ns) = iz(I);
        ix(I) = []; %remove point for future consideration
        iy(I) = [];
        iz(I) = [];
    end
    
    % Compute centerline distances
    pos.dS = sqrt( diff(pos.x).^2 + diff(pos.y).^2 + diff(pos.z).^2 );
    pos.S = cumsum([0 pos.dS]); %running distance from 1st endpoint
    branch(jbranch) = pos; %struct w/ location, distance info of branches
end
%% Cluster Analysis
% Uncomment the following code for cluster analysis, adjacency matrix,
% direcitonal adjacency matrix, and calculation of branchMap.
% Also, add any desired parameters as function outputs.


% clbin = cl;
% clbin(cl~=0) = 1;
% 
% % Identify connected centerline regions
% clgroups = bwlabeln(clbin); %assign labels to 26-pt connected regions 
% icluster = 1; %set primary cluster id for connected centerline regions
% jbcluster = brStartAdj; %initialize cluster list
% for jbranch = 1:Nbranch
%     if sum(sum(endpt(:,:,jbranch)))==0 % neglect branches w/ <2 endpoints
%         continue
%     end
%     py = endpt(1,2,jbranch);
%     px = endpt(2,2,jbranch);
%     pz = endpt(3,2,jbranch);
%     if clgroups(py,px,pz)==icluster
%         jbcluster(end+1) = jbranch;
%     end
% end
% jbcluster = unique(jbcluster,'stable');
% NbranchOfCluster = length(jbcluster);
% 
% %% Adjacency Matrix
% % For more info visit: en.wikipedia.org/wiki/Adjacency_matrix
% adj = zeros(NbranchOfCluster,NbranchOfCluster);
% for i = 1:NbranchOfCluster-1
%     pi1 = [endpt(1,1,jbcluster(i)) endpt(2,1,jbcluster(i)) endpt(3,1,jbcluster(i))];
%     pi2 = [endpt(1,2,jbcluster(i)) endpt(2,2,jbcluster(i)) endpt(3,2,jbcluster(i))];
%     for j = i+1:NbranchOfCluster
%         pj1 = [endpt(1,1,jbcluster(j)) endpt(2,1,jbcluster(j)) endpt(3,1,jbcluster(j))];
%         pj2 = [endpt(1,2,jbcluster(j)) endpt(2,2,jbcluster(j)) endpt(3,2,jbcluster(j))];
%         if sqrt(sum((pi1-pj1).^2)) < branchMinDist 
%             adj(i,j) = 1; %if endpoints w/in min. dist., add to adj matrix
%             adj(j,i) = 1;
%         end
%         if sqrt(sum((pi2-pj1).^2)) < branchMinDist
%             adj(i,j) = 1;
%             adj(j,i) = 1;
%         end
%         if sqrt(sum((pi1-pj2).^2)) < branchMinDist
%             adj(i,j) = 1;
%             adj(j,i) = 1;
%         end
%         if sqrt(sum((pi2-pj2).^2)) < branchMinDist
%             adj(i,j) = 1;
%             adj(j,i) = 1;
%         end
%     end
% end
% 
% %% Directional Adjacency Matrix
% initvec = @(i,N) full(sparse(1,i,1,1,N))';
% visited = zeros(NbranchOfCluster,1);
% toVisitDaughter = zeros(NbranchOfCluster,1);
% oldvisited = visited;
% adjDir = zeros(NbranchOfCluster,NbranchOfCluster);
% vec = initvec(brStartAdj,NbranchOfCluster);
% toVisitParent = max(0,adj*vec-visited);
% 
% visited(brStartAdj) = 1;
% adjDir(:,brStartAdj) = toVisitParent;
% 
% while sum(abs(visited-oldvisited))>0
%     oldvisited = visited;
%     toVisitParentIndx = find(toVisitParent);
%     nToVisitParent = length(toVisitParentIndx);
%     toVisitParent = zeros(NbranchOfCluster,1);
%     for i = 1:nToVisitParent
%         ib = toVisitParentIndx(i);
%         vec = initvec(ib,NbranchOfCluster);
%         toVisitDaughter = max(0,adj*vec-visited);
%         adjDir(:,ib) = toVisitDaughter;
%         visited(ib) = 1;
%         toVisitParent(find(toVisitDaughter)) = 1;
%     end
% end
% 
% %% Produce Tree Map
% A = adjDir;
% init0 = zeros(1,NbranchOfCluster); % initialize arrays and cells
% branchMap = cell(1,NbranchOfCluster);
% visited = zeros(1,NbranchOfCluster);
% init0(brStartAdj) = 1; % set staring place
% branchMap{brStartAdj} = brStartAdj;
% oldvisited = visited;
% visited(brStartAdj) = 1;
% branch(brStartAdj).ori = brStartOri;
% 
% %% Complete Branch Map
% while sum(abs(visited-oldvisited)) > 0
%     indx0 = find(init0);
%     init0 = zeros(1,NbranchOfCluster);
%     oldvisited = visited;
%     for i = 1:length(indx0)
%         ib = indx0(i);
%         init = initvec(ib,NbranchOfCluster);
%         PendParent = [branch(ib).x([1 end]); 
%             branch(ib).y([1 end]); 
%             branch(ib).z([1 end])];
%         parentOri = branch(ib).ori;
%         if parentOri==1
%             P0 = PendParent(:,2);
%             P1 = PendParent(:,1);
%         else
%             P0 = PendParent(:,1);
%             P1 = PendParent(:,2);
%         end
%         B = A*init;
%         B(find(visited)) = 0;
%         indx = find(B);
%         parent = indx0(i);
%         for j = 1:length(indx)
%             visited(indx(j)) = 1;
%             branchMap{indx(j)} = [branchMap{parent} indx(j)];
%             init0(indx(j)) = 1;
%             ibOri = indx(j);
%             Pend = [branch(ibOri).x([1 end]); 
%                 branch(ibOri).y([1 end]); 
%                 branch(ibOri).z([1 end])];
%             AA = [sum((P0-Pend(:,1)).^2) sum((P0-Pend(:,2)).^2) sum((P1-Pend(:,1)).^2) sum((P1-Pend(:,2)).^2)];
%             [~,I] = min(AA);
%             if mod(I,2)==1
%                 branch(ibOri).ori = 1;
%             else
%                 branch(ibOri).ori = -1;
%             end
%         end
%     end
% end


end