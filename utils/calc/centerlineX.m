function [CL, branchMat,junctionMat,branchTextList,junctionList,branchList] = centerlineX(Y, spurLength, sortingCriteria)
%CENTERLINEX: Produces centerline by labeling skeleton points as either:
% 1. end points = 1
% 2. middle/branch points = 2
% 3. junction points = 3 (if sortingCriteria=3)
%
%   The function initially labels all points according to list above, then 
%   iteratively removes short spurs until no changes  occur in skeleton.

%   INPUT
%   1. Y - skeleton in binary form
%   2. spurLength - length of vessel spurs to be removed
%   3. sortingCriteria - labels(=3) or doesnt label(=2) junctions
%
%   OUTPUT
%   1. CL - centerline with spurs removed and points classified
%   2. branchMat - branch indices & labels in matrix form
%   3. junctionMat - junction indices & labels in matrix form
%   4. branchTextList - accompaning text (number labels)
%   5. junctionList - junction indices & labels in list form
%
%       Erik Spaak, Umea University 2014
%       Used by: feature_extraction.m
%       Dependencies: NONE


%% Initial Variables
dim = size(Y);
modified = 1;
Niter = 0;
CL = 2*Y;

%% Big While Loop - Cut branches
while modified > 0  && Niter < 20 %do until convergence
    Niter = Niter + 1;
    
    % Deletion of branches
    if Niter > 1 %do after first iteration   
        modified = 0;
        uniqueBranchLabels = unique(branchList(:,4));
        for i = 1:length(uniqueBranchLabels)
            currentBranchLabel = uniqueBranchLabels(i);
            currentBranchIndices = find(branchList(:,4) == currentBranchLabel);
            currentBranchLength = length(currentBranchIndices);
            connectedToJunctions = 0;
            
            for j = currentBranchIndices'
                x0 = branchList(j,1); y0 = branchList(j,2); z0 = branchList(j,3);
                connectedToJunctions = [connectedToJunctions; unique(junctionMat(x0-1:x0+1, y0-1:y0+1, z0-1:z0+1))];
            end

            % Delete branch if too short and not b/w 2 separate junctions
            % OR deletes sections based only on length of the segmentations
            
            % This statement can remove all segments that are shorter than
            % desired pixel length by commenting the && section below
            if (currentBranchLength < spurLength) && (length(unique(connectedToJunctions)) < 3) %LOOK HERE
                for j = currentBranchIndices'
                    CL(branchList(j,1), branchList(j,2), branchList(j,3)) = 0;
                end
                modified = modified + 1;
            end
        end
    end

    %% Classify skeleton points
    CL = 2*logical(CL);
    CLindices = find(CL);
    for i = 1:length(CLindices)
        [x0, y0, z0] = ind2sub(dim, CLindices(i));
        % 26-neighborhood sum (cube around point of interest)
        neighSum = sum( logical(CL(x0-1:x0+1, y0-1:y0+1, z0-1:z0+1)), 'all');
        if neighSum > 3
            CL(CLindices(i)) = sortingCriteria; %mark as junction point
        end
    end

    %% Search for and label junctions
    if sortingCriteria==3
        junctionIndices = find(CL == 3); %find junction points (=3)
        [x0, y0, z0] = ind2sub(dim,junctionIndices);
        junctionMat = zeros(dim); %label matrix
        junctionList = [x0 y0 z0 zeros(length(x0),1)]; %label vector

        % Assign specific numeric label to each junction
        label = 0;
        for i = 1:length(x0)
            if junctionList(i,4) == 0
                label = label + 1;
                junctionMat(x0(i), y0(i), z0(i)) = label;
                junctionList(i,4) = label;

                investigatePointsList = [x0(i) y0(i) z0(i)];
                labeled = 1;
                while labeled > 0 %while still getting points w/ this label
                    labeled = 0;
                    newInvestigativePointsList = [];
                    for j = 1:length(investigatePointsList(:,1))
                        x1 = investigatePointsList(j,1);
                        y1 = investigatePointsList(j,2);
                        z1 = investigatePointsList(j,3);

                        % Collect 26-point neighborhoods
                        label26 = junctionMat(x1-1:x1+1, y1-1:y1+1, z1-1:z1+1);
                        antiLabel26 = imcomplement(imbinarize(label26));
                        CL26 = CL(x1-1:x1+1, y1-1:y1+1, z1-1:z1+1);
                        
                        % Find neighboring middle points not labeled
                        neigh = find(CL26.*antiLabel26 == 3);   
                        [x2, y2, z2] = ind2sub([3 3 3], neigh);
                        x3 = x1 + x2 - 2; %convert to a CL index
                        y3 = y1 + y2 - 2;
                        z3 = z1 + z2 - 2;

                        for k = 1:length(x3)
                            junctionMat(x3(k), y3(k), z3(k)) = label;
                            % Find the neighboring points in branchList
                            a = find(junctionList(:,1) == x3(k));
                            b = find(junctionList(:,2) == y3(k));
                            c = find(junctionList(:,3) == z3(k));
                            d = intersect(a,b); %finds common values
                            e = intersect(c,d); %locate idx in branchList

                            junctionList(e,4) = label;
                            labeled = labeled + 1; %count points collected
                        end
                        newInvestigativePointsList = [newInvestigativePointsList; x3 y3 z3];
                    end
                    investigatePointsList = newInvestigativePointsList;
                end
            end
        end
    end
    
    %% Search for and label middle points
    branchIndices = find(CL == 2);
    [x0, y0, z0] = ind2sub(dim, branchIndices);
    branchMat = zeros(dim); %label matrix
    branchList = [x0 y0 z0 zeros(length(x0), 2)]; %label vector

    % Label branches
    branchTextList = zeros(0,4);
    label = 0;
    for i = 1:length(x0)
        if branchList(i,4) == 0     
            label = label + 1;
            branchMat(x0(i), y0(i), z0(i)) = label;
            branchList(i,4) = label; 
            branchTextList = [branchTextList; x0(i) y0(i) z0(i) label]; % create a textlist 
            investigatePointsList = [x0(i) y0(i) z0(i)];
            labeled = 1;
            incrementer = 0;    
            while labeled > 0 %while still getting points under this label
                labeled = 0;
                newInvestigativePointsList = [];
                for j = 1:length(investigatePointsList(:,1))
                    x1 = investigatePointsList(j,1);
                    y1 = investigatePointsList(j,2);
                    z1 = investigatePointsList(j,3);

                    % Collect 26-neighborhoods
                    label26 = branchMat(x1-1:x1+1, y1-1:y1+1, z1-1:z1+1);
                    antiLabel26 = imcomplement(imbinarize(label26));
                    CL26 = CL(x1-1:x1+1, y1-1:y1+1, z1-1:z1+1);
                    
                    % Find neighboring branch points
                    neigh = find(CL26.*antiLabel26 == 2); 
                    [x2, y2, z2] = ind2sub([3 3 3], neigh);
                    x3 = x1 + x2 - 2;
                    y3 = y1 + y2 - 2;
                    z3 = z1 + z2 - 2;
                    incrementer = incrementer + 1;
                    
                    for k = 1:length(x3)
                        branchMat(x3(k), y3(k), z3(k)) = label;
                        % Find the neighboring points in branchList
                        a = find(branchList(:,1) == x3(k));
                        b = find(branchList(:,2) == y3(k));
                        c = find(branchList(:,3) == z3(k));
                        d = intersect(a,b); %finds common values
                        e = intersect(c,d); %locate idx in branchList

                        branchList(e,4) = label;
                        
                        % Sorting (find initial idx)
                        aa = find(branchList(:,1) == x1);
                        bb = find(branchList(:,2) == y1);
                        cc = find(branchList(:,3) == z1);
                        dd = intersect(aa,bb);
                        ee = intersect(cc,dd);
                        
                        % Start counting up along branch
                        if branchList(ee,5) == 0 && k == 1 %start count
                            branchList(e,5) = 1;                
                        elseif branchList(ee,5) == 0 && k == 2 %count back   
                            branchList(e,5) = -1;   
                        elseif branchList(ee,5) > 0 %keep counting                     
                            branchList(e,5) = branchList(ee,5) + 1;                         
                        elseif branchList(ee,5) < 0 %keep counting back                         
                            branchList(e,5) = branchList(ee,5) - 1;
                        end
                        labeled = labeled + 1; %count points collected 
                    end
                    newInvestigativePointsList = [newInvestigativePointsList; x3 y3 z3];
                end
                investigatePointsList = newInvestigativePointsList;
            end
        end 
    end
end