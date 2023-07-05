function [area_val, diam_val, segmentFull] = segment_cross_section_Jiantao(segments,...
    width,TOFsynthcross,res,r,InterpVals)
    area_val = zeros(segments,1);
    diam_val = zeros(segments,1);
    segmentFull = zeros([segments,(width).^2]);
    segmentFull(:,1) = 1;
    for n = 1:segments
        %%%%%% SLIDING THRESHOLD %%%%%%
        % Get Planes and normalize
        segment = reshape(TOFsynthcross(n,:),[(width),(width)]);
        segment = imbinarize(segment);
    
        areaThresh = round(sum(segment(:)).*0.05); %minimum area to keep
        conn = 6; %connectivity (i.e. 6-pt)
        segment = bwareaopen(segment,areaThresh,conn); %inverse fill holes
    
        % Remove all segments not closest to the center
        s = regionprops(segment,'centroid'); %centroids of unique lbls  
        CenterIm = [size(segment,1)/2,size(segment,2)/2]; %loc image center
        Centroids = reshape([s(:).Centroid],[2,length([s(:).Centroid])/2])';
        DisCen = sqrt(sum((Centroids - repmat(CenterIm,[size(Centroids,1),1])).^2,2));
        [~,CenIdx]  = min(DisCen); %find centroid closest to center
    
        % Fill in the holes and clean up
        [L,Num] = bwlabel(segment); %find centroid index
        LabUse = 1:Num;
        try
        segment = L==LabUse(CenIdx); %cut out other centroids
        segment = imopen(segment,ones(3,3)); %morphological opening
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
        catch
            area_val(n,1)=0;
            diam_val(n,1)=0;
        end
    end
end