function [area_val, diam_val, segmentFull] = segment_cross_section_thresh(segments,...
    width,timeMIPcrossection,vTimeFrameave,MAGcrossection,res,r,InterpVals,pixelSpace)
    area_val = zeros(segments,1);
    diam_val = zeros(segments,1);
    segmentFull = zeros([segments,(width).^2]);
    for n = 1:segments
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

        dArea = (res/10)*(pixelSpace(n)/10); %pixel size (cm^2)
        area_val(n,1) = sum(segment(:))*dArea*((2*r+1)/(2*r*InterpVals+1))^2;

        segmentFull(n,:) = segment(:);
        % New with ratios of areas. Ratio of smallest inner circle over
        % largest encompassing outer circle (assume circular area). Measure of
        % circularity of vessel (ratio =1 is circle,ratio<1 is irregular shape)
        D = bwdist(~segment); %euclidean distance transform, this uses an inverse segment (zeros in vessel, which then returns max distance to a corresponding vessel boundary.
        Rin = max(D(:)); %distance from center to closest non-zero entry
        [xLoc,yLoc] = find(bwperim(segment)); %get perimeter
        D = pdist2([xLoc,yLoc],[xLoc,yLoc]); %distance b/w perimeter points
        Rout = max(D(:))/2; %radius of largest outer circle
        diam_val(n,1) = Rin^2/Rout^2; %ratio of areas
        diam_val(diam_val==inf) = 0;
        %diam_val(n) = 2*sqrt(area_val(n)/pi); %equivalent diameter
    end
end