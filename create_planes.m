function [x_full,y_full,z_full,x,y,z,Tangent_V,Planes] = create_planes(branchList,r,matrix,InterpVals,width)
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
    %
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
    % quiver3(branchList(:,1),branchList(:,2),branchList(:,3),Tangent_V(:,1),Tangent_V(:,2),Tangent_V(:,3),0.5,'k.');
    % hold on
    % quiver3(branchList(:,1),branchList(:,2),branchList(:,3),V2(:,1),V2(:,2),V2(:,3),0.5,'r.');
    % hold on
    % quiver3(branchList(:,1),branchList(:,2),branchList(:,3),V3(:,1),V3(:,2),V3(:,3),0.5,'b.');
    % hold on
    
    % Interpolate
    % Get the full tangent plane for all the points
    %r = 10; %size of plane to select from non interpolated data is r*2+1. 10 is 10 pixel values assumed isotropic distance from the center.
    %InterpVals = 4; %choose the interpolation between pixel values
    Side = r*InterpVals; %creates correct number of points for interpolation
    %width = Side.*2+1; %width of plane in pixels
    Mid = zeros(length(branchList),1);
    % Find x values on line
    temp = repmat(V2(:,1)./InterpVals,[1 Side]);
    temp = cumsum(temp,2); %runs from 0 to +(r*interpVals) by unit dist/interp
    temp2 = -fliplr(temp); %runs from -(r*interpVals) to 0 by unit dist/interp
    x_val = [temp2 Mid temp]; %combine temps--size = N x (r*interpVals*2)+1 so x ranges from -10 pixels to +10 pixels of magnitude range.
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
    %
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
    x = 1:matrix(1);
    y = 1:matrix(2);
    z = 1:matrix(3);
end