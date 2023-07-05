function [crossSections] = interp_vol_to_planes(Vol,x,y,z,x_full,y_full,z_full,width,segments)
    CD_int = interp3(y,x,z,Vol,y_full(:),x_full(:),z_full(:),'linear',0);
    crossSections = reshape(CD_int,[segments,(width).^2]);
end