function [V,dcminfo] = shuffleDCM(path,PWD,flip)
cd(path)
DIR=dir;
DIR(1:2)=[];
filename=DIR(1).name;
slice = dicomread(filename);
[a,b]=size(slice);
dcminfo = dicominfo(filename);
numphase=dcminfo.CardiacNumberOfImages;
slices=length(DIR)/numphase;
V=zeros([a,b,slices,numphase],'single');
if flip==0
    pos=0;
else
    pos=slices+1;
end
for i=1:length(DIR)
    filename=DIR(i).name;
    slice = dicomread(filename);
    phase=rem(i,numphase);
    if phase == 0
        phase = numphase;
    end
    if phase ==1
        if flip==0
            pos = pos+1;
        else
            pos = pos-1;
        end
    end
    V(:,:,pos,phase)=slice(:,:);
end
cd(PWD)
