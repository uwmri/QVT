function [V,dcminfo] = shuffleDCM(path,PWD,flip)
cd(path)
DIR=dir(path);
Exp={'i\d*.MRDC.(\d*)$';'.*-(\d*).dcm$'};
expID=0;
flag=0;
for i=1:length(DIR)
    Name=DIR(i).name;
    if expID==0
        for j=1:length(Exp)
            [match] = regexp(Name,Exp{j},'tokens');
            if length(match) == 1
                expID=j;
                flag=1;
            end
        end
    end
    if flag==1
        [match] = regexp(Name,Exp{expID},'tokens');
        if length(match) == 1
            idx=str2double(match{1});
            SortedDir{idx,1}=Name;
        end
    end
end

filename=SortedDir{1};
slice = dicomread(filename);
[a,b]=size(slice);
dcminfo = dicominfo(filename);
numphase=dcminfo.CardiacNumberOfImages;
slices=length(SortedDir)/numphase;
V=zeros([a,b,slices,numphase],'single');
if flip==0
    pos=0;
else
    pos=slices+1;
end
for i=1:length(SortedDir)
    filename=SortedDir{i};
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
