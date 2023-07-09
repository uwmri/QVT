function BNums = readLabels(path)
fid = fopen(fullfile(path,'LabelledBranchesQVT.csv'));
%% This function just parses the filled out labelledbranches form for all branch numbers
%% It reads line bny line then outputs the branch numbers and names into a cell
tline = fgetl(fid);
BNums=cell(size(9,2));
for i=1:9
tline = fgetl(fid);
A = strsplit(tline,',');
LABLE=A{2};
for j=2:length(A)-1
    LABLE=append(LABLE,',',A{j+1});
end
B = strsplit(LABLE,'"');
BNums(i,1)=A(1);
BNums(i,2)={str2num(B{2})};
end
fclose(fid);