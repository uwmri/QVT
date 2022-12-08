basedir = 'I:\LIFE\patients\Visit1';
cd(basedir);

d=dir('life*');
pcvipr = zeros(1,length(d));
for i=1:length(d)
    disp(d(i).name)
    cd(d(i).name)
    life = findstr(d(i).name,'life');
    lifeid{i} = d(i).name(life(2):life(2)+8);
    dd = dir('*PCVIPR_Fast');
    if isempty(dd)
        pcvpir(i) = 0;
    else
        pcvipr(i) = 1;
        cd(dd(1).name)
        cd('processed_data/DICOM/')
        ddd = dir('*dcm');
        if isempty(ddd)
            istgz = dir('*.bz2');
            if ~isempty(istgz)
                system(['"C:\Program Files\7-Zip\7z.exe" x ' istgz(1).name ' -aoa']);
                ddd = dir('*dcm');
            end 
        end 
        info = dicominfo(ddd(1).name);
        dt = datetime([info.AcquisitionDate ' ' info.AcquisitionTime],'InputFormat','yyyyMMdd HHmmss');
        scan_date{i} = datestr(dt);
        mri_loc{i} = info.InstitutionName;
        coil{i} = info.ReceiveCoilName;
    
        scan_date = scan_date';
        mri_loc = mri_loc';
        coil = coil';
        lifeid = lifeid';

%         cd(dd(1).name)
%         try
%             cd('processed_data/DAT/')
%         catch
%             cd('processed_data/H5/')
%         end 
%         fid = fopen([pwd filesep 'pcvipr_header.txt'], 'r');
%         delimiter = ' ';
%         formatSpec = '%s%s%[^\n\r]'; %read 2 strings(%s%s),end line(^\n),new row(r)
%         % Info from headers are placed in dataArray, 1x2 cell array.
%         dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, ...
%             'MultipleDelimsAsOne', true, 'ReturnOnError', false);
%         fclose(fid);
%         dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
%         pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);
%         timeres(i) = pcviprHeader.timeres; %temporal resolution (ms)
    end 

    cd(basedir)
end 