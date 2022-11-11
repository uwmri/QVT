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
    
        lifeid = lifeid';
        scan_date = scan_date';
        mri_loc = mri_loc';
        coil = coil';
    end 

    cd(basedir)
end 