basedir = pwd;
now = datestr(now,'mm/dd/yyyy');

vesselNames = {'Left ICA Cavernous';'Left ICA Cervical';'Right ICA Cavernous'; ...
    'Right ICA Cervical';'Basilar';'Left MCA';'Right MCA';'Left PCA';...
    'Right PCA';'SS sinus';'Straight sinus';'Left Transverse';...
    'Right Transverse';'Left VA';'Right VA';'Left ACA';'Right ACA'};
params = {'FLOW';'PULSATILITY';'RESISTIVITY';'AREA';'MEAN_V';'BRANCH';'CENTERLINE'};
DATA = NaN(1,length(vesselNames),length(params));

saveOrder = {'LCER';'RCER';'LCAV';'RCAV';'BA';'LVA';'RVA';'LMCA';...
'RMCA';'LACA';'RACA';'LPCA';'RPCA';'SSS';'SS';'LTS';'RTS'};

warning('on')

dataDir = dir('*pcviprData*');
dataDir = dataDir([dataDir.isdir]);
cd(dataDir(1).name)
[num,txt,raw] = xlsread('SummaryParamTool.xls');

flow = cell2mat(raw(2:18,5));
flow = matchParams(flow);
DATA(1,:,1) = flow';

% PULSATILITY
pulsatility = cell2mat(raw(2:18,6));
pulsatility = matchParams(pulsatility);
DATA(1,:,2) = pulsatility';

% RESISTIVITY
warning('on');
resistivity = getParams(vesselNames,8);
DATA(1,:,3) = resistivity';
warning('off');

% VESSEL AREA
area = getParams(vesselNames,2);
DATA(1,:,4) = area';

% MEAN VELOCITY
mean_v = getParams(vesselNames,5);
DATA(1,:,5) = mean_v';

FLOW = squeeze(DATA(:,:,1))*60;
TCBF = FLOW(1)+FLOW(2)+FLOW(5);
PULSATILITY = squeeze(DATA(:,:,2));
RESISTIVITY = squeeze(DATA(:,:,3));
AREA = squeeze(DATA(:,:,4));
MEAN_V = squeeze(DATA(:,:,5));

cd(basedir)
% cd ../DICOM
% ddd = dir('*dcm');
% if isempty(ddd)
%     istgz = dir('*.bz2');
%     if ~isempty(istgz)
%         system(['"C:\Program Files\7-Zip\7z.exe" x ' istgz(1).name ' -aoa']);
%         ddd = dir('*dcm');
%     end 
% end 
%info = dicominfo(ddd(1).name);
%scan_date = info.FileModDate;
%mri_loc = info.InstitutionName;
%coil = info.ReceiveCoilName;
    
%mris = {now scan_date mri_loc};
flows =  num2cell([TCBF FLOW PULSATILITY RESISTIVITY AREA MEAN_V]);
params = flows;
%params = {mris{:} flows{:}};
cd(basedir);

clear area AREA basedir coil DATA dataDir ddd flow FLOW info mean_v MEAN_V
clear mri_loc now num pulsatility PULSATILITY raw resistivity RESISTIVITY 
clear saveOrder scan_date TCBF txt vesselNames rrInt mris flows


function param = getParams(vesselNames,col)
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{2}   '_T_averaged'],'');
        param(1) = t_avg(end-1,col);
    catch
        param(1) = NaN;
        warning(['NO ' vesselNames{2} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{4}   '_T_averaged'],'');
        param(2) = t_avg(end-1,col);
    catch
        param(2) = NaN;
        warning(['NO ' vesselNames{4} ' FOUND']);
    end 
    
    try
        t_avg  = xlsread('SummaryParamTool.xls',[vesselNames{1}   '_T_averaged'],'');
        param(3) = t_avg(end-1,col);
    catch
        param(3) = NaN;
        warning(['NO ' vesselNames{1} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{3}   '_T_averaged'],'');
        param(4) = t_avg(end-1,col);
    catch
        param(4) = NaN;
        warning(['NO ' vesselNames{3} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{5}   '_T_averaged'],'');
        param(5) = t_avg(end-1,col);
    catch
        param(5) = NaN;
        warning(['NO ' vesselNames{5} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{14}  '_T_averaged'],'');
        param(6) = t_avg(end-1,col);
    catch
        param(6) = NaN;
        warning(['NO ' vesselNames{14} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{15}  '_T_averaged'],'');
        param(7) = t_avg(end-1,col);
    catch
        param(7) = NaN;
        warning(['NO ' vesselNames{15} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{6}   '_T_averaged'],'');
        param(8) = t_avg(end-1,col);
    catch
        param(8) = NaN;
        warning(['NO ' vesselNames{6} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{7}   '_T_averaged'],'');
        param(9) = t_avg(end-1,col);
    catch
        param(9) = NaN;
        warning(['NO ' vesselNames{7} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{16} '_T_averaged'],'');
        param(10) = t_avg(end-1,col);
    catch
        param(10) = NaN;
        warning(['NO ' vesselNames{16} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{17} '_T_averaged'],'');
        param(11) = t_avg(end-1,col);
    catch
        param(11) = NaN;
        warning(['NO ' vesselNames{17} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{8}  '_T_averaged'],'');
        param(12) = t_avg(end-1,col);
    catch
        param(12) = NaN;
        warning(['NO ' vesselNames{8} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{9}  '_T_averaged'],'');
        param(13) = t_avg(end-1,col);
    catch
        param(13) = NaN;
        warning(['NO ' vesselNames{9} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{10} '_T_averaged'],'');
        param(14) = t_avg(end-1,col);
    catch
        param(14) = NaN;
        warning(['NO ' vesselNames{10} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{11} '_T_averaged'],'');
        param(15) = t_avg(end-1,col);
    catch
        param(15) = NaN;
        warning(['NO ' vesselNames{11} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{12} '_T_averaged'],'');
        param(16) = t_avg(end-1,col);
    catch
        param(16) = NaN;
        warning(['NO ' vesselNames{12} ' FOUND']);
    end 
    
    try
        t_avg = xlsread('SummaryParamTool.xls',[vesselNames{13} '_T_averaged'],'');
        param(17) = t_avg(end-1,col);
    catch
        param(17) = NaN;
        warning(['NO ' vesselNames{13} ' FOUND']);
    end 
end

function newParam = matchParams(param)
    newParam(1) = param(2);
    newParam(2) = param(4);
    newParam(3) = param(1);
    newParam(4) = param(3);
    newParam(5) = param(5);
    newParam(6) = param(14);
    newParam(7) = param(15);
    newParam(8) = param(6);
    newParam(9) = param(7);
    newParam(10) = param(16);
    newParam(11) = param(17);
    newParam(12) = param(8);
    newParam(13) = param(9);
    newParam(14) = param(10);
    newParam(15) = param(11);
    newParam(16) = param(12);
    newParam(17) = param(13);
end 