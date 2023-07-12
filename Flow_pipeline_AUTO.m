clear;clc;
%% Initialization
path2bids='Z:\Sergio\HITH_Control';
subject='sub-001';
path2labels=fullfile(path2bids,'derivatives\QVT',subject); % just to the root folder, can only need this if not using bids
percentileCutoff=0.2;
plotflag=1; %0 = no plot, 1 = make plot
plotraw=1; %0 = no plot, 1 = make plot with raw curves plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Don't change below %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Or do
%% Load Necessities from labels and processed data
DataName = dir(fullfile(path2labels,'*.mat'));
load(fullfile(path2labels,DataName.name));
foldername=strcat('Cutoff_',num2str(percentileCutoff));
try mkdir(fullfile(path2labels,foldername))
end
Labels=readLabels(path2labels);
BranchNum=Labels(:,2);
TitleNames=Labels(:,1);
BranchList=data_struct.branchList;
Quality=data_struct.StdvFromMean;
Quality=[[1:length(Quality)]' Quality BranchList(:,5)];
Flow=data_struct.flowPulsatile_val;
[~,frames]=size(Flow);
time=[data_struct.timeres:data_struct.timeres:(data_struct.timeres*(frames))]./1000;
%% Grab HQ flows for each artery
HQflows=cell(size(9,4));
for i=1:9
    QA=[];
    Blist=BranchNum{i}; %grab the branch labels for artery
    for n=1:length(Blist)
        Q=Quality(BranchList(:,4)==Blist(n),:);
        QA=[QA;Q]; % add all quality indexes from all branch labels into 1
    end
    a=sortrows(QA,2,"descend"); %sort quality into high - low
    cutoff=round(length(a)*percentileCutoff); %compute \% cutoff value
    if cutoff<5
        cutoff =5; %in case small vessel
    end
    Idxs=a(1:cutoff,1); %get row indices (of raw total branch list)
    HQflows{i,1}=Flow(Idxs,:); %store raw HQ flow data, mean, std, and mean flows 
    HQflows{i,2}=mean(Flow(Idxs,:));
    HQflows{i,3}=std(Flow(Idxs,:));
    HQflows{i,4}=mean(Flow(Idxs,:)');
end
%% Compute BBF's and Deviation Estimates
TotalInletFlow=mean(HQflows{1,2})+mean(HQflows{2,2})+mean(HQflows{9,2});
TotalInletFlowMin=min(HQflows{1,2})+min(HQflows{2,2})+min(HQflows{9,2});
TotalInletFlowMax=max(HQflows{1,2})+max(HQflows{2,2})+max(HQflows{9,2});
TotalInletFlowStd=sqrt(std(HQflows{1,4}).^2+std(HQflows{2,4}).^2+std(HQflows{9,4}).^2);
divisionDevIn=(TotalInletFlowStd./TotalInletFlow).^2;
%fprintf('Total Inlet Flow = %3.1f +- %3.1f mL/s\n',TotalInletFlow,TotalInletFlowStd)
for i=3:8
    TotalOutletFlow=mean(HQflows{i,2});
    TotalOutletFlowMin=min(HQflows{i,2});
    TotalOutletFlowMax=max(HQflows{i,2});
    TotalOutletFlowStd=std(HQflows{i,4});
    divisionDevOut=(TotalOutletFlowStd./TotalOutletFlow).^2;
    BFF(i-2,1)=TotalOutletFlow/TotalInletFlow; %mean BFF
    BFF(i-2,2)=BFF(i-2,1).*sqrt(divisionDevIn+divisionDevOut); %std BFF
    BFF(i-2,3)=TotalOutletFlowMax/TotalInletFlow; %max BFF
    BFF(i-2,4)=TotalOutletFlowMin/TotalInletFlow; %min BFF
    %fprintf(strcat(TitleNames{i},' flow = %3.1f +- %3.1f mL/s\n'),TotalOutletFlow,TotalOutletFlowStd)
end
T = max(time); %period of average heartbeat length from 4Dflow
v_bc_name = {TitleNames{1};TitleNames{2};TitleNames{9};'Sum';'T'};
v_bc = [mean(HQflows{1,2});mean(HQflows{2,2});mean(HQflows{9,2});TotalInletFlow;T];
v_bcstd = [std(HQflows{1,4});std(HQflows{2,4});std(HQflows{9,4});TotalInletFlowStd;0];
BFF_name = {'MCA_L';'MCA_R';'ACA_L';'ACA_R';'PCA_L';'PCA_R'};
outputname = strcat(path2labels,'\',foldername,'\Flow_ModelConfig_auto.xlsx'); %File Name in location of subject
writecell(v_bc_name, outputname, 'Sheet', 'BC', 'Range', 'A1');
writematrix([v_bc v_bcstd], outputname, 'Sheet', 'BC', 'Range', 'B1');
writecell(BFF_name, outputname, 'Sheet', 'BFF', 'Range', 'A1');
writematrix(BFF, outputname, 'Sheet', 'BFF', 'Range', 'B1');
%% Fit Fourier Series to Data
% sheetNames = {'Left ICA Cavernous_T_resolved','Right ICA Cavernous_T_resolved',...
%     'Left MCA_T_resolved','Right MCA_T_resolved',...
%     'Left ACA_T_resolved','Right ACA_T_resolved',...
%     'Left PCA_T_resolved','Right PCA_T_resolved',...
%     'Basilar_T_resolved'}; %for saving to comply with current pipeline, but should change eventually
sheetNames = {'Left ICA','Right ICA',...
    'Left MCA','Right MCA',...
    'Left ACA','Right ACA',...
    'Left PCA','Right PCA',...
    'Basilar'}; %for saving to comply with current pipeline, but should change eventually
for i = 1:9
    y = HQflows{i,2}; x = time;
    % Define the parameter value outside the fittype
    % Create a custom function that includes both coefficients and parameters
    myFunction = @(a0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,a8,b8,a9,b9, x) myEquation(a0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,a8,b8,a9,b9, x, T);
    % Create a fittype object using the custom function
    myFitType = fittype(myFunction, 'independent', 'x', 'coefficients', {'a0','a1','b1','a2','b2','a3','b3','a4','b4','a5','b5','a6','b6','a7','b7','a8','b8','a9','b9'});
    V = fit(x',y',myFitType);
    coeffNames = coeffnames(V); % Get the coefficient names
    C = coeffvalues(V); % Get the coefficient values
    sheet=sheetNames{i};

    % Write the coefficient data to the excel file
    writecell(coeffNames, outputname, 'Sheet', sheet, 'Range', 'A1');
    writematrix(C', outputname, 'Sheet', sheet, 'Range', 'B1');
    
    rawsheetname=strcat(sheet,' raw');
    % Write in the raw flows, error, and time (for plotting or comparing the
    % model against)
    FlowNames={'flow_mean';'flow_std';'time'};
    FlowVals=[y;HQflows{i,3};time];
    writecell(FlowNames, outputname, 'Sheet', rawsheetname, 'Range', 'A1');
    writematrix(FlowVals, outputname, 'Sheet', rawsheetname, 'Range', 'B1');

    % write in interpolation data (saves running the coefficients later)
    FitTime=0:T/99:T;
    FitFlow(i,:)=myEquation(C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),C(12),C(13),C(14),C(15),C(16),C(17),C(18),C(19), FitTime, T);
    FlowNames={'InterpFlow';'InterpTime'};
    FlowVals=[FitFlow(i,:);FitTime];
    writecell(FlowNames, outputname, 'Sheet', rawsheetname, 'Range', 'A4');
    writematrix(FlowVals, outputname, 'Sheet', rawsheetname, 'Range', 'B4');
end

if plotflag ~= 0
    %% Compute Equal Left and Right Artery Limits (for plotting)
    Mins=zeros([1,9]);
    Maxs=zeros([1,9]);
    for i=1:9
        mn=HQflows{i,2};
        Mins(i)=min(HQflows{i,1}(:))-0.05*min(HQflows{i,1}(:));
        if rem(i,2) == 0
            Mins((i-1):i)=min([Mins(i-1) Mins(i)]);
        end
        Maxs(i)=max(HQflows{i,1}(:))+0.05*max(HQflows{i,1}(:));
        if rem(i,2) == 0
            Maxs((i-1):i)=max([Maxs(i-1) Maxs(i)]);
        end
    end
    %% plot HQ flow
    figure(1)
    tiledlayout(5,2);
    for i=1:9
        nexttile
        MEAN=HQflows{i,2}';
        STD=HQflows{i,3}';
        h=fill([time';flipud(time')],[MEAN-STD;flipud(MEAN+STD)],[0 0 0],'Linestyle','None');
        set(h,'facealpha',.3)
        hold on
        if plotraw == 1
            rflow=HQflows{i,1};
            for j=1:length(rflow(:,1))
                plot(time,rflow(j,:),'-','Color',[0.2 0.5 0.9 0.3])
            end
        end
        plot(time,MEAN,'k*','MarkerSize',2)
        plot(FitTime,FitFlow(i,:),'k-')
        title(TitleNames{i})
        ylim([Mins(i) Maxs(i)])
        xlim([0 max(time)])
    end
    lgd = legend('STD of Flow','Flow Mean','Continuous Fit');
    lgd.FontSize=14;
    lgd.Layout.Tile = 10;
    %% Plot Bond Graph Results - BLANK FOR NOW
end
%% Custom function that incorporates coefficients and parameters
function y = myEquation(a0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,a8,b8,a9,b9, x, T)
    y= a0 + a1*cos(2*pi*x/T) + b1*sin(2*pi*x/T)+ ...
    + a2*cos(2*pi*2*x/T) + b2*sin(2*pi*2*x/T) ...
    + a3*cos(2*pi*3*x/T) + b3*sin(2*pi*3*x/T) ...
    + a4*cos(2*pi*4*x/T) + b4*sin(2*pi*4*x/T) ...
    + a5*cos(2*pi*5*x/T) + b5*sin(2*pi*5*x/T) ...
    + a6*cos(2*pi*6*x/T) + b6*sin(2*pi*6*x/T) ...
    + a7*cos(2*pi*7*x/T) + b7*sin(2*pi*7*x/T) ...
    + a8*cos(2*pi*8*x/T) + b8*sin(2*pi*8*x/T) ...
    + a9*cos(2*pi*9*x/T) + b9*sin(2*pi*9*x/T);
end