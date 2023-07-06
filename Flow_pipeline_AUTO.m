clear;clc;
%% Initialization
Table=readtable("C:\Users\sdem348\Desktop\patient 19\LabelledBranchesQVT.csv");
load("C:\Users\sdem348\Desktop\patient 19\qvtData_ISOfix_06Jul23_1557_v1-2.mat")
percentileCutoff=.10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Don't change below %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Or do
%% Load Necessities
BranchNum=Table.BranchNumber;
TitleNames=Table.Artery;
close all
BranchList=data_struct.branchList;
Quality=data_struct.StdvFromMean;
Quality=[[1:length(Quality)]' Quality BranchList(:,5)];
Flow=data_struct.flowPulsatile_val;
[~,frames]=size(Flow);
time=[data_struct.timeres:data_struct.timeres:(data_struct.timeres*(frames))]./1000;
%% Grab HQ flows for each artery

HQflows=cell(size(9,3));
for i=1:9
    QA=Quality(BranchList(:,4)==BranchNum(i),:);
    QA(QA(:,2)<1,:)=[];
    a=sortrows(QA,2,"descend");
    cutoff=round(length(a)*percentileCutoff);
    if cutoff<5
        cutoff =5;
    end
    Idxs=a(1:cutoff,1);
    HQflows{i,1}=Flow(Idxs,:);
    HQflows{i,2}=mean(Flow(Idxs,:));
    HQflows{i,3}=std(Flow(Idxs,:));
end
%% Fit Fourier Series to data
for i = 1:9
    y = HQflows{i,2}; x = time; T=max(time);
    % Define the parameter value outside the fittype
    % Create a custom function that includes both coefficients and parameters
    myFunction = @(a0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,a8,b8,a9,b9, x) myEquation(a0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,a8,b8,a9,b9, x, T);
    % Create a fittype object using the custom function
    myFitType = fittype(myFunction, 'independent', 'x', 'coefficients', {'a0','a1','b1','a2','b2','a3','b3','a4','b4','a5','b5','a6','b6','a7','b7','a8','b8','a9','b9'});
    V = fit(x',y',myFitType);
    coeffNames = coeffnames(V); % Get the coefficient names
    C = coeffvalues(V); % Get the coefficient values
    FitTime=0:0.01:max(time);
    FitFlow(i,:)=myEquation(C(1),C(2),C(3),C(4),C(5),C(6),C(7),C(8),C(9),C(10),C(11),C(12),C(13),C(14),C(15),C(16),C(17),C(18),C(19), FitTime, T);
end
%% Compute equal Left and Right Limits
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
    set(h,'facealpha',.2)
    hold on
    plot(time,MEAN,'k*','MarkerSize',2)
    plot(FitTime,FitFlow(i,:),'b-')
    title(TitleNames{i})
    ylim([Mins(i) Maxs(i)])
    xlim([0 max(time)])
end
lgd = legend('STD of Flow','Flow Mean','Continuous Fit');
lgd.FontSize=14;
lgd.Layout.Tile = 10;
%% Plot Bond Graph Results - BLANK FOR NOW

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