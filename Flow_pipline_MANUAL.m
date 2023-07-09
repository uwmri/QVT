clear all;clc;
%% Initialization
%currently follows bids, so only the path 2 bids is needed and subject
%number,
path2bids='Z:\Sergio\ExampleBIDS_Study';
subject='sub-001';
path2qvtderiv=fullfile(path2bids,'derivatives','QVT',subject);
path2flows=uigetdir(path2qvtderiv); %Pick the summary param FOLDER you want to process
%if you don't want to follow bids, you only need to input this full filename
% to the folder instead for initialization, but it's less elegant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Don't change below %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Or do
%% Setup params
filename= path2flows +"\SummaryParamTool.xls"; 
sheetNames = {'Left ICA Cavernous_T_resolved', 'Left MCA_T_resolved', 'Left PCA_T_resolved',...
    'Left ACA_T_resolved','Right ICA Cavernous_T_resolved','Right MCA_T_resolved',...
    'Right PCA_T_resolved','Right ACA_T_resolved','Basilar_T_resolved'};
rowNumber = 5; %This is the flow row of the branchpoint chosen
[~, Names] = xlsfinfo(filename); %names of sheets from QVT summary xls
%% Load Mean Flows and save mean BC's and BFF's
% Grab mean flow from the summary page of QVT for each artery
SumTable = readtable(filename);
% order ICA_L ICA_R Basilar MCA_L MCA_R ACA_L ACA_R PCA_L PCA_R
MeanFlows=SumTable.MeanFlowMl_s([2 4 6 7 8 17 18 9 10]-1);
MeanFlows(isempty(MeanFlows))=0; %in case flow wasn't measured in some arter
total_flow = sum(MeanFlows(1:3));
time = table2array(readtable(filename, 'Sheet', "Left ICA Cavernous_T_resolved", 'Range', 'B1:U1')); %read a raw data table for times
T = time(20); %period of average heartbeat length from 4Dflow
v_bc_name = {'ICA_L';'ICA_R';'Basilar';'Sum';'T'};
v_bc = [MeanFlows(1:3);total_flow;T];
BFF_name = {'MCA_L';'MCA_R';'ACA_L';'ACA_R';'PCA_L';'PCA_R'};
BFF = MeanFlows(4:end)./total_flow';
outputname = path2flows + "\Flow_BFFs_and_BCs_manual.xlsx"; %File Name in location of subject
writecell(v_bc_name, outputname, 'Sheet', 'BC', 'Range', 'A1');
writematrix(v_bc, outputname, 'Sheet', 'BC', 'Range', 'B1');
writecell(BFF_name, outputname, 'Sheet', 'BFF', 'Range', 'A1');
writematrix(BFF, outputname, 'Sheet', 'BFF', 'Range', 'B1');
%% Fit Fourier Series to data
for i = 1:numel(sheetNames)
    sheet = sheetNames{i};
    if ismember(sheet, Names)
        %SumTable = readtable(filename,'Sheet','Left ICA Cavernous_T_resolved');
        % Read the specific row of the specific sheet using readmatrix
        flow = readmatrix(filename, 'Sheet', sheet, 'Range', ['B' num2str(rowNumber) ':U' num2str(rowNumber)]);
        y = flow;
        x = time;
        % Define the parameter value outside the fittype
        % Create a custom function that includes both coefficients and parameters
        myFunction = @(a0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,a8,b8,a9,b9, x) myEquation(a0,a1,b1,a2,b2,a3,b3,a4,b4,a5,b5,a6,b6,a7,b7,a8,b8,a9,b9, x, T);
        % Create a fittype object using the custom function
        myFitType = fittype(myFunction, 'independent', 'x', 'coefficients', {'a0','a1','b1','a2','b2','a3','b3','a4','b4','a5','b5','a6','b6','a7','b7','a8','b8','a9','b9'});
        V = fit(x',y',myFitType);
        coeffNames = coeffnames(V); % Get the coefficient names
        C = coeffvalues(V); % Get the coefficient values
        outputname = path2flows + "\FlowCoeffs_manual.xlsx";
        % Write the data to the Excel file
        xlRange = 'A1';
        writecell(coeffNames, outputname, 'Sheet', sheet, 'Range', xlRange);
        xlRange = 'B1';
        writematrix(C', outputname, 'Sheet', sheet, 'Range', xlRange);
    end
end
% Custom function that incorporates coefficients and parameter
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