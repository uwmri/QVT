% script provided by James Rice 6/15/26, attributed to Tim Ruesink, maybe 
% used for his paper on the flow phantom? one edit to change filename
% SplineDistance.csv to Spine_Distance.csv, then discovered the script 
% expects that file to have four columns, and the file I have for the flow
% phantom only has two, so this script doesn't work? -- RVC 6/16/26


% mri_pwv
%
% SUMMARY performs a regional and local PWV assesment using 3 different
% transit time algorithms

clear all
close all

%% Default format

prop = get(0, 'default');
propname = fieldnames(prop);
for iprop = 1:length(propname)
   set(0, propname{iprop}, 'remove');
end

reset(gca)
reset(gcf)
set(0,'DefaultFigureColor','w','DefaultAxesFontSize',12,...
    'DefaultTextColor','k','DefaultAxesXColor','k',...
    'DefaultAxesYColor','k','DefaultAxesZColor','k',...
    'DefaultLineMarkerSize',5,'DefaultLineLineWidth',1,...
    'DefaultTextFontSize',12,...
    'defaultTextFontName','Times New Roman',...
    'defaultAxesFontName','Times New Roman');

[~,STUDY,~] = fileparts(pwd);

planeSelection = []; % empty (=[]) implies all available planes
if isempty(planeSelection),planeSelection = 1:100; end

tmp = xlsread('Waveform1.csv');
flow_time_init = tmp(:,1);
numFrames = length(flow_time_init); % arbitrary number of time frames

%% Read in csv files exported from EnSight

flow_init = [];
cnt = 1;
for i = planeSelection
    fn = ['Waveform',num2str(i),'.csv'];
    if ~exist(fn)
        planeSelection(i:end) = [];
        break;
    end
        
    fprintf('Reading in waveform from plane %i\n',i);
    tmp2 = xlsread(fn);
    flow_init(:,cnt) = tmp2(:,2);
    cnt = cnt + 1;
end

numPlanes = size(flow_init,2);

% edited here, RVC 6/16/26
tmp = xlsread('Spline_Distance.csv');
dist = tmp(3:end,4);

set(0,'DefaultAxesColorOrder',hot(numPlanes+2));

flag_interpolate = 1;

%% shift the waveforms so the max is in the middle of the series

[~,i_tmp] = max(flow_init(:,1));
for i = 1:numPlanes
    flow_init(:,i) = circshift(flow_init(:,i),round(length(flow_time_init)/2)-i_tmp);
end

%% filter signal

% perform a Savitzky-Golay filtering function

% need to choose based on time step
if rem(length(flow_time_init), 2) == 0
    framelen = length(flow_time_init)-1;
else framelen = length(flow_time_init)-2;
end

order = round(length(flow_time_init)/3); % need to choose based on time step
if rem(order,2) == 0
    order = order-1;
else order = order;
end

flow_filtered = zeros(size(flow_init));
for i = 1:numPlanes
    flow_filtered(:,i) = sgolayfilt(flow_init(:,i),order,framelen);
end

%% Visualize the original waveforms and filtered waveforms

fig1 = figure(1);
cmap = autumn(numPlanes);
spacing = round(numPlanes/5);
orig_sel = flow_init(:,1:spacing:end);
filtered_sel = flow_filtered(:,1:spacing:end);
hold on

for i = 1:5
    plt1(i) = plot(flow_time_init, orig_sel(:,i),'--b');
    plt2(i) = plot(flow_time_init, filtered_sel(:,i),'r');
end
legend([plt1(1) plt2(1)],{'Raw','Filtered'});
title('Raw vs. filtered waveforms');

hold off

%% Smooth and upsample the filtered waveforms

step = 1; % interpolant time step is 1 ms
t = flow_time_init(1):step:flow_time_init(end); % time vector for upsampling

for i = 1:numPlanes
    flow(:,i) = interp1(flow_time_init,flow_filtered(:,i),t,'spline');
end

%% Visualize upsampled and smoothed data

fig2 = figure(2);
cmap = autumn(numPlanes);
spacing = round(numPlanes/5);
filtered_sel = flow_filtered(:,1:spacing:end);
final_sel = flow(:,1:spacing:end);
hold on

for i = 1:5
    plt1(i) = plot(flow_time_init, filtered_sel(:,i),'--b');
    plt2(i) = plot(t, final_sel(:,i),'r');
end
legend([plt1(1) plt2(1)],{'Filtered','Smoothed'});
title('Filtered vs Smooth waveforms');
hold off

%% Initialize variables

N = numPlanes;
t = transpose(t);
t_s = -500:step:2000; % larger time vector for linear approximations

xcor = zeros(N-1,1); % transit time using cross correlation
ttf = zeros(N-1,1); % transit time using time to foot
ttu = zeros(N-1,1); % transit time using tim to upstroke

[F_max,i_max] = max(flow); % maximum of flow for all planes
[F_min,i_min] = min(flow(1:i_max,:)); % minimum of flow for all planes
dFlow = zeros(size(flow));
i_dmax = zeros(N,1);

k_20 = zeros(N,1); % point of 20% upstroke
k_80 = zeros(N,1); % point of 80% upstroke
ss = zeros(N,1); % index that start the linear approximation
ee = zeros(N,1); % index that ends the linear approximation
Lpwv_xcor = zeros(N,1);
Lpwv_ttf = zeros(N,1);
Lpwv_ttu = zeros(N,1);
Lpwv_xcor_filt = zeros(N,1);
Lpwv_ttf_filt = zeros(N,1);
Lpwv_ttu_filt = zeros(N,1);
F_line = zeros(length(t),N);

%%
figure;
plot(t,flow(:,1),'r')
title('Select start of upstroke')
[a,~] = ginput(1);
L = find(t <= (a + 0.5) & t >= (a - 0.5));

%% Analyze the first waveform
% find the indicies where the flow is between 19% and 21% of the maximum
% value.
k = find(flow(L:i_max(1),1) > (0.19*(F_max(1)-F_min(1))+F_min(1))...
    & flow(L:i_max(1),1) < (0.21*(F_max(1)-F_min(1))+F_min(1)));
k_20(1) = k(round(length(k)/2)); % index of 20% of maximum flow
k_20(1) = k_20(1) + L;

% find the indicies where the flow is between 79% and 81% of the maximum
% value.
k = find(flow(L:i_max(1),1) > (0.79*(F_max(1)-F_min(1))+F_min(1))...
    & flow(L:i_max(1),1) < (0.81*(F_max(1)-F_min(1))+F_min(1)));
k_80(1) = k(round(length(k)/2)); % index of 80% of maximum flow
k_80(1) = k_80(1) + L;

% Now we find the linear regression describing the points between 20% and
% 80% of the maximum flow
p = polyfit(t(k_20(1):k_80(1)),flow(k_20(1):k_80(1),1),1);
F_line(:,1) = p(1).*t + p(2);
throw1 = find(F_line(:,1) < flow(L,1)); % find all indicies that are less than 0
throw2 = find(F_line(:,1) > F_max(1)); % find all indicies that are more than the maximum flow
ss(1) = max(throw1);
ee(1) = min(throw2);

% Now we look at the 1st derivative of flow
dFlow(:,1) = gradient(flow(:,1));
[~,i_dmax(1)] = max(dFlow(L:i_max(1),1));
i_dmax(1) = i_dmax(1) + L;

%% Visualize the first waveform and its linear approximation

figure;
plot(t,flow(:,1),'r',t(ss(1):ee(1)),F_line(ss(1):ee(1),1),'b')
title('Visualize first waveform')

for i = 2:N
    
    % CROSS CORRELATION

    [C,lag] = xcorr(mat2gray(flow(:,i)),mat2gray(flow(:,1)));
    [~,I] = max(abs(C));
    xcor(i) = step*lag(I);

    % TIME TO FOOT

    % find the indicies where the flow is between 19% and 21% of the maximum
    % value.
    k = find(flow(L:i_max(i),i) > (0.15*(F_max(i)-F_min(i))+F_min(i))...
        & flow(L:i_max(i),i) < (0.25*(F_max(i)-F_min(i))+F_min(i)));
    k = k + L;
    k_20(i) = k(round(length(k)/2)); % index of 20% of maximum flow

    % find the indicies where the flow is between 79% and 81% of the maximum
    % value.
    k = find(flow(L:i_max(i),i) > (0.79*(F_max(i)-F_min(i))+F_min(i))...
        & flow(L:i_max(i),i) < (0.81*(F_max(i)-F_min(i))+F_min(i)));
    k = k + L;
    k_80(i) = k(round(length(k)/2)); % index of 80% of maximum flow

    % Now we find the linear regression describing the points between 20% and
    % 80% of the maximum flow
    p = polyfit(t(k_20(i):k_80(i)),flow(k_20(i):k_80(i),i),1);
    F_line(:,i) = p(1).*t + p(2);
    throw1 = find(F_line(:,i) < flow(L,i)); % find all indicies that are less than 0
    throw2 = find(F_line(:,i) > F_max(i)); % find all indicies that are more than the maximum flow
    ss(i) = max(throw1);
    ee(i) = min(throw2);
    ttf(i) = step*(ss(i)-ss(1));

    % TIME TO MAXIMUM UPSTROKE

    dFlow(:,i) = gradient(flow(:,i));
    [~,i_dmax(i)] = max(dFlow(L:i_max(i),i));
    i_dmax(i) = i_dmax(i) + L;
    ttu(i) = step*(i_dmax(i) - i_dmax(1));

end

%% savitzky golay filter for transit time vs distance data

% need to choose based on time step
if rem(numPlanes, 2) == 0
    framelen = numPlanes-1;
else framelen = numPlanes-2;
end

order = round(numPlanes/4); % need to choose based on time step
if rem(order,2) == 0
    order = order-1;
else order = order;
end

%% Calculate PWV
% full aorta PWV
P = zeros(3,2);
P(1,:) = polyfit(dist,xcor,1); pwv_xcor = 1/P(1,1)
P(2,:) = polyfit(dist,ttf,1); pwv_ttf = 1/P(2,1)
P(3,:) = polyfit(dist,ttu,1); pwv_ttu = 1/P(3,1)

PWV = [pwv_xcor;pwv_ttf;pwv_ttu];

xcor_filtered = sgolayfilt(xcor,order,framelen); 
ttf_filtered = sgolayfilt(ttf,order,framelen);
ttu_filtered = sgolayfilt(ttu,order,framelen);

% figure;
% hold on
% plot(dist,xcor,'ok');
for n = 9:N-8
    % region 1
    Q = zeros(3,2);
    Q_filt = zeros(3,2);
    Q(1,:) = polyfit(dist(n-8:n+7),xcor(n-8:n+7),1); Lpwv_xcor(n) = 1/Q(1,1);
    Q(2,:) = polyfit(dist(n-8:n+7),ttf(n-8:n+7),1); Lpwv_ttf(n) = 1/Q(2,1);
    Q(3,:) = polyfit(dist(n-8:n+7),ttu(n-8:n+7),1); Lpwv_ttu(n) = 1/Q(3,1);
    Q_filt(1,:) = polyfit(dist(n-8:n+7),xcor_filtered(n-8:n+7),1); Lpwv_xcor_filt(n) = 1/Q_filt(1,1);
    Q_filt(2,:) = polyfit(dist(n-8:n+7),ttf_filtered(n-8:n+7),1); Lpwv_ttf_filt(n) = 1/Q_filt(2,1);
    Q_filt(3,:) = polyfit(dist(n-8:n+7),ttu_filtered(n-8:n+7),1); Lpwv_ttu_filt(n) = 1/Q_filt(3,1);
%     plt = plot(dist(n-8:n+7),dist(n-8:n+7)*Q(1,1)+Q(1,2),'r');
%     pause(.1)
%     delete(plt)
end

max_xcor = max(Lpwv_xcor);
max_ttf = max(Lpwv_ttf);
max_ttu = max(Lpwv_ttu);

dist_step = dist(2)-dist(1);

%% Local PWV figure

figure;
subplot(3,1,1);
plot(dist(9:end-8),Lpwv_xcor(9:end-8),'r',dist(9:end-8),Lpwv_xcor_filt(9:end-8),'b')
legend('raw','filtered')
title('Local XCOR')
subplot(3,1,2);
plot(dist(9:end-8),Lpwv_ttf(9:end-8),'r',dist(9:end-8),Lpwv_ttf_filt(9:end-8),'b')
legend('raw','filtered')
title('Local TTF')
subplot(3,1,3);
plot(dist(9:end-8),Lpwv_ttu(9:end-8),'r',dist(9:end-8),Lpwv_ttu_filt(9:end-8),'b')
legend('raw','filtered')
title('Local TTU')
xlabel('Distance from aortic root [mm]')

Lpwv_xcor = transpose(Lpwv_xcor);
Lpwv_ttf = transpose(Lpwv_ttf);
Lpwv_ttu = transpose(Lpwv_ttu);
LPWV = [Lpwv_xcor;Lpwv_ttf;Lpwv_ttu];
Lpwv_xcor_filt = transpose(Lpwv_xcor_filt);
Lpwv_ttf_filt = transpose(Lpwv_ttf_filt);
Lpwv_ttu_filt = transpose(Lpwv_ttu_filt);
LPWV_filt = [Lpwv_xcor_filt;Lpwv_ttf_filt;Lpwv_ttu_filt];

%% Regional PWV figure

figure;
subplot(3,1,1)
plot(dist,xcor,'o',dist,dist*P(1,1)+P(1,2),'r',dist,xcor_filtered,'b')
xlabel('Distance [mm]'); ylabel('Transit-Time [ms]');
legend('actual data','linear regression','filtered transit-time');

subplot(3,1,2)
plot(dist,ttf,'o',dist,dist*P(2,1)+P(2,2),'r',dist,ttf_filtered,'b')
xlabel('Distance [mm]'); ylabel('Transit-Time [ms]');
legend('actual data','linear regression','filtered transit-time');

subplot(3,1,3)
plot(dist,ttu,'o',dist,dist*P(3,1)+P(3,2),'r',dist,ttu_filtered,'b')
xlabel('Distance [mm]'); ylabel('Transit-Time [ms]');
legend('actual data','linear regression','filtered transit-time');

%% Figures to analyze all linear regressions

% for n = 1:N
%     figure;
%     plot(t,flow(:,n),'r',t(ss(n):ee(n)),F_line(ss(n):ee(n),n),'b')
% end

%% Figures to analyze all maximum upstroke points
% 
% for n = 1:N
%     figure(n+4)
%     plot(t,flow(:,n),'r',t(i_dmax(n)),flow(i_dmax(n),n),'o')
% end

%% Write out PWV data to a file
LPWV = [transpose(dist);LPWV];
LPWV_filt = [transpose(dist);LPWV_filt];

fileID = fopen([STUDY,'_pwv_data.txt'],'w');
fprintf(fileID,'Regional PWV\n');
fprintf(fileID,'%6s, %6s, %6s\n','XCOR','TTF','TTU');
fprintf(fileID,'%6.3f, %6.3f, %6.3f\n',PWV);
fprintf(fileID,'\nRaw Local PWV\n');
fprintf(fileID,'%8s, %8s, %8s, %8s\n','Distance', 'XCOR','TTF', 'TTU');
fprintf(fileID,'%8.3f, %8.3f, %8.3f, %8.3f\r\n',LPWV);
fprintf(fileID,'\nFiltered Local PWV\n');
fprintf(fileID,'%8s, %8s, %8s, %8s\n','Distance', 'XCOR','TTF', 'TTU');
fprintf(fileID,'%8.3f, %8.3f, %8.3f, %8.3f\r\n',LPWV_filt);
fclose(fileID);
