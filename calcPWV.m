% script provided by James Rice 6/15/26, attributed to Tim Ruesink, used
% for his paper on the flow phantom.


set(0,'DefaultFigureColor','w','DefaultAxesFontSize',14,...
    'DefaultTextColor','k','DefaultAxesXColor','k',...
    'DefaultAxesYColor','k','DefaultAxesZColor','k',...
    'DefaultLineMarkerSize',30,'DefaultLineLineWidth',2,...
    'DefaultTextFontSize',14);

% % Black
% set(0,'DefaultFigureColor','k','DefaultAxesFontSize',14,...
%     'DefaultTextColor','w','DefaultAxesXColor','w',...
%     'DefaultAxesYColor','w','DefaultAxesZColor','w',...
%     'DefaultLineMarkerSize',30,'DefaultLineLineWidth',2,...
%     'DefaultLineColor','w','DefaultTextFontSize',14);

[~,STUDY,~] = fileparts(pwd);

planeSelection = []; % empty (=[]) implies all available planes
if isempty(planeSelection),planeSelection = 1:100; end

tmp = xlsread('Waveform1.csv');
flow_time_init = tmp(:,1);
numFrames = length(flow_time_init); % arbitrary number of time frames

%% READ IN WAVEFORM FILES FROM ENSIGHT
flow_init = [];
cnt = 1;
for i = planeSelection
    fn = ['Waveform',num2str(i),'.csv'];
    if ~exist(fn)
        planeSelection(i:end) = [];
        break;
    end
        
    fprintf('Reading in waveform from plane %i\n',i);
    flow_init(:,cnt) = xlsread(fn,['B2:B',num2str(numFrames+1)]);
    cnt = cnt + 1;
end

numPlanes = size(flow_init,2);
try
    tmp = xlsread('Spline_Distance.csv');
catch
    tmp = xlsread('Spline Distance.csv');
end
dist = tmp(planeSelection,1)';

set(0,'DefaultAxesColorOrder',hot(numPlanes+2));

flag_interpolate = 1;

%% INTERPOLATE
if flag_interpolate
    
    time_interpolant = (flow_time_init(end)-flow_time_init(1))/399;
    flow_time = flow_time_init(1):time_interpolant:flow_time_init(end);

    % Perform interpolation
    flow = [];
    for i = 1:numPlanes
        flow_init(:,i) = smooth(flow_init(:,i),'moving');
        flow(:,i) = (interp1(flow_time_init,flow_init(:,i),flow_time,'spline'));
    end    
    
    cmap = hot(numPlanes+1);
    leg_str = {};
    figure('units','inches','pos',[0.5,2,7,6]);
    for i = 1:numPlanes
        plot(flow_time,flow(:,i),'Color',cmap(i,:),'LineWidth',3);hold on;
        leg_str = [leg_str,{['Plane ',num2str(i)]}];
    end
    set(gca,'Box','off');
    axis tight;xlabel('Time [ms]');ylabel('Flow [m/s]');
    legend(leg_str);title(STUDY, 'Interpreter','none');

else
    flow_time = flow_time_init;
    flow = flow_init;
end
numPoints = length(flow_time);

%% TTP (Time-to-Peak)
[~,I] = max(flow,[],1);
ttp = flow_time(I);

%% TTU (Time-to-Upstroke) point of maximum acceleration of upslope
% Find where derivative is maximal on upslope (maximum acceleration)
[~,I] = max(gradient(flow'),[],2);
ttu = flow_time(I);

%% TTF (Time-to-Foot)
[maxVals,I] = max(flow,[],1);
normflow = [];for i = 1:numPlanes,normflow(:,i) = mat2gray(flow(:,i));end
ind80 = nan(1,numPlanes);
ind20 = nan(1,numPlanes);
P = [];upslope_fit = [];
for i = 1:numPlanes
    
    % Find 20% and 80% of peak
    [maxVal,~] = max(normflow(:,i));
    for j = I(i):-1:1
        if normflow(j,i) < 0.8*maxVal && isnan(ind80(i))
            ind80(i) = j;
        end
        if normflow(j,i) < 0.2*maxVal
            ind20(i) = j;
            break;
        end
    end
    
    if isnan(ind20(i))
        ind20(i) = 1; 
    end
    
    %  y = P(1)*x + P(2)   --->   0 = P(1)*x + P(2) ---->  x = -P(2)/P(1)
    P(:,i) = polyfit([flow_time(ind20(i)),flow_time(ind80(i))],[flow(ind20(i),i),flow(ind80(i),i)],1);
    upslope_fit(:,i) = polyval(P(:,i),flow_time);
    
end

% Find when the upslope hits the zero line (TTF)  x = -P(2)/P(1) from polyfit() output
ttf = -P(2,:)./P(1,:);

% figure;
% plot(flow_time,flow,'LineWidth',2);axis tight;hold on;
% plot(flow_time,upslope_fit,'--');
% xlabel('Time (ms)');ylabel('Flow (m/s)');title('Time-to-Foot (TTF)');
% axis([flow_time(1),flow_time(max(I)+30),min(flow(:)),max(flow(:))]);

% ttf = flow_time(ind20);
% figure;
% plot(flow_time,flow,'LineWidth',2);axis tight;hold on;
% plot(flow_time(ind20),flow(sub2ind(size(flow),ind20,1:numPlanes)),'o');
% xlabel('Time (ms)');ylabel('Flow (m/s)');title('Time-to-Foot (TTF)');
% axis([flow_time(1),flow_time(max(I)+30),min(flow(:)),max(flow(:))]);


%% XCorr (Cross correlation)
[~,I] = max(mean(flow,2));
firstDerivFlow = gradient(mean(flow,2)); 
iLow = 1;iHigh = numPoints;
for k = (I-2):-1:1
    if firstDerivFlow(k) < 0
        iLow = k;
        break;
    end
end
for k = (I+2):numPoints
    if firstDerivFlow(k) > 0
        iHigh = k;
        break;
    end
end
dt = flow_time(2)-flow_time(1);

xcor = [];
% figure;
for i = 1:numPlanes
    [c, lags] = xcorr(mat2gray(flow(iLow:iHigh,i)),mat2gray(flow(iLow:iHigh,1)),50);
    [~,I] = max(c);
    xcor(i) = dt*lags(I);
%     plot(flow_time(iLow:iHigh),mat2gray(flow(iLow:iHigh,1)),'k');hold on;plot(flow_time(iLow:iHigh)-dt*lags(I),mat2gray(flow(iLow:iHigh,i)),'r');hold off;pause;
end

%% Plot time shifts vs. distance
%  y = P(1)*x + P(2)
p = [];pfit = [];
p(:,1) = polyfit(dist,ttp,1);pfit(:,1) = polyval(p(:,1),dist);
p(:,2) = polyfit(dist,ttu,1);pfit(:,2) = polyval(p(:,2),dist);
p(:,3) = polyfit(dist,ttf,1);pfit(:,3) = polyval(p(:,3),dist);
p(:,4) = polyfit(dist,xcor,1);pfit(:,4) = polyval(p(:,4),dist);
pwv = 1./p(1,:)

%solo TTP
figure;plot(dist,ttp,'.');hold on;plot(dist,pfit(:,1),'r');
text('Units','normalized','position',[0.02,0.85,1],'String',['PWV_{TTP} = ',sprintf('%0.2f',(pwv(1))),'m/s']);
xlim([min(dist),max(dist)]);
ylabel('TTP [ms]');
title(STUDY,'Interpreter','none');

figure('Units','inches','position',[8,1,4,8]);
if sum(get(0,'DefaultFigureColor'))==0,set(gcf,'DefaultAxesColorOrder',[1,1,1]);end;

subplot(4,1,1);plot(dist,ttp,'.');hold on;plot(dist,pfit(:,1),'r');
text('Units','normalized','position',[0.02,0.85,1],'String',['PWV_{TTP} = ',sprintf('%0.2f',(pwv(1))),'m/s']);
xlim([min(dist),max(dist)]);
ylabel('TTP [ms]');
title(STUDY,'Interpreter','none');

subplot(4,1,2);plot(dist,ttu,'.');hold on;plot(dist,pfit(:,2),'r');
text('Units','normalized','position',[0.02,0.85,1],'String',['PWV_{TTU} = ',sprintf('%0.2f',(pwv(2))),'m/s']);
xlim([min(dist),max(dist)]);
ylabel('TTU [ms]');

subplot(4,1,3);plot(dist,ttf,'.');hold on;plot(dist,pfit(:,3),'r');
text('Units','normalized','position',[0.02,0.85,1],'String',['PWV_{TTF} = ',sprintf('%0.2f',(pwv(3))),'m/s']);
xlim([min(dist),max(dist)]);
ylabel('TTF [ms]');

subplot(4,1,4);plot(dist,xcor,'.');hold on;plot(dist,pfit(:,4),'r');
text('Units','normalized','position',[0.02,0.85,1],'String',['PWV_{XCOR} = ',sprintf('%0.2f',(pwv(4))),'m/s']);
xlim([min(dist),max(dist)]);
ylabel('XCorr Delay [ms]');
xlabel('Slice Position [mm]');


% Bar Plot
% figure;bar(pwv);colormap gray;set(gca,'XTickLabel',{'TTP','TTU','TTF','XCor'},'box','off');ylabel('PWV [m/s]');

drawnow;