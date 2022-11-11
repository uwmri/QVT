%% Reads/displays retrospectively gated PCVIPR gating track information.

% Author: Grant Roberts
% Date: September 18 2019
% Updated: December 15 2021

% Fully functional for the following gating file types:
%   ecg_track_xxxxxx WITH SCAN_INFO FILE (oldest)
%   Gating_Track_xxxxxx.pcvipr_track (newer)
%   Gating_Track_xxxxxx.pcvipr_track + Gating_Track_xxxxxx.pcvipr_track.full
%   ScanArchive Gating Files (newest)
%       ECG2Data,ECG2Trig,ECG3Data,ECG3Trig,PPGData,PPGTrig,RESPData,
%       RESPTrig,Gating_Track,Gating_Track_Full
% May not be accurate for prospectively gated data.

%% Filter Gating Directory
gatingDir = uigetdir(pwd, 'Select the directory with the gating tracks');
cd(gatingDir) %move to this directory
%gatingDir = pwd;

fullDir = dir(); 
fullDir(1:2) = [];
keepFiles = ~(contains({fullDir.name},'.md5sum') | contains({fullDir.name},'.h5'))'; 
filtDir = fullDir(keepFiles); %filter out md5sum and h5 files
fileNames = {filtDir.name}'; %get filenames from filtered directory
filesLoaded = {}; %initialize cell array to add loaded gating tracks
ecgData.gatingDir = gatingDir; %add master dir to composite data structure
isEcgTrack = 0;
isGatingTrack = 0;
isScanArchive = 0;
dataLoaded = 0;

diary ecgInformation.txt % start writing command window output to text file


%% Gating Tracks
%%%% Load Gating Data %%%%%
if sum(contains(fileNames,'pcvipr_track'))>1
    isGatingTrack = 1;
    loadFullTrack = 1; %if both tracks present, use full (not done in recon)
    idx = find(contains(fileNames,'pcvipr_track.full'));
elseif sum(contains(fileNames,'pcvipr_track.full'))
    isGatingTrack = 1;
    loadFullTrack = 1; %if we only find full track
    idx = find(contains(fileNames,'pcvipr_track.full'));
elseif sum(contains(fileNames,'pcvipr_track'))==1
    isGatingTrack = 1;
    loadFullTrack = 0; %if we only have gating track
    idx = find(contains(fileNames,'pcvipr_track'));
end 

if isGatingTrack
    name = fileNames{idx}; %grab name from fileNames
    fid = fopen(name);
    gate = fread(fid,'int32','b');
    if loadFullTrack
        gate = reshape(gate,[numel(gate)/5 5]); %full track is Nx5 format
    else
        gate = reshape(gate,[numel(gate)/4 4]);
    end 
    gate = sortrows(gate,3); % organize time chronologically
    fclose(fid);

    % Create gatingTrack struc (See ResearchGating.e from PCVIPR PSD) 
    gatingTrack.ecg = gate(:,1); %ecg (sawtooth) waveform
    gatingTrack.resp = 4095-gate(:,2); %resp waveform (12bit max)
    gatingTrack.time = (gate(:,3)-gate(1,3))/1e3; %time stamp (us/ms)
    %gatingTrack.prep = gate(:,4); %bin by time from prep pulses
    %gatingTrack.acq = gate(:,5); %projection # tag for each encode

    % Match recon RR calculation
    within = gatingTrack.ecg<2000 & gatingTrack.ecg>0; %see gating_lib.cpp
    gatingTrack.recon_rr = 2*median(gatingTrack.ecg(within)); %average RR
    %method above robust to decreased sampling rates of ecg sawtooth fx
    %e.g. full gating tracks vs. regular gating tracks

    % Approximate time-resolved RR intervals
    [rrFull,peaksFull] = getRR(gatingTrack.ecg); 
        gatingTrack.rr = nan(size(gate,1),1); %make array of NaNs
        gatingTrack.rr(peaksFull) = rrFull; %fill in NaNs with RRs
        % account for discretization error (to better match recon rr est.)
        mins = peaksFull+1;
        discShift = round(mean(gatingTrack.ecg(mins)));
        gatingTrack.rr = gatingTrack.rr + discShift;
        gatingTrack.bpm = 60000./gatingTrack.rr;

    % Detect Missed HBs and Early Triggers
    [missedHBIdx,earlyTrigIdx] = findMissedHB(rrFull,peaksFull);
        missedHBs = nan(length(gatingTrack.rr),1); % expand vector to match gate lengths
        earlyTrigs = nan(length(gatingTrack.rr),1); % expand vector to match gate lengths
        missedHBs(missedHBIdx) = gatingTrack.rr(missedHBIdx); % fill in expanded vector with missed HBs
        earlyTrigs(earlyTrigIdx) = gatingTrack.rr(earlyTrigIdx); % fill in expanded vector early triggers
        gatingTrack.missed = missedHBs;
        gatingTrack.early = earlyTrigs;
        gatingTrack.clean = gatingTrack.rr;
        gatingTrack.clean(missedHBIdx) = nan;
        gatingTrack.clean(earlyTrigIdx) = nan; % data which has missed HB and early trigger events removed
    
    % Second Method to DetectOutliers
    outliers = isoutlier(gatingTrack.rr);
    gatingTrack.outliers = nan(size(gate,1),1);
    gatingTrack.outliers(outliers) = gatingTrack.rr(outliers);
    gatingTrack.clean = gatingTrack.rr;
    gatingTrack.clean(outliers) = nan; %remove missed HBs

    fileNames(idx) = []; %clear name from filName for next gating track
    filesLoaded = [filesLoaded, name]; %add file to loaded files running cell
    ecgData.gatingTrack = gatingTrack; %master structure
    
%%%% Plotting and Output %%%%
    gatingTrack.timeres = mean(diff(gatingTrack.time)); %find time resolution by finding most common value
    time = gatingTrack.time/1000; %convert to seconds

    % Respiratory Subfigure
    hFig = figure; 
    hFig.WindowState = 'maximized'; %full-screen
    sgtitle('Gating Track Information'); %create figure title
    subplot(2,1,1); hold on; 
    title('Respiratory Waveforms'); %make subplot and subplot title
    plot(time,gatingTrack.resp,'Color',[0, 0.5, 0.19],'LineWidth',2);
    xlabel('Time (s)'); 
    ylabel('Respiratory Amplitude (a.u.)');
    legend('Respiratory Signal','Location','eastoutside');

    % PG/ECG Subfigure
    avg = nanmean(gatingTrack.rr);
    stdv = nanstd(gatingTrack.rr);
    upperLim = ones(length(gatingTrack.ecg),1).*(avg + 2*stdv); %get max clean RR interval
    meanRR = ones(length(gatingTrack.ecg),1).*avg; %get mean RR
    lowerLim = ones(length(gatingTrack.ecg),1).*(avg - 2*stdv); %get min RR
    reconRR = ones(length(gatingTrack.ecg),1).*gatingTrack.recon_rr; %recon-calculated RR

    subplot(2,1,2); hold on; 
    title('ECG/PG Triggers') %make second subplot
    plot(time,upperLim,'--','Color','black'); %plot dashed line of max RR
    plot(time,lowerLim,'--','Color','black','HandleVisibility','off'); %plot dashed line of min RR
    plot(time,meanRR,'-','Color','black'); %plot dashed line of mean RR
    plot(time,reconRR,'-','Color','cyan'); %plot recon RR
    scatter(time,gatingTrack.rr,35,'blue'); %scatter plot of rr intervals over time
    scatter(time,gatingTrack.outliers,35,'filled','blue'); %scatter plot of rr intervals over time
    scatter(time,gatingTrack.missed,35,'filled','red'); %fill in points with red if missed HB
    scatter(time,gatingTrack.early,35,'filled','yellow'); %fill in points with yellow if early trig
    xlabel('Time (s)'); 
    ylabel('RR Interval (ms)');
    if upperLim(1)<2200 && lowerLim(1)>500
        ylim([500 2200]);
    end 
    legend('2\sigma limits','Mean RR','Recon RR','ECG/PG RR-intervals','Outliers', ...
        'Likely Missed HBs','Likely Early Triggers','Location','eastoutside');
    savefig(hFig,'gatingTrack_resp_ecg_waveforms'); % save this figure

    % Histogram (BPM)
    histo = figure; 
    %histo.WindowState = 'maximized'; %full-screen
    subplot(2,1,1);     
    histogram(gatingTrack.bpm); % plot histogram of heart rate (not clean)
    xlabel('BPM'); 
    ylabel('Frequency (counts)'); 
    title('Histogram of Heart Rates');

    subplot(2,1,2); hold on;
    boxplot(gatingTrack.bpm,'Orientation','horizontal','Notch','on', ...
        'Widths',0.8,'OutlierSize',8,'Symbol','*r');
    jitter = (1+(rand(size(gatingTrack.bpm))-0.5)/10);
    scatter(gatingTrack.bpm,ones(size(gatingTrack.bpm)).*jitter,10,'+','k');
    savefig(histo,'histogram_boxplot'); % save this figure

%%%% Output/Save Results %%%%
    outputResults(gatingTrack,gatingDir,filesLoaded)
    dataLoaded = 1;
end 


%% ScanArchives -  ECG, PPG, Respiratory Waveforms/Triggers
%%%% WAVEFORMS %%%%
% Second ECG Lead - Full Waveform. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG2Data')); % does this file exist in fileNames?
if ~isempty(idx) % if we found something..
    name = fileNames{idx}; % grab name from fileNames
    waveform.ecg2 = importdata(name); % put data into waveform structure
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.ecg2 = waveform.ecg2;
end 

% Third ECG Lead - Full Waveform. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG3Data'));
if ~isempty(idx)
    name = fileNames{idx};
    waveform.ecg3 = importdata(name);
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.ecg3 = waveform.ecg3;
end 

% Peripheral Gating - Full Waveform. % If ECG gated, this will be noise.
idx = find(contains(fileNames,'PPGData')); % Obtained via pulse oximetry.
if ~isempty(idx)
    name = fileNames{idx};
    waveform.ppg = importdata(name);
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.ppg = waveform.ppg;
end 

% Respiratory Gating - Full Waveform. Taken from respiratory bellows.
idx = find(contains(fileNames,'RESPData'));
if ~isempty(idx)
    name = fileNames{idx};
    waveform.resp = importdata(name);
    waveform.resp = 4095-waveform.resp;
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.resp = waveform.resp;
end 

%%%% TRIGGERS %%%%
% QRS Triggers from Second ECG Lead. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG2Trig'));  % does this file exist in fileNames?
if ~isempty(idx) % if we found something..
    name = fileNames{idx}; % grab name from fileNames
    trigger.ecg2 = importdata(name)+1; % put data into trigger structure
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.ecg2 = trigger.ecg2;
end 

% QRS Triggers from Third ECG Lead. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG3Trig'));
if ~isempty(idx)
    name = fileNames{idx};
    trigger.ecg3 = importdata(name)+1;
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.ecg3 = trigger.ecg3;
end 

% Peak Signal Trigger from Peripheral Gating.
idx = find(contains(fileNames,'PPGTrig'));
if ~isempty(idx)
    name = fileNames{idx};
    trigger.ppg = importdata(name)+1;
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.ppg = trigger.ppg;
end 

% Triggers from Respiratory Bellows
idx = find(contains(fileNames,'RESPTrig'));
if ~isempty(idx)
    name = fileNames{idx};
    trigger.resp = importdata(name)+1;
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.resp = trigger.resp;
end 
        

%% Plotting (with complete waveforms)
if sum(contains(filesLoaded,'PPG'))==2 % if we have 2 ppg waveform files
    figure; plot(waveform.ppg); % unsure what the time resolution is..
    title('Raw PG Waveform with Triggers');
    hold on; trigPoints = nan(length(waveform.ppg),1); % expand trigger point vector to match ppg vector length
    trigPoints(trigger.ppg) = waveform.ppg(trigger.ppg); % find where the trigger points occur in waveform.ppg
    scatter(1:length(trigPoints),trigPoints,'filled','green'); % plot points of triggers on ppg waveform
    legend('Raw PPG Waveform','Triggers');
    fprintf('\nRaw Waveform Signal:\n');
    SNRppg = snr(waveform.ppg); % get SNR measure from power spectrum analysis
    SNRthresh = -8;
    if SNRppg>SNRthresh % this is an arbitray threshold value, can be changed
        fprintf('   PPG: SIGNAL DETECTED\n');
    else
        fprintf('   PPG: LIKELY NOISE\n');
    end 
    savefig('Raw_PG_waveform_triggers.fig')
end 

if sum(contains(filesLoaded,'RESP'))==2 % if we have 2 resp waveform files
    figure; plot(waveform.resp);
    title('Raw Respiratory Waveform with Triggers');
    hold on; trigPoints = nan(length(waveform.resp),1);
    trigPoints(trigger.resp) = waveform.resp(trigger.resp);
    scatter(1:length(trigPoints),trigPoints,'filled','green');
    legend('Raw Respiratory Waveform','Bellow Triggers');
    SNRresp = snr(waveform.resp);
    SNRthresh = -8;
    if SNRresp>SNRthresh % this is an arbitray threshold value, can be changed
        fprintf('   RESP: SIGNAL DETECTED\n');
    else
        fprintf('   RESP: LIKELY NOISE\n');
    end 
    savefig('Raw_resp_waveform_triggers.fig')
end 

if sum(contains(filesLoaded,'ECG'))==4 % if we have 4 ecg waveform files
    figure; sgtitle('Raw ECG Waveform with Triggers');
    subplot(2,1,1); plot(waveform.ecg2); title('ECG Lead 2');
    hold on; trigPoints = nan(length(waveform.ecg2),1);
    trigPoints(trigger.ecg2) = waveform.ecg2(trigger.ecg2);
    scatter(1:length(trigPoints),trigPoints,'filled','green');
    legend('ECG2','ECG2 Triggers'); hold off
    subplot(2,1,2); plot(waveform.ecg3); title('ECG Lead 3');
    hold on; trigPoints = nan(length(waveform.ecg3),1);
    trigPoints(trigger.ecg3) = waveform.ecg3(trigger.ecg3);
    scatter(1:length(trigPoints),trigPoints,'filled','green');
    legend('ECG3','ECG3 Triggers'); hold off
    SNRecg2 = snr(waveform.ecg2);
    SNRecg3 = snr(waveform.ecg3);
    SNRthresh = -8;
    if SNRecg2>SNRthresh % this is an arbitray threshold value, can be changed
        fprintf('   ECG2: SIGNAL DETECTED\n');
    else
        fprintf('   ECG2: LIKELY NOISE\n');
    end 
    if SNRecg3>SNRthresh % this is an arbitray threshold value, can be changed
        fprintf('   ECG3: SIGNAL DETECTED\n');
    else
        fprintf('   ECG3: LIKELY NOISE\n');
    end
    savefig('Raw_ecg_waveforms_triggers.fig')
end 


%% ECG Track (old gating)
%%%% Load ECG Gating Files %%%% 
ecgIdx = find(contains(fileNames,'ecg_track')); %does this file exist?
if ~isempty(ecgIdx) && ~dataLoaded %if we found something/haven't already loaded
    ecgName = fileNames{ecgIdx}; %grab name from fileNames
    fprintf('  \t A "scan_info.txt" file is required to load this ecg_track file\n');
    [infoFile,infoPath] = uigetfile('scan_info.txt','Find scan_info.txt File');
    infoName = fullfile(infoPath,infoFile);
    fid = fopen(infoName);  % open scan_info.txt

    % Search for Parameters to Rearrange ecg_track Info
    needTR = 1; % flag if we still need to find TR
    needNPROJ = 1; % flag if we still need to find # projections
    needINTER = 1; % flag if we still need to find interleave #
    needORDER = 1; % flag if we still need to find projection order
    looking = needTR + needNPROJ + needINTER + needORDER; %if >0, keep looking
    count = 0;
    while looking                                     
        tline = fgetl(fid); % grab single lines from scan_info                             
        if ischar(tline)                       
            if needTR
                param1 = contains(tline, 'TR'); % find "TR" line  
                if param1             
                    needTR = 0; % turn flag off to stop searching
                    %tr=sscanf(tline, '%*s %*s %s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s');
                    tr = sscanf(tline, '%*s %*s %f');
                    %tr=str2double(tr); % turn from string into number
                end 
            end 
            if needNPROJ
                param2 = contains(tline,'rhuser8') | contains(tline,'rhuser08'); 
                if param2                
                    needNPROJ = 0; 
                    %nproj=sscanf(tline, '%*s %*s %s');
                    nproj=sscanf(tline, '%*s %*s %d');
                    %nproj=str2double(nproj);
                end 
            end 
            if needINTER
                param3 = contains(tline, 'rhuser10'); 
                if param3            
                    needINTER = 0; 
                    inter = sscanf(tline, '%*s %*s %d');
                end 
            end 
            if needORDER
                param4 = strfind(tline, 'PR Order');
                if isfinite(param4)                
                    needORDER = 0; 
                    PR_ORDER = sscanf(tline, ' %*s %*s %s');
                    PR_ORDER = str2double(PR_ORDER);
                end 
            end   
        end
        looking = needTR + needNPROJ + needINTER + needORDER; % update while loop flag
        if count>3000 %break while loop if stuck
            break
        end 
        count = count+1;
        %PR_ORDER = 2;
    end 
    fclose(fid); %from scan info

    % Load ECG Track
    subproj = nproj / inter; % calculate # projections in interleave
    fid = fopen(ecgName, 'rb');
    ecg = fread(fid,'int', 'b');
    fclose(fid) %from ecg track
    
    % Determine projection order, changes depending on PR_ORDER
    acquisition_order = zeros(nproj,1); % array for all projections     
    sub_order = bitreverse(subproj); % bit reverse for the subframes

    if PR_ORDER == 2 
        frame_order = 0:inter-1; % subframes are ordered sequentially
    elseif PR_ORDER == 4
        frame_order = bitrevorder(inter); % subframes are bit reversed
    end

    for pos = 0:nproj-1   
        inter_index = mod(pos,inter); % which interleave am I in?
        inter_pos = floor(pos/inter); % index position w/in interleave
        sub_pos = inter_index + inter_pos; % not sure what this is
        if PR_ORDER == 2
            acquisition_order(pos+1) = subproj*frame_order(floor(pos/subproj)+1) + sub_order(mod(pos,subproj)+1);
        end

        if  PR_ORDER == 4
            acquisition_order(pos+1) = sub_order(mod(sub_pos,subproj)+1) + subproj*frame_order(inter_index+1);
        end
    end

    ecg_sorted = ecg(acquisition_order+1); %amplitude
    TR = tr*4; % for 4point encode
    time = (TR:TR:TR*length(ecg_sorted))'; %time in ms

    % Put data into ecgTrack structure
    ecgTrack.ecg = ecg_sorted; % ecg (sawtooth) waveform
    ecgTrack.time = time; % time at each data point (ms)
    ecgTrack.acq = (1:length(ecg_sorted))'; % data acquisition number
    
    % Match recon RR calculation
    within = ecgTrack.ecg<2000 & ecgTrack.ecg>0; %see gating_lib.cpp
    ecgTrack.recon_rr = 2*median(ecgTrack.ecg(within)); %average RR
    %method above robust to decreased sampling rates of ecg sawtooth fx
    %e.g. full gating tracks vs. regular gating tracks
    
    [rrECG, peaksECG] = getRR(ecgTrack.ecg); % get rr intervals from ecg data
        ecgTrack.rr = nan(length(ecgTrack.ecg),1); 
        ecgTrack.rr(peaksECG) = rrECG; % fill in NaNs with rr intervals
        % account for discretization error (to better match recon rr est.)
        mins = peaksECG+1;
        discShift = round(mean(ecgTrack.ecg(mins)));
        ecgTrack.rr = ecgTrack.rr + discShift;
        ecgTrack.bpm = 60000./ecgTrack.rr;

    % Detect Missed HBs and Early Triggers
    [missedHBIdx,earlyTrigIdx] = findMissedHB(rrECG,peaksECG);
        missedHBs = nan(length(ecgTrack.rr),1); % expand vector to match gate lengths
        earlyTrigs = nan(length(ecgTrack.rr),1); % expand vector to match gate lengths
        missedHBs(missedHBIdx) = ecgTrack.rr(missedHBIdx); % fill in expanded vector with missed HBs
        earlyTrigs(earlyTrigIdx) = ecgTrack.rr(earlyTrigIdx); % fill in expanded vector early triggers
        ecgTrack.missed = missedHBs;
        ecgTrack.early = earlyTrigs;
        ecgTrack.clean = ecgTrack.rr;
        ecgTrack.clean(missedHBIdx) = nan;
        ecgTrack.clean(earlyTrigIdx) = nan; % data which has missed HB and early trigger events removed
    
    % Second Method to DetectOutliers
    outliers = isoutlier(gatingTrack.rr);
    gatingTrack.outliers = nan(size(gate,1),1);
    gatingTrack.outliers(outliers) = gatingTrack.rr(outliers);
    gatingTrack.clean = gatingTrack.rr;
    gatingTrack.clean(outliers) = nan; %remove missed HBs
    
%%%% Plotting and Output %%%%
    ecgTrack.timeres = mean(diff(ecgTrack.time)); %find time resolution by finding most common value
    time = ecgTrack.time/1000; %convert to seconds
    
    % PG/ECG Subfigure
    avg = nanmean(ecgTrack.rr);
    stdv = nanstd(ecgTrack.rr);
    upperLim = ones(length(gatingTrack.ecg),1).*(avg + 2*stdv); %get max clean RR interval
    meanRR = ones(length(gatingTrack.ecg),1).*avg; %get mean RR
    lowerLim = ones(length(gatingTrack.ecg),1).*(avg - 2*stdv); %get min RR
    reconRR = ones(length(gatingTrack.ecg),1).*ecgTrack.recon_rr; %recon-calculated RR

    hFig = figure; hold on; 
    hFig.WindowState = 'maximized'; %full-screen
    title('ECG/PG Triggers') %make second subplot
    plot(time,upperLim,'--','Color','black'); %plot dashed line of max RR
    plot(time,lowerLim,'--','Color','black','HandleVisibility','off'); %plot dashed line of min RR
    plot(time,meanRR,'--','Color','blue'); %plot dashed line of mean RR
    plot(time,reconRR,'--','Color','cyan'); %plot recon RR
    scatter(time,ecgTrack.rr,35,'blue'); %scatter plot of rr intervals over time
    scatter(time,ecgTrack.outliers,35,'filled','blue'); %scatter plot of rr intervals over time
    scatter(time,ecgTrack.missed,35,'filled','red'); %fill in points with red if missed HB
    scatter(time,ecgTrack.early,35,'filled','yellow'); %fill in points with yellow if early trig
    xlabel('Time (s)'); 
    ylabel('RR Interval (ms)');
    legend('2\sigma limits','Mean RR','Recon RR','ECG/PG RR-intervals','Outliers',...
        'Likely Missed HBs','Likely Early Triggers','Location','eastoutside');
    savefig(hFig,'ecgTrack_ecg_waveforms'); % save this figure

    % Histogram (BPM)
    histo = figure; 
    %histo.WindowState = 'maximized'; %full-screen
    subplot(2,1,1);     
    histogram(ecgTrack.bpm); % plot histogram of heart rate (not clean)
    xlabel('BPM'); 
    ylabel('Frequency (counts)'); 
    title('Histogram of Heart Rates');

    subplot(2,1,2); hold on;
    boxplot(ecgTrack.bpm,'Orientation','horizontal','Notch','on', ...
        'Widths',0.8,'OutlierSize',8,'Symbol','*r');
    jitter = (1+(rand(size(ecgTrack.bpm))-0.5)/10);
    scatter(ecgTrack.bpm,ones(size(ecgTrack.bpm)).*jitter,10,'+','k');
    savefig(histo,'histogram_boxplot'); % save this figure

%%%% Output/Save Results %%%%
    outputResults(ecgTrack,gatingDir,filesLoaded)
    dataLoaded = 1;
end


%% Catch
if isempty(filesLoaded) % if we didn't find anything to load..
    disp('No gating tracks were found');
    return % kick us out
end 

%% Save Master Data Structure
save ecgData.mat ecgData

%% Clear unnecessary variables
clear ans acquisition_order count dataLoaded earlyTrigs ecg ecg_sorted
clear ecgIdx ecgName ecgTrack fid fil* frame_order fullDir gatingDir
clear idx info* is* inter* keepFiles looking maxRR meanRR minRR missedHBs
clear need* nproj numPoints outliers param* peaksECG pos PR_ORDER reconRR
clear rrECG sub* time tline tr TR within discShift mins avg stdv
clear gat* loadFullTrack name peaksFull rrFull SNR* trig* waveform
clear earlyTrigIdx missedHBIdx lowerLim upperLim jitter

diary off % turn off diary, stop writing to ecgInformation.txt

%% Ancillary functions
function outputResults(datastruct,gatingDir,filesLoaded)
    %%%%% Command Line Output %%%%%
    fprintf('ECG/PG Gating Information:\n');
    fprintf('    Time: %s\n',datetime('now'));
    fprintf('    Directory: %s\n',gatingDir);
    for i=1:length(filesLoaded)
        fprintf('    Files loaded: %s\n',filesLoaded{i});
    end
    A1 = {'ECG/PG Gating Information';'Time';'Directory'};
    B1 = {'';datestr(datetime('now'));gatingDir};
    for i=1:length(filesLoaded)
        A1{end+1} = 'Files loaded';
        B1{end+1} = filesLoaded{i};
    end
    
    fprintf('\n');
    fprintf('Temporal\n');
    fprintf('    ECG/Resp. time resolution (TR): %7.2f ms.\n',datastruct.timeres);
    fprintf('    Scan duration: %7.2f s.\n\n',datastruct.time(end)/1000);
    A2 = {'';'Temporal'; 'ECG/Resp. time resolution (TR, ms)'; 'Scan duration (s)'};
    B2 = {'';''; datastruct.timeres; datastruct.time(end)/1000};
        
    fprintf('RR statistics\n');
    fprintf('    Recon RR interval: %7.2f ms.\n',datastruct.recon_rr);
    fprintf('    Mean RR interval: %7.2f ms.\n',nanmean(datastruct.rr));
    fprintf('    Median RR interval: %7.2f ms.\n',nanmedian(datastruct.rr));
    fprintf('    Range: [%d - %d] ms.\n\n',min(datastruct.rr),max(datastruct.rr));
    A3 = {'';'RR statistics';'Recon RR interval (ms)';'Mean RR interval (ms)'; ...
        'Median RR interval (ms)';'Max RR (ms)';'Min RR (ms)'};
    B3 = {'';'';datastruct.recon_rr;nanmean(datastruct.rr); ...
        nanmedian(datastruct.rr);min(datastruct.rr);max(datastruct.rr)};
    
    fprintf('Heart Rate Statistics\n');
    fprintf('    Mean heart rate: %7.2f bpm.\n',nanmean(datastruct.bpm));
    fprintf('    Median heart rate: %7.2f bpm.\n',nanmedian(datastruct.bpm));
    fprintf('    Standard deviation: %7.2f bpm.\n',nanstd(datastruct.bpm));
    fprintf('    Coefficient variation: %7.2f%%.\n\n',(nanstd(datastruct.bpm)/nanmean(datastruct.bpm))*100);
    A4 = {'';'Heart Rate Statistics';'Mean heart rate (bpm)';'Median heart rate (bpm)'; ...
        'Standard deviation (bpm)';'Coefficient variation'};
    B4 = {'';'';nanmean(datastruct.bpm);nanmedian(datastruct.bpm); ...
        nanstd(datastruct.bpm);(nanstd(datastruct.bpm)/nanmean(datastruct.bpm))*100};
    
    numMissedFull = sum(~isnan(datastruct.missed));
    numEarlyFull = sum(~isnan(datastruct.early));
    numOutliers = sum(~isnan(datastruct.outliers));
    totalRR = sum(~isnan(datastruct.rr));
    meanHR = nanmean(datastruct.bpm);
    upperHR = meanHR + 5;
    lowerHR = meanHR - 5;
    withinLims = nansum(datastruct.bpm>lowerHR & datastruct.bpm<upperHR);
    RRs = datastruct.rr(~isnan(datastruct.rr));
    dRRs = diff(RRs);
    fprintf('Error Estimation and Heart Rate Variability\n');
    fprintf('    Total recorded heartbeats: %d beats.\n',totalRR);
    fprintf('    Estimated missed heartbeats: %d beats (%2.2f%%).\n',numMissedFull,100*(numMissedFull/totalRR));
    fprintf('    Estimated early triggers: %d beats (%2.2f%%).\n',numEarlyFull,100*(numEarlyFull/totalRR));
    fprintf('    Total error: %2.2f%%.\n',100*((numEarlyFull+numMissedFull)/totalRR));
    fprintf('    Total outliers: %2.2f%%.\n',100*(numOutliers/totalRR));
    fprintf('    Percent data within mean HR +/- 5 BPM: %2.2f%%.\n',100*(withinLims/totalRR));
    fprintf('    SDNN: %2.1f.\n',std(RRs));
    fprintf('    RMSSD: %2.1f.\n',sqrt(mean(dRRs.^2)));
    fprintf('    pNN50: %2.2f%%.\n',sum(abs(dRRs)>50)/length(dRRs)*100);
    fprintf('   For info on SDNN, RMSSD, pNN50: pubmed.ncbi.nlm.nih.gov/29034226 \n\n');
    A5 = {'';'Error Estimation';'Total recorded heartbeats';'Estimated missed heartbeats'; ...
        'Estimated early triggers';'Total error (%)';'Total outliers (%)'; ...
        'Data within mean HR +/- 5 BPM (%)';'SDNN';'RMSSD';'pNN50'};
    B5 = {'';'';totalRR;numMissedFull;numEarlyFull;100*((numEarlyFull+numMissedFull)/totalRR); ...
        100*(numOutliers/totalRR); 100*(withinLims/totalRR); nanstd(datastruct.rr); ...
        sqrt(mean(dRRs.^2)); sum(abs(dRRs)>50)/length(dRRs)*100};
    
    within_rr = datastruct.ecg < datastruct.recon_rr;
    ecg_filtered = datastruct.ecg(within_rr);
    sum_within = size(ecg_filtered,1);
    sum_total = size(datastruct.ecg,1);
    lowBPM = datastruct.recon_rr > 2000;
    highBPM = datastruct.recon_rr < 500;
    fprintf('PCVIPR Recon Output (may not match perfectly due to TR rounding)\n') 
    fprintf('    Median RR is %d ms.\n',datastruct.recon_rr);
    fprintf('    Expected heart rate is %7.4f bpm.\n',60000/datastruct.recon_rr);
    fprintf('    Values within expected RR = %7.4f %%.\n\n',100*(sum_within/sum_total));

    A6 = {'';'PCVIPR Recon Output (may not match exactly due to TR rounding)'; ...
        'Median RR (ms)'; 'Expected heart rate (bpm)'; 'Values within expected RR (%)'};
    B6 = {'';'';datastruct.recon_rr;60000/datastruct.recon_rr; ...
        100*(sum_within/sum_total)};
    
    lm = fitlm(datastruct.time/60000,datastruct.bpm);
    coefs = lm.Coefficients(2,:);
    fprintf('Regression and Warnings\n');
    fprintf('    RR linear fit slope (BPM/min) = %4.4f\n',coefs.Estimate(1));
    fprintf('    p-value = %0.5f\n',coefs.pValue(1));
    if lowBPM
        fprintf('WARNING -- MEDIAN RR > 2000 (HR < 30 BPM)\n');
    elseif highBPM
        fprintf('WARNING -- MEDIAN RR < 500 (HR > 120 BPM)\n');
    else
        fprintf('NO WARNINGS');
    end 
    A7 = {''; 'RR linear fit slope (BPM/time)'; 'p-value'; ...
        'HR < 30 bpm?'; 'HR > 120 bpm?'};
    B7 = {''; coefs.Estimate(1); coefs.pValue(1); lowBPM; highBPM};
    %%%%% Excel Saving %%%%%
    Parameter = [A1;A2;A3;A4;A5;A6;A7];
    Value = [B1;B2;B3;B4;B5;B6;B7];
    writetable(table(Parameter,Value),'gatingStats.xlsx');
end 


function [RRs,peaks] = getRR(ecg)
    df = diff(ecg); % differentiate ecg with respect to time
    dt = mode(df); % get time resolution of ecg by finding most common value in derivative
    peaks = find(df<-12*dt); % find negative peaks, 12*dt is an arbitray value
    
    %%%% OPTION 1 - Take Peak ECG Trigger Time %%%%
    RRs = ecg(peaks); % get the actual rr value at the peak
    % ignores trigger delays (good for regular gating tracks)
    
    %%%% OPTION 2 - Account for ECG Start Time Offset %%%%
    %starts = peaks + 1; %next index is start of new RR
    %peaks(1) = []; %kill first trigger (don't know start)
    %starts(end) = []; %kill last start of RR (don't know when it ends)
    %for i=1:length(peaks)
    %    RRs(i)=ecg(peaks(i))-ecg(starts(i)); %e.g. 1289-19 
    %end
    
    %%%% OPTION 3 - Use Acquisition Times to Compute RR %%%%
    %starts = peaks + 1; %next index is start of new RR
    %peaks(1) = []; %kill first trigger (don't know start)
    %starts(end) = []; %kill last start of RR (don't know when it ends)
    %for i=1:length(peaks)
    %    RRs(i)=time(peaks(i))-time(starts(i));
    %end  
end 


function [missedHBidx,earlyTrigIdx] = findMissedHB(rr,peaks)
    % Find missed heartbeats
    missedHBidx = []; % create missed hearbeat array
    earlyTrigIdx = []; % create early trigger array
    if length(rr)>100 % if we have more than 100 rr intervals
        for L = 1:length(rr)
            if L<25 % starting on the left side
                RR_window = linspace(1,50,50)'; % create sliding window over first 50 points
            elseif L>length(rr)-26 % in the middle, move the sliding window
                RR_window = linspace(length(rr)-50,length(rr)-1,50)'; % create sliding window of width 50
            else % if we're at the end of the rr vector
                RR_window = linspace(L-24,L+25,50)'; % create sliding window over last 50 points
            end
            sliding_RR = median(rr(RR_window)); % grab median rr from sliding window
            scale = rr(L)/sliding_RR; % find how the current rr interval compares to the median
            if scale > 1.5 % if scale is greater than 1.66 times the median rr, it is likely a missed heartbeat. 1.66 is arbitrary
                missedHBidx = [missedHBidx,peaks(L)]; % find index in the peak array and place in running array
            end 
            if scale <0.50 % if scale is less than 0.66 times the median rr, it is likely an early trigger. 0.66 is arbitrary
                earlyTrigIdx = [earlyTrigIdx,peaks(L)]; % find index in the peak array and place in running array
            end 
        end
    elseif length(rr)>40 % if we have less than 100 rr intervals but more than 40
        for L = 1:length(rr)
            if L<10
                RR_window = linspace(1,20,20)';
            elseif L>length(rr)-11
                RR_window = linspace(length(rr)-20,length(rr)-1,20)';
            else
                RR_window = linspace(L-9,L+10,20)';
            end
            sliding_RR = median(rr(RR_window));
            scale = rr(L)/sliding_RR;
            if scale > 1.5
                missedHBidx = [missedHBidx,peaks(L)];
            end 
            if scale <0.5
                earlyTrigIdx = [earlyTrigIdx,peaks(L)];
            end 
        end
    else % if we have a very small sample size (unreliable)
        for L = 1:length(rr)
            scale = rr(L)/median(rr);
            if scale > 1.5
                missedHBidx = [missedHBidx,peaks(L)];
            end 
            if scale <0.5
                earlyTrigIdx = [earlyTrigIdx,peaks(L)];
            end 
        end 
    end 
end 


function B = bitreverse(N)
    N_bits = ceil(log2(N));
    B_tmp = zeros(1, N);
    for j = 0 : N-1
         tmp = j;
         for i = 0 : N_bits
             B_tmp(j+1) = B_tmp(j+1) + 2^(N_bits - i - 1) * mod(tmp,2);
             tmp = floor(tmp/2);
         end
    end

    B(1) = B_tmp(1);
    curr_high = 0;

    for j = 2 : N
         tmp = N+1;
         diff_elt = 1;
         for i = 2 : N
             if B_tmp(i) > curr_high
                 diff = B_tmp(i) - curr_high;
                 if diff < tmp
                    diff_elt = i;
                    tmp = diff;
                 end
             end
         end
         curr_high = B_tmp(diff_elt);
         B(diff_elt) = j-1;
    end
end 