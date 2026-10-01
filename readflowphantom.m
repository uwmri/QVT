function readflowphantom(dirname)
% function [ waveformScaled, waveformExtremes, fitlimits ] = readflowphantom(dirname)
% function [ waveformScaled, waveformExtremes, fitlimits, idxminmax, treedist ] = readflowphantom(dirname)

nfiles = 40; % same as number of centerline points
% nfiles = 7; % first few for testing purpoes
ntimes = 40;
% tqfirst = nfiles/4;
% thalf = 2*tqfirst;
% tqlast = 3*tqfirst;
tres = 25;
wavecells = cell(1,nfiles);
for jf = 1:nfiles
    fn = strcat(dirname,filesep,sprintf('Waveform%d.csv',jf));
    if ~exist(fn,'file')
        fprintf('No file %s\n',fn)
        return
    end
    wavecells{jf} = readtable(fn);
    dim = size(wavecells{jf});
    if dim(1) ~= ntimes && dim(2) ~= 2
        fprintf('Incorrect dimensions in %s\n',fn)
    end
    % for k = 1:ntimes
    %     if abs(wavecells{jf}.Time(k) - k*tres) > 1e-5
    %         fprintf('Incorrect time resolution in %s\n',fn)
    %     end
    % end
end

fn = strcat(dirname,filesep,'Spline_Distance.csv');
if ~exist(fn,'file')
    fprintf('no file %s\n',fn)
    return
end

Tdist = readtable(fn);
treedist = Tdist.Distance;

waveforms = zeros(nfiles,ntimes);
for k = 1:nfiles
    a = table2array(wavecells{k});
    waveforms(k,:) = a(:,2)';
end

processFlow(waveforms, treedist(1:nfiles), tres);

% [ par,~,~ ] = minandsavestats(waveforms,treedist,ntimes,tres);
% 
% fprintf('esitmated pulse wave velocity is %.3f m/s\n',par(end))
% 
% % need to smooth the waveforms
% % that's not here yet
% 
% waveformScaled = zeros(nfiles,ntimes);
% waveformExtremes = zeros(nfiles,2);
% idxminmax = zeros(nfiles,2);
% fitlimits = zeros(nfiles,2);
% wrpflag = 0;
% for kk = 1:nfiles
%     [ wfmax, idxmax ] = max(waveforms(kk,:));
%     waveformExtremes(kk,2) = wfmax;
%     idxminmax(kk,2) = idxmax;
%     [ wfmin, idxmin ] = min(waveforms(kk,:));
%     waveformExtremes(kk,1) = wfmin;
%     idxminmax(kk,1) = idxmin;
%     waveformScaled(kk,:) = (waveforms(kk,:) - wfmin)/(wfmax - wfmin);
% 
%     if mod(wrpflag,2) == 0
%         if idxmax > tqlast
%             wrpflag = wrpflag + 1;
%         end
%     else
%         if idxmax > tqfirst && idxmax < thalf
%             wrpflag = wrpflag + 1;
%         end
%     end
% 
%     % anything I do here will skew the slope?
%     idxveclow = find(waveformScaled(kk,:) < 0.25);
%     idxvechi = find(waveformScaled(kk,:) > 0.75);
%     % start is largest idxveclow less than idxmax
%     % end is smallest idxvechigh greater than idmin
%     v = find(idxveclow < idxmax);
%     if isempty(v)
%         idxstart = idxveclow(end);
%     else
%         idxstart = v(end);
%     end
%     v = find(idxvechi > idxmin);
%     if isempty(v)
%         idxend = idxvechi(1);
%     else
%         idxend = v(1);
%     end 
%     fitlimits(kk,1) = idxstart;
%     fitlimits(kk,2) = idxend;
%     % if idxstart < idxend
%     %     idxrange = idxstart:idxend;
%     %     x = tres*idxrange;
%     % else
%     %     idxwraprange = idxstart:nfiles;
%     %     idxrange = [ idxwraprange 1:idxend ];
%     %     irngtmp = [ (idxwraprange - nfiles) 1:idxend ];
%     %     x = tres*irngtmp;
%     % end
%     % x = x + tres*nfiles*floor(wrpflag/2 + 0.001);
%     % y = waveformScaled(kk,idxrange);
%     % mdl = fitlm(y,x);
%     % mdl.Coefficients
% end
% 
