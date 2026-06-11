function [ params,exitflag,output ] = minandsavestats(treeflow,treedist,...
    nframes,timeres)

nptstree = size(treeflow,1);
treeflownorm = treeflow;
flowmeans = zeros(nptstree,1);
flowstddevs = flowmeans;

% need to normalize flow
for ii = 1:nptstree
    v = treeflow(ii,:);
    [sig, mu] = std(v);
    treeflownorm(ii,:) = (v - mu)/sig;
    flowmeans(ii) = mu;
    flowstddevs(ii) = sig;
end

% initguess = [ zeros(1,nframes), 5 ]; % initial is zero waveform, pwv = 5 m/s
initguess = [ zeros(1,6), 5 ];
wavemodel = @PWVsinecosinemodel;
fitwts = ones(size(treedist));

% need to "recalculate" costfun on each call because the values, not
% the variables, become part of costfun ???
costfun = @(inParams)PWVchisq(inParams,treedist,treeflownorm,...
    timeres,fitwts,wavemodel);
% timeres in ms; mm/ms should be the same as m/s
options = optimset('Display','iter', 'TolCon', 1e-7, 'TolX', 1e-7, 'TolFun', 1e-7,'DiffMinChange', 1e-3);
[params,exitflag,output] = fminunc(costfun, initguess, options);

d.treeflow = treeflow;
d.treedist = treedist;
d.treeflownorm = treeflownorm;
d.lentree = nptstree;
d.flowmeans = flowmeans;
d.flowdevs = flowstddevs;
d.par = params;
d.exflag = exitflag;
d.out = output;
d.model = wavemodel;
d.ntimes = size(treeflow,2);
d.T_heart = timeres*d.ntimes;
d.timeres = timeres;
d.wts = fitwts;

fnmat = strcat('H:\LinuxHome\sandbox\PWVfits\fitdata', ...
    string(datetime('now','Format','yyyyMMdd_HHmmssSSS')),'.mat');
save(fnmat,'d')