function plotgenerator(d,ysep,nrow,ncol,idxstart)

% d is a structure

colorcels = { 'b-', 'r-', 'g-', 'm-', 'k-' };
ncolors = length(colorcels);

modeltimes = 200;
waveModel = d.model(d.par,d.treedist,d.T_heart,modeltimes);
xvals = (0:(d.ntimes-1))*d.T_heart/d.ntimes;
xvalmodel = (0:(modeltimes-1))*d.T_heart/modeltimes;
idxlast = idxstart + nrow*ncol - 1;
npts = length(d.flowmeans);
if idxlast > npts
    idxlast = npts;
end
tiledlayout("horizontal")
nexttile
kmax = idxlast - idxstart;
for k = 0:kmax
    idxcol = 1 + mod(k,ncolors); % this is the color to use
    ii = k + idxstart;
    idxrow = 1 + mod(k,nrow);
    ydel = ysep*(nrow/2-idxrow);
    plot(xvalmodel,waveModel(ii,:)*d.flowdevs(ii) + d.flowmeans(ii) + ydel,...
        strcat(colorcels{idxcol},'-'),...
        xvals,d.treeflow(ii,:) + ydel,colorcels{idxcol}, ...
        'LineWidth',2)
    if idxrow == nrow
        hold off
        if k < kmax
            nexttile
        end
    else
        hold on
    end
end