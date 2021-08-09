%% plot predictive densities for figures and charts

clear
close all
fclose all;

%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% setup
resultsdir = '~/jam/lager/var2021-matfiles/baseline';

startSam   = [];
datalabel  = 'fredMD16-2021-04';
modellabel = strcat(datalabel, '-censoredYields');
p          = 12; % number of lags

% modellabel = strcat(modellabel, '-covidRobustPrior');

fontsize = 14;
wrap = [];
initwrap

doConnectDots = false;

datestrfmt = 'mmmm yyyy';

%#ok<*UNRCH>

%% select plots


corelist = [1 5 6 12];

% ylist = corelist; % ylist(~ismember(ylist,corelist));

ylist = 6;

ylist    = 1 : 16;

%% load matfiles (oos runs)

do2020 = '';


prettylabelCONST  = 'CONST';
mfilename         = sprintf('%s-CONST%s-p%d.mat', modellabel, do2020, p);
matCONST          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelCONSTsansCOVID  = 'CONST ex COVID';
mfilename                  = sprintf('%s-CONSTsansCOVID%s-p%d.mat', modellabel, do2020, p);
matCONSTsansCOVID          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSV  = 'SV';
mfilename      = sprintf('%s-SV%s-p%d.mat', modellabel, do2020, p);
matSV          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

% prettylabelSVar1  = 'SVar1';
% mfilename      = sprintf('%s-SVar1%s-p%d.mat', modellabel, do2020, p);
% matSVar1          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVtdof = 'SV-t(5)';
mfilename         = sprintf('%s-SVt-dof5%s-p%d.mat', modellabel, do2020, p);
matSVtdof         = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVt = 'SV-t';
mfilename      = sprintf('%s-SVt%s-p%d.mat', modellabel, do2020, p);
matSVt         = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVO   = 'SVO';
mfilename        = sprintf('%s-SVOmax20%s-p%d.mat', modellabel, do2020, p);
matSVO           = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVoutMiss  = 'SV-OutMiss';
Ofactor               = 5;
mfilename             = sprintf('%s-SVnanO%d%s-p%d.mat', modellabel, Ofactor, do2020, p);
matSVoutMiss          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVdummy  = 'SV-Dummies';
dummyPrecision      = 0;
mfilename           = sprintf('%s-SV-COVIDEACH-precision%d-p%d.mat', modellabel, dummyPrecision, p);
matSVdummy          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVOtdof     = 'SVO-t(9)';
mfilename           = sprintf('%s-SVOtmax20-dof9%s-p%d.mat', modellabel, do2020, p);
matSVOtdof             = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVOt     = 'SVO-t';
mfilename           = sprintf('%s-SVOtmax20%s-p%d.mat', modellabel, do2020, p);
matSVOt             = matfile(fullfile(resultsdir, mfilename), 'writable', false);

%% grab a few things
N              = matCONST.N;


ydates         = matCONST.ydates;
fcstNhorizons  = matCONST.fcstNhorizons;

Tjumpoffs      = matCONST.Tjumpoffs;


np             = 12;
ncode          = matCONST.ncode;
tcode          = matCONST.tcode;

setQuantiles   = matCONST.setQuantiles;

Ylabels                             = fredMDprettylabel(ncode);
Ylabelslong                         = Ylabels;
Ylabelslong(ismember(tcode, [2 5])) = strcat(Ylabelslong(ismember(tcode, [2 5])), ' (APR growth)');

%% patch in data
data = matCONST.data;
if ~isequal(data, matSVO.data)
    error houston
end
if ~isequal(data, matSVt.data)
    error houston
end

if ~isequal(data, matSVoutMiss.data)
    error houston
end

if ~isequal(data, matSV.data)
    error houston
end

% if ~isequal(data, matSVar1.data)
%     error houston
% end

if ~isequal(data, matSVdummy.data)
    error houston
end

if ~isequal(data, matSVOt.data)
    error houston
end

Tdata = length(data);

%% construct fcstdates
yearStart = year(ydates(1));
yearStop  = year(ydates(end)) + 9;

if matSV.doQuarterly
    fcstdates = genrQdates(yearStart, yearStop);
else
    fcstdates = genrMdates(yearStart, yearStop, 1);
end
ndx = find(fcstdates == ydates(1));
if isempty(ndx)
    error houston
end
fcstdates = fcstdates(ndx:end);

ndxCI = find(ismember(matCONST.setQuantiles, normcdf([-1 1]) * 100));


%% collect densities
thismat = matCONST;
fcstTailsCONST  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidCONST    = thismat.fcstYmedian;

thismat      = matSV;
fcstTailsSV  = NaN(size(fcstTailsCONST));
fcstMidSV    = NaN(size(fcstMidCONST));
ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSV(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSV(:,:,ndx)      = thismat.fcstYmedian;

% thismat         = matSVar1;
% fcstTailsSVAR1  = NaN(size(fcstTailsCONST));
% fcstMidSVAR1    = NaN(size(fcstMidCONST));
% ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
% fcstTailsSVAR1(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
% fcstMidSVAR1(:,:,ndx)      = thismat.fcstYmedian;

% thismat = matSVO;
% fcstTailsSVO  = NaN(size(fcstTailsCONST));
% fcstMidSVO    = NaN(size(fcstMidCONST));
% ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
% fcstTailsSVO(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
% fcstMidSVO(:,:,ndx)      = thismat.fcstYmedian;
%
%
% thismat = matSVt;
% fcstTailsSVt  = NaN(size(fcstTailsCONST));
% fcstMidSVt    = NaN(size(fcstMidCONST));
% ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
% fcstTailsSVt(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
% fcstMidSVt(:,:,ndx)      = thismat.fcstYmedian;
%
% thismat = matSVtdof;
% fcstTailsSVtdof  = NaN(size(fcstTailsCONST));
% fcstMidSVtdof    = NaN(size(fcstMidCONST));
% ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
% fcstTailsSVtdof(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
% fcstMidSVtdof(:,:,ndx)      = thismat.fcstYmedian;
%
% thismat = matSVOt;
% fcstTailsSVOt  = NaN(size(fcstTailsCONST));
% fcstMidSVOt    = NaN(size(fcstMidCONST));
% ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
% fcstTailsSVOt(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
% fcstMidSVOt(:,:,ndx)      = thismat.fcstYmedian;
%
% thismat = matSVOtdof;
% fcstTailsSVOtdof  = NaN(size(fcstTailsCONST));
% fcstMidSVOtdof    = NaN(size(fcstMidCONST));
% ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
% fcstTailsSVOtdof(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
% fcstMidSVOtdof(:,:,ndx)      = thismat.fcstYmedian;
%
%
% thismat = matSVoutMiss;
% fcstTailsSVoutMiss  = NaN(size(fcstTailsCONST));
% fcstMidSVoutMiss    = NaN(size(fcstMidCONST));
% ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
% fcstTailsSVoutMiss(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
% fcstMidSVoutMiss(:,:,ndx)      = thismat.fcstYmedian;
%
% outliermissYmid   = NaN(length(ydates), N, length(Tjumpoffs));
% outliermissYtails = NaN(length(ydates), N, length(ndxCI), length(Tjumpoffs));
% outliermissYmid(:,:,ndx)     = matSVoutMiss.drawsYmid;
% outliermissYtails(:,:,:,ndx) = matSVoutMiss.drawsYtails(:,:,ndxCI,:);
%
%

thismat                      = matSVdummy;
fcstTailsSVdummy             = NaN(size(fcstTailsCONST));
fcstMidSVdummy               = NaN(size(fcstMidCONST));
ndx                          = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSVdummy(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSVdummy(:,:,ndx)      = thismat.fcstYmedian;

thismat                             = matCONSTsansCOVID;
fcstTailsCONSTsansCOVID             = NaN(size(fcstTailsCONST));
fcstMidCONSTsansCOVID               = NaN(size(fcstMidCONST));
ndx                                 = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsCONSTsansCOVID(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidCONSTsansCOVID(:,:,ndx)      = thismat.fcstYmedian;



%% select dates to be shown
ndxT  = find(ydates(Tjumpoffs) == datenum(2020,4,1));


darkgreen = [0 .75 0];
black     = [0 0 0];


%% define color green
green = [0 .5 0];

%% plot CONST vs CONSTexCOVID 
close all
thisPlotLabel = 'predictiveDensityChartCONSTvsCONSTexCOVID';
for n = ylist
    
    
    jumpoff   = Tjumpoffs(ndxT);
    thisdateT = ydates(jumpoff);
    
    % compare predictive densities (jumpoff)
    jj                = find(fcstdates == ydates(jumpoff));
    %     allthesedates     = fcstdates(jj-np:jj+fcstNhorizons);
    theseFcstdates    = fcstdates(jj+(1:fcstNhorizons));
    %     theseJumpoffdates = fcstdates(jj+(-np:0));
    %     if ~isequal(theseJumpoffdates, ydates(jj+(-np:0)))
    %         error houston
    %     end
    
    theseRealized = NaN(fcstNhorizons, 1);
    ndx = jj+(1:fcstNhorizons);
    ndx = ndx(ndx <= Tdata);
    theseRealized(1 : length(ndx)) = data(ndx,n);
    %     theseFcstdatesCumJumpoff = [theseJumpoffdates(end); theseFcstdates];
    
    
    hanni = NaN(2,1);
    
    
    % plot
    thisfig   = figure;
    set(thisfig,'defaultLegendAutoUpdate','off');
    ax = gca;
    set(ax, 'fontsize', fontsize)
    hold on
    hanni(1)   = plotCI(fcstMidCONST(n,:,ndxT), squeeze(fcstTailsCONST(n,:,:,ndxT)), ...
        theseFcstdates, [], 'k-', 'linewidth', 3);
    
    hanni(2)   = plot(theseFcstdates,fcstMidCONSTsansCOVID(n,:,ndxT), '-', 'color', green, 'linewidth', 3);
    plot(theseFcstdates,squeeze(fcstTailsCONSTsansCOVID(n,:,:,ndxT)), ':', 'color', green, 'linewidth', 2);
    
    %         hanni(3)   = plot(theseFcstdates,fcstMidSVO(n,:,ndxT), 'b-.', 'linewidth', 3);
    %         plot(theseFcstdates,squeeze(fcstTailsSVO(n,:,:,ndxT)), 'b-.', 'linewidth', 2);
    
    %     plot(theseJumpoffdates, data(jj+(-np:0),n), '-', 'color', darkgreen, 'linewidth', 2);
    %     if doConnectDots
    %         plot([theseJumpoffdates(end) theseFcstdates(1)], [data(jj,n) fcstMidCONST(n,1,ndxT)], ':', 'color', black, 'linewidth', 2);
    %     end
    
    if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
    
    xlim(theseFcstdates([1 end]))
    xticks = theseFcstdates(1 : 6 : end);
    datetick('x', 'yyyy:mm', 'keepticks', 'keeplimits')
    ax.XAxis.MinorTick       = 'on';
    
    wrapthisfigure(thisfig, sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
    
    ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
    wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLE', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
    
    %     hl = legend(hanni(1:2), prettylabelCONST, prettylabelCONSTsansCOVID, ...
    %         'location', 'best');
    hl = legend(hanni(1:2), 'standard BVAR per April 2020', 'parameters estimated prior March 2020', ...
        'location', 'south');
    wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
    
    delete(ht)
    wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
    
end

%% plot CONST vs CONSTexCOVID vs SV
close all
thisPlotLabel = 'predictiveDensityChartCONSTvsCONSTexCOVIDvsSV';
for n = ylist
    
    
    jumpoff   = Tjumpoffs(ndxT);
    thisdateT = ydates(jumpoff);
    
    % compare predictive densities (jumpoff)
    jj                = find(fcstdates == ydates(jumpoff));
    %     allthesedates     = fcstdates(jj-np:jj+fcstNhorizons);
    theseFcstdates    = fcstdates(jj+(1:fcstNhorizons));
    %     theseJumpoffdates = fcstdates(jj+(-np:0));
    %     if ~isequal(theseJumpoffdates, ydates(jj+(-np:0)))
    %         error houston
    %     end
    
    theseRealized = NaN(fcstNhorizons, 1);
    ndx = jj+(1:fcstNhorizons);
    ndx = ndx(ndx <= Tdata);
    theseRealized(1 : length(ndx)) = data(ndx,n);
    %     theseFcstdatesCumJumpoff = [theseJumpoffdates(end); theseFcstdates];
    
    
    hanni = NaN(3,1);
    
    
    % plot
    thisfig   = figure;
    set(thisfig,'defaultLegendAutoUpdate','off');
    ax = gca;
    set(ax, 'fontsize', fontsize)
    hold on
    hanni(1)   = plotCI(fcstMidCONST(n,:,ndxT), squeeze(fcstTailsCONST(n,:,:,ndxT)), ...
        theseFcstdates, [], 'k-', 'linewidth', 3);
    
    hanni(2)   = plot(theseFcstdates,fcstMidCONSTsansCOVID(n,:,ndxT), '-', 'color', green, 'linewidth', 3);
    plot(theseFcstdates,squeeze(fcstTailsCONSTsansCOVID(n,:,:,ndxT)), ':', 'color', green, 'linewidth', 2);
    
    hanni(3)   = plot(theseFcstdates,fcstMidSV(n,:,ndxT), 'r-', 'linewidth', 3);
    plot(theseFcstdates,squeeze(fcstTailsSV(n,:,:,ndxT)), 'r--', 'linewidth', 3);
    
%     hanni(4)   = plot(theseFcstdates,fcstMidSVAR1(n,:,ndxT), 'b-.', 'linewidth', 3);
%     plot(theseFcstdates,squeeze(fcstTailsSVAR1(n,:,:,ndxT)), 'b-.', 'linewidth', 3);
     
    %     plot(theseJumpoffdates, data(jj+(-np:0),n), '-', 'color', darkgreen, 'linewidth', 2);
    %     if doConnectDots
    %         plot([theseJumpoffdates(end) theseFcstdates(1)], [data(jj,n) fcstMidCONST(n,1,ndxT)], ':', 'color', black, 'linewidth', 2);
    %     end
    
    if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
    
    xlim(theseFcstdates([1 end]))
    xticks = theseFcstdates(1 : 6 : end);
    datetick('x', 'yyyy:mm', 'keepticks', 'keeplimits')
    ax.XAxis.MinorTick       = 'on';
    
    wrapthisfigure(thisfig, sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
    
    ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
    wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLE', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
    
    %     hl = legend(hanni(1:2), prettylabelCONST, prettylabelCONSTsansCOVID, ...
    %         'location', 'best');
    hl = legend(hanni, 'standard BVAR per April 2020', 'parameters estimated prior March 2020', 'SV-RW per April 2020', ...
        'location', 'south');
    wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], false)
    
    delete(ht)
    wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
    
    
end





%% finish
dockAllFigures
finishwrap
finishscript

function tabulatePanel(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, prettylabel1, prettylabel2, prettylabel3, fcstMid1, fcstMid2, fcstMid3, fcstTails1, fcstTails2, fcstTails3, ncode, ndxCI, setQuantiles) %#ok<DEFNU>

if isempty(wrap)
    return
end

tabname = sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28));

datalabels = {sprintf('%s median', prettylabel1), ...
    sprintf('%s %4.2f%% quantile', prettylabel1, setQuantiles(1,ndxCI(1))), ...
    sprintf('%s %4.2f%% quantile', prettylabel1, setQuantiles(1,ndxCI(2))), ...
    sprintf('%s median', prettylabel2), ...
    sprintf('%s %4.2f%% quantile', prettylabel2, setQuantiles(1,ndxCI(1))), ...
    sprintf('%s %4.2f%% quantile', prettylabel2, setQuantiles(1,ndxCI(2))), ...
    sprintf('%s median', prettylabel3), ...
    sprintf('%s %4.2f%% quantile', prettylabel3, setQuantiles(1,ndxCI(1))), ...
    sprintf('%s %4.2f%% quantile', prettylabel3, setQuantiles(1,ndxCI(2)))};

data4table = [transpose(fcstMid1(n,:,ndxT)), squeeze(fcstTails1(n,:,:,ndxT)), ...
    transpose(fcstMid2(n,:,ndxT)), squeeze(fcstTails2(n,:,:,ndxT)), ...
    transpose(fcstMid3(n,:,ndxT)), squeeze(fcstTails3(n,:,:,ndxT))];


% write into csv
writedatatable(wrap, tabname, theseFcstdates, data4table, datalabels, 'yyyy:mm');

% tabulate in tex
tabdir = wrap.dir;
fid = fopen(fullfile(tabdir, strcat(tabname, '.tex')), 'wt');
fprintf(fid, '\\begin{small}\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, size(data4table,2)));
fprintf(fid, '\\toprule\n');
fprintf(fid, ' & \\multicolumn{3}{c}{%s}', prettylabel1);
fprintf(fid, ' & \\multicolumn{3}{c}{%s}', prettylabel2);
fprintf(fid, ' & \\multicolumn{3}{c}{%s}', prettylabel3);
fprintf(fid, '\\\\ \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10} \n');
fprintf(fid, 'Dates ');
fprintf(fid, ' & 50.00\\%% &  %4.2f\\%% &  %4.2f\\%%', setQuantiles(1,ndxCI(1)), setQuantiles(1,ndxCI(2)));
fprintf(fid, ' & 50.00\\%% &  %4.2f\\%% &  %4.2f\\%%', setQuantiles(1,ndxCI(1)), setQuantiles(1,ndxCI(2)));
fprintf(fid, ' & 50.00\\%% &  %4.2f\\%% &  %4.2f\\%%', setQuantiles(1,ndxCI(1)), setQuantiles(1,ndxCI(2)));
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');

for t = 1 : size(data4table,1)
    fprintf(fid, '%s ', datestr(theseFcstdates(t), 'yyyy:mm'));
    fprintf(fid, '& %6.2f ', data4table(t,:));
    fprintf(fid, '\\\\\n');
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');
fprintf(fid, '\\end{small}\n');
fclose(fid);
type(fullfile(tabdir, strcat(tabname, '.tex')))
latexwrapper(wrap, 'add', 'sidetab', strcat(tabname, '.tex'), [])
end

function tabulatePanel2(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, prettylabel1, prettylabel2, fcstMid1, fcstMid2, fcstTails1, fcstTails2, ncode, ndxCI, setQuantiles) %#ok<DEFNU>

if isempty(wrap)
    return
end

tabname = sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28));

datalabels = {sprintf('%s median', prettylabel1), ...
    sprintf('%s %4.2f%% quantile', prettylabel1, setQuantiles(1,ndxCI(1))), ...
    sprintf('%s %4.2f%% quantile', prettylabel1, setQuantiles(1,ndxCI(2))), ...
    sprintf('%s median', prettylabel2), ...
    sprintf('%s %4.2f%% quantile', prettylabel2, setQuantiles(1,ndxCI(1))), ...
    sprintf('%s %4.2f%% quantile', prettylabel2, setQuantiles(1,ndxCI(2))), ...
    };

data4table = [transpose(fcstMid1(n,:,ndxT)), squeeze(fcstTails1(n,:,:,ndxT)), ...
    transpose(fcstMid2(n,:,ndxT)), squeeze(fcstTails2(n,:,:,ndxT)), ...
    ];

% write into csv
writedatatable(wrap, tabname, theseFcstdates, data4table, datalabels, 'yyyy:mm');

% tabulate in tex
tabdir = wrap.dir;
fid = fopen(fullfile(tabdir, strcat(tabname, '.tex')), 'wt');
fprintf(fid, '\\begin{small}\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, size(data4table,2)));
fprintf(fid, '\\toprule\n');
fprintf(fid, ' & \\multicolumn{3}{c}{%s}', prettylabel1);
fprintf(fid, ' & \\multicolumn{3}{c}{%s}', prettylabel2);
fprintf(fid, '\\\\ \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n');
fprintf(fid, 'Dates ');
fprintf(fid, ' & 50.00\\%% &  %4.2f\\%% &  %4.2f\\%%', setQuantiles(1,ndxCI(1)), setQuantiles(1,ndxCI(2)));
fprintf(fid, ' & 50.00\\%% &  %4.2f\\%% &  %4.2f\\%%', setQuantiles(1,ndxCI(1)), setQuantiles(1,ndxCI(2)));
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');

for t = 1 : size(data4table,1)
    fprintf(fid, '%s ', datestr(theseFcstdates(t), 'yyyy:mm'));
    fprintf(fid, '& %6.2f ', data4table(t,:));
    fprintf(fid, '\\\\\n');
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');
fprintf(fid, '\\end{small}\n');
fclose(fid);
type(fullfile(tabdir, strcat(tabname, '.tex')))
latexwrapper(wrap, 'add', 'sidetab', strcat(tabname, '.tex'), [])
end