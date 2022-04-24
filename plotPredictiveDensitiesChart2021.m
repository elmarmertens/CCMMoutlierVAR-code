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

fontsize = 18;
wrap = [];
initwrap

doConnectDots = false;

datestrfmt = 'mmmm yyyy';

%#ok<*UNRCH>

%% select plots
doCONST = false; % const only
doSV    = false;  % const vs sv

doRow1  = true;
doRow2  = true;
doRow3  = true;

doSVobar = false;


doCompareSVt  = false;
doCompareSVOt = false;

doShowOnlyLatest = false;

corelist = [1 5 6 12];

% ylist = corelist;
% ylist    = 1 : 16;
% ylist    = ylist(~ismember(ylist,corelist));

ylist = 6;

%% load matfiles (oos runs)

do2020 = '';


prettylabelCONST  = 'CONST';
mfilename         = sprintf('%s-CONST%s-p%d.mat', modellabel, do2020, p);
matCONST          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSV  = 'SV';
mfilename      = sprintf('%s-SV%s-p%d.mat', modellabel, do2020, p);
matSV          = matfile(fullfile(resultsdir, mfilename), 'writable', false);


prettylabelSVt = 'SV-t';
mfilename      = sprintf('%s-SVt%s-p%d.mat', modellabel, do2020, p);
matSVt         = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVO   = 'SVO';
mfilename        = sprintf('%s-SVOmax20%s-p%d.mat', modellabel, do2020, p);
matSVO           = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVobar   = 'SVo';
mfilename        = sprintf('%s-SVobarmax20%s-p%d.mat', modellabel, do2020, p);
matSVobar           = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVoutMiss  = 'SV-OutMiss';
Ofactor               = 5;
mfilename             = sprintf('%s-SVnanO%d%s-p%d.mat', modellabel, Ofactor, do2020, p);
matSVoutMiss          = matfile(fullfile(resultsdir, mfilename), 'writable', false);

prettylabelSVdummy  = 'SV-Dummies';
dummyPrecision      = 0;
mfilename           = sprintf('%s-SV-COVIDEACH-precision%d-p%d.mat', modellabel, dummyPrecision, p);
matSVdummy          = matfile(fullfile(resultsdir, mfilename), 'writable', false);


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

thismat = matSVO;
fcstTailsSVO  = NaN(size(fcstTailsCONST));
fcstMidSVO    = NaN(size(fcstMidCONST));
ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSVO(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSVO(:,:,ndx)      = thismat.fcstYmedian;


thismat = matSVobar;
fcstTailsSVobar  = NaN(size(fcstTailsCONST));
fcstMidSVobar    = NaN(size(fcstMidCONST));
ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSVobar(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSVobar(:,:,ndx)      = thismat.fcstYmedian;


thismat = matSVt;
fcstTailsSVt  = NaN(size(fcstTailsCONST));
fcstMidSVt    = NaN(size(fcstMidCONST));
ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSVt(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSVt(:,:,ndx)      = thismat.fcstYmedian;


thismat = matSVOt;
fcstTailsSVOt  = NaN(size(fcstTailsCONST));
fcstMidSVOt    = NaN(size(fcstMidCONST));
ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSVOt(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSVOt(:,:,ndx)      = thismat.fcstYmedian;



thismat = matSVoutMiss;
fcstTailsSVoutMiss  = NaN(size(fcstTailsCONST));
fcstMidSVoutMiss    = NaN(size(fcstMidCONST));
ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSVoutMiss(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSVoutMiss(:,:,ndx)      = thismat.fcstYmedian;

outliermissYmid   = NaN(length(ydates), N, length(Tjumpoffs));
outliermissYtails = NaN(length(ydates), N, length(ndxCI), length(Tjumpoffs));
outliermissYmid(:,:,ndx)     = matSVoutMiss.drawsYmid;
outliermissYtails(:,:,:,ndx) = matSVoutMiss.drawsYtails(:,:,ndxCI,:);


thismat = matSVdummy;
fcstTailsSVdummy = NaN(size(fcstTailsCONST));
fcstMidSVdummy    = NaN(size(fcstMidCONST));
ndx                     = find(ismember(Tjumpoffs, thismat.Tjumpoffs));
fcstTailsSVdummy(:,:,:,ndx)  = squeeze(thismat.fcstYquantiles(:,:,ndxCI,:));
fcstMidSVdummy(:,:,ndx)      = thismat.fcstYmedian;


%% list ylabels and ncodes
if ~isempty(wrap)
    tabname = sprintf('datalist-%s.tex', datalabel);
    filename = fullfile(wrap.dir, tabname);
    fid = fopen(filename, 'wt');
    
    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{tabular}{lll}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, 'Variable & FRED-MD code & transformation ');
    fprintf(fid, '\\\\\n');
    fprintf(fid, '\\midrule\n');
    for n = 1 : N
        fprintf(fid, '%s ', Ylabels{n});
        fprintf(fid, '& %s ', ncode{n});
        
        switch tcode(n)
            case 1
                % % no transformation
            case 2
                fprintf(fid, '& %s\n', '\ensuremath{\Delta x_t}');
            case 4
                fprintf(fid, '& %s\n', '\ensuremath{\log(x_t)}');
            case 5
                fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 1200}');
            otherwise
                error houston
        end
        
        fprintf(fid, '\\\\\n');
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{center}\n');
    fprintf(fid, '\n');
    %     fprintf(fid, 'Note: ');
    %     fprintf(fid, 'Data obtained from the %s vintage of FRED-MD. ', strtok(csv_in, '.'));
    %     fprintf(fid, 'Monthly observations ');
    %     fprintf(fid, 'from %s to %s. ', datestryymm(ydates(1)), datestryymm(ydates(end)));
    % fprintf(fid, '%s \n', fredMDtcodeNote);
    fclose(fid);
    type(filename)
    latexwrapper(wrap, 'add', 'sidetab', tabname, [])
end

%% select dates to be shown
ndxRange = find(ydates(Tjumpoffs) >= datenum(2020,1,1));

ndxRange = ndxRange(:)'; % ensure it is a row vector

if doShowOnlyLatest
    ndxRange = ndxRange(end);
end
    
darkgreen = [0 .75 0];
black     = [0 0 0];

tickdates = [datenum(2019,[6 12],1),datenum(2020,[6 12],1), datenum(2021,[6 12],1), datenum(2022,[6 12],1), datenum(2023,[6 12],1)];

%% plot CONST vs SV vs SVOt
if doRow1
    close all
    thisPlotLabel = 'predictiveDensityChartCONSTvsSVvsSVOt';
    for n = ylist
        
        for ndxT = ndxRange
            
            jumpoff   = Tjumpoffs(ndxT);
            thisdateT = ydates(jumpoff);
            
            %% compare predictive densities (jumpoff)
            jj                = find(fcstdates == ydates(jumpoff));
            %             allthesedates     = fcstdates(jj:jj+fcstNhorizons);
            theseFcstdates    = fcstdates(jj+(1:fcstNhorizons));
            %             theseJumpoffdates = fcstdates(jj+(-np:0));
            %             if ~isequal(theseJumpoffdates, ydates(jj+(-np:0)))
            %                 error houston
            %             end
            
            theseRealized = NaN(fcstNhorizons, 1);
            ndx = jj+(1:fcstNhorizons);
            ndx = ndx(ndx <= Tdata);
            theseRealized(1 : length(ndx)) = data(ndx,n);
            %             theseFcstdatesCumJumpoff = [theseJumpoffdates(end); theseFcstdates];
            
            
            hanni = NaN(3,1);
            
            
            %% plot
            thisfig   = figure;
            set(thisfig,'defaultLegendAutoUpdate','off');
            ax = gca;
            set(ax, 'fontsize', fontsize)
            hold on
            hanni(1)   = plotCI(fcstMidCONST(n,:,ndxT), squeeze(fcstTailsCONST(n,:,:,ndxT)), ...
                theseFcstdates, [], 'k-', 'linewidth', 3);
            
            hanni(2)   = plot(theseFcstdates,fcstMidSV(n,:,ndxT), 'r-', 'linewidth', 3);
            plot(theseFcstdates,squeeze(fcstTailsSV(n,:,:,ndxT)), 'r--', 'linewidth', 2);
            
            hanni(3)   = plot(theseFcstdates,fcstMidSVOt(n,:,ndxT), 'b-.', 'linewidth', 3);
            plot(theseFcstdates,squeeze(fcstTailsSVOt(n,:,:,ndxT)), 'b-.', 'linewidth', 2);
            
            %             plot(theseJumpoffdates, data(jj+(-np:0),n), '-', 'color', darkgreen, 'linewidth', 2);
            %             if doConnectDots
            %                 plot([theseJumpoffdates(end) theseFcstdates(1)], [data(jj,n) fcstMidCONST(n,1,ndxT)], ':', 'color', black, 'linewidth', 2);
            %             end
            
            if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
            
            xlim(theseFcstdates([1 end]))
            xticks(tickdates)
            datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
            ax.XAxis.MinorTick       = 'on';
            
            wrapthisfigure(thisfig, sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLE', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend(hanni, prettylabelCONST, prettylabelSV, prettylabelSVOt, ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(ht)
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(hl);
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            hd = plot(theseFcstdates, theseRealized, 'd', 'color', darkgreen, 'linewidth', 2);
            plot(theseFcstdates, theseRealized, ':', 'color', darkgreen, 'linewidth', 2);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLEDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], false)
            
            delete(ht);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend([hanni; hd], prettylabelCONST, prettylabelSV, prettylabelSVOt, 'realized', ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGENDDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            
        end
    end
    
end

%% plot SVO vs SVt
if doRow2
    close all
    thisPlotLabel = 'predictiveDensityChartSVOt';
    for n = ylist
        
        for ndxT = ndxRange
            
            jumpoff   = Tjumpoffs(ndxT);
            thisdateT = ydates(jumpoff);
            
            
            %% compare predictive densities (jumpoff)
            jj                = find(fcstdates == ydates(jumpoff));
            % allthesedates     = fcstdates(jj-np:jj+fcstNhorizons);
            theseFcstdates    = fcstdates(jj+(1:fcstNhorizons));
            % theseJumpoffdates = fcstdates(jj+(-np:0));
            % if ~isequal(theseJumpoffdates, ydates(jj+(-np:0)))
            %     error houston
            % end
            
            theseRealized = NaN(fcstNhorizons, 1);
            ndx = jj+(1:fcstNhorizons);
            ndx = ndx(ndx <= Tdata);
            theseRealized(1 : length(ndx)) = data(ndx,n);
            % theseFcstdatesCumJumpoff = [theseJumpoffdates(end); theseFcstdates];
            
            
            hanni = NaN(3,1);
            
            
            %% plot
            thisfig   = figure;
            set(thisfig,'defaultLegendAutoUpdate','off');
            ax = gca;
            set(ax, 'fontsize', fontsize)
            hold on
            hanni(1)   = plotCIaltcolor(fcstMidSVOt(n,:,ndxT), squeeze(fcstTailsSVOt(n,:,:,ndxT)), ...
                theseFcstdates, [], 'b-', 'linewidth', 3);
            
            hanni(2)   = plot(theseFcstdates,fcstMidSVO(n,:,ndxT), 'k-.', 'linewidth', 3);
            plot(theseFcstdates,squeeze(fcstTailsSVO(n,:,:,ndxT)), 'k-.', 'linewidth', 2);
            
            hanni(3)   = plot(theseFcstdates,fcstMidSVt(n,:,ndxT), 'r--', 'linewidth', 3);
            hh = plot(theseFcstdates,squeeze(fcstTailsSVt(n,:,:,ndxT)), 'r--', 'linewidth', 2);
            hanni(3) = hh(1);
            
            %         hanni(4)   = plot(theseFcstdates,fcstMidSVOt(n,:,ndxT), 'c:', 'linewidth', 3);
            %         plot(theseFcstdates,squeeze(fcstTailsSVOt(n,:,:,ndxT)), 'c:', 'linewidth', 2);
            
            
            % plot(theseJumpoffdates, data(jj+(-np:0),n), '-', 'color', darkgreen, 'linewidth', 2);
            % if doConnectDots
            %     plot([theseJumpoffdates(end) theseFcstdates(1)], [data(jj,n) fcstMidSVt(n,1,ndxT)], ':', 'color', black, 'linewidth', 2);
            % end
            
            if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
            
            xlim(theseFcstdates([1 end]))
            xticks(tickdates)
            datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
            ax.XAxis.MinorTick       = 'on';
            
            wrapthisfigure(thisfig, sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLE', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend(hanni, prettylabelSVOt, prettylabelSVO, prettylabelSVt, ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(ht)
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(hl);
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            hd = plot(theseFcstdates, theseRealized, 'd', 'color', darkgreen, 'linewidth', 2);
            plot(theseFcstdates, theseRealized, ':', 'color', darkgreen, 'linewidth', 2);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLEDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], false)
            
            delete(ht);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend([hanni; hd], prettylabelSVOt, prettylabelSVO, prettylabelSVt, 'realized', ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGENDDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            %         tabulatePanel(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, ...
            %             prettylabelSV, prettylabelSVO, prettylabelSVt, ...
            %             fcstMidSV, fcstMidSVO, fcstMidSVt, ...
            %             fcstTailsSV, fcstTailsSVO, fcstTailsSVt, ...
            %             ncode, ndxCI, setQuantiles)
        end
    end
end

%% plot SVobar vs SVt
if doSVobar
    close all
    thisPlotLabel = 'predictiveDensityChartSVobar';
    for n = ylist
        
        for ndxT = ndxRange
            
            jumpoff   = Tjumpoffs(ndxT);
            thisdateT = ydates(jumpoff);
            
            
            %% compare predictive densities (jumpoff)
            jj                = find(fcstdates == ydates(jumpoff));
            %             allthesedates     = fcstdates(jj-np:jj+fcstNhorizons);
            theseFcstdates    = fcstdates(jj+(1:fcstNhorizons));
            %             theseJumpoffdates = fcstdates(jj+(-np:0));
            %             if ~isequal(theseJumpoffdates, ydates(jj+(-np:0)))
            %                 error houston
            %             end
            
            theseRealized = NaN(fcstNhorizons, 1);
            ndx = jj+(1:fcstNhorizons);
            ndx = ndx(ndx <= Tdata);
            theseRealized(1 : length(ndx)) = data(ndx,n);
            %             theseFcstdatesCumJumpoff = [theseJumpoffdates(end); theseFcstdates];
            
            
            hanni = NaN(3,1);
            
            
            %% plot
            thisfig   = figure;
            set(thisfig,'defaultLegendAutoUpdate','off');
            ax = gca;
            set(ax, 'fontsize', fontsize)
            hold on
            hanni(1)   = plotCIaltcolor(fcstMidSVOt(n,:,ndxT), squeeze(fcstTailsSVOt(n,:,:,ndxT)), ...
                theseFcstdates, [], 'b-', 'linewidth', 3);
            
            hanni(2)   = plot(theseFcstdates,fcstMidSVO(n,:,ndxT), 'k-.', 'linewidth', 3);
            plot(theseFcstdates,squeeze(fcstTailsSVO(n,:,:,ndxT)), 'k-.', 'linewidth', 2);
            
            hanni(3)   = plot(theseFcstdates,fcstMidSVobar(n,:,ndxT), 'r--', 'linewidth', 3);
            hh = plot(theseFcstdates,squeeze(fcstTailsSVobar(n,:,:,ndxT)), 'r--', 'linewidth', 2);
            hanni(3) = hh(1);
            
            %         hanni(4)   = plot(theseFcstdates,fcstMidSVOt(n,:,ndxT), 'c:', 'linewidth', 3);
            %         plot(theseFcstdates,squeeze(fcstTailsSVOt(n,:,:,ndxT)), 'c:', 'linewidth', 2);
            
            
            %             plot(theseJumpoffdates, data(jj+(-np:0),n), '-', 'color', darkgreen, 'linewidth', 2);
            %             if doConnectDots
            %                 plot([theseJumpoffdates(end) theseFcstdates(1)], [data(jj,n) fcstMidSVt(n,1,ndxT)], ':', 'color', black, 'linewidth', 2);
            %             end
            
            if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
            
            xlim(theseFcstdates([1 end]))
            xticks(tickdates)
            datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
            ax.XAxis.MinorTick       = 'on';
            
            wrapthisfigure(thisfig, sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLE', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend(hanni, prettylabelSVOt, prettylabelSVO, prettylabelSVobar, ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(ht)
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(hl);
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            hd = plot(theseFcstdates, theseRealized, 'd', 'color', darkgreen, 'linewidth', 2);
            plot(theseFcstdates, theseRealized, ':', 'color', darkgreen, 'linewidth', 2);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLEDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(ht);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend([hanni; hd], prettylabelSVOt, prettylabelSVO, prettylabelSVobar, 'realized', ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGENDDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGENDDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], false)
            
            %             tabulatePanel(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, ...
            %                 prettylabelSVOt, prettylabelSVO, prettylabelSVobar, ...
            %                 fcstMidSVOt, fcstMidSVO, fcstMidSVobar, ...
            %                 fcstTailsSVOt, fcstTailsSVO, fcstTailsSVobar, ...
            %                 ncode, ndxCI, setQuantiles)
        end
    end
end

%% plot SVO vs SVdummy vs SVoutmiss

if doRow3
    ndxRange = find(ydates(Tjumpoffs) >= datenum(2020,3,1));
    ndxRange = ndxRange(:)';
    
    
    if doShowOnlyLatest
        ndxRange = ndxRange(end);
    end
    
    close all
    thisPlotLabel = 'predictiveDensityChartOutmiss';
    for n = ylist
        
        for ndxT = ndxRange
            
            jumpoff   = Tjumpoffs(ndxT);
            thisdateT = ydates(jumpoff);
            
            %% identify outlier data per sample end
            dataO           = data(1:jumpoff,:);
            dev             = abs(dataO - median(dataO));
            iqr             = range(prctile(dataO, [25 75]));
            dataOnan        = dev > Ofactor * iqr;
            outlier         = dataO;
            outlier(~dataOnan) = NaN;
            dataO(dataOnan) = NaN;
            
            
            
            %% compare predictive densities (jumpoff)
            jj                = find(fcstdates == ydates(jumpoff));
            %             allthesedates     = fcstdates(jj-np:jj+fcstNhorizons);
            theseFcstdates    = fcstdates(jj+(1:fcstNhorizons));
            %             theseJumpoffdates = fcstdates(jj+(-np:0));
            %             if ~isequal(theseJumpoffdates, ydates(jj+(-np:0)))
            %                 error houston
            %             end
            
            theseRealized = NaN(fcstNhorizons, 1);
            ndx = jj+(1:fcstNhorizons);
            ndx = ndx(ndx <= Tdata);
            theseRealized(1 : length(ndx)) = data(ndx,n);
            %             theseFcstdatesCumJumpoff = [theseJumpoffdates(end); theseFcstdates];
            
            
            hanni = NaN(3,1);
            
            
            %% plot
            thisfig   = figure;
            set(thisfig,'defaultLegendAutoUpdate','off');
            ax = gca;
            set(ax, 'fontsize', fontsize)
            hold on
            hanni(1)   = plotCI(fcstMidSVoutMiss(n,:,ndxT), squeeze(fcstTailsSVoutMiss(n,:,:,ndxT)), ...
                theseFcstdates, [], 'k--', 'linewidth', 3);
            
            hanni(2)   = plot(theseFcstdates,fcstMidSVOt(n,:,ndxT), 'b-.', 'linewidth', 3);
            plot(theseFcstdates,squeeze(fcstTailsSVOt(n,:,:,ndxT)), 'b-.', 'linewidth', 2);
            
            hanni(3)   = plot(theseFcstdates,fcstMidSVdummy(n,:,ndxT), 'm-', 'linewidth', 3); % 'color', [0 .5 0],
            plot(theseFcstdates,squeeze(fcstTailsSVdummy(n,:,:,ndxT)), 'm-', 'linewidth', 1);
            
            
            % plot(theseJumpoffdates, outliermissYmid(jj+(-np:0),n, ndxT), ':', 'color', [0 0 0], ...
            %     'markersize', 8, 'linewidth', 2);
            % plot(theseJumpoffdates, squeeze(outliermissYtails(jj+(-np:0),n, :, ndxT)), ':', 'color', [0 0 0], ...
            %     'markersize', 8, 'linewidth', 2);
            %             plotCI(outliermissYmid(jj+(-np:0),n, ndxT), squeeze(outliermissYtails(jj+(-np:0),n, :, ndxT)), ...
            %                 theseJumpoffdates, [], 'k:', 'linewidth', 2);
            %
            %             plot(theseJumpoffdates, data(jj+(-np:0),n), '-', 'color', darkgreen, 'linewidth', 2);
            %             plot(theseJumpoffdates, outlier(jj+(-np:0),n), 'o', 'color', [0 0 0], ...
            %                 'markersize', 8, 'linewidth', 2);
            %             if doConnectDots
            %                 plot([theseJumpoffdates(end) theseFcstdates(1)], [data(jj,n) fcstMidSVoutMiss(n,1,ndxT)], ':', 'color', black, 'linewidth', 2);
            %             end
            
            if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
            
            xlim(theseFcstdates([1 end]))
            xticks(tickdates)
            datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
            ax.XAxis.MinorTick       = 'on';
            
            wrapthisfigure(thisfig, sprintf('%s-%s-%s', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLE', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend(hanni, prettylabelSVoutMiss, prettylabelSVOt, prettylabelSVdummy,  ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(ht)
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGEND', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(hl);
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            hd = plot(theseFcstdates, theseRealized, 'd', 'color', darkgreen, 'linewidth', 2);
            plot(theseFcstdates, theseRealized, ':', 'color', darkgreen, 'linewidth', 2);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLEDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            delete(ht);
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            hl = legend([hanni; hd], prettylabelSVoutMiss, prettylabelSVOt, prettylabelSVdummy, 'realized',  ...
                'location', 'best');
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHLEGENDDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], true)
            
            ht = title(sprintf('%s', datestr(thisdateT, datestrfmt)));
            wrapthisfigure(thisfig, sprintf('%s-%s-%s-WITHTITLELEGENDDATA', thisPlotLabel, ncode{n}, datestr(thisdateT, 28)), wrap, [], [], [], [], false)
            
            %             tabulatePanel(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, ...
            %                 prettylabelSVoutMiss, prettylabelSVOt, prettylabelSVdummy, ...
            %                 fcstMidSVoutMiss, fcstMidSVOt, fcstMidSVdummy, ...
            %                 fcstTailsSVoutMiss, fcstTailsSVOt, fcstTailsSVdummy, ...
            %                 ncode, ndxCI, setQuantiles)
        end
    end
end



%% finish
dockAllFigures
finishwrap
finishscript

function tabulatePanel(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, prettylabel1, prettylabel2, prettylabel3, fcstMid1, fcstMid2, fcstMid3, fcstTails1, fcstTails2, fcstTails3, ncode, ndxCI, setQuantiles)

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

function tabulatePanel2(wrap, thisPlotLabel, n, ndxT, thisdateT, theseFcstdates, prettylabel1, prettylabel2, fcstMid1, fcstMid2, fcstTails1, fcstTails2, ncode, ndxCI, setQuantiles)

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