%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% Initial operations
clear; close all; clc;

datadir = '~/jam/lager/var2021-matfiles/baseline';

matSV         = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SV-p12.mat'));
matSVO        = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVOmax20-p12.mat'));
matSVoutmiss  = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVnanO5-p12.mat'));
matSVt        = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVt-p12.mat'));
matSVOt       = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVOtmax20-p12.mat'));
% matSVOt2      = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVOtmax20-hipriorprob-p12.mat'));
matSVobar     = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVobarmax20-p12.mat'));

wrap = [];
initwrap


%#ok<*UNRCH>

%% pull some objects
ydates = matSV.ydates;
ncode  = matSV.ncode;
N      = length(ncode);
Tjumpoffs = matSV.Tjumpoffs;

p = 12;

Ylabels = fredMDprettylabel(ncode);

fontsize = 18;

corelist = [1 5 6 12];

% ylist = corelist;
ylist = 1 : N;

orange = [255,125,0] / 255;
green  = [0 .5 0];

tickdates = datenum(1960:10:2020,1,1);
datelim   = [datenum(1960,1,1) ydates(end)];

%% plot SV w/ and w/o outliers for each model
for ndxT = find(ydates(Tjumpoffs) == datenum(2020,1,1)) 
    
    % [find(ydates(Tjumpoffs) == datenum(2021,3,1)), find(ydates(Tjumpoffs) == datenum(2020,1,1))]
    
    thisT     = Tjumpoffs(ndxT);
    
    for n = ylist
        
        svfig = figure;
        axSV = gca;
        set(axSV, 'fontsize', fontsize)
        hold on
        plot(ydates, matSV.drawsSVmid(n,:,ndxT), '-', 'color', green, 'linewidth', 2)
        xtickdates(datelim)
        xticks(tickdates)
        ylimSV = ylim;
        
        svofig = figure;
        axSVO = gca;
        set(axSVO, 'fontsize', fontsize)
        hold on
        plot(ydates, matSVO.drawsSVmid(n,:,ndxT), ':', 'color', 'blue', 'linewidth', 3);
        plot(ydates, matSVO.drawsLambdaSVmid(n,:,ndxT), 'k-', 'linewidth', 2)
        xtickdates(datelim)
        xticks(tickdates)
        ylimSVO = ylim;
        
        svotfig = figure;
        axSVOt = gca;
        set(axSVOt, 'fontsize', fontsize)
        hold on
        plot(ydates, matSVOt.drawsSVmid(n,:,ndxT), 'r:', 'linewidth', 3);
        plot(ydates, matSVOt.drawsLambdaSVmid(n,:,ndxT), 'k-', 'linewidth', 2)
        xtickdates(datelim)
        xticks(tickdates)
        ylimSVOt = ylim;
        
        %         svotfig2 = figure;
        %         axSVOt2 = gca;
        %         set(axSVOt2, 'fontsize', fontsize)
        %         hold on
        %         plot(ydates, matSVOt2.drawsSVmid(n,:,ndxT), 'm:', 'linewidth', 3);
        %         plot(ydates, matSVOt2.drawsLambdaSVmid(n,:,ndxT), 'k-', 'linewidth', 2)
        %         xtickdates(datelim)
        %         xticks(tickdates)
        %         ylimSVOt2 = ylim;
        
        svtfig = figure;
        axSVt = gca;
        set(axSVt, 'fontsize', fontsize)
        hold on
        plot(ydates, matSVt.drawsSVmid(n,:,ndxT), ':', 'color', [125,175,255] / 255, 'linewidth', 3);
        plot(ydates, matSVt.drawsLambdaSVmid(n,:,ndxT), 'k-', 'linewidth', 3)
        xtickdates(datelim)
        xticks(tickdates)
        ylimSVt = ylim;
        
        svobarfig = figure;
        axSVobar = gca;
        set(axSVobar, 'fontsize', fontsize)
        hold on
        plot(ydates, matSVobar.drawsSVmid(n,:,ndxT), ':', 'color', 'magenta', 'linewidth', 3);
        plot(ydates, matSVobar.drawsLambdaSVmid(n,:,ndxT), 'k-', 'linewidth', 2)
        xtickdates(datelim)
        xticks(tickdates)
        ylimSVobar = ylim;
        
        
        lims = [0 max([ylimSV ylimSVO, ylimSVOt  ylimSVt ylimSVobar])];
        %         ylim(axSV, lims)
        ylim(axSVO, lims)
        ylim(axSVOt, lims)
        %         ylim(axSVOt2, lims)
        ylim(axSVt, lims)
        %         ylim(axSVobar, lims)
        
        wrapthisfigure(svfig, sprintf('stochvol-SV-%s-%s', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], true);
        wrapthisfigure(svtfig, sprintf('stochvol-SVt-%s-%s', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], true);
        wrapthisfigure(svotfig, sprintf('stochvol-SVOt-%s-%s', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], true);
        %         wrapthisfigure(svotfig2, sprintf('stochvol-SVOt2-%s-%s', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], true);
        wrapthisfigure(svofig, sprintf('stochvol-SVO-%s-%s', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], true);
        wrapthisfigure(svobarfig, sprintf('stochvol-SVobar-%s-%s', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], true);
        
        title(axSV, sprintf('SV \n%s', datestr(ydates(thisT), 'mmm yyyy')))
        title(axSVt, sprintf('SV-t \n%s', datestr(ydates(thisT), 'mmm yyyy')))
        ht = title(axSVOt, sprintf('SVO-t \n%s', datestr(ydates(thisT), 'mmm yyyy')));
        %         title(axSVOt2, sprintf('SVO-t ALT \n%s', datestr(ydates(thisT), 'mmm yyyy')))
        title(axSVO, sprintf('SVO \n%s', datestr(ydates(thisT), 'mmm yyyy')))
        title(axSVobar, sprintf('SV-o \n%s', datestr(ydates(thisT), 'mmm yyyy')))
        wrapthisfigure(svfig, sprintf('stochvol-SV-%s-%s-WITHTITLE', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], false);
        wrapthisfigure(svtfig, sprintf('stochvol-SVt-%s-%s-WITHTITLE', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], false);
        wrapthisfigure(svotfig, sprintf('stochvol-SVOt-%s-%s-WITHTITLE', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], false);
        %         wrapthisfigure(svotfig2, sprintf('stochvol-SVOt2-%s-%s-WITHTITLE', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], false);
        wrapthisfigure(svofig, sprintf('stochvol-SVO-%s-%s-WITHTITLE', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], false);
        wrapthisfigure(svobarfig, sprintf('stochvol-SVobar-%s-%s-WITHTITLE', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], false);
        
        delete(ht)
        legend(axSVOt, 'Total', 'Persistent SV', 'location', 'best')
        wrapthisfigure(svotfig, sprintf('stochvol-SVOt-%s-%s-WITHLEGEND', ncode{n}, datestr(ydates(thisT), 'mmmyyyy')), wrap, [], [], [], [], true);

    end
    
end

%% finish
dockAllFigures
finishwrap
finishscript