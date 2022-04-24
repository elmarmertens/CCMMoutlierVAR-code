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

doSVtdof = false;

matSV         = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SV-p12.mat'));
matSVO        = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVOmax20-p12.mat'));
matSVoutmiss  = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVnanO5-p12.mat'));
if doSVtdof
    figlabel      = 'SVOtdof2021'; 
    matSVt        = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVt-dof5-p12.mat'));
    matSVOt       = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVOtmax20-dof9-p12.mat'));
else
    figlabel      = 'SVOt2021'; 
    matSVt        = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVt-p12.mat'));
    matSVOt       = matfile(fullfile(datadir, 'fredMD16-2021-04-censoredYields-SVOtmax20-p12.mat'));
end

wrap = [];
initwrap

do2 = true;
 %#ok<*UNRCH>
 
%% pull some objects
ydates = matSV.ydates;
ncode  = matSV.ncode;
N      = length(ncode);
Tjumpoffs = matSV.Tjumpoffs; % skip Jan/Feb
ndx = ismember(ydates(Tjumpoffs), [datenum(2020, [3 4 9 12], 1) datenum(2021, 3, 1)]);
Tjumpoffs = Tjumpoffs(ndx);
Njumpoffs = length(Tjumpoffs);

Ylabels = fredMDprettylabel(ncode);

fontsize = 12;

% lineTypes = {'-d', '--o', 'd-.', ':o'};
% lineTypes = {'-d', '--o', 'd-.', ':o', '-d', '--o', 'd-.', ':o'};
lineTypes = {'-', '-d', 'd-.', '--', ':'};


ylist = 6;
%% plot COVIDSV - 4 columns

close all

ndx = find(ydates >= datenum(2020,1,1)); % logical does not work with matfile objects
thesedates = ydates(ndx);
tickdates = [datenum(2020,1 , 1) datenum(2021,1, 1)];

hanni = NaN(Njumpoffs,1);
for n = ylist
    
    % multi-panel version
    thisfig = figure;
    
    subplot(2,4,1)
    axSV = gca;
    set(axSV, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSV.Tjumpoffs, Tjumpoffs))';
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        hanni(ndxT) = plot(thesedates, matSV.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2);
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVylim = ylim;
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title('(a) SV')
    
    
    % SVO
    subplot(2,4,2)
    axSVO = gca;
    set(axSVO, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSVO.Tjumpoffs, Tjumpoffs))';
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        plot(thesedates, matSVO.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVOylim = ylim;
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title('(b) SVO')
    
    % SVOX
    subplot(2,4,6)
    axSVOX = gca;
    set(axSVOX, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSVO.Tjumpoffs, Tjumpoffs))';
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        plot(thesedates, matSVO.drawsLambdaSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVOXylim = ylim;
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title(sprintf('(f) SVO \n ex outlier'))
    
    % SVt
    subplot(2,4,3)
    axSVt = gca;
    set(axSVt, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSVt.Tjumpoffs, Tjumpoffs))';
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        plot(thesedates, matSVt.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVtylim = ylim;
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title('(c) SV-t')
    
    % SVtX
    subplot(2,4,7)
    axSVtX = gca;
    set(axSVtX, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSVt.Tjumpoffs, Tjumpoffs))';
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        plot(thesedates, matSVt.drawsLambdaSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVtXylim = ylim;
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title(sprintf('(g) SV-t \n ex outlier'))
    
    % SVOt
    subplot(2,4,4)
    axSVOt = gca;
    set(axSVOt, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSVOt.Tjumpoffs, Tjumpoffs))';
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        plot(thesedates, matSVOt.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVOtylim = ylim;
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title('(d) SVO-t')
    
    % SVOtX
    subplot(2,4,8)
    axSVOtX = gca;
    set(axSVOtX, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSVOt.Tjumpoffs, Tjumpoffs))';
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        plot(thesedates, matSVOt.drawsLambdaSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVOtXylim = ylim;
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title(sprintf('(h) SVO-t \n ex outlier'))
    
    % SVoutmiss
    hanni = NaN(Njumpoffs, 1);
    subplot(2,4,5)
    axSVoutmiss = gca;
    set(axSVoutmiss, 'fontsize', fontsize)
    hold on
    theseJumpoffs = find(ismember(matSVoutmiss.Tjumpoffs, Tjumpoffs));
    iter = 0;
    for i = 1 : length(theseJumpoffs)
        iter = iter + 1;
        ndxT = theseJumpoffs(i);
        hanni(i) = plot(thesedates, matSVoutmiss.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2);
        if iter == length(lineTypes)
            iter = 0;
        end
    end
    SVoutmissylim = ylim;
    lgd = legend(hanni, datestr(ydates(Tjumpoffs), 'yyyy:mm'), ...
        'Orientation', 'horizontal', ...
        'location', 'layout', 'box', 'off'); 
    lgd.Position = [0 0 1 .05];
    lgd.FontSize = 11;
    % lgd.Layout.Tile = 'south';
    
    xticks(tickdates)
    datetick('x', 'yyyy', 'keepticks')
    xlim(thesedates([1 end]))
    title('(e) SV-OutMiss')
    
    
    lims = [0 max([SVOylim,SVOXylim, SVylim, SVtylim, SVtXylim, SVOtylim, SVOtXylim, SVoutmissylim])];
    upperrowlims = [0 max([SVOylim,SVylim, SVtylim, SVOtylim])];
    lowerrowlims = [0 max([SVOXylim, SVtXylim, SVOtXylim, SVoutmissylim])];
    ylim(axSV, upperrowlims)
    ylim(axSVO, upperrowlims)
    ylim(axSVt, upperrowlims)
    ylim(axSVOt, upperrowlims)
    
    ylim(axSVOX, upperrowlims)
    ylim(axSVtX, upperrowlims)
    ylim(axSVOtX, upperrowlims)
    ylim(axSVoutmiss, upperrowlims)
    
    orient landscape
    wrapthisfigure(thisfig, sprintf('covid%s-%s-all', figlabel, ncode{n}), wrap, [], [], [], [], true);
    
    ylim(axSVOX, lowerrowlims)
    ylim(axSVtX, lowerrowlims)
    ylim(axSVOtX, lowerrowlims)
    ylim(axSVoutmiss, lowerrowlims)
    wrapthisfigure(thisfig, sprintf('covid%s-%s-all-seplims', figlabel, ncode{n}), wrap, [], [], [], [], false);
    
    % wrapthisfigure(thisfig, sprintf('covidSV2020-%s', ncode{n}), wrap, [], [], [], [], true);
    %     sgtitle(Ylabels{n})
    %     wrapthisfigure(thisfig, sprintf('covidSV2020-%s-WITHTITLE', ncode{n}), wrap, [], [], [], [], false);
    
    
    
    
end

%% plot COVIDSV - 2 columns

if do2
    close all
    
    
    ndx = find(ydates >= datenum(2020,1,1)); % logical does not work with matfile objects
    thesedates = ydates(ndx);
    tickdates = [datenum(2020, [1 7], 1) datenum(2021,1, 1)];

    hanni = NaN(Njumpoffs,1);
    for n =  ylist
        
        % multi-panel version
        thisfig = figure;
        
        subplot(2,2,1)
        axSV = gca;
        set(axSV, 'fontsize', fontsize)
        hold on
        theseJumpoffs = find(ismember(matSV.Tjumpoffs, Tjumpoffs))';
        iter = 0;
        for i = 1 : length(theseJumpoffs)
            iter = iter + 1;
            ndxT = theseJumpoffs(i);
            hanni(ndxT) = plot(thesedates, matSV.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2);
            if iter == length(lineTypes)
                iter = 0;
            end
        end
        SVylim = ylim;
        xticks(tickdates)
        datetick('x', 'yyyy:mm', 'keepticks')
        xlim(thesedates([1 end]))
        title('(a) SV')
        
        
        
        % SVOt
        subplot(2,2,2)
        axSVOt = gca;
        set(axSVOt, 'fontsize', fontsize)
        hold on
        theseJumpoffs = find(ismember(matSVOt.Tjumpoffs, Tjumpoffs))';
        iter = 0;
        for i = 1 : length(theseJumpoffs)
            iter = iter + 1;
            ndxT = theseJumpoffs(i);
            plot(thesedates, matSVOt.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
            if iter == length(lineTypes)
                iter = 0;
            end
        end
        SVOtylim = ylim;
        xticks(tickdates)
        datetick('x', 'yyyy:mm', 'keepticks')
        xlim(thesedates([1 end]))
        title('(b) SVO-t')
        
        % SVOtX
        subplot(2,2,4)
        axSVOtX = gca;
        set(axSVOtX, 'fontsize', fontsize)
        hold on
        theseJumpoffs = find(ismember(matSVOt.Tjumpoffs, Tjumpoffs))';
        iter = 0;
        for i = 1 : length(theseJumpoffs)
            iter = iter + 1;
            ndxT = theseJumpoffs(i);
            plot(thesedates, matSVOt.drawsLambdaSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
            if iter == length(lineTypes)
                iter = 0;
            end
        end
        SVOtXylim = ylim;
        xticks(tickdates)
        datetick('x', 'yyyy:mm', 'keepticks')
        xlim(thesedates([1 end]))
        title('(d) SVO-t ex outlier')
        
        % SVoutmiss
        hanni = NaN(Njumpoffs, 1);
        subplot(2,2,3)
        axSVoutmiss = gca;
        set(axSVoutmiss, 'fontsize', fontsize)
        hold on
        theseJumpoffs = find(ismember(matSVoutmiss.Tjumpoffs, Tjumpoffs));
        iter = 0;
        for i = 1 : length(theseJumpoffs)
            iter = iter + 1;
            ndxT = theseJumpoffs(i);
            hanni(i) = plot(thesedates, matSVoutmiss.drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2);
            if iter == length(lineTypes)
                iter = 0;
            end
        end
        SVoutmissylim = ylim;
        lgd = legend(hanni, datestr(ydates(Tjumpoffs), 'yyyy:mm'), ...
            'Orientation', 'horizontal', ...
            'location', 'layout', 'box', 'off');
        lgd.Position = [0 0 1 .05];
        lgd.FontSize = 11;
        % lgd.Layout.Tile = 'south';
        
        xticks(tickdates)
        datetick('x', 'yyyy:mm', 'keepticks')
        xlim(thesedates([1 end]))
        title('(c) SV-OutMiss')
        
        
        lims = [0 max([SVylim, SVOtylim, SVOtXylim, SVoutmissylim])];
        upperrowlims = [0 max([SVylim, SVOtylim])];
        lowerrowlims = [0 max([SVOtXylim, SVoutmissylim])];
        ylim(axSV, upperrowlims)
        ylim(axSVOt, upperrowlims)
        
        ylim(axSVOtX, upperrowlims)
        ylim(axSVoutmiss, upperrowlims)
        
        orient landscape
        wrapthisfigure(thisfig, sprintf('covid%s-%s', figlabel, ncode{n}), wrap, [], [], [], [], true);
        
        ylim(axSVOtX, lowerrowlims)
        ylim(axSVoutmiss, lowerrowlims)
        wrapthisfigure(thisfig, sprintf('covid%s-%s-seplims', figlabel, ncode{n}), wrap, [], [], [], [], false);
        
        % wrapthisfigure(thisfig, sprintf('covidSV2020-%s', ncode{n}), wrap, [], [], [], [], true);
        %     sgtitle(Ylabels{n})
        %     wrapthisfigure(thisfig, sprintf('covidSV2020-%s-WITHTITLE', ncode{n}), wrap, [], [], [], [], false);
        
        
        
        
    end
    
end
%% finish
dockAllFigures
finishwrap
finishscript