%% collect and tabulate predictive scores for various models with AR(1) in SV

clear
close all
fclose all;
clc

%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% setup

resultsdir     = pwd;

wrap = [];
initwrap

%% define models

% note: for SVO/SVOt set useRB to true to use the RB scores (for other models, the standard scores use RB anyway)

% SV
m = 1; % m + 1;
models(m).datalabel    = 'fredMD16-2021-04';
models(m).resultlabel  = 'censoredYields-SVar1-p12';
models(m).prettylabel  = 'SV-AR(1)';
models(m).useRB        = false;

% SVO
m = m + 1;
models(m).datalabel    = 'fredMD16-2021-04';
models(m).resultlabel  = 'censoredYields-SVOar1max20-p12';
models(m).prettylabel  = 'SVO-AR(1)';
models(m).useRB        = true;

% SV-t
m = m + 1;
models(m).datalabel    = 'fredMD16-2021-04';
models(m).resultlabel  = 'censoredYields-SVtar1-p12';
models(m).prettylabel  = 'SV-t-AR(1)';
models(m).useRB        = false;

% SVO-t
m = m + 1;
models(m).datalabel    = 'fredMD16-2021-04';
models(m).resultlabel  = 'censoredYields-SVOtar1max20-p12';
models(m).prettylabel  = 'SVO-t-AR(1)';
models(m).useRB        = true;

% SV-OutMiss
m = m + 1;
models(m).datalabel    = 'fredMD16-2021-04';
models(m).resultlabel  = 'censoredYields-SVar1nanO5-p12';
models(m).prettylabel  = 'SV-AR(1)-OutMiss';
models(m).useRB        = false;

models   = models([4 2 3 5 1]); % order benchmark last

%% load models
oos = cell(0);

varlist0 = {'ydates', 'Tjumpoffs', 'N', 'tcode', 'ncode', 'cumcode', ...
    'MCMCdraws', 'fcstNdraws', ...
    'fcstNhorizons', 'fcstYmvlogscore', ...
    };

varlistRB = {'ydates', 'Tjumpoffs', 'N', 'tcode', 'ncode', 'cumcode', ...
    'MCMCdraws', 'fcstNdraws', ...
    'fcstNhorizons', 'fcstYmvlogscore', ...
    'fcstYmvlogscoreRB' ...
    };

for m = 1 : length(models)
    
    %#ok<*UNRCH>
    %% load data
    
    if models(m).useRB
        varlist = varlistRB;
    else
        varlist = varlist0;
    end
    
    matfilename = sprintf('%s-%s.mat', models(m).datalabel, models(m).resultlabel);
    oos{m} = load(fullfile(resultsdir, matfilename), ...
        varlist{:});
    if models(m).useRB
        oos{m}.fcstYmvlogscore  = oos{m}.fcstYmvlogscoreRB;
    end
    
    if m > 1
        if oos{m}.MCMCdraws ~= oos{1}.MCMCdraws
            warning('unequal numbers of MCMCdraws, 0 has %d, 1 has %d', oos{m}.MCMCdraws, oos{1}.MCMCdraws)
        end
        
        if oos{m}.ydates(1) ~= oos{1}.ydates(1)
            error('oos estimates based on different sample starts: %s vs %s', datestr(oos{m}.ydates(1)), datestr(oos{1}.ydates(1)))
        end
        if oos{m}.ydates(end) ~= oos{1}.ydates(end)
            error('oos estimates based on different sample starts: %s vs %s', datestr(oos{m}.ydates(1)), datestr(oos{1}.ydates(1)))
        end
        
        if ~isequal(oos{m}.Tjumpoffs, oos{1}.Tjumpoffs)
            error('oos jumpoffs differ')
        end
    end
end

%% collect and merge scores
SCORES    = cell2mat(cellfun(@(x) x.fcstYmvlogscore', oos, 'UniformOutput', false));
TJUMPOFFS = oos{1}.Tjumpoffs;
YDATES    = oos{1}.ydates;

% SCORES2 = cell2mat(cellfun(@(x) x.fcstYmvlogscore', oos2, 'UniformOutput', false));
%
% SCORES    = cat(1, SCORES2, SCORES1); % mind the order
% TJUMPOFFS = cat(1, oos2{1}.Tjumpoffs, oos{1}.Tjumpoffs);
% YDATES    = cat(1, oos2{1}.ydates, oos{1}.ydates);
%
% % check
% if ~all(diff(TJUMPOFFS) == 1)
%     error('merge error TJUMPOFFS')
% end
%
%
% % patch in SCORES3
% SCORES3 = cell2mat(cellfun(@(x) x.fcstYmvlogscore', oos3, 'UniformOutput', false));
% ndx = ismember(TJUMPOFFS, oos3{1}.Tjumpoffs);
%
% SCORES(ndx,:) = SCORES3;

%% collect some global parameters

Ylabels = fredMDprettylabel(oos{1}.ncode);
N       = length(Ylabels);

Ylabels = strrep(Ylabels, '_', '');


%% loop over samples
s = 0;

s= s + 1;
sam(s).evalStop    = datenum(2021,2,1); % datenum(2017,12,1);
sam(s).evalStart   = datenum(1975,1,1);
sam(s).label       = '1975:01-2021:02';
sam(s).prettylabel  = 'Full sample';

% s= s + 1;
% sam(s).evalStop    = datenum(2019,12,1); % datenum(2017,12,1);
% sam(s).evalStart   = datenum(1975,1,1);
% sam(s).label       = '1975--2019';
% sam(s).prettylabel  = 'Pre 2020';

s= s + 1;
sam(s).evalStop  = datenum(1984,12,1); % datenum(2017,12,1);
sam(s).evalStart = datenum(1975,1,1);
sam(s).label     = '1975-1984';
sam(s).prettylabel  = 'G Inflation';

s= s + 1;
sam(s).evalStop  = datenum(2007,12,1); % datenum(2017,12,1);
sam(s).evalStart = datenum(1985,1,1);
sam(s).label     = '1985-2007';
sam(s).prettylabel  = 'G Moderation';

s= s + 1;
sam(s).evalStart = datenum(2008,1,1);
sam(s).evalStop  = datenum(2014,12,1); % datenum(2017,12,1);
sam(s).label     = '2008-2014';
sam(s).prettylabel  = 'GFC';

s = s + 1;
sam(s).evalStop  = datenum(2021,2,1); % datenum(2017,12,1);
sam(s).evalStart = datenum(2020,3,1);
sam(s).label     = '2020:03-2021:02';
sam(s).prettylabel  = 'COVID-19';

s = s + 1;
sam(s).evalStop  = datenum(2021,2,1); % datenum(2017,12,1);
sam(s).evalStart = datenum(2020,7,1);
sam(s).label     = '2020:07-2021:02';
sam(s).prettylabel  = 'COVID-19 since Jul 2020';

%% prepare results table across  samples
scoresTABLE       = NaN(length(models) - 1, length(sam));
scoresdmstatTABLE = NaN(length(models) - 1, length(sam));

%% loop over samples
for s = 1 : length(sam)
    
    %% reset Tjumpoff and ydates
    Tjumpoffs = TJUMPOFFS;
    ydates    = YDATES;
    
    %% eval window
    evalStart = sam(s).evalStart;
    evalStop  = sam(s).evalStop;
    
    samLabel = sprintf('%s-%s', datestr(evalStart, 'yyyymm'), datestr(evalStop, 'yyyymm'));
    
    
    
    %% cut eval sample if desired
    if isempty(evalStop)
        evalStop = ydates(Tjumpoffs(end));
    end
    if isempty(evalStart)
        evalStart = ydates(Tjumpoffs(1));
    end
    
    ndxJumpoff = ismember(Tjumpoffs, find((ydates >= evalStart) & (ydates <= evalStop)));
    
    Tjumpoffs   = Tjumpoffs(ndxJumpoff);
    
    
    dates    = ydates(Tjumpoffs);
    
    comparisonNote = sprintf('Evaluation window from %s through %s.', ...
        datestryymm(dates(1)), datestryymm(dates(end)));
    
    %% collect stuff
    
    % cut sample as need
    scores = SCORES(ndxJumpoff,:);
    
    
    avgscores   = mean(scores);
    sumscores   = sum(scores);
    dscores     = scores(:,1:end-1) - scores(:,end);
    sumdscores  = sum(dscores);
    
    %% compute dm
    dmTstats = zeros(length(models)-1, 1);
    
    nwLag = 2;
    
    for m = 1 : length(models) - 1
        [~,dmTstats(m)] = dmtest(scores(:,end), scores(:,m), nwLag);
    end
    
    
    
    %% collect  scores into table
    scoresTABLE(:,s)       = sumdscores;
    scoresdmstatTABLE(:,s) = dmTstats;
    
    %% plot differences in scores
    if s == 1
        cumDscores = cumsum(dscores);
        theselabel = {models(1:end-1).prettylabel};
        
        if range(year(dates)) > 40
            tickyears = datenum(1980:10:2020,1,1);
        elseif range(year(dates)) > 10
            tickyears = datenum(1975:5:2025,1,1);
        else
            tickyears = [];

        end
        
      
        thisfig = figure;
        hold on
        set(gca, 'fontsize', 18)
        hh = plot(dates, cumDscores, 'linewidth', 2);
        nbershades(dates)
        xlim(dates([1 end]))
        if ~isempty(tickyears)
            xticks(tickyears)
        end
        datetick('x', 'keepticks', 'keeplimits')
        plotOrigin('k:', [], [], 2)
        legend(hh, theselabel, 'location', 'best', 'box', 'on')
        wrapthisfigure(thisfig, sprintf('cumlogscores-%s', samLabel), wrap, [], [], [], [], true)
        title(sprintf('cumulative scores relative to %s', models(end).prettylabel))
        wrapthisfigure(thisfig, sprintf('cumlogscores-%s-WITHTITLE', samLabel), wrap)
        
    end
    
    
    %% finish sam
    dockAllFigures
end % sam


%% tabulate scores (PAPER VERSION)
tabname    = 'MVlogscoresTableSV-AR1.tex';
tabcaption = sprintf('Cumulative differences in predictive log-scores (vs. %s)', models(end).prettylabel);

if isempty(wrap)
    tabdir = pwd;
else
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
end

% write table
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
% fprintf(fid, '\\small\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.3', 1, length(models)-1));
fprintf(fid, '\\toprule\n');
fprintf(fid, '& \\multicolumn{%d}{c}{\\bf Models} ', length(models) - 1);
fprintf(fid, '\\\\\\cmidrule{2-%d}', 1 + length(models) - 1);
fprintf(fid, '\n');

fprintf(fid, '{\\bf Samples} ');
for m = 1 : length(models) - 1
    fprintf(fid, '& \\multicolumn{1}{c}{\\bf %s} ', models(m).prettylabel);
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');

for s = 1 : length(sam)
    
    [~,maxndx]= max(scoresTABLE(:,s));
    
    samLabel = sam(s).label;
    fprintf(fid, '{\\bf %s} ', sam(s).prettylabel);
    fprintf(fid, '\\\\\n');
    fprintf(fid, '{\\quad %s} ', samLabel);
    for m = 1 : length(models) - 1
        fprintf(fid, '& %6.2f ', scoresTABLE(m,s));
    end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');


% fprintf(fid, 'Significance of differences relative to model %s assessed by Diebold-Mariano test using Newey-West standard errors with $2$ lags.\n', ...
%     models(end).prettylabel);
fclose(fid);
type(fullfile(tabdir, tabname))





%% finish script
finishwrap
finishscript
