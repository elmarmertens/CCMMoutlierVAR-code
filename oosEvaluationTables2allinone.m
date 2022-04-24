%% compare OOS results from two estimates
% load quantico*.mat files and assess OOS performance

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

doCharts = false;

%% one big wrapper

titlename = 'oosPairwiseEvaluationTables';
initwrap
    
%% eval window

for sam = 1 : 4
    
    switch sam
        case 1
            % baseline
            evalStart  = datenum(1975,1,1);
            evalStop   = datenum(2017,12,1);
        case 2
            % GFC
            evalStart = datenum(2007,1,1);
            evalStop  = datenum(2014,12,1);
        case 3
            % 2020
            evalStart = datenum(2020,3,1);
            evalStop  = datenum(2021,3,1);

    end
    
    
    evaltxt = sprintf('evalStart%sevalEnd%s', datestr(evalStart, 'yyyymm'), datestr(evalStop, 'yyyymm'));
    
    
    %% define models
    % CONST
    m = 1;
    models(m).datalabel    = 'fredMD16-2021-04'; %#ok<*SAGROW>
    models(m).resultlabel  = 'censoredYields-CONST-p12';
    models(m).prettylabel  = 'CONST';
    models(m).shortlabel   = 'CONST';
    
    % SV
    m = 2;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SV-p12';
    models(m).prettylabel  = 'SV';
    models(m).shortlabel   = 'SV';
    
    % SVO
    m = 3;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVOmax20-p12';
    models(m).prettylabel  = 'SVO';
    models(m).shortlabel   = 'SVOar1';
    
    % SV-t
    m = 4;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVt-p12';
    models(m).prettylabel  = 'SV-t';
    models(m).shortlabel   = 'SVt';
    
    % SVO-t
    m = 5;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVOtmax20-p12';
    models(m).prettylabel  = 'SVO-t';
    models(m).shortlabel   = 'SVOt';
    
    % SV-nanOutlier
    m = 6;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVnanO5-p12';
    models(m).prettylabel  = 'SV-OutMiss';
    models(m).shortlabel   = 'SVoutmiss';
    
    % Obar
    m = 7;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVobarmax20-p12';
    models(m).prettylabel  = 'SV-common outlier';
    models(m).shortlabel   = 'SVobar';
    

    
    
    %% define model sets
    
    MODELS  = {[2 7] [3 7] [2 3]};
    MLABELS = {'SV vs SVo', 'SVO vs SVo', 'SV vs SVO'};
    
    
    
    
    %% loop over model sets
    
    for m = 1 : length(MODELS) 
        
        
        
        m0 = MODELS{m}(1);
        m1 = MODELS{m}(2);
        
        
        doQuarterly = false;
        
        %#ok<*UNRCH>
        %% load data
        clear oos0 oos1 oos2 ydates Tjumpoffs
        
        varlist = {'ydates', 'Tjumpoffs', 'N', 'tcode', 'ncode', 'cumcode', 'MCMCdraws', ...
            'fcstNhorizons', 'fcstYrealized', 'fcstYhaterror', 'fcstYmederror', 'fcstCRPS', ...
            'fcstYhat', 'fcstYmedian', ...
            };
        
        oos0 = load(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)), ...
            varlist{:});
        oos1 = load(fullfile(resultsdir, sprintf('%s-%s.mat', models(m1).datalabel, models(m1).resultlabel)), ...
            varlist{:});
        
        if oos0.MCMCdraws ~= oos1.MCMCdraws
            warning('unequal numbers of MCMCdraws, 0 has %d, 1 has %d', oos0.MCMCdraws, oos1.MCMCdraws)
        end
        
        
        %% check for identical samples
        if oos0.ydates(1) ~= oos1.ydates(1)
            error('oos estimates based on different sample starts: %s vs %s', datestr(oos0.ydates(1)), datestr(oos1.ydates(1)))
        end
        ydates    = oos0.ydates;
        
        if ~isequal(oos0.Tjumpoffs, oos1.Tjumpoffs)
            error('oos jumpoffs differ')
        end
        
        Tjumpoffs = oos0.Tjumpoffs;
        
        %% cut eval sample if desired
        ndxJumpoff = ismember(Tjumpoffs, find((ydates >= evalStart) & (ydates <= evalStop)));
        Tjumpoffs   = Tjumpoffs(ndxJumpoff);
        
        
        dates    = ydates(Tjumpoffs);
        
        %         comparisonNote = sprintf('Evaluation window from %s through %s.', ...
        %             datestryymm(dates(1)), datestryymm(dates(end)));
        comparisonNote = sprintf('Evaluation window from %s through %s.', ...
            datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));
        shortcomparisonNote = sprintf('from %s through %s.', ...
            datestr(dates(1), 'yyyy:mm'), datestr(dates(end), 'yyyy:mm'));
        
        
        %% some parameters
        Nhorizons  = oos0.fcstNhorizons;
        
        if doCharts
            Ylabels = fredMDshortlabel(oos0.ncode);
        else
            Ylabels = fredMDprettylabel(oos0.ncode);
        end
        N       = length(Ylabels);
        
        Ylabels = strrep(Ylabels, '_', '');
        
        %% setup monthly tables
        if sam >= 3
            theseHorizons = [1 3 6 9];
        else
            theseHorizons = [1 3 12 24];
        end
        
        
        doBold      = @(x) round(x * 100) < 100;
        
        %% RMSE
        
        mseloss0 = oos0.fcstYhaterror.^2;
        mseloss1 = oos1.fcstYhaterror.^2;
        
        % match samples
        mseloss0 = mseloss0(:,:,ndxJumpoff);
        mseloss1 = mseloss1(:,:,ndxJumpoff);
        
        RMSE0 = sqrt(nanmean(mseloss0,3));
        RMSE1 = sqrt(nanmean(mseloss1,3));
        relativeRMSE01 =  RMSE1 ./ RMSE0; % here: RMSE
        
        
        %% MAE
        maeloss0     = abs(oos0.fcstYmederror);
        maeloss1     = abs(oos1.fcstYmederror);
        
        % match samples
        maeloss0 = maeloss0(:,:,ndxJumpoff);
        maeloss1 = maeloss1(:,:,ndxJumpoff);
        
        mae0      = nanmean(maeloss0,3);
        mae1      = nanmean(maeloss1,3);
        relativeMAD01 = mae1 ./ mae0; % here: RMAE
        
        
        
        %% CRPS
        crpsloss0     = oos0.fcstCRPS;
        crpsloss1     = oos1.fcstCRPS;
        % match samples
        crpsloss0 = crpsloss0(:,:,ndxJumpoff);
        crpsloss1 = crpsloss1(:,:,ndxJumpoff);
        
        
        relativeCRPS01 = nanmean(crpsloss1,3) ./ nanmean(crpsloss0,3); % here: Relative CRPS
        
        
        %% compare all
        statlabels = {'RMSE', 'MAE', 'CRPS'};
        if doCharts
            tabname = sprintf('allinoneChart-%s-%s-vs-%s-%s.tex', models(m0).datalabel, ...
                models(m0).shortlabel, models(m1).shortlabel, evaltxt);
        else
            tabname = sprintf('allinone-%s-%s-vs-%s-%s.tex', models(m0).datalabel, ...
                models(m0).shortlabel, models(m1).shortlabel, evaltxt);
        end
        tabcaption =  sprintf('%s %s', MLABELS{m}, shortcomparisonNote);
        compareAllinone(tabname, wrap, doCharts, ...
            mseloss0, mseloss1, relativeRMSE01, ...
            maeloss0, maeloss1, relativeMAD01, ...
            crpsloss0, crpsloss1, relativeCRPS01, ...
            models(m0).prettylabel, models(m1).prettylabel, ...
            Ylabels, theseHorizons, tabcaption, statlabels, doBold, comparisonNote)
        
        
    end
end

%% finish script
finishwrap
finishscript


function compareAllinone(tabname, wrap, doCharts, ...
    mseloss0, mseloss1, relativeRMSE01, ...
    maeloss0, maeloss1, relativeMAE01, ...
    crpsloss0, crpsloss1, relativeCRPS01, ...
    prettylabel0, prettylabel1, ...
    Ylabels, theseHorizons, tabcaption, statlabels, doBold, comparisonNote)


%% parse inputs
N = length(Ylabels);
Nhorizons = length(theseHorizons);

%% DM tests

[relativeRMSE01, dmRMSEtstat] = dodm(mseloss0, mseloss1, relativeRMSE01, theseHorizons);

[relativeMAE01, dmMADtstat] = dodm(maeloss0, maeloss1, relativeMAE01, theseHorizons);

[relativeCRPS01, dmCRPStstat] = dodm(crpsloss0, crpsloss1, relativeCRPS01, theseHorizons);


%% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
end

%% tabulate
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
if doCharts
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.3', 1, 3 * Nhorizons));
    fprintf(fid, ' & \\multicolumn{%d}{c}{\\bf %s}  & \\multicolumn{%d}{c}{\\bf %s}  & \\multicolumn{%d}{c}{\\bf %s} \\\\ \\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \n', ...
        Nhorizons, statlabels{1}, Nhorizons, statlabels{2}, Nhorizons, statlabels{3}, ...
        1+1,1+Nhorizons,1+Nhorizons+1,1+2*Nhorizons,1+2*Nhorizons+1, 1+3*Nhorizons);
    fprintf(fid, 'Var. / Hor. ');
    
else
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, 3 * Nhorizons));
    fprintf(fid, ' & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s} \\\\ \\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \n', ...
        Nhorizons, statlabels{1}, Nhorizons, statlabels{2}, Nhorizons, statlabels{3}, ...
        1+1,1+Nhorizons,1+Nhorizons+1,1+2*Nhorizons,1+2*Nhorizons+1, 1+3*Nhorizons);
    fprintf(fid, 'Variable / Horizon ');
end
for h = 1 : Nhorizons
    fprintf(fid, '& \\multicolumn{1}{c}{$%d$} ', theseHorizons(h));
end
for h = 1 : Nhorizons
    fprintf(fid, '& \\multicolumn{1}{c}{$%d$} ', theseHorizons(h));
end
for h = 1 : Nhorizons
    fprintf(fid, '& \\multicolumn{1}{c}{$%d$} ', theseHorizons(h));
end
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for n = 1 : N
    fprintf(fid, '%s ', Ylabels{n});
    for h = 1 : Nhorizons
        if isfinite(relativeRMSE01(n,h))
            if doBold(relativeRMSE01(n,h))
                fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeRMSE01(n,h), Zstar(dmRMSEtstat(n,h)))));
                %             elseif doRed(relativeRMSE01(n,h))
                %                 fprintf(fid, '& %s ', dcolred(sprintf('%6.2f%s ', relativeRMSE01(n,h), Zstar(dmRMSEtstat(n,h)))));
            else
                fprintf(fid, '& %6.2f%s ', relativeRMSE01(n,h), Zstar(dmRMSEtstat(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeMAE01(n,h))
            if doBold(relativeMAE01(n,h))
                fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)))));
            else
                fprintf(fid, '& %6.2f%s ', relativeMAE01(n,h), Zstar(dmMADtstat(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    for h = 1 : Nhorizons
        if isfinite(relativeCRPS01(n,h))
            if doBold(relativeCRPS01(n,h))
                fprintf(fid, '& %s ', dcolgreen(sprintf('%6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)))));
            else
                fprintf(fid, '& %6.2f%s ', relativeCRPS01(n,h), Zstar(dmCRPStstat(n,h)));
            end
        else
            fprintf(fid, '& -.- ');
        end
    end
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');

if ~doCharts
    sigone = cat(1, abs(dmRMSEtstat) > norminv(0.95, 0, 1) & (round(relativeRMSE01,2) == 1), ...
        abs(dmMADtstat) > norminv(0.95, 0, 1) & (round(relativeMAE01,2) == 1), ...
        abs(dmCRPStstat) > norminv(0.95, 0, 1) & (round(relativeCRPS01,2) == 1));
    
    fprintf(fid, 'Note: Comparison of ``%s'''' (baseline, in denominator) against ``%s.'''' Values below 1 indicate improvement over baseline. \n', ...
        prettylabel0, prettylabel1);
    fprintf(fid, '%s \n', comparisonNote);
    
    fprintf(fid, 'Significance assessed by Diebold-Mariano-West test using Newey-West standard errors with $h + 1$ lags.\n');
    
    if any(sigone, 'all')
        if sum(sigone(:)) > 1
            fprintf(fid, 'Due to the close behavior of some of the models compared, and rounding of the reported values, a few comparisons show significant ratios  of 1.00.\n');
            fprintf(fid, 'These cases arise from persistent differences in performance that are, however, too small to be relevant after rounding.\n');
        else
            fprintf(fid, 'Due to the close behavior of some of the models compared, and rounding of the reported values, one of the comparisons shows a significant ratio of 1.00.\n');
            fprintf(fid, 'This case arises from persistent differences in performance that are, however, too small to be relevant after rounding.\n');
        end
    end
    
    
    if ~all(isfinite(relativeMAE01(:)))
        fprintf(fid, 'In some cases, due to strong performance of the baseline model, relative MAD may involve divisions by zero. These cases are reported as blank entries.');
    end
end
fclose(fid);

type(fullfile(tabdir, tabname))

end

function [deltaLoss, tstat] = dodm(loss0, loss1, deltaLoss, theseHorizons)

loss0 = loss0(:,theseHorizons,:);
loss1 = loss1(:,theseHorizons,:);
deltaLoss = deltaLoss(:,theseHorizons);

[N, Nhorizons,~] = size(loss0);

tstat = NaN(N,Nhorizons);

for h = 1 : Nhorizons
    nwLag = theseHorizons(h) + 1;
    for n = 1 : N
        thisloss0 = squeeze(loss0(n,h,:));
        thisloss1 = squeeze(loss1(n,h,:));
        
        if isequaln(thisloss0, thisloss1) || any(isinf(thisloss0)) || any(isinf(thisloss1))
            % do noting
        else
            [~,tstat(n,h)] = dmtest(thisloss0,thisloss1, nwLag);
        end
    end
end

end