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

resultsdir     = pwd;

%% eval window

for sam =  1 : 4
    switch sam
        case 1
            % baseline
            evalStart  = [];
            evalStop   = datenum(2017,12,1);
        case 2
            % GFC
            evalStart = datenum(2007,1,1);
            evalStop  = datenum(2014,12,1);
        case 3
            % 2020
            evalStart = datenum(2020,3,1);
            evalStop   = [];
        case 4 % 1985-2017
            % baseline
            evalStart  = datenum(1985,1,1);
            evalStop   = datenum(2017,12,1);
    end
    titlename = 'oosEvaluationTablesSVAR1';
    if ~isempty(evalStart)
        titlename = strcat(titlename, '-from', datestr(evalStart, 'yyyymmm'));
    end
    if ~isempty(evalStop)
        titlename = strcat(titlename, '-to', datestr(evalStop, 'yyyymmm'));
    end
    initwrap
    
    
    %% list of models
    % CONST
    m = 1;
    models(m).datalabel    = 'fredMD16-2021-04'; %#ok<*SAGROW>
    models(m).resultlabel  = 'censoredYields-CONST-p12';
    models(m).prettylabel  = 'CONST';
    
    % SV
    m = 2;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVar1-p12';
    models(m).prettylabel  = 'SV-AR(1)';
    
    % SVO
    m = 3;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVOar1max20-p12';
    models(m).prettylabel  = 'SV-AR(1)-O';
    
    % SV-t
    m = 4;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVtar1-p12';
    models(m).prettylabel  = 'SV-AR(1)-t';
    
    % SVO-t
    m = 5;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVOtar1max20-p12';
    models(m).prettylabel  = 'SV-AR(1)-Ot';
    
    % SV-nanOutlier
    m = 6;
    models(m).datalabel    = 'fredMD16-2021-04';
    models(m).resultlabel  = 'censoredYields-SVar1nanO5-p12';
    models(m).prettylabel  = 'SV-AR(1)-OutMiss';

    
	%% define model sets
    
    MODELS  = {[2 3 4 5] [2 3 5 6]};
    MLABELS = {'outlier-augmented SV', 'SV-OutMiss'};


    doLogscore  = false;
    
    %% loop over MODELS sets
    for M = 1 : length(MODELS)
        
        m0 = MODELS{M}(1);
        m1 = MODELS{M}(2);
        m2 = MODELS{M}(3);
        m3 = MODELS{M}(4);
        
        
        
        doQuarterly = contains(models(1).datalabel, 'quarterly');
        
        %#ok<*UNRCH>
        %% load data
        clear oos0 oos1 oos2 oos3 ydates Tjumpoffs
        
        varlist = {'ydates', 'Tjumpoffs', 'N', 'tcode', 'ncode', 'cumcode', 'MCMCdraws', ...
            'fcstNhorizons', 'fcstYrealized', 'fcstYhaterror', 'fcstYmederror', 'fcstCRPS', 'fcstLogscore', ...
            'fcstYquantiles', 'setQuantiles', ...
            'fcstYhat', 'fcstYhatRB', 'fcstYmedian', ...
            };
        
        
        oos0 = load(fullfile(resultsdir, sprintf('%s-%s.mat', models(m0).datalabel, models(m0).resultlabel)), ...
            varlist{:});
        oos1 = load(fullfile(resultsdir, sprintf('%s-%s.mat', models(m1).datalabel, models(m1).resultlabel)), ...
            varlist{:});
        oos2 = load(fullfile(resultsdir, sprintf('%s-%s.mat', models(m2).datalabel, models(m2).resultlabel)), ...
            varlist{:});
        oos3 = load(fullfile(resultsdir, sprintf('%s-%s.mat', models(m3).datalabel, models(m3).resultlabel)), ...
            varlist{:});
        
       
        
        %% check inputs
        if oos0.MCMCdraws ~= oos1.MCMCdraws
            warning('unequal numbers of MCMCdraws, 0 has %d, 1 has %d', oos0.MCMCdraws, oos1.MCMCdraws)
        end
        if oos0.MCMCdraws ~= oos2.MCMCdraws
            warning('unequal numbers of MCMCdraws, 0 has %d, 1 has %d', oos0.MCMCdraws, oos2.MCMCdraws)
        end
        
        if oos0.setQuantiles ~= oos1.setQuantiles
            error('unequal setQuantiles')
        end
        if oos0.setQuantiles ~= oos2.setQuantiles
            error('unequal setQuantiles')
        end
        if oos0.setQuantiles ~= oos3.setQuantiles
            error('unequal setQuantiles')
        end
        setQuantiles = oos0.setQuantiles;
        
        %% check for identical samples
        if oos0.ydates(1) ~= oos1.ydates(1)
            error('oos estimates based on different sample starts: %s vs %s', datestr(oos0.ydates(1)), datestr(oos1.ydates(1)))
        end
        if oos0.ydates(1) ~= oos2.ydates(1)
            error('oos estimates based on different sample starts: %s vs %s', datestr(oos0.ydates(1)), datestr(oos2.ydates(1)))
        end
        if oos0.ydates(1) ~= oos3.ydates(1)
            error('oos estimates based on different sample starts: %s vs %s', datestr(oos0.ydates(1)), datestr(oos3.ydates(1)))
        end
        YDATES    = oos0.ydates;
        
        if ~isequal(oos0.Tjumpoffs, oos1.Tjumpoffs)
            error('oos jumpoffs differ')
        end
        if ~isequal(oos0.Tjumpoffs, oos2.Tjumpoffs)
            error('oos jumpoffs differ')
        end
        if ~isequal(oos0.Tjumpoffs, oos3.Tjumpoffs)
            error('oos jumpoffs differ')
        end
        
        TJUMPOFFS = oos0.Tjumpoffs;
        
        
        %% cut eval sample if desired
        if isempty(evalStop)
            evalStop = YDATES(TJUMPOFFS(end));
        end
        if isempty(evalStart)
            evalStart = YDATES(TJUMPOFFS(1));
        end
        
        ndxJumpoff = ismember(TJUMPOFFS, find((YDATES >= evalStart) & (YDATES <= evalStop)));
        
        TJUMPOFFS   = TJUMPOFFS(ndxJumpoff);
        
        dates    = YDATES(TJUMPOFFS);
        
        samtxt = sprintf('from%s-to%s', datestr(dates(1), 'yyyymm'), datestr(dates(end), 'yyyymm'));
        
        comparisonNote = sprintf('Evaluation window from %s through %s.', ...
            datestryymm(dates(1)), datestryymm(dates(end)));
        
        
        %% some parameters
        Nhorizons  = oos0.fcstNhorizons;
        
        Ylabels = fredMDprettylabel(oos0.ncode);
        N       = length(Ylabels);
        
        Ylabels = strrep(Ylabels, '_', '');
        
        %% setup monthly tables
        if dates(1) > datenum(2019,1,1)
            theseHorizons = [1 3 6 9];
        else
            theseHorizons = [1 3 12 24];
        end
        
        
        %% RMSE
        thisStat = 'RMSE';
        tabname = sprintf('relative%s-%s-%s-SVAR1set%d.tex', thisStat, models(m0).datalabel, samtxt, M);
        
        loss0 = oos0.fcstYhaterror.^2;
        loss1 = oos1.fcstYhaterror.^2;
        loss2 = oos2.fcstYhaterror.^2;
        loss3 = oos3.fcstYhaterror.^2;
        
        % match samples
        loss0 = loss0(:,:,ndxJumpoff);
        loss1 = loss1(:,:,ndxJumpoff);
        loss2 = loss2(:,:,ndxJumpoff);
        loss3 = loss3(:,:,ndxJumpoff);
        
        RMSE0 = sqrt(nanmean(loss0,3));
        RMSE1 = sqrt(nanmean(loss1,3));
        RMSE2 = sqrt(nanmean(loss2,3));
        RMSE3 = sqrt(nanmean(loss3,3));
        deltaLoss01 =  RMSE1 ./ RMSE0; % here: RMSE
        deltaLoss02 =  RMSE2 ./ RMSE0; % here: RMSE
        deltaLoss03 =  RMSE3 ./ RMSE0; % here: RMSE
        tabcaption = sprintf('Relative %s', thisStat);
        if ~isempty(MLABELS{M})
            tabcaption = strcat(tabcaption, sprintf(' (%s)', MLABELS{M}));
        end
        compareLosses4(tabname, wrap, dates, loss0, loss1, loss2, loss3, deltaLoss01, deltaLoss02, deltaLoss03, ...
            models(m0).prettylabel, models(m1).prettylabel, models(m2).prettylabel, models(m3).prettylabel, ...
            Ylabels, theseHorizons, thisStat, tabcaption, comparisonNote)
        
        
        %% MAE
        thisStat = 'MAE';
        tabname = sprintf('relative%s-%s-%s-SVAR1set%d.tex', thisStat, models(m0).datalabel, samtxt, M);
        
        loss0     = abs(oos0.fcstYmederror);
        loss1     = abs(oos1.fcstYmederror);
        loss2     = abs(oos2.fcstYmederror);
        loss3     = abs(oos3.fcstYmederror);
        
        % match samples
        loss0 = loss0(:,:,ndxJumpoff);
        loss1 = loss1(:,:,ndxJumpoff);
        loss2 = loss2(:,:,ndxJumpoff);
        loss3 = loss3(:,:,ndxJumpoff);
        
        mid0      = nanmean(loss0,3);
        mid1      = nanmean(loss1,3);
        mid2      = nanmean(loss2,3);
        mid3      = nanmean(loss3,3);
        deltaLoss01 = mid1 ./ mid0; % here: RMAE
        deltaLoss02 = mid2 ./ mid0; % here: RMAE
        deltaLoss03 = mid3 ./ mid0; % here: RMAE
        tabcaption = sprintf('Relative %s', thisStat);
        if ~isempty(MLABELS{M})
            tabcaption = strcat(tabcaption, sprintf(' (%s)', MLABELS{M}));
        end
        
        compareLosses4(tabname, wrap, dates, loss0, loss1, loss2, loss3, deltaLoss01, deltaLoss02, deltaLoss03, ...
            models(m0).prettylabel, models(m1).prettylabel, models(m2).prettylabel, models(m3).prettylabel, ...
            Ylabels, theseHorizons, thisStat, tabcaption, comparisonNote)
        
        %% CRPS
        thisStat = 'CRPS';
        tabname = sprintf('relative%s-%s-%s-SVAR1set%d.tex', thisStat, models(m0).datalabel, samtxt, M);
        
        loss0     = oos0.fcstCRPS;
        loss1     = oos1.fcstCRPS;
        loss2     = oos2.fcstCRPS;
        loss3     = oos3.fcstCRPS;
        % match samples
        loss0 = loss0(:,:,ndxJumpoff);
        loss1 = loss1(:,:,ndxJumpoff);
        loss2 = loss2(:,:,ndxJumpoff);
        loss3 = loss3(:,:,ndxJumpoff);
        
        crps0       = nanmean(loss0,3);
        deltaLoss01 = nanmean(loss1,3) ./ crps0; % here: Relative CRPS
        deltaLoss02 = nanmean(loss2,3) ./ crps0; % here: Relative CRPS
        deltaLoss03 = nanmean(loss3,3) ./ crps0; % here: Relative CRPS
        
        tabcaption = sprintf('Relative %s', thisStat);
        if ~isempty(MLABELS{M})
            tabcaption = strcat(tabcaption, sprintf(' (%s)', MLABELS{M}));
        end
        
        compareLosses4(tabname, wrap, dates, loss0, loss1, loss2, loss3, deltaLoss01, deltaLoss02, deltaLoss03, ...
            models(m0).prettylabel, models(m1).prettylabel, models(m2).prettylabel, models(m3).prettylabel, ...
            Ylabels, theseHorizons, thisStat, tabcaption, comparisonNote)
        
        %% Logscore
        if doLogscore
            
            thisStat = 'Logscore';
            tabname = sprintf('relative%s-%s-%s-SVAR1set%d.tex', thisStat, models(m0).datalabel, samtxt, M);
            
            loss0     = oos0.fcstLogscore;
            loss1     = oos1.fcstLogscore;
            loss2     = oos2.fcstLogscore;
            loss3     = oos3.fcstLogscore;
            % match samples
            loss0 = loss0(:,:,ndxJumpoff);
            loss1 = loss1(:,:,ndxJumpoff);
            loss2 = loss2(:,:,ndxJumpoff);
            loss3 = loss3(:,:,ndxJumpoff);
            
            
            score0 = nanmean(loss0,3);
            deltaLoss01 = nanmean(loss1,3) - score0;
            deltaLoss02 = nanmean(loss2,3) - score0;
            deltaLoss03 = nanmean(loss2,3) - score0;
            
            tabcaption = 'Log-score differences';
            if ~isempty(MLABELS{M})
                tabcaption = strcat(tabcaption, sprintf(' (%s)', MLABELS{M}));
            end
            
            compareLosses4(tabname, wrap, dates, loss0, loss1, loss2, loss3, deltaLoss01, deltaLoss02, deltaLoss03, ...
                models(m0).prettylabel, models(m1).prettylabel, models(m2).prettylabel, models(m3).prettylabel, ...
                Ylabels, theseHorizons, thisStat, tabcaption, comparisonNote)
            
        end
    end
    
    %% finish script
    dockAllFigures
    finishwrap
end
finishscript

function compareLosses4(tabname, wrap, dates, loss0, loss1, loss2, loss3, deltaLoss01, deltaLoss02, deltaLoss03, prettylabel0, prettylabel1, prettylabel2, prettylabel3, Ylabels, theseHorizons, thisStat, tabcaption, comparisonNote) %#ok<INUSL>

%% DM test
loss0 = loss0(:,theseHorizons,:);
loss1 = loss1(:,theseHorizons,:);
loss2 = loss2(:,theseHorizons,:);
loss3 = loss3(:,theseHorizons,:);
deltaLoss01 = deltaLoss01(:,theseHorizons);
deltaLoss02 = deltaLoss02(:,theseHorizons);
deltaLoss03 = deltaLoss03(:,theseHorizons);

[N, Nhorizons,Nobs] = size(loss0);

doDMW = Nobs > 25;

if doDMW
    dmTstat01 = NaN(N,Nhorizons);
    for h = 1 : Nhorizons
        nwLag = theseHorizons(h) + 1;
        for n = 1 : N
            thisloss0 = squeeze(loss0(n,h,:));
            thisloss1 = squeeze(loss1(n,h,:));
            
            if any(isinf(thisloss0)) || any(isinf(thisloss1))
                % do noting
            else
                [~,dmTstat01(n,h)] = dmtest(thisloss0,thisloss1, nwLag);
            end
        end
    end
    
    dmTstat02 = NaN(N,Nhorizons);
    for h = 1 : Nhorizons
        nwLag = theseHorizons(h) + 1;
        for n = 1 : N
            thisloss0 = squeeze(loss0(n,h,:));
            thisloss2 = squeeze(loss2(n,h,:));
            
            if any(isinf(thisloss0)) || any(isinf(thisloss2))
                % do noting
            else
                [~,dmTstat02(n,h)] = dmtest(thisloss0,thisloss2, nwLag);
            end
            
            %             if abs(dmTstat02(n,h)) > 2 && round(deltaLoss02(n,h),2) == 1
            %                 thisfig = figure;
            %                 subplot(2,1,1)
            %                 hold on
            %                 plot(dates, thisloss0, '-')
            %                 plot(dates, thisloss2, '-.')
            %                 legend('baseline', 'alternative')
            %                 datetick('x')
            %                 xlim(dates([1 end]))
            %                 subplot(2,1,2)
            %                 hold on
            %                 h1 = plot(dates, thisloss0 - thisloss2);
            %                 h2 = plot(dates, cumsum(thisloss0 - thisloss2) ./ (1:length(dates))');
            %                 plothorzline(0, [], 'k:')
            %                 legend([h1 h2], 'difference', 'cumulated avg diff')
            %                 datetick('x')
            %                 xlim(dates([1 end]))
            %                 sgtitle(sprintf('%s (h=%d) \n %s' , Ylabels{n}, theseHorizons(h), tabname))
            %                 wrapthisfigure(thisfig, sprintf('var%d-h%d-%s' , n, theseHorizons(h), tabname), wrap)
            %             end
            
        end
    end
    dmTstat03 = NaN(N,Nhorizons);
    for h = 1 : Nhorizons
        nwLag = theseHorizons(h) + 1;
        for n = 1 : N
            thisloss0 = squeeze(loss0(n,h,:));
            thisloss3 = squeeze(loss3(n,h,:));
            
            if any(isinf(thisloss0)) || any(isinf(thisloss3))
                % do noting
            else
                [~,dmTstat03(n,h)] = dmtest(thisloss0,thisloss3, nwLag);
            end
        end
    end
else
    dmTstat01 = zeros(N,Nhorizons);
    dmTstat02 = zeros(N,Nhorizons);
    dmTstat03 = zeros(N,Nhorizons);
end

%% set up tab
if ~isempty(wrap)
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'sidetab', tabname, tabcaption)
end

%% tabulate
fid = fopen(fullfile(tabdir, tabname), 'wt');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.4', 1, 3 * Nhorizons));
fprintf(fid, '\\toprule\n');
fprintf(fid, ' & \\multicolumn{%d}{c}{%s}  & \\multicolumn{%d}{c}{%s} & \\multicolumn{%d}{c}{%s} \\\\ \\cmidrule(lr){2-%d}\\cmidrule(lr){%d-%d}\\cmidrule(lr){%d-%d} \n', ...
    Nhorizons, prettylabel1, Nhorizons, prettylabel2, Nhorizons, prettylabel3, 1+Nhorizons, 1+Nhorizons+1, 1+2*Nhorizons, 1+2*Nhorizons+1, 1+3*Nhorizons);
fprintf(fid, 'Variable / Horizon');
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
        if isnan(deltaLoss01(n,h))
            fprintf(fid, '& \\ccol{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaLoss01(n,h), Zstar(dmTstat01(n,h)));
        end
    end
    for h = 1 : Nhorizons
        if isnan(deltaLoss02(n,h))
            fprintf(fid, '& \\ccol{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaLoss02(n,h), Zstar(dmTstat02(n,h)));
        end
    end
    for h = 1 : Nhorizons
        if isnan(deltaLoss03(n,h))
            fprintf(fid, '& \\ccol{--} ');
        else
            fprintf(fid, '& %6.2f%s ', deltaLoss03(n,h), Zstar(dmTstat03(n,h)));
        end
    end
    
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');

sigone = cat(1, abs(dmTstat01) > norminv(0.95, 0, 1) & (round(deltaLoss01,2) == 1), ...
    abs(dmTstat02) > norminv(0.95, 0, 1) & (round(deltaLoss02,2) == 1));

if contains(lower(tabcaption), 'log-score')
    fprintf(fid, 'Note: Comparison of ``%s'''' (baseline, subtracted) against ``%s,'''' ``%s,'''' and ``%s.'''' Positive values indicate improvement over baseline. \n', ...
        prettylabel0, prettylabel1, prettylabel2, prettylabel3);
else
    fprintf(fid, 'Note: Comparison of ``%s'''' (baseline, in denominator) against ``%s,'''' ``%s,'''' and ``%s.'''' Values below 1 indicate improvement over baseline. \n', ...
        prettylabel0, prettylabel1, prettylabel2, prettylabel3);
end
fprintf(fid, '%s \n', comparisonNote);
if doDMW
    fprintf(fid, 'Significance assessed by Diebold-Mariano-West test using Newey-West standard errors with $h + 1$ lags.\n');
else
    fprintf(fid, 'Due to the low number of observations in the evaluation window, significance tests have not been performed.\n');
end
if ~contains(lower(tabcaption), 'log-score')
    if any(sigone, 'all')
        if sum(sigone(:)) > 1
            fprintf(fid, 'Due to the close behavior of some of the models compared, and rounding of the report values, a few comparisons show a significant relative %s of 1.00.\n', thisStat);
            fprintf(fid, 'These cases arise from persistent differences in performance that are, however, too small to be relevant after rounding.\n');
        else
            fprintf(fid, 'Due to the close behavior of some of the models compared, and rounding of the report values, one of the comparisons shows a significant relative %s of 1.00.\n', thisStat);
            fprintf(fid, 'This case arises from persistent differences in performance that are, however, too small to be relevant after rounding.\n');
        end
    end
end
fclose(fid);
type(fullfile(tabdir, tabname))

end

