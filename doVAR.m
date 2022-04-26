%% estimate a single VAR models for a given sample
% select modeltype in third cell below
% Carriero, Clark, Marcellino and Mertens (forthcoming, REStat)

%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>

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

Nstreams    = max(1,getparpoolsize);
rndStreams  = initRandStreams(Nstreams, [], 1);


%% set parameters for VAR and MCMC

modeltype = 'SVO';

jumpoffDate = datenum(2021,03,1);
% jumpoffDate = datenum(2020,04,1);
% jumpoffDate = datenum(1975,01,1);
% jumpoffDate = datenum(2019,12,1);

SVOmaxscale = 20;

datalabel           = 'fredMD16-2021-04';
% datalabel           = 'fredMD16levels-2021-04';
doQuarterly         = false;
MCMCdraws           = 1e3;               % Final number of MCMC draws after burn in
fcstNdraws          = 100 * MCMCdraws;    % draws sampled from predictive density
doCensorYields      = true;

samStart            = []; % datenum(1988,12,1);                 % truncate start of sample if desired (leave empty if otherwise)

% SED-PARAMETERS-HERE
MCMCdraws           = 1e2;               % Final number of MCMC draws after burn in
fcstNdraws          = MCMCdraws;    % draws sampled from predictive density


doStoreXL           = false; %#ok<*NASGU>
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)
Compute_diagnostics = false;              % compute Inefficiency Factors and Potential

doRobustPrior       = false;

if doCensorYields
    ELBbound            = .25;
else
    ELBbound            = [];
end

if doQuarterly
    p  = 4;
    np = 4; % number of periods per year, used for calibrating priors
    datalabel = strcat(datalabel, '-quarterly');
else
    p  = 12;
    np = 12;
end


%% load data

% load CSV file
dum=importdata(sprintf('%s.csv', datalabel),',');


ydates=dum.data(3:end,1);
% Variable names
ncode=dum.textdata(1,2:end);
% Transformation codes (data are already transformed)
tcode  =dum.data(1,2:end);
cumcode=logical(dum.data(2,2:end));
% Data
data=dum.data(3:end,2:end);

% define index of yields that need to obey ELB (out of sample)
if doCensorYields
    ndxYIELDS = find(ismember(ncode, {'FEDFUNDS', 'GS1', 'GS5', 'GS10'}));
else
    ndxOTHERYIELDS = [];
end

Tdata = length(ydates);

Ylabels = fredMDprettylabel(ncode);
Ylabelslong = Ylabels;
Ylabelslong(ismember(tcode, [2 5])) = strcat(Ylabelslong(ismember(tcode, [2 5])), ' (APR growth)');


%% process settings
N = size(data,2);

doTightPrior = false; % N >= 10;

Kbvar = N * p + 1; % number of regressors per equation
K = Kbvar;


% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

% define oos jump offs

thisT = find(ydates == jumpoffDate);

% other settings
fractiles = [5, normcdf(-1) * 100, 100 - normcdf(-1) * 100, 95]; % legacy, used for some plots

setQuantiles = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
Nquantiles    = length(setQuantiles);
ndxCI         = ismember(setQuantiles, fractiles);

%% mean for Minnesota prior: zero (diff) or RW (level)

if contains(lower(datalabel), 'levels')
    
    minnesotaPriorMean = ones(N,1);
    
else
    
    minnesotaPriorMean = NaN(N,1);
    
    for n = 1 : N
        switch ncode{n}
            case {'CUMFNS', 'UNRATE', ...
                    'WPSFD49207',  'PPICMM', 'PCEPI', ...
                    'FEDFUNDS', 'HOUST', 'GS1', 'GS5', 'GS10', 'BAAFFM', 'WUXIASHADOWRATE',...
                    'GDPCTPI', 'PCECTPI', 'GPDICTPI'}
                minnesotaPriorMean(n) = 1;
            otherwise
                minnesotaPriorMean(n) = 0;
        end
    end
end

%% allocate memory for out-of-sample forecasts
fcstNhorizons     = 24;  % number of steps forecasted (1:fcstNhorizon)

yearStart = year(ydates(1));
yearStop  = year(ydates(end)) + 9;

if doQuarterly
    fcstdates = genrQdates(yearStart, yearStop);
else
    fcstdates = genrMdates(yearStart, yearStop, 1);
end
ndx = find(fcstdates == ydates(1));
if isempty(ndx)
    error houston
end
fcstdates = fcstdates(ndx:end);


%% start latexwrapper to collect results
titlename=sprintf('%s-p%d', datalabel, p);
if ~isempty(samStart)
    titlename = strcat(titlename, '-', datestr(samStart, 'yyyymmm'));
end

titlename = strcat(titlename, '-', modeltype);
if contains(modeltype, 'SVO')
    titlename=strcat(titlename, sprintf('maxscale%d', SVOmaxscale));
end

if doQuarterly
    titlename = strcat(titlename, '-', datestr(ydates(thisT), 'yyyyQQ'), '-since',  datestr(ydates(1), 'yyyyQQ'));
else
    titlename = strcat(titlename, '-', datestr(ydates(thisT), 'yyyymmm'), '-since',  datestr(ydates(1), 'yyyymmm'));
end

wrap = [];
% initwrap

%% plot input data

% thesedates = ydates(1:thisT);
% Ofactor = 5;
% 
% for n = 1 : N
% 
%     dev = abs(data(1:thisT,n) - median(data(1:thisT,n)));
%     iqr = range(prctile(data(1:thisT,n), [25 75]));
%     outndx = dev > Ofactor * iqr;
%     
%     this = figure;
%     hold on
%     plot(thesedates, data(1:thisT,n))
%     xtickdates(thesedates)
% 
%     plot(thesedates(outndx), data(outndx,n), 'rd', 'linewidth', 2)
% 
%     %     title('full sample')
%     
%     %     subplot(2,1,2)
%     %     bar(ydates(thisT-23:thisT), data(thisT-23:thisT,n))
%     %     xtickdates([ydates(thisT-24:thisT); ydates(thisT+1) + mean(diff(ydates))])
%     %     datetick('x', 28, 'keeplimits')
%     %
%     %     title('most recent data')
%     
%     sgtitle(Ylabelslong{n})
%     
%     wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
% end

%% setup forecast origin

TID   = parid;
T     = thisT - p;

% collect realized values
yrealized = NaN(N, fcstNhorizons);
for h = 1 : fcstNhorizons
    if thisT + h < Tdata
        yrealized(:,h) = data(thisT+h,:)';
    end
end

% set Funds Rate equal to ELB when at ELB
% (Note: ELB may be set higher than actual funds rate readings, e.g.
% 25p)
if ~isempty(ELBbound)
    yieldsrealized         = yrealized(ndxYIELDS,:);
    ndx                    = yieldsrealized < ELBbound;
    yieldsrealized(ndx)    = ELBbound;
    yrealized(ndxYIELDS,:) = yieldsrealized;
end

% cumulate realizations and predictions if necessary
yrealized(cumcode,:) = cumsum(yrealized(cumcode,:),2);

%% MCMC sampler

tic

switch modeltype
    
    case 'CONST'
        [PAI_all, SIGMA_all, ~, ...
            ydraws, yhat, yhatRB] ...
            = mcmcVARconstvcv(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'CONSTnanO'
        
        Ofactor = 5;
        [PAI_all, SIGMA_all, ~, Ydatadraws, Ynan, ...
            ydraws, yhat, yhatRB] ...
            = mcmcVARconstvcvNaNoutlier(thisT, MCMCdraws, Ofactor, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'CONSTsansCOVID'
        [PAI_all, SIGMA_all, ~, ...
            ydraws, yhat, yhatRB] ...
            = mcmcVARconstvcvSansCOVID(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'CONSTdummy'
        
        dummyPrecision = 0;
        [PAI_all, SIGMA_all, ~, ...
            ydraws, yhat, yhatRB] ...
            = mcmcVARconstvcvDummy(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, dummyPrecision, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
        K = size(PAI_all,1);
        
    case 'SVdummy'
        
        dummyPrecision = 0;
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            ydraws, yhat, yhatRB] ...
            = mcmcVARSVdummy(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, dummyPrecision, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'SV'
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            ydraws, yhat, yhatRB, ~, logscoredraws] ...
            = mcmcVARSV(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            yrealized,...
            fcstNdraws, fcstNhorizons, rndStreams{TID}, true);
        
    case 'SVar1'
        [PAI_all, PHI_all, invA_all, SVrho_all, sqrtht_all, ...
            ydraws, yhat, yhatRB] ...
            = mcmcVARSVar1(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID}, true);
        
    case 'SVnanO'
        
        Ofactor = 5;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, Ydatadraws, Ynan, ...
            ydraws, yhat, yhatRB, ~, ydraws2, yhat2, logscoredraws] ...
            = mcmcVARSVnanOutlier(thisT, MCMCdraws, Ofactor, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            yrealized, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID}, true);
        
    case 'SVt'
        tdofGrid  = 3 : 40;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            SVtscalelog2_all, SVtdof_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws, logscoredraws] ...
            = mcmcVARSVt(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            tdofGrid, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            yrealized, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'SVtar1'
        tdofGrid  = 3 : 40;
        
        [PAI_all, PHI_all, invA_all, SVrho_all, sqrtht_all, ...
            SVtscalelog2_all, SVtdof_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws] ...
            = mcmcVARSVtar1(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            tdofGrid, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'SVOt'
        
        
        SVOpriorobs = 10 * np;
        SVOalpha    = 1 / (10 * np) * SVOpriorobs; % 10 years of data with 1 outlier every 10 years
        SVObeta     = SVOpriorobs - SVOalpha;
        
        tdofGrid  = 3 : 40;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            SVOprob_all, SVOscale_all, ...
            SVtscalelog2_all, SVtdof_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws, fcstSVoutliers] ...
            = mcmcVARSVOt(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVOalpha, SVObeta, SVOmaxscale, ...
            tdofGrid, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
            
    case 'SVOtar1'
        
        
        SVOpriorobs = 10 * np;
        SVOalpha    = 1 / (10 * np) * SVOpriorobs; % 10 years of data with 1 outlier every 10 years
        SVObeta     = SVOpriorobs - SVOalpha;
        
        tdofGrid  = 3 : 40;
        
        [PAI_all, PHI_all, invA_all, SVrho_all, sqrtht_all, ...
            SVOprob_all, SVOscale_all, ...
            SVtscalelog2_all, SVtdof_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws, fcstSVoutliers] ...
            = mcmcVARSVOtar1(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVOalpha, SVObeta, SVOmaxscale, ...
            tdofGrid, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'SVobar'
        SVobarpriorobs = 10 * np;
        SVobaralpha    = 1 / (4 * np) * SVobarpriorobs; % 10 years of data with 1 outlier every 4 years
        SVobarbeta     = SVobarpriorobs - SVobaralpha;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            SVobarprob_all, SVobarscale_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws] ...
            = mcmcVARSVObar(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVobaralpha, SVobarbeta, SVOmaxscale, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'SVOobar'
        SVOpriorobs = 10 * np;
        SVOalpha    = 1 / (4 * np) * SVOpriorobs; % 10 years of data with 1 outlier every 4 years
        SVObeta     = SVOpriorobs - SVOalpha;
        
        SVobarpriorobs = 10 * np;
        SVobaralpha    = 1 / (4 * np) * SVobarpriorobs; % 10 years of data with 1 outlier every 4 years
        SVobarbeta     = SVobarpriorobs - SVobaralpha;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            SVOprob_all, SVOscale_all, ...
            SVobarprob_all, SVobarscale_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws] ...
            = mcmcVARSVOobar(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVOalpha, SVObeta, SVOmaxscale, ...
            SVobaralpha, SVobarbeta, SVOmaxscale, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
    case 'SVO'
        SVOpriorobs = 10 * np;
        SVOalpha    = 1 / (4 * np) * SVOpriorobs; % 10 years of data with 1 outlier every 4 years
        SVOalpha    = repmat(SVOalpha, N, 1);
        SVObeta     = SVOpriorobs - SVOalpha;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            SVOprob_all, SVOscale_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws, fcstSVoutliers, logscoredraws] ...
            = mcmcVARSVO(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVOalpha, SVObeta, SVOmaxscale, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            yrealized, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID},true);
        
    case 'SVOar1'
        SVOpriorobs = 10 * np;
        SVOalpha    = 1 / (4 * np) * SVOpriorobs; % 10 years of data with 1 outlier every 4 years
        SVOalpha    = repmat(SVOalpha, N, 1);
        SVObeta     = SVOpriorobs - SVOalpha;
        
        [PAI_all, PHI_all, invA_all, SVrho_all, sqrtht_all, ...
            SVOprob_all, SVOscale_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws, fcstSVoutliers] ...
            = mcmcVARSVOar1(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVOalpha, SVObeta, SVOmaxscale, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID},true);
        
    otherwise
        error('modeltype <<%s>> not reognized')
end

toc

%% LOGSCORE

offsetlogscoredraw = max(logscoredraws);

figure
subplot(2,1,1)
histogram(exp(logscoredraws - offsetlogscoredraw))
title('exp(log score draws)')
subplot(2,1,2)
histogram(logscoredraws - offsetlogscoredraw)
title('log score draws')

yLogscore       = log(mean(exp(logscoredraws - offsetlogscoredraw))) + offsetlogscoredraw;
fprintf('Logscore: %6.4f\n', yLogscore)

%MV logscore
thesedraws  = squeeze(ydraws(:,1,:));
MU          = mean(thesedraws, 2);
Sigma       = cov(thesedraws', 1); % normalize variance by N
sqrtSigma   = chol(Sigma)';
logdetSigma = 2 * sum(log(diag(sqrtSigma)));
dev         = sqrtSigma \ (yrealized(:,1) - MU);
SSR         = dev' * dev;
yMVlogscore = -.5 * (N * log(2 * pi) + logdetSigma + SSR);

fprintf('Mean-variance Logscore: %6.4f\n', yMVlogscore)


%% Diagnostics

if Compute_diagnostics
    if contains(lower(modeltype), 'const')
        DiagnosticsCONST(PAI_all,SIGMA_all,N,K,MCMCdraws);
    else
        Diagnostics(sqrtht_all,invA_all,PAI_all,PHI_all,N,K,MCMCdraws);
    end
end

%% compute out-of-sample forecasts

ydraws(cumcode,:,:)  = cumsum(ydraws(cumcode,:,:),2);
yhat(cumcode,:)      = cumsum(yhat(cumcode,:),2);
yhatRB(cumcode,:)      = cumsum(yhatRB(cumcode,:),2);


% compute median
ymed = median(ydraws,3);

%% collect (missing) data draws
switch modeltype
    case {'SVnanO', 'CONSTnanO'}
        Ymid   = median(Ydatadraws, 3);
        Ytails = prctile(Ydatadraws, fractiles, 3);
        
        ydraws2(cumcode,:,:)  = cumsum(ydraws2(cumcode,:,:),2);
        yhat2(cumcode,:)      = cumsum(yhat2(cumcode,:),2);
        ymed2 = median(ydraws2,3);
        
        for n = 1 : N
            
            thisfig = figure;
            hold on
            %     plot(ydates(p+1:thisT),Ymid(:,n), 'k-', 'linewidth', 2);
            %     plot(ydates(p+1:thisT),squeeze(Ytails(:,n, [2 3])), 'k--', 'linewidth', 1);
            %     plot(ydates(p+1:thisT),squeeze(Ytails(:,n, [1 4])), 'k:', 'linewidth', 1);
            
            plotCI(Ymid(:,n), squeeze(Ytails(:,n, :)), ydates(p+1:thisT),[],  'k-', 'linewidth', 2);
            xtickdates(ydates)
            
            nanmid = Ymid(:,n);
            nanmid(~Ynan(:,n)) = NaN;
            plot(ydates(p+1:thisT), nanmid, 'rd', 'linewidth', 3)
            
            wrapthisfigure(thisfig, sprintf('dataX-%s', ncode{n}), wrap)
        end
        
end

%% plot outmiss predictive densities per Jumpoff for Outmiss
switch modeltype
    case {'SVnanO', 'CONSTnanO'}
        jj = find(fcstdates == ydates(thisT));
        if jj ~= thisT
            error houston
        end
        allthesedates = fcstdates(jj-np:jj+fcstNhorizons);
        theseFcstdates = fcstdates(jj+(1:fcstNhorizons));
        theseJumpoffdates = fcstdates(jj+(-np:0));
        for n = 1 : N
            
            hanni = NaN(7,1);
            
            thesedraws = squeeze(ydraws(n,:,:));
            
            thisfig = figure;
            
            fcstMid   = mean(thesedraws, 2);
            fcstMed   = median(thesedraws, 2);
            fcstTails = prctile(thesedraws, fractiles, 2);
            
            thesedraws2 = squeeze(ydraws2(n,:,:));
            fcstMid2   = mean(thesedraws2, 2);
            fcstTails2 = prctile(thesedraws2, fractiles, 2);
            
            %     subplot(2,1,1)
            hold on
            hanni(1) = plot(theseFcstdates,fcstMid, 'k-', 'linewidth', 3);
            %     hanni(2) = plot(theseFcstdates,yhatRB(n,:), 'bd', 'linewidth', 3);
            %     hanni(3) = plot(theseFcstdates,fcstMed, 'r-.', 'linewidth', 3);
            hanni(4:5) = plot(theseFcstdates,fcstTails(:, [2 3]), 'k--', 'linewidth', 2);
            %     hanni(6:7) = plot(theseFcstdates,fcstTails(:, [1 4]), 'k:', 'linewidth', 2);
            
            hanni(9) = plot(theseFcstdates,fcstMid2, 'r-', 'linewidth', 3);
            plot(theseFcstdates,fcstTails2(:, [2 3]), 'r:', 'linewidth', 2);
            
            hanni(8)  = plot(theseJumpoffdates, data(jj+(-np:0),n), 'g-', 'linewidth', 2);
            if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
            xtickdates(allthesedates)
            
            title(sprintf('%s \n per %s', Ylabelslong{n}, datestr(ydates(thisT))))
            
            
            
            
            
            legend(hanni([8 1 9]), 'data', 'mean', 'YhatMiss', ...
                'location', 'best')
            %     subplot(2,1,2)
            %     hold on
            %     hanni(1) = plot(theseFcstdates,fcstMid, 'k-', 'linewidth', 3);
            %     hanni(2) = plot(theseFcstdates,yhatRB(n,:), 'bd', 'linewidth', 3);
            %     hanni(3) = plot(theseFcstdates,fcstMed, 'r-.', 'linewidth', 3);
            %     plot(theseJumpoffdates, data(jj+(-np:0),n), 'g-', 'linewidth', 2);
            %     if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
            %     xtickdates(allthesedates)
            
            wrapthisfigure(thisfig, sprintf('outmisspredictiveDensity-%s', ncode{n}), wrap)
        end
end

%% plot predictive densities per Jumpoff
jj = find(fcstdates == ydates(thisT));
if jj ~= thisT
    error houston
end
allthesedates = fcstdates(jj-np:jj+fcstNhorizons);
theseFcstdates = fcstdates(jj+(1:fcstNhorizons));
theseJumpoffdates = fcstdates(jj+(-np:0));
for n = 1 : N
    
    hanni = NaN(7,1);
    
    thesedraws = squeeze(ydraws(n,:,:));
    
    thisfig = figure;
    
    fcstMid   = mean(thesedraws, 2);
    fcstMed   = median(thesedraws, 2);
    fcstTails = prctile(thesedraws, fractiles, 2);
    
    %     subplot(2,1,1)
    hold on
    hanni(1) = plot(theseFcstdates,fcstMid, 'k-', 'linewidth', 3);
    %     hanni(2) = plot(theseFcstdates,yhatRB(n,:), 'bd', 'linewidth', 3);
    hanni(3) = plot(theseFcstdates,fcstMed, 'r-.', 'linewidth', 3);
    hanni(4:5) = plot(theseFcstdates,fcstTails(:, [2 3]), 'k--', 'linewidth', 2);
    hanni(6:7) = plot(theseFcstdates,fcstTails(:, [1 4]), 'k:', 'linewidth', 2);
    hanni(8)  = plot(theseJumpoffdates, data(jj+(-np:0),n), 'g-', 'linewidth', 2);
    if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
    xtickdates(allthesedates)
    
    title(sprintf('%s \n per %s', Ylabelslong{n}, datestr(ydates(thisT))))
    
    
    legend(hanni([8 1 3 4 6]), 'data', 'mean', 'Median', '68% band', '90% band', ...
        'location', 'best')
    %     subplot(2,1,2)
    %     hold on
    %     hanni(1) = plot(theseFcstdates,fcstMid, 'k-', 'linewidth', 3);
    %     hanni(2) = plot(theseFcstdates,yhatRB(n,:), 'bd', 'linewidth', 3);
    %     hanni(3) = plot(theseFcstdates,fcstMed, 'r-.', 'linewidth', 3);
    %     plot(theseJumpoffdates, data(jj+(-np:0),n), 'g-', 'linewidth', 2);
    %     if min(ylim) < 0 && max(ylim) > 0, plotOrigin, end
    %     xtickdates(allthesedates)
    
    wrapthisfigure(thisfig, sprintf('predictiveDensity-%s', ncode{n}), wrap)
end


%% report Ainv
if exist('invA_all', 'var')
    
    medAinv = median(invA_all,3);
    meanAinv = mean(invA_all,3);
    
    ubound = 5;
    thisfig = figure;
    bar3(abs(medAinv))
    hold on
    [row, col] = find(abs(medAinv) > 5);
    plot3(col, row, repmat(ubound, length(row), 1), 'rd', 'linewidth', 5)
    xlim([0 N+1])
    ylim([0 N+1])
    zlim([0 5])
    xticks(2 : 2 : N);
    yticks(2 : 2 : N);
    xticks('manual')
    yticks('manual')
    title(sprintf('%s model: median(|A^{-1}|)', modeltype))
    wrapthisfigure(thisfig, sprintf('medAinv-%s', modeltype), wrap)
    
    thisfig = figure;
    bar3(abs(meanAinv))
    hold on
    [row, col] = find(abs(meanAinv) > 5);
    plot3(col, row, repmat(ubound, length(row), 1), 'rd', 'linewidth', 5)
    xlim([0 N+1])
    ylim([0 N+1])
    zlim([0 5])
    xticks(2 : 2 : N);
    yticks(2 : 2 : N);
    xticks('manual')
    yticks('manual')
    title(sprintf('%s model: mean(|A^{-1}|)', modeltype))
    wrapthisfigure(thisfig, sprintf('midAinv-%s', modeltype), wrap)
    
end



%% SV Graphs
switch modeltype
    
    case {'SV', 'SVdummy', 'SVar1', 'SVnanO'}
        SVdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = sqrtht_all(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVdraws(:,:,m) = thisSVdraw;
        end
        SVmid   = mean(SVdraws,3);
        SVtails = prctile(SVdraws, [5 95], 3);
        
        
        for n = 1 : N
            thisfig = figure;
            
            hold on
            h1 = plot(ydates(1:thisT), SVmid(n,:), 'k-', 'linewidth', 2);
            plot(ydates(1:thisT), squeeze(SVtails(n,:,:)), 'k-', 'linewidth', 1)
            title(sprintf('%s', Ylabels{n}))
            xtickdates(ydates(1:thisT))
            wrapthisfigure(thisfig, sprintf('SV%s', ncode{n}), wrap)
            
        end
        
    case {'SVt', 'SVtdof', 'SVtar1'}
        
        SVtdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = sqrtht_all(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVtdraws(:,:,m) = thisSVdraw;
        end
        SVtmid   = mean(SVtdraws,3);
        SVttails = prctile(SVtdraws, [5 95], 3);
        
        
        
        % compute SV sans outliers
        stochvol = sqrtht_all .* exp(-.5 * SVtscalelog2_all);
        SVdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = stochvol(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVdraws(:,:,m) = thisSVdraw;
        end
        SVmid   = mean(SVdraws,3);
        SVtails = prctile(SVdraws, [5 95], 3);
        
        
        for n = 1 : N
            thisfig = figure;
            
            hold on
            h1 = plot(ydates(1:thisT), SVtmid(n,:), 'r-', 'linewidth', 2);
            %             plot(ydates(1:thisT), squeeze(SVOtails(n,:,:)), 'r-', 'linewidth', 1)
            h2 = plot(ydates(1:thisT), SVmid(n,:), 'k-.', 'linewidth', 2);
            title(sprintf('%s', Ylabels{n}))
            xtickdates(ydates(1:thisT))
            legend([h1 h2], 'SV (incl. t outlier)', 'SV (proper)', 'location', 'best')
            wrapthisfigure(thisfig, sprintf('SV-%s-%s', modeltype, ncode{n}), wrap)
            
        end
        
        
        
        
        % plot density of SVt scale
        
        oscales = permute(exp(SVtscalelog2_all / 2), [3 2 1]);
        
        ogrid = 0 : .1 : 20;
        outlierpdf = NaN(length(ogrid),T,N);
        parfor n = 1 : N
            for t = 1 : T
                outlierpdf(:,t,n) = ksdensity(oscales(:,t,n), ogrid);
            end
        end
        
        for n = 1 : N
            thisfig = figure;
            surf(ydates(p+1:thisT), ogrid, outlierpdf(:,:,n));
            shading interp
            datetick('x')
            title(sprintf('%s -- SV-t scale', Ylabels{n}))
            view(88,33)
            wrapthisfigure(thisfig, sprintf('SVtscale-%s', ncode{n}),wrap)
        end
        
        
        
    case {'SVO', 'SVOar1'}
        
        SVOdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = sqrtht_all(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVOdraws(:,:,m) = thisSVdraw;
        end
        SVOmid   = mean(SVOdraws,3);
        SVOtails = prctile(SVOdraws, [5 95], 3);
        
        
        
        % compute SV sans outliers
        stochvol = sqrtht_all ./ SVOscale_all;
        SVdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = stochvol(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVdraws(:,:,m) = thisSVdraw;
        end
        SVmid   = mean(SVdraws,3);
        SVtails = prctile(SVdraws, [5 95], 3);
        
        
        for n = 1 : N
            thisfig = figure;
            
            hold on
            h1 = plot(ydates(1:thisT), SVOmid(n,:), 'r-', 'linewidth', 2);
            %             plot(ydates(1:thisT), squeeze(SVOtails(n,:,:)), 'r-', 'linewidth', 1)
            h2 = plot(ydates(1:thisT), SVmid(n,:), 'k-.', 'linewidth', 2);
            title(sprintf('%s', Ylabels{n}))
            xtickdates(ydates(1:thisT))
            legend([h1 h2], 'SV', 'SV w/o outlier', 'location', 'best')
            wrapthisfigure(thisfig, sprintf('SV%s', ncode{n}), wrap)
            
        end
        
        % plot density of SVO scale (for comparison against SVt)
        
        oscales = permute(SVOscale_all, [3 2 1]);
        
        ogrid = 0 : .1 : 20;
        outlierpdf = NaN(length(ogrid),T,N);
        parfor n = 1 : N
            for t = 1 : T
                outlierpdf(:,t,n) = ksdensity(oscales(:,t,n), ogrid);
            end
        end
        
        for n = 1 : N
            thisfig = figure;
            surf(ydates(p+1:thisT), ogrid, outlierpdf(:,:,n));
            shading interp
            datetick('x')
            title(sprintf('%s -- SVO scale', Ylabels{n}))
            view(88,33)
            wrapthisfigure(thisfig, sprintf('SVOfullscale-%s', ncode{n}),wrap)
        end
        
        
        % plot density of SVO scale
        if strcmpi(modeltype, 'SVO')
            SVOgrid   = 1 : SVOmaxscale;
        else
            SVOgrid   = [1, 2, 3, 4, 5, 6, 8, 10, 15, 20];
        end
        outlierpdf = NaN(SVOmaxscale,T,N);
        
        
        for t = 1 : T
            for n = 1 : N
                outlierpdf(:,t,n) = histcounts(squeeze(SVOscale_all(n,t,:)), 1:SVOmaxscale+1, 'Normalization', 'probability');
            end
        end
        
        for n = 1 : N
            thisfig = figure;
            subplot(1,2,1)
            surf(ydates(p+1:thisT), 2:SVOmaxscale, outlierpdf(2:end,:,n));
            shading interp
            if length(SVOgrid) > 10
                yticks(SVOgrid(2:2:end))
            else
                yticks(SVOgrid)
            end
            datetick('x')
            thisZlim = zlim;
            zlim([0 max(thisZlim(2), .01)])
            title('SV-O scale')
            
            priordraws = betadraw(SVOalpha(n), SVObeta(n), 1e4);
            [priorpdf, priorx] = ksdensity(priordraws, 'Support', 'positive');
            
            subplot(1,2,2)
            [posteriorpdf, posteriorx] = ksdensity(SVOprob_all(n,:), 'Support', 'positive');
            hold on
            h1 = plot(priorx, priorpdf, 'k-.');
            h2 = plot(posteriorx, posteriorpdf, 'r-');
            legend([h1 h2], 'prior', 'posterior')
            title('SV-O prob')
            wrapthisfigure(thisfig, sprintf('SVOscale-%s', ncode{n}),wrap)
            sgtitle(Ylabels{n})
            wrapthisfigure(thisfig, sprintf('SVOscale-%s-WITHTITLE', ncode{n}),wrap)
            
        end
        
        %% outlier prob, standalone
        
        
        for n = 1 : N
            priordraws = betadraw(SVOalpha(n), SVObeta(n), 1e4);
            [priorpdf, priorx] = ksdensity(priordraws, 'Support', 'positive');
            
            thisfig = figure;
            set(gca, 'fontsize', 20)
            [posteriorpdf, posteriorx] = ksdensity(SVOprob_all(n,:), 0 : .001 : .1, 'Support', 'positive');
            hold on
            h1 = plot(priorx, priorpdf, 'k-.', 'linewidth', 3);
            h2 = plot(posteriorx, posteriorpdf, 'r-', 'linewidth', 3);
            xlim([0 .1])
            wrapthisfigure(thisfig, sprintf('SVOprob-%s', ncode{n}),wrap, [], [], [], [], true);
            legend([h1 h2], 'prior', 'posterior')
            wrapthisfigure(thisfig, sprintf('SVOprob-%s-WITHLEGEND', ncode{n}),wrap, [], [], [], [], true);
            title(Ylabels{n})
            wrapthisfigure(thisfig, sprintf('SVOprob-%s-WITHTITLE', ncode{n}),wrap)
        end
        
        
    case {'SVobar'}
        
        SVOdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = sqrtht_all(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVOdraws(:,:,m) = thisSVdraw;
        end
        SVOmid   = mean(SVOdraws,3);
        SVOtails = prctile(SVOdraws, [5 95], 3);
        
        
        
        % compute SV sans outliers
        stochvol = sqrtht_all ./ SVobarscale_all;
        SVdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = stochvol(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVdraws(:,:,m) = thisSVdraw;
        end
        SVmid   = mean(SVdraws,3);
        SVtails = prctile(SVdraws, [5 95], 3);
        
        
        %% SV
        for n = 1 : N
            
            thisfig = figure;
            set(gca, 'fontsize', 20)
            
            hold on
            h1 = plot(ydates(1:thisT), SVOmid(n,:), 'r:', 'linewidth', 2);
            %             plot(ydates(1:thisT), squeeze(SVOtails(n,:,:)), 'r-', 'linewidth', 1)
            h2 = plot(ydates(1:thisT), SVmid(n,:), 'k-', 'linewidth', 2);
            ht = title(sprintf('%s', Ylabels{n}));
            xticks(datenum(1960:10:2030,1,1))
            xtickdates(ydates(1:thisT), 'keepticks')
            hl = legend([h1 h2], 'SV', 'SV w/o outlier', 'location', 'best');
            wrapthisfigure(thisfig, sprintf('SV%s-%s-WITHTITLELEGEND', ncode{n}, modeltype), wrap)
            delete(ht)
            wrapthisfigure(thisfig, sprintf('SV%s-%s-WITHLEGEND', ncode{n}, modeltype), wrap)
            delete(hl)
            wrapthisfigure(thisfig, sprintf('SV%s-%s', ncode{n}, modeltype), wrap)
            
        end
        
        %% outlier prob, standalone
        priordraws = betadraw(SVobaralpha, SVobarbeta, 1e4);
        [priorpdf, priorx] = ksdensity(priordraws, 'Support', 'positive');
        
        
        thisfig = figure;
        set(gca, 'fontsize', 20)
        [posteriorpdf, posteriorx] = ksdensity(SVobarprob_all(1,:), 0 : .001 : 1, 'Support', 'positive');
        hold on
        h1 = plot(priorx, priorpdf, 'k-.', 'linewidth', 3);
        h2 = plot(posteriorx, posteriorpdf, 'r-', 'linewidth', 3);
        xlim([0 1])
        wrapthisfigure(thisfig, 'SVobarprob',wrap, [], [], [], [], true);
        legend([h1 h2], 'prior', 'posterior')
        wrapthisfigure(thisfig, sprintf('SVobarprob-WITHLEGEND'),wrap, [], [], [], [], false);
        
    case {'SVOobar'}
        
        SVOdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = sqrtht_all(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVOdraws(:,:,m) = thisSVdraw;
        end
        SVOmid   = mean(SVOdraws,3);
        SVOtails = prctile(SVOdraws, [5 95], 3);
        
        
        
        % compute SV sans outliers
        stochvol = sqrtht_all ./ SVobarscale_all ./ SVOscale_all;
        stochvol2 = sqrtht_all ./ SVobarscale_all;
        SVdraws  = NaN(N,thisT,MCMCdraws);
        SVdraws2 = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = stochvol(:,:,m);
            thisSVdraw = NaN(N,thisT);
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            SVdraws(:,:,m) = thisSVdraw;
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = stochvol2(:,:,m);
            thisSVdraw = NaN(N,thisT);
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            SVdraws2(:,:,m) = thisSVdraw;
            
            
            
        end
        SVmid   = mean(SVdraws,3);
        SVtails = prctile(SVdraws, [5 95], 3);
        SVmid2   = mean(SVdraws2,3);
        SVtails2 = prctile(SVdraws2, [5 95], 3);
        
        %% plot SV
        for n = 1 : N
            thisfig = figure;
            subplot(2,1,1)
            hold on
            h1 = plot(ydates(1:thisT), SVOmid(n,:), 'r-', 'linewidth', 2);
            %             plot(ydates(1:thisT), squeeze(SVOtails(n,:,:)), 'r-', 'linewidth', 1)
            h2 = plot(ydates(1:thisT), SVmid(n,:), 'k-.', 'linewidth', 2);
            %             h3 = plot(ydates(1:thisT), SVmid2(n,:), 'b-.', 'linewidth', 2);
            %             title(sprintf('%s', Ylabels{n}))
            xtickdates(ydates(1:thisT))
            legend([h1 h2], 'SV', 'SV w/o outlier', 'location', 'best')
            
            subplot(2,1,2)
            hold on
            h1 = plot(ydates(1:thisT), SVOmid(n,:), 'r-', 'linewidth', 2);
            %             plot(ydates(1:thisT), squeeze(SVOtails(n,:,:)), 'r-', 'linewidth', 1)
            h3 = plot(ydates(1:thisT), SVmid2(n,:), 'b-.', 'linewidth', 2);
            xtickdates(ydates(1:thisT))
            legend([h1  h3], 'SV', 'SV w/o common outlier', 'location', 'best')
            
            sgtitle(sprintf('%s', Ylabels{n}))
            
            wrapthisfigure(thisfig, sprintf('SV%s-%', ncode{n}, modeltype), wrap)
            
        end
        
        %% outlier prob, standalone
        priordraws = betadraw(SVobaralpha, SVobarbeta, 1e4);
        [priorpdf, priorx] = ksdensity(priordraws, 'Support', 'positive');
        
        
        thisfig = figure;
        set(gca, 'fontsize', 20)
        [posteriorpdf, posteriorx] = ksdensity(SVobarprob_all(1,:), 0 : .001 : 1, 'Support', 'positive');
        hold on
        h1 = plot(priorx, priorpdf, 'k-.', 'linewidth', 3);
        h2 = plot(posteriorx, posteriorpdf, 'r-', 'linewidth', 3);
        xlim([0 1])
        wrapthisfigure(thisfig, 'SVobarprob',wrap, [], [], [], [], true);
        legend([h1 h2], 'prior', 'posterior')
        wrapthisfigure(thisfig, sprintf('SVobarprob-WITHLEGEND'),wrap, [], [], [], [], false);
        
        %% outlier prob, standalone
        priordraws = betadraw(SVOalpha, SVObeta, 1e4);
        [priorpdf, priorx] = ksdensity(priordraws, 'Support', 'positive');
        
        
        for n = 1 : N
            thisfig = figure;
            set(gca, 'fontsize', 20)
            [posteriorpdf, posteriorx] = ksdensity(SVOprob_all(n,:), 0 : .001 : .1, 'Support', 'positive');
            hold on
            h1 = plot(priorx, priorpdf, 'k-.', 'linewidth', 3);
            h2 = plot(posteriorx, posteriorpdf, 'r-', 'linewidth', 3);
            xlim([0 .1])
            wrapthisfigure(thisfig, sprintf('SVOprob-%s', ncode{n}),wrap, [], [], [], [], true);
            legend([h1 h2], 'prior', 'posterior')
            wrapthisfigure(thisfig, sprintf('SVOprob-%s-WITHLEGEND', ncode{n}),wrap, [], [], [], [], true);
            title(Ylabels{n})
            wrapthisfigure(thisfig, sprintf('SVOprob-%s-WITHTITLE', ncode{n}),wrap)
        end
        
    case {'SVOt', 'SVOtar1', 'SVOtdof'}
        
        SVdraws = NaN(N,thisT,MCMCdraws);
        parfor m=1:MCMCdraws
            invA = invA_all(:,:,m);
            
            thisSV = NaN(N,thisT);
            thisSV(:,p+1:end) = sqrtht_all(:,:,m);
            
            thisSVdraw = NaN(N,thisT);
            
            for t=p+1:thisT
                thisSVdraw(:,t) = sqrt(diag(invA*diag(thisSV(:,t).^2)*invA'));
            end
            
            SVdraws(:,:,m) = thisSVdraw;
        end
        SVmid   = mean(SVdraws,3);
        SVtails = prctile(SVdraws, [5 95], 3);
        
        
        for n = 1 : N
            thisfig = figure;
            
            hold on
            h1 = plot(ydates(1:thisT), SVmid(n,:), 'k-', 'linewidth', 2);
            plot(ydates(1:thisT), squeeze(SVtails(n,:,:)), 'k-', 'linewidth', 1)
            title(sprintf('%s', Ylabels{n}))
            xtickdates(ydates(1:thisT))
            wrapthisfigure(thisfig, sprintf('SV%s', ncode{n}), wrap)
            
        end
        
        %% outlier prob, standalone
        priordraws = betadraw(SVOalpha, SVObeta, 1e4);
        [priorpdf, priorx] = ksdensity(priordraws, 'Support', 'positive');
        
        
        for n = 1 : N
            thisfig = figure;
            set(gca, 'fontsize', 20)
            [posteriorpdf, posteriorx] = ksdensity(SVOprob_all(n,:), 0 : .001 : .1, 'Support', 'positive');
            hold on
            h1 = plot(priorx, priorpdf, 'k-.', 'linewidth', 3);
            h2 = plot(posteriorx, posteriorpdf, 'r-', 'linewidth', 3);
            xlim([0 .1])
            wrapthisfigure(thisfig, sprintf('SVOtprob-%s', ncode{n}),wrap, [], [], [], [], true);
            legend([h1 h2], 'prior', 'posterior')
            wrapthisfigure(thisfig, sprintf('SVOtprob-%s-WITHLEGEND', ncode{n}),wrap, [], [], [], [], true);
            title(Ylabels{n})
            wrapthisfigure(thisfig, sprintf('SVOtprob-%s-WITHTITLE', ncode{n}),wrap)
        end
end

 

%% plot t dof prob
if any(ismember({'SVt', 'SVOt'}, modeltype))
    
    close all
    tdofposterior = NaN(N, length(tdofGrid));
    for n = 1 : N
        % tposterior(n,:) = ksdensity(SVtdof_all(n,:), tdofGrid);
        tdofposterior(n,:) = histcounts(SVtdof_all(n,:), [tdofGrid tdofGrid(end) + 1], 'Normalization', 'probability');
    end
    
    priorprob = repmat(1 / length(tdofGrid), 1, length(tdofGrid));
    for n = 1 : N
        thisfig = figure;
        hold on
        hpost = bar(tdofGrid, tdofposterior(n,:), 1, 'facecolor', 'r', 'edgecolor', 'w');
        hprior = plot(tdofGrid, priorprob, 'k--', 'linewidth', 2);
        set(gca, 'fontsize', 20)
        
        xlim([tdofGrid(1) - 1, tdofGrid(end) + 1])
        xticks([3, 5 : 5 : tdofGrid(end)])
        wrapthisfigure(thisfig, sprintf('%s-dofposterior-%s', modeltype, ncode{n}),wrap)
        legend([hprior, hpost], 'prior', 'posterior')
        wrapthisfigure(thisfig, sprintf('%s-dofposterior-%s-WITHLEGEND', modeltype, ncode{n}),wrap)
        
        title(Ylabels{n})
        wrapthisfigure(thisfig, sprintf('%s-dofposterior-%s-WITHTITLELEGEND', modeltype, ncode{n}),wrap)
    end
    
end

%% plot SV-AR1 rho
if contains(lower(modeltype), 'ar1')
    
    priordraws = .8 + .2 * randn(MCMCdraws,1);
    for n = 1 : N
        thisfig = figure;
        hold on
        
        
        plotpriorposteriordraws(SVrho_all(n,:)', priordraws, 0:.01:1)
        
        title(Ylabels{n})
        wrapthisfigure(thisfig, sprintf('%s-SVar1-%s-WITHTITLELEGEND', modeltype, ncode{n}),wrap)
    end
    
end


%% wrap up
dockAllFigures
finishwrap





