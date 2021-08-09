%% ************************************************************************************
% Triangular algorithm, Carriero Clark and Marcellino (2015), Large Vector Autoregressions
% with stochastic volatility and flexible priors. The code takes approx. 5
% mins on an average desktop.
% ************************************************************************************
% Model is:
%
%     Y(t) = Pi(L)Y(t-1) + v(t); Y(t) is Nx1;
%     v(t) = inv(A)*(LAMBDA(t)^0.5)*e(t); e(t) ~ N(0,I);
%                   _                                         _
%                  |    1          0       0       ...      0  |
%                  |  A(2,1)       1       0       ...      0  |
%         A =      |  A(3,1)     A(3,2)    1       ...      0  |
%                  |   ...        ...     ...      ...     ... |
%                  |_ A(N,1)      ...     ...   A(N,N-1)    1 _|
%
%    Lambda(t)^0.5 = diag[sqrt_h(1,t)  , .... , sqrt_h(N,t)];
%
%    ht=exp(Vol_states)
%    sqrtht=sqrt(ht)=sqrt(exp(Vol_states))=exp(Vol_states/2)
%    Vol_states=2*ln(sqrtht)
%
%    Vol_states(i,t)   = Vol_states(i,t-1) + eta(i,t),
%    eta(t) ~ N(0,PHI); with PHI full
% ------------------------------------------------------------------------------------

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



%% ctd.
ndxSHADOWRATE = find(strcmp(ncode, 'FEDFUNDS'));

% define index of yields that need to obey ELB (at least out of sample)
if doCensorYields
    ndxOTHERYIELDS = find(ismember(ncode, {'GS1', 'GS5', 'GS10'}));
else
    ndxOTHERYIELDS = [];
end

ndxYIELDS = union(ndxSHADOWRATE, ndxOTHERYIELDS);
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

minnesotaPriorMean = NaN(N,1);

if contains(datalabel, 'levels')
    minnesotaPriorMean(:) = 1;
else
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



%% setup estimation windo

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
        
    case 'SVtdof'
        SVtdof  = 5;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            SVtscalelog2_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws] ...
            = mcmcVARSVtdof(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVtdof, ...
            ndxYIELDS, ELBbound, ...
            check_stationarity, ...
            fcstNdraws, fcstNhorizons, rndStreams{TID});
        
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
        
    case 'SVOtdof'
        
        
        SVOpriorobs = 10 * np;
        SVOalpha    = 1 / (10 * np) * SVOpriorobs; % 10 years of data with 1 outlier every 10 years
        SVObeta     = SVOpriorobs - SVOalpha;
        
        SVtdof = 9;
        
        [PAI_all, PHI_all, invA_all, sqrtht_all, ...
            SVOprob_all, SVOscale_all, ...
            SVtscalelog2_all, ...
            ydraws, yhat, yhatRB, fcstSVdraws, fcstSVoutliers] ...
            = mcmcVARSVOtdof(thisT, MCMCdraws, p, np, data, ydates, ...
            minnesotaPriorMean, doTightPrior, doRobustPrior, ...
            SVOalpha, SVObeta, SVOmaxscale, ...
            SVtdof, ...
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
            yrealized, ...
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

% offsetlogscoredraw = max(logscoredraws);
% 
% figure
% subplot(2,1,1)
% histogram(exp(logscoredraws - offsetlogscoredraw))
% title('exp(log score draws)')
% subplot(2,1,2)
% histogram(logscoredraws - offsetlogscoredraw)
% title('log score draws')
% 
% yLogscore       = log(mean(exp(logscoredraws - offsetlogscoredraw))) + offsetlogscoredraw;
% fprintf('Logscore: %6.4f\n', yLogscore)
% 
% %MV logscore
% thesedraws  = squeeze(ydraws(:,1,:));
% MU          = mean(thesedraws, 2);
% Sigma       = cov(thesedraws', 1); % normalize variance by N
% sqrtSigma   = chol(Sigma)';
% logdetSigma = 2 * sum(log(diag(sqrtSigma)));
% dev         = sqrtSigma \ (yrealized(:,1) - MU);
% SSR         = dev' * dev;
% yMVlogscore = -.5 * (N * log(2 * pi) + logdetSigma + SSR);
% 
% fprintf('Mean-variance Logscore: %6.4f\n', yMVlogscore)


%% Diagnostics

if Compute_diagnostics
    if contains(lower(modeltype), 'const')
        DiagnosticsCONST(PAI_all,SIGMA_all,N,K,MCMCdraws);
    else
        Diagnostics(sqrtht_all,invA_all,PAI_all,PHI_all,N,K,MCMCdraws);
    end
end

%% compute out-of-sample forecasts
% ydraws(cumcode,:,:)  = cumsum(ydraws(cumcode,:,:),2);
% yhat(cumcode,:)      = cumsum(yhat(cumcode,:),2);
% yhatRB(cumcode,:)    = cumsum(yhatRB(cumcode,:),2);



        
%% store
matfilename = sprintf('%s-%s-p%d-%s-alldraws', datalabel, datestr(ydates(thisT), 'yyyymm'), p, modeltype);
varlist = {'ydates', 'p', 'N', ...
    'data', 'np',  ...
    'ncode', 'tcode', 'cumcode', ...
    'fcst*', 'fcstNhorizons', ...
    'ndxYIELDS', 'ELBbound', ...
    '*draws*', ...
    '*_all', ...
    'datalabel', ...
    'setQuantiles', ...
    'MCMCdraws'};

save(matfilename, varlist{:}, '-v7.3');

%% wrap up
dockAllFigures
finishwrap





