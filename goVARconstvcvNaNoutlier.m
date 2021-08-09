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
rndStreams  = initRandStreams(Nstreams, [], 0);


%% set parameters for VAR and MCMC

datalabel           = 'fredMD16-2021-04'; 
doQuarterly         = false;
MCMCdraws           = 1e3;               % Final number of MCMC draws after burn in
fcstNdraws          = 100 * MCMCdraws;    % draws sampled from predictive density
doCensorYields      = true;
Ofactor             = 10;
do2020 = false;

% SED-PARAMETERS-HERE

doStoreXL           = false; %#ok<*NASGU>
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)
Compute_diagnostics = false;              % compute Inefficiency Factors and Potential

doLoMem             = true; % do not store memory intensive stuff, just do oos forecasts

doPlotData          = false;

samStart            = []; % datenum(1988,12,1);                 % truncate start of sample if desired (leave empty if otherwise)
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

if doLoMem
    doStoreXL = false;
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
    ndxYIELDS = [];
end

Tdata = length(ydates);

Ylabels = fredMDprettylabel(ncode);

%% process settings
N = size(data,2);

doTightPrior = false; % N >= 10;

Kbvar = N * p + 1; % number of regressors per equation
K = Kbvar;

if isempty(ELBbound)
    labelSampling = 'NOshadowrate';
else
    labelSampling = 'censoredYields';
end

if doTightPrior
    labelSampling = strcat(labelSampling, '-tightBVARshrinkage');
end

labelSampling = strcat(labelSampling, sprintf('-CONSTnanO%d', Ofactor));

% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

% define oos jump offs
if isempty(samStart)
    Tjumpoffs = find(ydates >= datenum(1985,1,1));
else
    Tjumpoffs = find(ydates >= datenum(2000,1,1));
end

if do2020
    Tjumpoffs = find(ydates >= datenum(2020,1,1));
    labelSampling = strcat(labelSampling, '-2020');
end

Njumpoffs = length(Tjumpoffs);

% ELB settings

% other settings

setQuantiles   = [.5, 2.5, 5, normcdf(-1) * 100, 25 , 75,  (1 - normcdf(-1)) * 100, 95, 97.5, 99.5];
Nquantiles    = length(setQuantiles);
ndxCI         = ismember(setQuantiles, [5, normcdf(-1) * 100, 100 - normcdf(-1) * 100, 95]);
%% mean for Minnesota prior: zero (diff) or RW (level)

minnesotaPriorMean = NaN(N,1);

for n = 1 : N
    switch ncode{n}
        case {'CUMFNS', 'UNRATE', ...
                'WPSFD49207',  'PPICMM', 'PCEPI', ...
                'FEDFUNDS', 'HOUST', 'GS5', 'GS10', 'BAAFFM', 'WUXIASHADOWRATE'}
            minnesotaPriorMean(n) = 1;
        otherwise
            minnesotaPriorMean(n) = 0;
    end
end


%% allocate memory for out-of-sample forecasts

fcstNhorizons     = 24;  % number of steps forecasted (1:fcstNhorizon)

% fcstYdraws        = NaN(N,fcstNhorizons,fcstNdraws,Njumpoffs);
fcstYrealized     = NaN(N,fcstNhorizons,Njumpoffs);
fcstYhat          = NaN(N,fcstNhorizons,Njumpoffs); % predictive mean
fcstYhatRB        = NaN(N,fcstNhorizons,Njumpoffs); % predictive mean (linear RB)
fcstYmedian       = NaN(N,fcstNhorizons,Njumpoffs); % predictive median
fcstYhaterror     = NaN(N,fcstNhorizons,Njumpoffs);
fcstYhatRBerror   = NaN(N,fcstNhorizons,Njumpoffs);
fcstYmederror     = NaN(N,fcstNhorizons,Njumpoffs);
fcstLogscore      = NaN(N,fcstNhorizons,Njumpoffs);
fcstCRPS          = NaN(N,fcstNhorizons,Njumpoffs);
fcstYquantiles    = NaN(N,fcstNhorizons,Nquantiles, Njumpoffs);

fcstYmvlogscore   = NaN(1,Njumpoffs); % one-step ahead only



%% allocate memory for MCMC output (ex forecast)
if ~doLoMem
    drawsPAI       = NaN(K, N, MCMCdraws, Njumpoffs);
    drawsSIGMA  = NaN(N, N, MCMCdraws, Njumpoffs);
    
    drawsY      = NaN(Tdata,N,MCMCdraws, Njumpoffs);
end

drawsYmid        = NaN(Tdata, N, Njumpoffs);
drawsYtails      = NaN(Tdata, N, Nquantiles, Njumpoffs);

drawsMaxVARroot = NaN(MCMCdraws, Njumpoffs);

%% start latexwrapper to collect results
titlename=sprintf('%s-%s-p%d', datalabel, labelSampling, p);
if ~isempty(samStart)
    titlename = strcat(titlename, '-', datestr(samStart, 'yyyymmm'));
end
initwrap
% wrap = [];

%% plot input data
if doPlotData
    for n = 1 : N
        this = figure;
        plot(ydates, data(:,n))
        xtickdates(ydates)
        wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
    end
end

%% loop over QRT estimates

% progressbar(0)
parfor ndxT = 1 : Njumpoffs % parfor
    
    TID   = parid;
    thisT = Tjumpoffs(ndxT);
    T     = thisT - p;
    
    fprintf('loop %d, thisT %d, with TID %d\n', ndxT, thisT, TID)
    
    %% MCMC sampler
    
    OK = false;
    while ~OK
        try
            
            [PAI_all, SIGMA_all, ~, Ydatadraws_all, Ynan,...
                ydraws, yhat, yhatRB] ...
                = mcmcVARconstvcvNaNoutlier(thisT, MCMCdraws, Ofactor, p, np, data, ydates, ...
                minnesotaPriorMean, doTightPrior, ...
                ndxYIELDS, ELBbound, ...
                check_stationarity, ...
                fcstNdraws, fcstNhorizons, rndStreams{TID}); %#ok<PFBNS>
            OK = true;
            
        catch ME
            
            warning('trying again with %d due to \n << %s >>',ndxT, ME.message)
            
        end
    end
    
    %% Convergence diagnostics
    if Compute_diagnostics
        % display('computing convergence diagnostics..')
        DiagnosticsCONST(PAI_all,SIGMA_all,N,K,MCMCdraws);
    end
    
    %% compute out-of-sample forecasts
    
    % a word on parfor strategy:
    % to make matlab better see the intended use of sliced variabes, use
    % local temp variables and then copy those into the slices at end of
    % loop
    
    % collect realized values
    yrealized = NaN(N, fcstNhorizons);
    for h = 1 : fcstNhorizons
        if thisT + h <= Tdata
            yrealized(:,h) = data(thisT+h,:)';
        end
    end
    
    
    
    % cumulate realizations and predictions if necessary
    yrealized(cumcode,:) = cumsum(yrealized(cumcode,:),2);
    ydraws(cumcode,:,:)  = cumsum(ydraws(cumcode,:,:),2);
    yhat(cumcode,:)      = cumsum(yhat(cumcode,:),2);
    yhatRB(cumcode,:)    = cumsum(yhatRB(cumcode,:),2);
    
    % compute median
    ymed = median(ydraws,3);
    
    % mv logscore one step ahead (ignoring ELB)
    thesedraws  = squeeze(ydraws(:,1,:));
    MU          = mean(thesedraws, 2);
    Sigma       = cov(thesedraws', 1); % normalize variance by N
    sqrtSigma   = chol(Sigma)';
    logdetSigma = 2 * sum(log(diag(sqrtSigma)));
    dev         = sqrtSigma \ (yrealized(:,1) - MU);
    SSR         = dev' * dev;
    yMVlogscore = -.5 * (N * log(2 * pi) + logdetSigma + SSR);
    
    % logscore
    yLogscore = NaN(N,fcstNhorizons);
    if isempty(ELBbound)
        ndxBoundedSupport = false(N,1);
    else
        ndxBoundedSupport = ismember(1:N, ndxYIELDS)';
    end
    ndxUnboundedSupport = ~ndxBoundedSupport;
    % a) compute logscore for variables with unbounded support via gaussian approximation
    thesedraws = ydraws(ndxUnboundedSupport,:,:);
    mu     = mean(thesedraws, 3);
    sigma2 = var(thesedraws, 1, 3); % normalize variance by N (rather than N-1)
    yLogscore(ndxUnboundedSupport,:) = -.5 * (log(2 * pi) + log(sigma2) + ((yrealized(ndxUnboundedSupport,:) - mu).^2 ./ sigma2));
    % b) compute logscore for variables with bounded support vis kernel density of truncated normal
    if ~isempty(ndxBoundedSupport)
        ndx = find(ndxBoundedSupport);
        for h = 1 : fcstNhorizons
            for n = 1 : length(ndx) % loop over elements of Y
                if ~isnan(yrealized(ndx(n),h))
                    % note: ksdensity around NaN returns 0, log(ksdensity) is then also NaN
                    thesedraws = squeeze(ydraws(ndx(n),h,:));
                    adjust4ELB = ELBbound - eps; % eps to ensure *positive* draws
                    yLogscore(ndx(n),h) = log(ksdensity(thesedraws - adjust4ELB, yrealized(ndx(n),h) - adjust4ELB, 'support', 'positive'));
                end
            end
        end
    end
    
    % CRPS
    yCRPS = NaN(N,fcstNhorizons);
    for h = 1 : fcstNhorizons
        for n = 1 : N % loop over elements of Y
            yCRPS(n,h) = crpsDraws(yrealized(n,h), ydraws(n,h,:));
        end
    end
    
    
 
    
    %% compute maxVARroot
    theseMaxVARroots = NaN(MCMCdraws, 1); % placed before doLoMem to avoif parfor warning
    % setup companion form matrix
    comp                        = zeros(N * p);
    comp(N + 1 : end,1:N*(p-1)) = eye(N*(p-1));
    for m = 1 : MCMCdraws
        thisPAI = PAI_all(:,:,m);
        comp(1:N,:) = thisPAI(2:Kbvar,:)';
        % compute maxLambda
        theseMaxVARroots(m) = max(abs(eig(comp)));
    end
    
    
    %% copy results into sliced variables
    
    fcstYmvlogscore(:,ndxT) = yMVlogscore;
    
    
    % forecast
    fcstYhat(:,:,ndxT)      = yhat; % mean(ydraws,3) or analytically
    fcstYhatRB(:,:,ndxT)    = yhatRB;
    fcstYmedian(:,:,ndxT)   = ymed; % median(ydraws,3);
    fcstYhaterror(:,:,ndxT) = yrealized - yhat;
    fcstYmederror(:,:,ndxT) = yrealized - ymed;
    % fcstYdraws(:,:,:,ndxT)  = ydraws;
    fcstYrealized(:,:,ndxT) = yrealized;
    fcstCRPS(:,:,ndxT)      = yCRPS;
    fcstLogscore(:,:,ndxT)  = yLogscore;
    fcstYquantiles(:,:,:,ndxT)  = prctile(ydraws, setQuantiles, 3);
    
    % copy mcmc output
    drawsMaxVARroot(:,ndxT)     = theseMaxVARroots;
    
    dummy                      = NaN(Tdata, N, MCMCdraws);
    dummy(p+1:thisT,:,:)       = Ydatadraws_all;
    drawsYmid(:,:,ndxT)     = median(dummy, 3);
    drawsYtails(:,:,:,ndxT) = prctile(dummy, setQuantiles, 3);
    
    if ~doLoMem
        drawsY(:, :, :, ndxT)      = dummy;
        drawsPAI(:,:,:,ndxT)  = PAI_all;
        % drawsPHI(:,:,ndxT)     = PHI_all;
        % drawsINVA(:,:,:,ndxT)  = invA_all;
        drawsSIGMA(:,:,:,ndxT) = SIGMA_all;
        % prepare dummy to make parfor work
    end
end

%% plot evolution of predictive densities
theseHorizons = [1 8 16 24];
for n = 1 : N
    
    thisfig = figure;
    
    for ii = 1 : length(theseHorizons)
        h = theseHorizons(ii);
        
        subplot(2,2,ii)
        fcstMid    = squeeze(fcstYhat(n,h,:));
        theseTails = squeeze(fcstYquantiles(n,h,ndxCI,:))';
        
        
        hold on
        plotCI(fcstMid, theseTails, ydates(Tjumpoffs));
        if any(n == ndxYIELDS)
            plot(ydates(Tjumpoffs),squeeze(fcstYmedian(n,h,:)), 'r-.', 'linewidth', 3)
        end
        plot(ydates(Tjumpoffs),squeeze(fcstYhatRB(n,h,:)), 'b--', 'linewidth', 3)
        
        title(sprintf('h=%d', h))
        sgtitle(sprintf('%s', Ylabels{n}))
        xtickdates(ydates(Tjumpoffs))
        
    end
    wrapthisfigure(thisfig, sprintf('predictiveDensity-%s', ncode{n}), wrap)
end



%% plot evolution of forecast-error variances

if ~doLoMem
    % collect FEV
    drawsFEV = NaN(MCMCdraws, Njumpoffs, N);
    for n = 1 : N;
        drawsFEV(:,:,n) = squeeze(drawsSIGMA(n,n,:,:));
    end
    FEVmid   = squeeze(mean(drawsFEV,1));
    FEVtails = permute(prctile(drawsFEV,[5 95],1), [2 1 3]);
    
    for n = 1 : N
        thisfig = figure;
        
        hold on
        plot(ydates(Tjumpoffs), FEVmid(:,n), 'k-', 'linewidth', 2)
        plot(ydates(Tjumpoffs), FEVtails(:,:,n), 'k-', 'linewidth', 1)
        title(sprintf('%s', Ylabels{n}))
        xtickdates(ydates(Tjumpoffs))
        wrapthisfigure(thisfig, sprintf('fev%s', ncode{n}), wrap)
        
    end
end

%% plot companion maxLambda
meanMaxVARroot = mean(drawsMaxVARroot,1);
medMaxVARroot = median(drawsMaxVARroot,1);
tailsMaxVARroot = prctile(drawsMaxVARroot, [5 95], 1);

this = figure;
hold on
plot(ydates(Tjumpoffs), meanMaxVARroot, 'k-', 'linewidth', 2)
plot(ydates(Tjumpoffs), medMaxVARroot, 'r--', 'linewidth', 2)
plot(ydates(Tjumpoffs), tailsMaxVARroot', 'k-', 'linewidth', 1)
xtickdates(ydates(Tjumpoffs))
wrapthisfigure(this, 'maxVARroot', wrap)

%% store qrt summary
matfilename = sprintf('%s-%s-p%d', datalabel, labelSampling, p);
if ~isempty(samStart)
    matfilename = strcat(matfilename, '-', datestr(samStart, 'yyyymmm'));
end
varlist = {'ydates', 'p', 'Tjumpoffs', 'N', ...
    'data', 'np',  ...
    'ncode', 'tcode', 'cumcode', ...
    'fcst*', 'fcstNhorizons', ...
    'ndxYIELDS', 'ELBbound', ...
    'drawsY*', ...
    'meanMaxVARroot', 'medMaxVARroot', 'tailsMaxVARroot', ...
    'datalabel', 'labelSampling', ...
    'doQuarterly', ...
    'setQuantiles', ...
    'MCMCdraws'};

if doStoreXL
    matfilename = sprintf('%s-%s-p%d-draws', datalabel, labelSampling, p);
    save(matfilename, varlist{:}, 'draws*', '-v7.3');
end

clear *_all
clear drawsMaxVARroot 
clear fcstYmvlogscore*Draws
save(matfilename, varlist{:}, '-v7.3');

%% wrap up
dockAllFigures
finishwrap

