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
doRobustPrior = false;

% SED-PARAMETERS-HERE

doStoreXL           = false; %#ok<*NASGU>
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)
Compute_diagnostics = false;              % compute Inefficiency Factors and Potential

dummyType           = 'COVIDEACH'; 
dummyPrecision      = 0;

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
    doStoreXL = false; % set to true to store draws (in separate matfile)
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

Ndummy = sum(ydates >= datenum(2020, 3,1));
Kbvar = N * p + 1; % number of VAR regressors per equation
K = Kbvar + Ndummy; % number of all regressors

if isempty(ELBbound)
    labelSampling = 'NOshadowrate';
else
    labelSampling = 'censoredYields';
end

if doTightPrior
    labelSampling = strcat(labelSampling, '-tightBVARshrinkage');
end
if doRobustPrior
    labelSampling = strcat(labelSampling, '-covidrobustPrior');
end

labelSampling = strcat(labelSampling, '-SV');
labelSampling = strcat(labelSampling, sprintf('-%s-precision%d', dummyType, floor(dummyPrecision)));

% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

% define oos jump offs
Tjumpoffs = find(ydates >= datenum(2020,3,1));

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

[PAImedian, PAImean, PAIstdev] = deal(NaN(K, N, Njumpoffs));
PAIquantiles                   = NaN(K, N, Nquantiles, Njumpoffs);
    

drawsSVmid        = NaN(N, Tdata, Njumpoffs);
drawsSVtails      = NaN(N, Tdata, Nquantiles, Njumpoffs);


%% allocate memory for MCMC output (ex forecast)
if ~doLoMem
    drawsPAI       = NaN(K, N, MCMCdraws, Njumpoffs);
    drawsPHI       = NaN(N*(N-1)/2+N, MCMCdraws, Njumpoffs);
    drawsINVA      = NaN(N, N, MCMCdraws, Njumpoffs);
    drawsSQRTHT    = NaN(N, Tdata, MCMCdraws, Njumpoffs);
end

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
    [PAI_all, PHI_all, invA_all, sqrtht_all, ...
        ydraws, yhat, yhatRB] ...
        = mcmcVARSVdummy(thisT, MCMCdraws, p, np, data, ydates, ...
        minnesotaPriorMean, doTightPrior, doRobustPrior, dummyPrecision, ...
        ndxYIELDS, ELBbound, ...
        check_stationarity, ...
        fcstNdraws, fcstNhorizons, rndStreams{TID}); %#ok<PFBNS>
    
    
    %% Convergence diagnostics
    if Compute_diagnostics
        % display('computing convergence diagnostics..')
        thisK = size(PAI_all,1);
        Diagnostics(sqrtht_all,invA_all,PAI_all,PHI_all,N,thisK,MCMCdraws);
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
    theseMaxVARroots = NaN(MCMCdraws, 1); % placed before doLoMem to avoid parfor warning
    % setup companion form matrix
    comp                        = zeros(N * p);
    comp(N + 1 : end,1:N*(p-1)) = eye(N*(p-1));
    for m = 1 : MCMCdraws
        thisPAI = PAI_all(:,:,m);
        comp(1:N,:) = thisPAI(2:Kbvar,:)';
        % compute maxLambda
        theseMaxVARroots(m) = max(abs(eig(comp)));
    end
    
    %% collect PAI moments
	thisK     = size(PAI_all,1);
    
    thisDummy = NaN(K,N);
    thisDummy(1:thisK,:)  = median(PAI_all,3); 
    PAImedian(:,:,ndxT)   = thisDummy;
    
    thisDummy = NaN(K,N);
    thisDummy(1:thisK,:)      = mean(PAI_all,3); 
    PAImean(:,:,ndxT)         = thisDummy;
    
    thisDummy = NaN(K,N);
    thisDummy(1:thisK,:)      = std(PAI_all,1,3); 
    PAIstdev(:,:,ndxT)        = thisDummy;
    
    thisDummy                 = NaN(K,N,Nquantiles);
    thisDummy(1:thisK,:,:)    = prctile(PAI_all,setQuantiles,3);
    PAIquantiles(:,:,:,ndxT)  = thisDummy;
    
    %% compute SV
    SVdraws = NaN(N,Tdata,MCMCdraws);
    stochvol                      = NaN(N, Tdata, MCMCdraws);
    stochvol(:,p+1:thisT, :)      = sqrtht_all;
    for m=1:MCMCdraws
        invA     = invA_all(:,:,m);
        for t=1:Tdata
            SVdraws(:,t,m) = sqrt(diag(invA*diag(stochvol(:,t,m).^2)*invA'));
        end
    end
    
    drawsSVmid(:,:,ndxT)     = median(SVdraws, 3);
    drawsSVtails(:,:,:,ndxT) = prctile(SVdraws, setQuantiles, 3);
    
    
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
    
    if ~doLoMem
        thisK = Kbvar + sum(ydates(1:thisT) >= datenum(2020, 3,1));
        thisPAI                     = NaN(K, N, MCMCdraws);
        thisPAI(1:thisK,:,:)        = PAI_all;
        drawsPAI(:,:,:,ndxT)        = thisPAI;
        drawsPHI(:,:,ndxT)    = PHI_all;
        drawsINVA(:,:,:,ndxT) = invA_all;
        % prepare dummy to make parfor work
        dummy                      = NaN(N, Tdata, MCMCdraws);
        dummy(:,p+1:thisT, :)      = sqrtht_all;
        drawsSQRTHT(:, :, :, ndxT) = dummy;
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



%% plot COVIDSV
ndx = ydates >= datenum(2020,1,1);
thesedates = ydates(ndx);
theseJumpoffs = Njumpoffs - 3 : Njumpoffs;

for n = 1 : N
    thisfig = figure;
    ax = gca;
    set(ax, 'fontsize', 16)
    hold on
    
    lineTypes = {'-d', '--o', 'd-.', ':o'};
    iter = 0;
    for ndxT = theseJumpoffs
        iter = iter + 1;
        plot(thesedates, drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
    end
    
    
    ylim([0 max(ylim)])
    xtickdates(thesedates)
    legend(datestr(ydates(Tjumpoffs(theseJumpoffs)), 'yyyy:mm'), 'location', 'northwest', 'box', 'off')
    
    
    wrapthisfigure(thisfig, sprintf('covidSV%s', ncode{n}), wrap)
    
end


%% plot companion maxLambda
meanMaxVARroot   = mean(drawsMaxVARroot,1);
medMaxVARroot   = median(drawsMaxVARroot,1);
tailsMaxVARroot = prctile(drawsMaxVARroot, [5 95], 1);

this = figure;
hold on
plot(ydates(Tjumpoffs), meanMaxVARroot, 'k-', 'linewidth', 2)
plot(ydates(Tjumpoffs), medMaxVARroot, 'r--', 'linewidth', 2)
plot(ydates(Tjumpoffs), tailsMaxVARroot', 'k-', 'linewidth', 1)
xtickdates(ydates(Tjumpoffs))
wrapthisfigure(this, 'maxVARroot', wrap)

%% plot all coviddummies

if ~doLoMem 
    
    drawsDummy = drawsPAI(end-Ndummy+1:end,:,:,:);
    
    bardates = ydates(Tjumpoffs);
    bardates = cat(1, bardates(1) - 30, bardates, bardates(end) + 30);
    for nn = 1 : Ndummy
        for n = 1 : N
            
            thisMid = squeeze(median(drawsDummy(nn,n,:,:), 3));
            thisCI  = squeeze(prctile(drawsDummy(nn,n,:,:), [5 95], 3));
            
            thisfig = figure;
            
            hold on
            bar(ydates(Tjumpoffs), thisMid)
            % plot(ydates(Tjumpoffs), thisMid, 'r-d', 'linewidth', 2)
            plot(ydates(Tjumpoffs), thisCI, 'r--d')
            title(sprintf('%s', Ylabels{n}))
            plotOrigin
            
            % xlim(ydates(Tjumpoffs([1 end])))
            xlim(bardates([1 end]))
            datetick('x', 'yyyy:mm', 'keeplimits')
            
            title(sprintf('Dummy %d for %s', nn, Ylabels{n}))
            
            wrapthisfigure(thisfig, sprintf('covidDummy%d-%s', nn, ncode{n}), wrap)
            
        end
    end
end

%% plot coviddummies for last vintage

if ~doLoMem && ~strcmpi(dummyType, 'CONSTsansCOVID')
    
    drawsDummy = drawsPAI(end-Ndummy+1:end,:,:,end); % final vintage only
    
    bardates = ydates(Tjumpoffs);
    bardates = cat(1, bardates(1) - 30, bardates, bardates(end) + 30);
    for n = 1 : N
        
        thisMid = squeeze(median(drawsDummy(:,n,:), 3));
        thisCI  = squeeze(prctile(drawsDummy(:,n,:), [5 95], 3));
        
        thisfig = figure;
        
        hold on
        bar(ydates(Tjumpoffs), thisMid)
        % plot(ydates(Tjumpoffs), thisMid, 'r-d', 'linewidth', 2)
        h1 = plot(ydates(Tjumpoffs), thisCI, 'b--d');
        % title(sprintf('%s', Ylabels{n}))
        plotOrigin
        
        xticks(ydates(Tjumpoffs))
        xlim(bardates([1 end]))
        datetick('x', 'yyyy:mm', 'keeplimits', 'keepticks')
        
        % title(sprintf('Dummy %d for %s', nn, Ylabels{n}))
        
        YLIM = ylim;
        wrapthisfigure(thisfig, sprintf('FINALcovidDummies-%s', ncode{n}), wrap)
        delete(h1) 
        ylim(YLIM);
        wrapthisfigure(thisfig, sprintf('FINALcovidDummies-%s-WITHOUTCI', ncode{n}), wrap)
        
    end
end

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
    'PAI*', ...
    'meanMaxVARroot', 'medMaxVARroot', 'tailsMaxVARroot', ...
    'datalabel', 'labelSampling', ...
    'doQuarterly', ...
    'setQuantiles', ...
    'MCMCdraws'};

if doStoreXL
    matfilename = sprintf('%s-%s-p%d-draws', datalabel, labelSampling, p);
    save(matfilename, varlist{:}, 'draws*', '-v7.3');
end
clear drawsSVtdof drawsSVtscale* *_all
clear drawsSVtails drawsMaxVARroot 
clear fcstYmvlogscore*Draws
save(matfilename, varlist{:}, '-v7.3');

%% wrap up
dockAllFigures
finishwrap

