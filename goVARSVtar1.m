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
do2020              = false;               
do1975              = false;               
doRobustPrior       = false;
tdofGrid            = 3 : 40; % grid for uniform prior

% SED-PARAMETERS-HERE

doStoreXL           = false; %#ok<*NASGU>
check_stationarity  = 0;                  % Truncate nonstationary draws? (1=yes)
Compute_diagnostics = false;              % compute Inefficiency Factors and Potential

doLoMem             = ~do2020; % do not store memory intensive stuff, just do oos forecasts

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
if doRobustPrior
    labelSampling = strcat(labelSampling, '-covidrobustPrior');
end

labelSampling = strcat(labelSampling, '-SVtar1');

% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

% define oos jump offs
if isempty(samStart)
    Tjumpoffs = find(ydates >= datenum(1975,1,1));
else
    Tjumpoffs = find(ydates >= datenum(2000,1,1));
end

if do2020
    Tjumpoffs = find(ydates >= datenum(2020,1,1));
    labelSampling = strcat(labelSampling, '-2020');
end

if do1975
    Tjumpoffs = find(ydates >= datenum(1975,1,1) & ydates < datenum(1985,1,1));
    labelSampling = strcat(labelSampling, '-1975');
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

fcstYmvlogscoreDraws = NaN(fcstNdraws,Njumpoffs); % one-step ahead only
fcstYmvlogscore   = NaN(1,Njumpoffs); % one-step ahead only
fcstYmvlogscore2  = NaN(1,Njumpoffs); % one-step ahead only

[PAImedian, PAImean, PAIstdev] = deal(NaN(K, N, Njumpoffs));
PAIquantiles                   = NaN(K, N, Nquantiles, Njumpoffs);
    

drawsSVmid        = NaN(N, Tdata, Njumpoffs);
drawsSVtails      = NaN(N, Tdata, Nquantiles, Njumpoffs);

drawsLambdaSVmid    = NaN(N, Tdata, Njumpoffs);
drawsLambdaSVtails  = NaN(N, Tdata, Nquantiles, Njumpoffs);

%% allocate memory for MCMC output (ex forecast)
if ~doLoMem
    drawsPAI       = NaN(K, N, MCMCdraws, Njumpoffs);
    drawsPHI       = NaN(N*(N-1)/2+N, MCMCdraws, Njumpoffs);
    drawsINVA      = NaN(N, N, MCMCdraws, Njumpoffs);
    drawsSQRTHT    = NaN(N, Tdata, MCMCdraws, Njumpoffs);
    
    drawsSVtscalelog2   = NaN(N, Tdata, MCMCdraws, Njumpoffs);
    
end

drawsSVtdof    = NaN(N, MCMCdraws, Njumpoffs);
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
    
    thisdata = data; % parfor
    
    %% prepare realized values
    yrealized = NaN(N, fcstNhorizons);
    for h = 1 : fcstNhorizons
        if thisT + h <= Tdata
            yrealized(:,h) = thisdata(thisT+h,:)';
        end
    end
    yrealized(cumcode,:) = cumsum(yrealized(cumcode,:),2);
    
    %% MCMC sampler
    [PAI_all, PHI_all, invA_all, SVrho_all, sqrtht_all, ...
        SVtscalelog2_all, SVtdof_all, ...
        ydraws, yhat, yhatRB, fcstSVdraws, logscoredraws] ...
        = mcmcVARSVtar1(thisT, MCMCdraws, p, np, thisdata, ydates, ...
        minnesotaPriorMean, doTightPrior, doRobustPrior, ...
        tdofGrid, ...
        ndxYIELDS, ELBbound, ...
        check_stationarity, ...
        yrealized,...
        fcstNdraws, fcstNhorizons, rndStreams{TID}); %#ok<PFBNS>
    
    
    %% Convergence diagnostics
    if Compute_diagnostics
        % display('computing convergence diagnostics..')
        Diagnostics(sqrtht_all,invA_all,PAI_all,PHI_all,N,K,MCMCdraws);
    end
    
    %% compute out-of-sample forecasts
    
    % a word on parfor strategy:
    % to make matlab better see the intended use of sliced variabes, use
    % local temp variables and then copy those into the slices at end of
    % loop
    
    
    % cumulate realizations and predictions if necessary
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
    
    %% collect PAI moments
    PAImedian(:,:,ndxT)       = median(PAI_all,3);
    PAImean(:,:,ndxT)         = mean(PAI_all,3);
    PAIstdev(:,:,ndxT)        = std(PAI_all,1,3);
    PAIquantiles(:,:,:,ndxT)  = prctile(PAI_all,setQuantiles,3);
    
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
    
    
    % pure SV part
    SVdraws                       = NaN(N,Tdata,MCMCdraws);
    stochvol                      = NaN(N, Tdata, MCMCdraws);
    stochvol(:,p+1:thisT, :)      = sqrtht_all .* exp(-.5 * SVtscalelog2_all);
    for m=1:MCMCdraws
        invA     = invA_all(:,:,m);
        for t=1:Tdata
            SVdraws(:,t,m) = sqrt(diag(invA*diag(stochvol(:,t,m).^2)*invA'));
        end
    end
    
    drawsLambdaSVmid(:,:,ndxT)     = median(SVdraws, 3);
    drawsLambdaSVtails(:,:,:,ndxT) = prctile(SVdraws, setQuantiles, 3);
    
    
    
    
    
    %% copy results into sliced variables
    
    fcstYmvlogscoreDraws(:,ndxT)  = logscoredraws;
    maxlogscoredraw               = max(logscoredraws);
    fcstYmvlogscore(:,ndxT)       = log(mean(exp(logscoredraws - maxlogscoredraw))) + maxlogscoredraw;
    fcstYmvlogscore2(:,ndxT)      = yMVlogscore;
    
    
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
    drawsSVtdof(:,:,ndxT)      = SVtdof_all;
    
    if ~doLoMem
        drawsPAI(:,:,:,ndxT)  = PAI_all;
        drawsPHI(:,:,ndxT)    = PHI_all;
        drawsINVA(:,:,:,ndxT) = invA_all;
        % prepare dummy to make parfor work
        dummy                      = NaN(N, Tdata, MCMCdraws);
        dummy(:,p+1:thisT, :)      = sqrtht_all;
        drawsSQRTHT(:, :, :, ndxT) = dummy;
        
        
        % prepare dummy to make parfor work
        dummy                      = NaN(N, Tdata, MCMCdraws);
        dummy(:, p+1:thisT,:)      = SVtscalelog2_all;
        drawsSVtscalelog2(:,:,:,ndxT)  = dummy;
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



%% SVtscale
close all
if ~doLoMem
    % collect QRT SVtscales
    thisSVtscaleQRT       = NaN(N,length(ydates));
    thisSVtscaleQRTtails  = NaN(N,2,length(ydates));
    
    drawsSVtscale = exp(0.5 * drawsSVtscalelog2);
    theseSVtscaleVintages = squeeze(median(drawsSVtscale, 3)); % returns N x Tdata x Njumpoff
    theseSVtscaleVintagesTails = prctile(drawsSVtscale, [normcdf(1) 90], 3); % returns N x Tdata x Njumpoff
    for v = 1 : Njumpoffs
        ndx = find(~isnan(theseSVtscaleVintages(1,:,v)),1, 'last'); % note: sufficient to check only for first element
        if isnan(theseSVtscaleVintages(1,ndx,v))
            error houston
        end
        thisSVtscaleQRT(:,ndx)          = theseSVtscaleVintages(:,ndx,v);
        thisSVtscaleQRTtails(:,:,ndx)   = squeeze(theseSVtscaleVintagesTails(:,ndx,:,v));
    end
    
    thisSVtscaleT   = theseSVtscaleVintages(:,:,end);
    
    % plot SVt scale QRT vs final
    bardates = ydates(Tjumpoffs);
    bardates = cat(1, bardates(1) - 15, bardates, bardates(end) + 15);
    for n = 1 : N
        
        thisfig = figure;
        hold on
        h1 = bar(ydates(Tjumpoffs), thisSVtscaleT(n,Tjumpoffs));
        h2 = plot(ydates(Tjumpoffs), thisSVtscaleQRT(n,Tjumpoffs), 'rd', 'linewidth', 5);
        h3 = plot(ydates(Tjumpoffs), squeeze(thisSVtscaleQRTtails(n,1,Tjumpoffs)), 'r--', 'linewidth', 2);
        h4 = plot(ydates(Tjumpoffs), squeeze(thisSVtscaleQRTtails(n,2,Tjumpoffs)), 'r:', 'linewidth', 2);
        xtickdates(bardates)
        legend([h1 h2 h3 h4], 'Final (median)', 'QRT (median)', 'QRT (upper 84)', 'QRT (upper 90)')
        title(sprintf('SV-t scale for %s', Ylabels{n}))
        
        wrapthisfigure(thisfig, sprintf('SVtscale%d', n),wrap)
        
    end
end


clear drawsSVtscale


%% plot COVIDSV
% ndx = ydates >= datenum(2020,1,1);
% thesedates = ydates(ndx);
% 
% % theseJumpoffs = [1 4 9 13];
% if sum(ndx) > 9
% theseJumpoffs = [1 4 9 sum(ndx)];
% else
%     theseJumpoffs = 3 : sum(ndx);
% end
% 
% for n = 1 : N
%     thisfig = figure;
%     subplot(1,2,1)
%     ax = gca;
%     set(ax, 'fontsize', 16)
%     hold on
%     
%     lineTypes = {'-d', '--o', 'd-.', ':o'};
%     iter = 0;
%     for ndxT = theseJumpoffs
%         iter = iter + 1;
%         plot(thesedates, drawsSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
%         if iter == length(lineTypes)
%             iter = 0;
%         end
%     end
%     ylim([0 max(ylim)])
%     xtickdates(thesedates)
%     title('including t scales')
%     
%     subplot(1,2,2)
%     axo = gca;
%     set(axo, 'fontsize', 16)
%     hold on
%     lineTypes = {'-d', '--o', 'd-.', ':o'};
%     iter = 0;
%     for ndxT = theseJumpoffs
%         iter = iter + 1;
%         plot(thesedates, drawsLambdaSVmid(n,ndx,ndxT), lineTypes{iter}, 'linewidth', 2)
%         if iter == length(lineTypes)
%             iter = 0;
%         end
%     end
%     ylim([0 max([ylim(axo), ylim(ax)])])
%     xtickdates(thesedates)
%     legend(datestr(ydates(Tjumpoffs(theseJumpoffs)), 'yyyy:mm'), 'location', 'best', 'box', 'off')
%     title('ex t scales')
%     
%     sgtitle(sprintf('\\bf %s FE-SV', Ylabels{n}))
%     
%     %         legend(datestr(ydates(Tjumpoffs), 'yyyy:mm'), 'location', 'northwest', 'box', 'off')
%     
%     
%     wrapthisfigure(thisfig, sprintf('covidSV%s', ncode{n}), wrap)
%     
% end

%% plot tdof
% generate histograms
% tdofposterior = NaN(length(tdofGrid), Njumpoffs, N);
% for n = 1 : N
%     for j = 1 : Njumpoffs
%         % tposterior(n,:) = ksdensity(SVtdof_all(n,:), tdofGrid);
%         tdofposterior(:,j,n) = histcounts(drawsSVtdof(n,:,j), [tdofGrid tdofGrid(end) + 1], 'Normalization', 'probability');
%     end
% end
% 
% % plot histograms (As surface)
% for n = 1 : N
%     thisfig = figure;
%     surf(ydates(Tjumpoffs), tdofGrid, tdofposterior(:,:,n));
%     shading interp
%     yticks([3, 5 : 5 : tdofGrid(end)])
%     datetick('x', 'keeplimits', 'keepticks')
%     title(Ylabels{n})
%     wrapthisfigure(thisfig, sprintf('SVt-dofposterior-%s', ncode{n}),wrap)
% end

% medians
tdofmed = squeeze(median(drawsSVtdof, 2));
tdofCI  = prctile(drawsSVtdof, [5 95], 2);
for n = 1 : N
    thisfig = figure;
    plotCI(tdofmed(n,:)', squeeze(tdofCI(n,:,:))', ydates(Tjumpoffs));
    yticks([3, 5 : 5 : tdofGrid(end)])
    ylim([2 40])
    datetick('x', 'keeplimits', 'keepticks')
    title(Ylabels{n})
    wrapthisfigure(thisfig, sprintf('SVt-dofposterior-med-%s', ncode{n}),wrap)
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

%% plot logscores
this = figure;
hold on
h1 = plot(ydates(Tjumpoffs), fcstYmvlogscore, 'r-', 'linewidth', 2);

% jack = mean(fcstYmvlogscoreDraws,1);
% plot(ydates(Tjumpoffs), jack, 'm-', 'linewidth', 2)

h2 = plot(ydates(Tjumpoffs), fcstYmvlogscore2 , 'k--', 'linewidth', 2);
legend([h1 h2], 'mixture approx.', 'Gaussian approx.', 'location', 'best')
xtickdates(ydates(Tjumpoffs))
wrapthisfigure(this, 'MVlogscore', wrap)

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
    'drawsSV*', 'drawsLambda*', ...
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
clear drawsLambdaSVtails drawsSVtails drawsMaxVARroot 
clear fcstYmvlogscore*Draws
save(matfilename, varlist{:}, '-v7.3');

%% wrap up
dockAllFigures
finishwrap

