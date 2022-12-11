function [PAI_all, PHI_all, invA_all, sqrtht_all, ...
    Ydraws_all, Ynan, ...
    fcstYdraws, fcstYhat, fcstYhatRB, fcstSVdraws, fcstYdraws2, fcstYhat2, ...
    fcstLogscoreDraws, fcstLogscoreNaNdraws] = mcmcVARSVnanOutlier(thisT, MCMCdraws, ...
    Ofactor, p, np, data0, ydates0, ...
    minnesotaPriorMean, doTightPrior, doRobustPrior, ...
    ndxYIELDS, ELBbound, check_stationarity, ...
    yrealized, ...
    fcstNdraws, fcstNhorizons, rndStream, doprogress, ndxFOREX)

% mcmc of BVAR-SV without shadowrate sampling (treating funds rate as
% regular data) predictive density for funds rate is truncated at given
% value of ELBbound

if nargin < 18
    doprogress = false;
end

if nargin < 19
    ndxFOREX = [];
end

if isempty(Ofactor)
    Ofactor = 5;
end

%% get TID
% used to provide context for warning messages
TID   = parid;

%% truncate sample
samEnd = ydates0(thisT);
ndx    = ydates0 <= samEnd;

data   = data0(ndx,:);
ydates = ydates0(ndx); % needed if doRobustPrior is true

%% --------------------------- OPTIONS ----------------------------------
if doTightPrior
    theta=[0.05 0.5 100 2];       % hyperparameters of Minnesota prior:
else
    theta=[0.1 0.5 100 2];       % hyperparameters of Minnesota prior:
end
% [lambda1 lambda2 int lambda3], int is the
% prior on the intercept. 

% lambda1 is the overall shrinkage, lambda2 the
% cross shrinkage and lambda 3 the lag decay
% (quadratic if =2). Note lambda2~=1 implies
% the prior becomes asymmetric across eqation,
% so this would not be implementable in the
% standard conjugate setup.

burnin            = 2 * ceil(0.1*MCMCdraws);    % burn in
MCMCreps           = MCMCdraws + burnin;    % total MCMC draws


%% -------------------------Create data matrices-------------------------
% preliminary since subject to missing data
% pointers
[Nobs,N]=size(data);
% matrix X - skeletons
lags=zeros(Nobs,N*p);
for l=1:p
    lags(p+1:Nobs,(N*(l-1)+1):N*l) = data(p+1-l:Nobs-l,1:N);
end
X = [ones(Nobs-p,1) lags(p+1:Nobs,:)];
% trim Y
Y = data(p+1:end,:); % outliers will be taken out in first step of MCMC
% update pointers
[T,K]=size(X);
Klagreg    = K - 1; % number of lag regressors (without intercept)
ndxKlagreg = 1 + (1 : Klagreg); % location of the Klag regressors in X

% prepare Xjumpoff in case there are no outliers to be drawn
Xjumpoff    = zeros(K,1);
Xjumpoff(1) = 1;
for l=1:p
    Xjumpoff(1+(l-1)*N+(1:N)) = data(Nobs-(l-1),1:N);
end

Xjumpoff2 = Xjumpoff;

logtwopi  = log(2 * pi);
Nlogtwopi = N * logtwopi;

%% prepare missing data



% determine SW outliers
dev             = abs(data - median(data));
iqr             = range(prctile(data, [25 75]));
outlierNdx      = dev > Ofactor * iqr;

% ignore outliers for FOREX pre 1980s
if any(ndxFOREX)
    if ydates(end) < datenum(1985,1,1)
        outlierNdx(:,ndxFOREX) = false;
    end
end

dataOnan        = outlierNdx;
dataOnan(1:p,:) = false; % initial conditions taken as data
Ynan            = dataOnan(p+1:end,:);

Tstar  = find(any(Ynan,2), 1, 'first') - 1; % number of last NaN free obs in the VAR

ndxOutlierRealized1 = abs(yrealized(:,1)' - median(data)) > Ofactor * iqr;

%% allocate memory for out-of-sample forecasts
Ndraws      = fcstNdraws / MCMCdraws;
if mod(fcstNdraws, MCMCdraws) ~= 0
    error('fcstNdraws must be multiple of MCMCdraws')
end
fcstYdraws         = NaN(N,fcstNhorizons, Ndraws, MCMCdraws); % see reshape near end of script
fcstSVdraws        = NaN(N,fcstNhorizons, Ndraws, MCMCdraws);
fcstLogscoreDraws     = NaN(Ndraws, MCMCdraws); % see reshape near end of script
fcstLogscoreNaNdraws  = NaN(Ndraws, MCMCdraws); % see reshape near end of script

fcstYdraws2 = NaN(N,fcstNhorizons, Ndraws, MCMCdraws); % see reshape near end of script

yhatdraws   = NaN(N,fcstNhorizons, MCMCdraws);

%% prepare state space for forecasting
fcstA                  = zeros(K,K);
fcstA(1,1)             = 1; % unit root for the constant
fcstA(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form


ndxfcstY          = 1+(1:N);
fcstB             = zeros(K,N);
fcstB(ndxfcstY,:) = eye(N);


%% -----------------Prior hyperparameters for bvar model

if doRobustPrior
    trainingT = sum(ydates < datenum(2020,3,1)) - p;
else
    trainingT = T;
end    

% Prior on conditional mean coefficients, use Minnesota setup
ARresid=NaN(trainingT-1,N);
ARrho  = NaN(2,N);
for i=1:N
    yt_0=[ones(trainingT-1,1) Y(1:trainingT-1,i)];
    yt_1=Y(2:trainingT,i);
    ARrho(:,i) = (yt_0' * yt_0) \ (yt_0' * yt_1);
    ARresid(:,i)=yt_1-yt_0*(yt_0\yt_1);
    
    % resid = yt_1 - yt_0 * ARrho(:,i);
    % checkdiff(resid, ARresid(:,i));
end
AR_s2= diag(diag(ARresid'*ARresid))./(trainingT-2);

% recompute ARreisd over full sample to create start values for SV
ARresid=NaN(T-1,N);
for i=1:N
    yt_0=[ones(T-1,1) Y(1:T-1,i)];
    yt_1=Y(2:T,i);
    ARresid(:,i) = yt_1 - yt_0 * ARrho(:,i);
end

Pi_pm=zeros(N * Klagreg,1); Pi_pv=eye(N * Klagreg); co=0;
% todo: Pi_pv could be encoded as vector, only using diagonal elements anyhow
sigma_const = NaN(1,N);
for i=1:N
    sigma_const(i)=AR_s2(i,i)*theta(3);  % this sets the prior variance on the intercept
    for l=1:p; %#ok<*NOSEL>
        for j=1:N
            co=co+1;
            if (i==j)
                if l==1
                    Pi_pm(co)=minnesotaPriorMean(i); % this sets the prior means for the first own lag coefficients.
                end
                Pi_pv(co,co)=theta(1)/(l^theta(4)); % prior variance, own lags
            else
                Pi_pv(co,co)=(AR_s2(i,i)/AR_s2(j,j)*theta(1)*theta(2)/(l^theta(4))); % prior variance, cross-lags
            end
        end
    end
end

% Pai~N(vec(MU_pai),OMEGA_pai)
diagOMEGA_pai = [sigma_const;reshape(diag(Pi_pv),Klagreg,N)];
OMEGA_pai     = diag(diagOMEGA_pai(:)); % prior variance of Pai
MU_pai        = [zeros(1,N);reshape(Pi_pm,Klagreg,N)];                   % prior mean of Pai

% A~N(MU_A,inv(OMEGA_A_inv))
MU_A = NaN(N-1,N);
OMEGA_A_inv = NaN(N-1,N-1,N);
for i = 2:N;
    MU_A(1:i-1,i) = zeros(i-1,1);             % prior mean of A
    OMEGA_A_inv(1:i-1,1:i-1,i) = 0*eye(i-1);  % prior precision of A
end;

% PHI~IW(s_PHI,d_PHI)
d_PHI = N+3;                 % prior dofs
s_PHI = d_PHI*(0.15*eye(N)) * 12 / np; % prior scale, where eye(N)=PHI_ 

% prior on initial states
Vol_0mean     = zeros(N,1);   % time 0 states prior mean
% Vol_0var      = 100*eye(N);   % time 0 states prior variance
Vol_0vcvsqrt  = 10*eye(N); % chol(Vol_0var)';

%% prepare state space for missing data sampling
% state space with K states, N observables and N shocks
% constant transition matrix miss.A, constant shock loadings miss.B (ex SV)
% and time-varying measurement loadings miss.C (to account for missing
% shadow-rate values)

if ~isempty(Tstar)
    miss.T              = T - Tstar;
    miss.Ynan           = Ynan(Tstar+1:end,:)'; % note the transpose, needed for call of sampler later
    
    % prepare state space matrices
    miss.A               = fcstA;
    miss.B               = fcstB;
    
    % prepare missC
    miss.C               = zeros(N,K,miss.T); % time-varying to account for missing values
    miss.C(:,ndxfcstY,:) = repmat(eye(N), [1 1 miss.T]);
    % zero out missing values
    for t = 1 : miss.T
        for n = 1 : N
            if miss.Ynan(n,t)
                miss.C(n,:,t) = 0;
            end
        end
    end
    
    % encode "prior" over initial conditions (which are fixed)
    miss.X0  = X(Tstar+1,:)';
end

%% PREPARE KSC sampler

[gridKSC, gridKSCt, logy2offset] = getKSC7values(T, N);

%% >>>>>>>>>>>>>>>>>>>>>>>>>> Gibbs sampler <<<<<<<<<<<<<<<<<<<<<<<<<<<

% Storage arrays for posterior draws
PAI_all     = NaN(K,N,MCMCdraws);
PHI_all     = NaN(N*(N-1)/2+N,MCMCdraws);
invA_all    = NaN(N,N,MCMCdraws);
sqrtht_all  = NaN(N,T,MCMCdraws);
Ydraws_all  = NaN(T,N,MCMCdraws);

% define some useful matrices prior to the MCMC loop
% PAI  = zeros(K,N);                                            % pre-allocate space for PAI
comp = [eye(N*(p-1)),zeros(N*(p-1),N)];                       % companion form
iV   = diag(1./diag(OMEGA_pai)); iVb_prior=iV*vec(MU_pai);    % inverses of prior matrices
EYEn = eye(N);

%% start of MCMC loop
% the algorithm is as described in page 15, but it starts from step 2b,
% this is the same as starting from step 1 but is more convenient as
% it requires less initializations (one can think of steps 2b to 2d in
% repetition 1 as an initialization).

if doprogress
    progressbar(0);
end
m = 0;
while m < MCMCreps % using while, not for loop to allow going back in MCMC chain
    
    if m == 0
        % initializations
        A_                  = eye(N);                                  % initialize A matrix
        PREVdraw.A_         = A_; 
        PREVdraw.invA_      = EYEn / A_;
        
        trainingT           = sum(ydates < datenum(2020,3,1)) - p;
        PREVdraw.PAI        = X(1:trainingT,:)\Y(1:trainingT,:);
        % PREVdraw.PAI        = MU_pai; % + sqrt(diagOMEGA_pai) .* randn(rndStream, K, N);
        
        PREVdraw.sqrtht     = sqrt([ARresid(1,:).^2; ARresid.^2]);     % Initialize sqrt_sqrtht
        PREVdraw.sqrtht(PREVdraw.sqrtht < 1e-3) = 1e-3;
        PREVdraw.sqrtht(PREVdraw.sqrtht > 25)   = 25;
        
        PREVdraw.Vol_states = 2*log(PREVdraw.sqrtht)';  %note the transpose               % Initialize states
        
        %         PREVdraw.PHI_       = 0.0001*eye(N);                           % Initialize PHI_, a draw from the covariance matrix W
        PREVdraw.sqrtPHI_   = sqrt(0.0001)*eye(N);                           % Initialize PHI_, a draw from the covariance matrix W
        
    end % m == 0
    m = m + 1;
    
    
    % init with previous draws values
    A_         = PREVdraw.A_;
    invA_      = PREVdraw.invA_;
    sqrtht     = PREVdraw.sqrtht;
    Vol_states = PREVdraw.Vol_states;
    %     PHI_       = PREVdraw.PHI_;
    PAI        = PREVdraw.PAI;
    sqrtPHI_   = PREVdraw.sqrtPHI_;
    
    % if mod(m,10) == 0; clc; disp(['percentage completed:' num2str(100*m/MCMCreps) '%']); toc; end
    
    %% STEP 1 draw missing values
    
    if ~isempty(Tstar)
        miss.Y               = Y(Tstar+1:end,:)';
        miss.Y(miss.Ynan)    = 0; % set missing values to zeros;
        
        miss.A(ndxfcstY, :)  = PAI';
        miss.B(ndxfcstY,:)   = invA_;
        
        % map SV into sqrtSigma
        miss.sqrtSigma = sqrtht(Tstar+1:end,:)';
        
        missStates = stateABnanDraw1(miss.A, miss.B, ndxfcstY, miss.Y, miss.Ynan, ...
            miss.X0, miss.sqrtSigma, rndStream);
        
        
        % check that draws of data match data
        Ydraws = missStates(ndxfcstY,:);
        delta  = bsxfun(@minus, Ydraws, miss.Y);
        delta(miss.Ynan) = 0;
        if max(abs(delta(:))) > 1e-4
            error('thisT = %d, MCMC = %d: missing data sampler does not replicate data, MAE=%e', thisT, m, max(abs(delta(:))))
        end
        
        % patch new data draws into Y and X
        Y(Tstar+1:end,:)   = Ydraws';
        
        % reconstruct X
        dataO              = data;
        dataO(dataOnan)    = Y(Ynan);
        %     checkdiff(Y, dataO(p+1:end,:));
        lags=zeros(Nobs,N*p);
        for l=1:p
            lags(p+1:Nobs,(N*(l-1)+1):N*l) = dataO(p+1-l:Nobs-l,1:N);
        end
        X = [ones(Nobs-p,1) lags(p+1:Nobs,:)];
        
        % generate state vector for forecast jumpoff
        Xjumpoff    = zeros(K,1);
        Xjumpoff(1) = 1;
        for l=1:p
            Xjumpoff(1+(l-1)*N+(1:N)) = dataO(Nobs-(l-1),1:N);
        end
    end
    
    %% STEP 2b: Draw from the conditional posterior of PAI
    stationary=0;
    while stationary==0;
        
        
        % CCM: This is the only new step (triangular algorithm).
        % PAI=triang(Y,X,N,K,T,invA_,sqrtht,iV,iVb_prior,rndStream);
 
        PAI=CTA(Y,X,N,K,T,A_,sqrtht,iV,iVb_prior,PAI,rndStream);

        
        if (check_stationarity==0 || max(abs(eig([PAI(ndxKlagreg,:)' ; comp]))) < 1); stationary = 1; end;
    end
    RESID = Y - X*PAI; % compute the new residuals
    
    %% STEP 2c: Draw the covariances
    for ii = 2:N
        % weighted regression to get Z'Z and Z'z (in Cogley-Sargent 2005 notation)
        y_spread_adj=RESID(:,ii)./sqrtht(:,ii);
        %         X_spread_adj=[]; for vv=1:ii-1;  X_spread_adj=[X_spread_adj RESID(:,vv)./sqrtht(:,ii)]; end  %#ok<AGROW>
        X_spread_adj = RESID(:,1 : ii - 1) ./ sqrtht(:,ii); % note: use of implicit vector expansion
        ZZ=X_spread_adj'*X_spread_adj; Zz=X_spread_adj'*y_spread_adj;
        % computing posteriors moments
        Valpha_post = (ZZ + OMEGA_A_inv(1:ii-1,1:ii-1,ii))\eye(ii-1);
        alpha_post  = Valpha_post*(Zz + OMEGA_A_inv(1:ii-1,1:ii-1,ii)*MU_A(1:ii-1,ii));
        % draw and store
        alphadraw   = alpha_post+chol(Valpha_post,'lower')*randn(rndStream,ii-1,1);
        A_(ii,1:ii-1)= -1*alphadraw';
        % EM: Note: need to init A_ at least once to set A_(:,end)
    end
    invA_=A_\EYEn; % compute implied draw from A^-1, needed in step 2b.
    
    %% STEP 2d and STEP 1: Draw mixture states and then volatility states
    
    
    logy2 = log((RESID*A_').^2 + logy2offset);
    
    
    [Vol_states, ~, eta] = StochVolKSCcorrsqrt(logy2', Vol_states, sqrtPHI_, Vol_0mean, Vol_0vcvsqrt, gridKSC, gridKSCt, N, T, rndStream);
    
    sqrtht = exp(Vol_states/2)'; %compute sqrtht^0.5 from volatility states, needed in step 2b
    
    %% STEP 2a: Draw volatility variances
    
    Zdraw     = randn(rndStream, N, T + d_PHI);
    
    sqrtPHIpost = chol(s_PHI + eta*eta', 'lower');
    sqrtZZ      = chol(Zdraw * Zdraw'); % note: right uppper choleski
    sqrtPHI_    = sqrtPHIpost / sqrtZZ; % just a square root, not choleski
    PHI_        = sqrtPHI_ * sqrtPHI_'; % derive posterior draw from PHI, equation 12.
    
    
    %% post burnin: store draws and draw from oos-predictive density
    if m > burnin;
        
        thisdraw = m-burnin;
        % STORE DRAWS
        PAI_all(:,:,thisdraw)      = PAI;
        PHI_all(:,thisdraw)        = PHI_((tril(PHI_))~=0);
        invA_all(:,:,thisdraw)     = invA_;
        sqrtht_all(:,:,thisdraw)   = sqrtht'; % note the transpose
        Ydraws_all(:,:,thisdraw)   = Y;
        
        
        %% compute OOS draws
        % 
        
        % draw and scale SV shocks
        logSV0      = Vol_states(:,end); % Note: Vol_states record logs of *variances*
        logSVshocks = sqrtPHI_ * randn(rndStream, N, fcstNhorizons * Ndraws);
        logSVshocks = reshape(logSVshocks, N, fcstNhorizons, Ndraws);
        
        % draw random numbers
        zdraws  = randn(rndStream, N, fcstNhorizons, Ndraws);
        
        for nn = 1 : Ndraws
            
            logSV       = bsxfun(@plus, logSV0, cumsum(logSVshocks(:,:,nn),2));
            
            
            fcstSVdraws(:,:,nn,thisdraw)    = exp(logSV * 0.5);
            nushocks            = zeros(N, fcstNhorizons+1); % padding with a line of zeros for use with ltitr
            nushocks(:,1:end-1) = invA_ * (fcstSVdraws(:,:,nn,thisdraw) .* zdraws(:,:,nn));
            
            %% c) update VAR companion form and iterate
            fcstA(ndxfcstY, :)       = PAI';
            fcstX0 = Xjumpoff;
            fcstX2 = Xjumpoff2;
			
			%% d) logscore (anticipating no future outlier)
            
            % compute for full vector of realized values
            muX          = fcstA * fcstX0;
            muY          = muX(ndxfcstY);
            sqrtOmegaY   = invA_ * diag(fcstSVdraws(:,1,nn,thisdraw));
            logdetOmegaY = sum(logSV(:,1)); % logSV stores log variances!
            
            Ydev = sqrtOmegaY \ (yrealized(:,1) - muY);
            fcstLogscoreDraws(nn,thisdraw) = -.5 * (Nlogtwopi + logdetOmegaY + Ydev' * Ydev);
            
            % compute for non-outlier values of the realizations
            if any(ndxOutlierRealized1)
                thismuY        = muY(~ndxOutlierRealized1);
                OmegaY     = sqrtOmegaY(~ndxOutlierRealized1,:) * sqrtOmegaY(~ndxOutlierRealized1,:)';
                thissqrtOmegaY = chol(OmegaY, 'lower');
                logdetOmegaY = 2 * sum(log(diag(thissqrtOmegaY)));
                % checkdiff(logdetOmegaY, log(det(OmegaY)));
                
                Ydev = thissqrtOmegaY \ (yrealized(~ndxOutlierRealized1,1) - thismuY);
                
                thisN = length(Ydev);
                fcstLogscoreNaNdraws(nn,thisdraw) = -.5 * (thisN * logtwopi + logdetOmegaY + Ydev' * Ydev);

            else
                fcstLogscoreNaNdraws(nn,thisdraw) = fcstLogscoreDraws(nn,thisdraw);
            end
            % check for ELB
            if ~isempty(ELBbound)
                
                if any(yrealized(ndxYIELDS,1) <= ELBbound) 
                    warning('realized Yields were at ELB, but logscores are not yet adjusted')
                end
                
            end

            %% e) density draws
            if isempty(ELBbound)
                fcstXdraws                  = ltitr(fcstA, fcstB, nushocks', fcstX0); % faster forecast simulation using ltitr
                fcstYdraws(:,:,nn,thisdraw) = fcstXdraws(2:end,ndxfcstY)';

                fcstXdraws                   = ltitr(fcstA, fcstB, nushocks', fcstX2); % faster forecast simulation using ltitr
                fcstYdraws2(:,:,nn,thisdraw) = fcstXdraws(2:end,ndxfcstY)';

            else
                % need sequential simulation
                for n = 1 : fcstNhorizons
                    
                    % regular foreacst
                    fcstXdraw    = fcstA * fcstX0 + fcstB * nushocks(:,n);
                    ydraw        = fcstXdraw(ndxfcstY);
                    these = ydraw(ndxYIELDS);
                    if any(these < ELBbound)
                        these(these < ELBbound) = ELBbound;
                        ydraw(ndxYIELDS)        = these;
                        fcstXdraw(ndxfcstY)     = ydraw;
                    end
                    % collect draw
                    fcstYdraws(:,n,nn,thisdraw) = ydraw;
                    % prepare next iteration
                    fcstX0                   = fcstXdraw;
                    
                    % forecst jumping off actual data
                    fcstXdraw    = fcstA * fcstX2 + fcstB * nushocks(:,n);
                    ydraw        = fcstXdraw(ndxfcstY);
                    these = ydraw(ndxYIELDS);
                    if any(these < ELBbound)
                        these(these < ELBbound) = ELBbound;
                        ydraw(ndxYIELDS)        = these;
                        fcstXdraw(ndxfcstY)     = ydraw;
                    end
                    % collect draw
                    fcstYdraws2(:,n,nn,thisdraw) = ydraw;
                    % prepare next iteration
                    fcstX2                   = fcstXdraw;
                    
                    
                end % n
            end % if ELBbound
            
            

        end % nn
        
        % RB moments: mean
        fcstX0                   = Xjumpoff;
        nushocks(:)              = 0;
        fcstXdraws               = ltitr(fcstA, fcstB, nushocks', fcstX0);
        yhatdraws(:,:,thisdraw)  = fcstXdraws(2:end,ndxfcstY)';
        
    end
    
    
    %% store current draw into PREVdraw
    PREVdraw.invA_      = invA_;
    PREVdraw.A_         = A_;
    PREVdraw.sqrtht     = sqrtht;
    PREVdraw.Vol_states = Vol_states;
    PREVdraw.PAI        = PAI;
    %     PREVdraw.PHI_       = PHI_;
    PREVdraw.sqrtPHI_   = sqrtPHI_; % only needed for simulating SV forecasts
    
    if doprogress
        progressbar(m / MCMCreps)
    end
    
end %end of the Gibbs sampler

fcstYdraws        = reshape(fcstYdraws, N, fcstNhorizons, fcstNdraws);
fcstLogscoreDraws    = reshape(fcstLogscoreDraws, fcstNdraws, 1);
fcstLogscoreNaNdraws = reshape(fcstLogscoreNaNdraws, fcstNdraws, 1);

fcstYhatRB     = mean(yhatdraws,3);
fcstYhat       = mean(fcstYdraws,3);

fcstYdraws2    = reshape(fcstYdraws2, N, fcstNhorizons, fcstNdraws);
fcstYhat2      = mean(fcstYdraws2,3);

fcstSVdraws    = reshape(fcstSVdraws, N, fcstNhorizons, fcstNdraws);



fprintf('DONE with thisT %d, TID %d \n', thisT, TID)

return

