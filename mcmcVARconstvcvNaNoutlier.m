function [PAI_all, SIGMA_all, sqrtht_all, Ydraws_all, Ynan, ...
    fcstYdraws, fcstYhat, fcstYhatRB] = mcmcVARconstvcvNaNoutlier(thisT, MCMCdraws, ...
    Ofactor, p, np, data0, ydates0, ...
    minnesotaPriorMean, doTightPrior, ...
    ndxYIELDS, ELBbound, check_stationarity, ...
    fcstNdraws, fcstNhorizons, rndStream) %#ok<INUSL>

% mcmc of BVAR-SV without shadowrate sampling (treating funds rate as
% regular data) predictive density for funds rate is truncated at given
% value of ELBbound

if isempty(Ofactor)
    Ofactor = 5;
end

doprogress = false;

%% get TID
% used to provide context for warning messages
TID   = parid;

%% truncate sample
samEnd = ydates0(thisT);
ndx    = ydates0 <= samEnd;

data   = data0(ndx,:);

%% --------------------------- OPTIONS ----------------------------------
if doTightPrior
    theta=[0.05 0.5 100 2];       % hyperparameters of Minnesota prior:
else
    theta=[0.1 0.5 100 2];       % hyperparameters of Minnesota prior:
end
% [lambda1 lambda2 int lambda3], int is the
% prior on the intercept. lambda1, lambda2
% and lambda3 are as in equation (42) with
% lambda1 the overall shrinkage, lambda2 the
% cross shrinkage and lambda 3 the lag decay
% (quadratic if =2). Note lambda2~=1 implies
% the prior becomes asymmetric across eqation,
% so this would not be implementable in the
% standard conjugate setup.
% Minn_pmean = 0;              % Prior mean of the 1-st own lag for each
% equation. For nonstationary variables, this
% is usually set to 1. For transformed
% stationary variables this is set to 0.

burnin            = 2 * ceil(0.1*MCMCdraws);    % burn in
MCMCreps           = MCMCdraws + burnin;        % total MCMC draws


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

%% prepare missing data



% determine SW outliers
dev             = abs(data - median(data));
iqr             = range(prctile(data, [25 75]));
dataOnan        = dev > Ofactor * iqr;
dataOnan(1:p,:) = false; % initial conditions taken as data
Ynan            = dataOnan(p+1:end,:);

Tstar  = find(any(Ynan,2), 1, 'first') - 1; % number of last NaN free obs in the VAR

%% allocate memory for out-of-sample forecasts
Ndraws      = fcstNdraws / MCMCdraws;
if mod(fcstNdraws, MCMCdraws) ~= 0
    error('fcstNdraws must be multiple of MCMCdraws')
end
fcstYdraws  = NaN(N,fcstNhorizons, Ndraws, MCMCdraws); % see reshape near end of script


yhatdraws   = NaN(N,fcstNhorizons, MCMCdraws);

%% prepare state space for forecasting
fcstA                  = zeros(K,K);
fcstA(1,1)             = 1; % unit root for the constant
fcstA(1+N+1:end,2:end) = [eye(N*(p-1)),zeros(N*(p-1),N)]; % fill in lower part of companion form


ndxfcstY          = 1+(1:N);
fcstB             = zeros(K,N);
fcstB(ndxfcstY,:) = eye(N);



%% -----------------Prior hyperparameters for bvar model

% Prior on conditional mean coefficients, use Minnesota setup
% constructed from actual data, note that initial Y was constructed from data
ARresid=NaN(T-1,N);
for i=1:N
    yt_0=[ones(T-1,1) Y(1:end-1,i)];
    yt_1=Y(2:end,i);
    ARresid(:,i)=yt_1-yt_0*(yt_0\yt_1);
end
AR_s2= diag(diag(ARresid'*ARresid))./(T-2);

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

% Pai~N(vec(MU_pai),OMEGA_pai), equation 7.
diagOMEGA_pai = [sigma_const;reshape(diag(Pi_pv),Klagreg,N)];
OMEGA_pai     = diag(diagOMEGA_pai(:)); % prior variance of Pai
MU_pai        = [zeros(1,N);reshape(Pi_pm,Klagreg,N)];                   % prior mean of Pai

% A~N(MU_A,inv(OMEGA_A_inv)), equation 8.
MU_A = NaN(N-1,N);
OMEGA_A_inv = NaN(N-1,N-1,N);
for i = 2:N;
    MU_A(1:i-1,i) = zeros(i-1,1);             % prior mean of A
    OMEGA_A_inv(1:i-1,1:i-1,i) = 0*eye(i-1);  % prior precision of A
end;

% constVCV prior

% training sample prior 
trainingT = 40;
SigmaBar = ARresid(1:trainingT,:)' * ARresid(1:trainingT,:) / trainingT;
% SigmaBar is used for MCMC init (and prior, if desired)
 
% SigmaDof = N + 2;
% SigmaT   = SigmaBar * (SigmaDof - N - 1);

% uninformative prior
SigmaDof = 0;
SigmaT   = zeros(N);

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
%% >>>>>>>>>>>>>>>>>>>>>>>>>> Gibbs sampler <<<<<<<<<<<<<<<<<<<<<<<<<<<

% Storage arrays for posterior draws
PAI_all     = NaN(K,N,MCMCdraws); 
 
sqrtht_all  = NaN(N,T,MCMCdraws);

SIGMA_all   = NaN(N,N,MCMCdraws);
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

m = 0;
while m < MCMCreps % using while, not for loop to allow going back in MCMC chain
    
    if m == 0
        % initializations
        %         Sigma      = SigmaT / (SigmaDof - N - 1);
        Sigma      = SigmaBar;
        
        cholSigma  = chol(Sigma)';
        volSigma   = diag(cholSigma)';
        invA_      = bsxfun(@rdivide, cholSigma, volSigma);
		A_         = EYEn / invA_;
        sqrtht    = repmat(volSigma, T, 1);
        
        
        PREVdraw.A_         = A_;                                 
        trainingT           = sum(ydates < datenum(2020,3,1)) - p;
        PREVdraw.PAI        = X(1:trainingT,:)\Y(1:trainingT,:);

        PREVdraw.sqrtht     = sqrtht;
        % PREVdraw.cholSigma  = cholSigma;
        
		        % init PAI for missing data draws
        PREVdraw.PAI        = MU_pai; % + sqrt(diagOMEGA_pai) .* randn(rndStream, K, N);

    end % m == 0
    m = m + 1;
    
    
    % init with previous draws values
    % Sigma      = PREVdraw.Sigma; % not needed
    A_         = PREVdraw.A_;
    sqrtht    = PREVdraw.sqrtht;
    PAI        = PREVdraw.PAI;
    
    % if mod(m,10) == 0; clc; disp(['percentage completed:' num2str(100*m/MCMCreps) '%']); toc; end
    
    %% STEP 1 draw missing values
    
    if ~isempty(Tstar)
        miss.Y               = Y(Tstar+1:end,:)';
        miss.Y(miss.Ynan)    = 0; % set missing values to zeros;
        miss.A(ndxfcstY, :)  = PAI';
        
        miss.B(ndxfcstY,:)   = invA_;
        
        % map SV into sqrtSigma
        miss.sqrtSigma = sqrtht(Tstar+1:end,:)';
        
        missStates = stateABCdraws(miss.A, miss.B, miss.C, miss.Y, miss.Ynan, ...
            miss.X0, miss.sqrtSigma, 1, rndStream);
        
        % check that draws of data match data
        Ydraws = missStates(ndxfcstY,:);
        delta  = bsxfun(@minus, Ydraws, miss.Y);
        delta(miss.Ynan) = 0;
        if max(abs(delta(:))) > 1e-9
            warning('missing data sampler does not replicate data, MAE=%e', max(abs(delta(:))))
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
    
    %% STEP 2b: Draw from the conditional posterior of PAI, equation 10.
    stationary=0;
    while stationary==0;
        
        
        % CCM: This is the only new step (triangular algorithm).
        % PAI=triang(Y,X,N,K,T,invA_,sqrtht,iV,iVb_prior,rndStream);
 
        PAI=CTA(Y,X,N,K,T,A_,sqrtht,iV,iVb_prior,PAI,rndStream);

        
        if (check_stationarity==0 || max(abs(eig([PAI(ndxKlagreg,:)' ; comp]))) < 1); stationary = 1; end;
    end
    RESID = Y - X*PAI; % compute the new residuals
    
    %% STEP 2c: Draw the covariances, equation 11.
    Sigma = bayesVCVgibbsDraw1(SigmaT, SigmaDof, RESID, rndStream);
    
    cholSigma  = chol(Sigma)';
    volSigma   = diag(cholSigma)';
    invA_      = bsxfun(@rdivide, cholSigma, volSigma);
	A_         = EYEn / invA_;
    sqrtht     = repmat(volSigma, T, 1);
    
    
    
    %% post burnin: store draws and draw from oos-predictive density
    if m > burnin;
        
        thisdraw = m-burnin;
        % STORE DRAWS
        PAI_all(:,:,thisdraw)      = PAI;
        sqrtht_all(:,:,thisdraw)   = sqrtht'; % note the transpose
        SIGMA_all(:,:,thisdraw)    = Sigma;
        Ydraws_all(:,:,thisdraw)   = Y;
        
        
        %% compute OOS draws
        % 
        
        
        % draw random numbers
        zdraws  = randn(rndStream, N, fcstNhorizons, Ndraws);
        
        for nn = 1 : Ndraws
            nushocks            = zeros(N, fcstNhorizons+1); % padding with a line of zeros for use with ltitr
            nushocks(:,1:end-1) = cholSigma * zdraws(:,:,nn);
            
            % c) update VAR companion form and iterate
            fcstA(ndxfcstY, :)       = PAI';
            fcstX0 = Xjumpoff;
            
            if isempty(ELBbound)
                fcstXdraws                  = ltitr(fcstA, fcstB, nushocks', fcstX0); % faster forecast simulation using ltitr
                fcstYdraws(:,:,nn,thisdraw) = fcstXdraws(2:end,ndxfcstY)';
            else
                % need sequential simulation
                for n = 1 : fcstNhorizons
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
    PREVdraw.A_         = A_;
    PREVdraw.sqrtht    = sqrtht;
    % PREVdraw.cholSigma  = cholSigma;
    PREVdraw.PAI        = PAI;
    
    if doprogress
        progressbar(m / MCMCreps)
    end
    
end %end of the Gibbs sampler

fcstYdraws     = reshape(fcstYdraws, N, fcstNhorizons, fcstNdraws);

fcstYhatRB     = mean(yhatdraws,3);
fcstYhat       = mean(fcstYdraws,3);





fprintf('DONE with thisT %d, TID %d \n', thisT, TID)

return

