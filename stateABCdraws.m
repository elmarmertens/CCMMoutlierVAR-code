function Xdraws = stateABCdraws(A, B, C, Ydata, yNaNndx, X0, sqrtSigma, Ndraws, rndStream)
% STATEABCDRAW
% ....
% supposes X0 is deterministically given
%
%   Coded by  Elmar Mertens, em@elmarmertens.com


%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Ydata             = permute(Ydata, [1 3 2]); 
Nw                = size(B,2);

if nargin < 8 || isempty(Ndraws)
    Ndraws = 1;
end
if nargin < 9 || isempty(rndStream)
    rndStream = getDefaultStream;
end



%% init Variables and allocate memory
I                 = eye(Nx);

yDataNdx          = ~yNaNndx;

%% allocate memory
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, 1, T);
Yplustilde                  = zeros(Ny, Ndraws, T);
[XtT, Xttm1]                = deal(zeros(Nx, 1, T));
[Xplus, XplustT, Xplusttm1] = deal(zeros(Nx, Ndraws, T));


%% generate plus data

wplus   = randn(rndStream, Nw, Ndraws, T);
X0plus  = repmat(X0, [1 Ndraws]); % redundant given Matlab's implict matrix expansion since 2016a

%% Forward Loop: Kalman Forecasts
Sigmatt = zeros(Nx,Nx);
Xtt     = X0;
Xplustt = repmat(X0, 1, Ndraws);

disturbanceplus  = zeros(Nx, Ndraws, T);

for t = 1 : T
    
    % "plus" States and priors
    Bsv                     = B * diag(sqrtSigma(:,t));
    disturbanceplus(:,:,t)  = Bsv * wplus(:,:,t);
    BSigmaB                 = Bsv * Bsv';
    
    if t == 1
        % Xplus(:,:,t) = A * X0plus       + disturbanceplus(:,:,t);
        Xlagplus = X0plus;
    else
        % Xplus(:,:,t) = A * Xplus(:,:,t-1) + disturbanceplus(:,:,t);
        Xlagplus = Xplus(:,:,t-1);
    end
    
    Xplus(:,:,t) = A * Xlagplus + disturbanceplus(:,:,t);
    % priors
    Sigmattm1(:,:,t)        = A * Sigmatt * A' + BSigmaB;
    Xttm1(:,:,t)            = A * Xtt;
    Xplusttm1(:,:,t)        = A * Xplustt;
    
    
    % observed innovation
    Yplus            = C(:,:,t) * Xplus(:,:,t);
    SigmaYttm1       = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';
    
    
    Ytilde(:,:,t)         = Ydata(:,:,t)    - C(:,:,t) * Xttm1(:,:,t);
    Yplustilde(:,:,t)     = Yplus           - C(:,:,t) * Xplusttm1(:,:,t);
    
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) = eye(sum(yDataNdx(:,t))) / SigmaYttm1(yDataNdx(:,t),yDataNdx(:,t));
    
     
    % Kalman Gain
    K                       = (Sigmattm1(:,:,t) * C(:,:,t)') * invSigmaYttm1(:,:,t);
    ImKC(:,:,t)             = I - K * C(:,:,t);
    
    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t) * ImKC(:,:,t)'; % Joseph form for better numerical stability
    
    Xtt                     = Xttm1(:,:,t)     + K * Ytilde(:,:,t);
    Xplustt                 = Xplusttm1(:,:,t) + K * Yplustilde(:,:,t);
    
end

%% Backward Loop: Disturbance Smoother
XplustT(:,:,T)  = Xplustt;
XtT(:,T)        = Xtt;

StT             = C(:,:,T)' * (invSigmaYttm1(:,:,T) * Ytilde(:,T));
SplustT         = C(:,:,T)' * (invSigmaYttm1(:,:,T) * Yplustilde(:,:,T));


for t = (T-1) : -1 : 1
    Atilde      = A * ImKC(:,:,t);
    
    StT             = Atilde' * StT + C(:,:,t)' * (invSigmaYttm1(:,:,t) * Ytilde(:,:,t));
    XtT(:,:,t)      = Xttm1(:,:,t) + Sigmattm1(:,:,t) * StT;
    
    SplustT         = Atilde' * SplustT + C(:,:,t)' * (invSigmaYttm1(:,:,t) * Yplustilde(:,:,t));
    XplustT(:,:,t)  = Xplusttm1(:,:,t) + Sigmattm1(:,:,t) * SplustT;
end

%% sample everything together (and reorder output dimensions)
Xdraws  = bsxfun(@minus, XtT, XplustT) + Xplus;
Xdraws  = permute(Xdraws, [1 3 2]);



