function Xdraws = stateABCdraw1delta(A, B, C, Ydata, yNaNndx, X0, sqrtSigma, rndStream)
% STATEABCDRAW
% ....
% supposes X0 is deterministically given

% note: rows of C with yNaNndx == true must be zero (no checks!)

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Nw                = size(B,2);

if nargin < 8 || isempty(rndStream)
    rndStream = getDefaultStream;
end



%% init Variables and allocate memory
I                 = eye(Nx);

yDataNdx          = ~yNaNndx;

%% allocate memory
[Sigmattm1, ImKC]           = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, T);
[XtT, Xttm1, Xplus]         = deal(zeros(Nx, T));


%% generate plus data
wplus   = randn(rndStream, Nw, T);

%% Forward Loop: Kalman Forecasts
Sigmatt = zeros(Nx,Nx);
Xtt     = zeros(Nx,1); % mean of differences between actual and plus (and thus zeros)

disturbanceplus  = zeros(Nx, T);

for t = 1 : T
    
    % "plus" States and priors
    Bsv                     = B * diag(sqrtSigma(:,t));
    disturbanceplus(:,t)    = Bsv * wplus(:,t);
    BSigmaB                 = Bsv * Bsv';
    
    if t == 1
        Xlagplus = X0;
    else
        Xlagplus = Xplus(:,t-1);
    end
    
    Xplus(:,t)              = A * Xlagplus + disturbanceplus(:,t);
    
    % priors
    Sigmattm1(:,:,t)        = A * Sigmatt * A' + BSigmaB;
    Xttm1(:,t)              = A * Xtt;
    
    % observed innovation
    Yplus            = C(:,:,t) * Xplus(:,t);
    SigmaYttm1       = C(:,:,t) * Sigmattm1(:,:,t) * C(:,:,t)';
    
    
    Ytilde(:,t)      = Ydata(:,t) - Yplus - C(:,:,t) * Xttm1(:,t);
    
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(yDataNdx(:,t),yDataNdx(:,t),t) = eye(sum(yDataNdx(:,t))) / SigmaYttm1(yDataNdx(:,t),yDataNdx(:,t));
    
     
    % Kalman Gain
    K                       = (Sigmattm1(:,:,t) * C(:,:,t)') * invSigmaYttm1(:,:,t);
    ImKC(:,:,t)             = I - K * C(:,:,t);
    
    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t) * ImKC(:,:,t)'; % Joseph form for better numerical stability
    
    Xtt                     = Xttm1(:,t) + K * Ytilde(:,t);
    %     Xplustt                 = Xplusttm1(:,t) + K * Yplustilde(:,t);
    
end

%% Backward Loop: Disturbance Smoother
XtT(:,T)        = Xtt;

StT             = C(:,:,T)' * (invSigmaYttm1(:,:,T) * Ytilde(:,T));


for t = (T-1) : -1 : 1
    
    Atilde          = A * ImKC(:,:,t);
    
    StT             = Atilde' * StT + C(:,:,t)' * (invSigmaYttm1(:,:,t) * Ytilde(:,t));
    XtT(:,t)        = Xttm1(:,t) + Sigmattm1(:,:,t) * StT;
    
end

%% sample everything together (and reorder output dimensions)
% Xdraws  = bsxfun(@minus, XtT, XplustT) + Xplus;
Xdraws  = XtT + Xplus;
% Xdraws  = permute(Xdraws, [1 3 2]);



