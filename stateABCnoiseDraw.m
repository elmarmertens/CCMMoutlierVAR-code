function [Xdraw, noiseDraw, disturbanceDraw] = stateABCnoiseDraw(A, B, C, Ydata, X0, sqrtSigma, sqrtR, rndStream)
% STATEABCDRAW
% ....
% supposes X0 is deterministically given

%   Coded by  Elmar Mertens, em@elmarmertens.com


%% parse inputs
Nx                = size(A, 1);
[Ny, T]           = size(Ydata);
Nw                = size(B,2);

if nargin < 8
    sqrtR = [];
end

if nargin < 9 || isempty(rndStream)
    rndStream = getDefaultStream;
end



%% init Variables and allocate memory
if ~isempty(sqrtR) && ismatrix(sqrtR)
    sqrtR = repmat(sqrtR, [1 1 T]);
end
Ix  = eye(Nx);
Iy  = eye(Ny);

%% allocate memory
[Sigmattm1, ImKC, BSB]      = deal(zeros(Nx, Nx, T));
invSigmaYttm1               = zeros(Ny, Ny, T);
Ytilde                      = zeros(Ny, T);
[XtT, Xttm1, Xplus]         = deal(zeros(Nx, T));


%% generate plus data
wplus   = randn(rndStream, Nw, T);
if ~isempty(sqrtR)
    Nnoise  = size(sqrtR, 2);
    eplus   = randn(rndStream, Nnoise, T);
end

%% Forward Loop: Kalman Forecasts
Sigmatt = zeros(Nx,Nx);
Xtt     = zeros(Nx,1); % mean of differences between actual and plus (and thus zeros)

disturbanceplus  = zeros(Nx, T);
if isempty(sqrtR)
    noiseplus        = [];
else
    noiseplus        = zeros(Ny, T);
end

for t = 1 : T
    
    % "plus" States and priors
    Bsv                     = B * diag(sqrtSigma(:,t));
    disturbanceplus(:,t)    = Bsv * wplus(:,t);
    BSB(:,:,t)              = Bsv * Bsv';
    
    if t == 1
        Xlagplus = X0;
    else
        Xlagplus = Xplus(:,t-1);
    end
    
    Xplus(:,t)              = A * Xlagplus + disturbanceplus(:,t);
    
    % priors
    Sigmattm1(:,:,t)        = A * Sigmatt * A' + BSB(:,:,t);
    Xttm1(:,t)              = A * Xtt;
    %     Xplusttm1(:,t)          = A * Xplustt;
    
    
    % observed innovation
    Yplus            = C * Xplus(:,t);
    SigmaYttm1       = C * Sigmattm1(:,:,t) * C';
    
    if ~isempty(sqrtR)
        noiseplus(:,t)    = sqrtR(:,:,t) * eplus(:,t); %#ok<AGROW>
        Yplus             = Yplus + noiseplus(:,t);
        SigmaYttm1        = SigmaYttm1 + sqrtR(:,:,t) * sqrtR(:,:,t)';
    end
    
    Ytilde(:,t)      = Ydata(:,t) - Yplus - C * Xttm1(:,t);
    
    % Block Inverse of Y-VCV, accounting for missing obs (with zero VCV)
    invSigmaYttm1(:,:,t) = Iy / SigmaYttm1;
    
    % Kalman Gain
    K                       = (Sigmattm1(:,:,t) * C') * invSigmaYttm1(:,:,t);
    ImKC(:,:,t)             = Ix - K * C;
    
    % posteriors
    Sigmatt                 = ImKC(:,:,t) * Sigmattm1(:,:,t) * ImKC(:,:,t)'; % Joseph form for better numerical stability
    
    Xtt                     = Xttm1(:,t) + K * Ytilde(:,t);
    
end

%% Backward Loop: Disturbance Smoother
XtT(:,T)        = Xtt;

StT             = C' * (invSigmaYttm1(:,:,T) * Ytilde(:,T));


if nargout > 2
    disturbancetT        = zeros(Nx, T);
    disturbancetT(:,T)   = BSB(:,:,T) * StT;
else
    disturbancetT        = [];
end
    

for t = (T-1) : -1 : 1
    
    Atilde          = A * ImKC(:,:,t);
    
    StT             = Atilde' * StT + C' * (invSigmaYttm1(:,:,t) * Ytilde(:,t));
    XtT(:,t)        = Xttm1(:,t) + Sigmattm1(:,:,t) * StT;
 
    if ~isempty(disturbancetT)
        disturbancetT(:,t)   = BSB(:,:,t) * StT;
    end
    
end


% sample everything together (and reorder output dimensions)
Xdraw  = XtT + Xplus;


%% extra outputs

if nargout > 1 && ~isempty(sqrtR)
    noisetT      = zeros(Ny, T);
    for t = T : -1 : 1
        noisetT(:,t) = Ytilde(:,t) - C * (XtT(:,t) - Xttm1(:,t));
    end
    noiseDraw  = noiseplus + noisetT;
else
    noiseDraw  = [];
end

if nargout > 2
   disturbanceDraw = disturbanceplus + disturbancetT;
end


       
