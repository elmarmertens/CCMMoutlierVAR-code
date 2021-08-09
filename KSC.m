%% Draws mixture states and then volatilities using KSC algorithm
function [h, h0] = KSC(v_tilda,h,Qt,states_pmean,states_pvar, rndStream)
% Based on modifications of code by J. Chan / D. Korobilis
%
% further amendments by Elmar Mertens:
% - added rndStream argument
% - h0 output


if nargin < 6 || isempty(rndStream)
    rndStream = getDefaultStream;
end

% pointers
[T, N]= size(v_tilda);

% normal mixture moments
pi =   [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
mi =   [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;
si =   [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];

% 1. draw mixture states from a 7-point discrete distribution
S=zeros(T,N);
for i=1:N
    q = repmat(pi,T,1).*normpdf(repmat(v_tilda(:,i),1,7),repmat(h(:,i),1,7)+repmat(mi,T,1), repmat(sqrt(si),T,1));
    q = q./repmat(sum(q,2),1,7);
    S(:,i) = 7 - sum(repmat(rand(rndStream,T,1),1,7)<cumsum(q,2),2) +1 ;
end

% 2. draw volatilities conditional on mixture states, using CK (1994)
[h, h0] = CK(v_tilda-mi(S),si(S),Qt,N,T,states_pmean,states_pvar,rndStream);


%% Function for Carter and Kohn (1994) smoother
function [St_draw, S0_draw] = CK(y,Ht,Qt,N,T,S0,P0,rndStream)

% transpose input
y=y';

% prepare matrices for storage
St_collect = zeros(N,T);  Pt_collect = zeros(N,N,T); St_draw = zeros(N,T);

% forward recursions (Kalman filter)
St = S0; Pt = P0;
for t=1:T
    % prediction using transition equation (prior)
    St_1 = St; Pt_1 = Pt + Qt;
    % use observation equation to compute forecast error and Kalman gain
    vt = y(:,t) - St_1;              % conditional forecast error
    Varvt = Pt_1 + diag(Ht(t,:));    % variance of forecast error
    Kt=Pt_1/Varvt;                   % Kalman Gain
    % update
    St = St_1 + Kt*vt; Pt = Pt_1 - Kt*Pt_1;
    % store
    St_collect(:,t) = St; Pt_collect(:,:,t) = Pt;
end

% Backward simulation

% for better performance: draw all random variables at once
zdraws = randn(rndStream,N,T+1); % T+1 to include also time zero

St_draw(:,T) =St_collect(:,T) + chol(Pt_collect(:,:,T),'lower') * zdraws(:,T);
for t=T-1:-1:1
    % compute moments
    Kt    = Pt_collect(:,:,t)/(Pt_collect(:,:,t) + Qt);
    Smean = St_collect(:,t)   + Kt*(St_draw(:,t+1) - St_collect(:,t));
    Svar  = Pt_collect(:,:,t) - Kt*Pt_collect(:,:,t);
    % draw and store
    St_draw(:,t) = Smean + chol(Svar,'lower') * zdraws(:,t);
end

% add time 0 smoothed state
Kt    = P0 / (P0 + Qt);
Smean = S0   + Kt * (St_draw(:,1) - S0);
Svar  = P0 - Kt * P0;
S0_draw = (Smean + chol(Svar,'lower') * zdraws(:,T+1))';

% transpose output
St_draw = St_draw';

%% Function for Carter and Kohn (1994) smoother -- implemented via QR factorizations (EM)
% not used, since the QR is slower than the regular algorithm
% use only if the choleski's in the regular CK function fail
% function [St_draw, S0_draw] = CKqr(y,Ht,Qt,N,T,S0,P0,rndStream)
% 
% % transpose input
% y=y';
% sqrtP0 =  chol(P0)';
% 
% % prepare matrices for storage
% St_collect = zeros(N,T);  sqrtPt_collect = zeros(N,N,T); St_draw = zeros(N,T);
% 
% % forward recursions (Kalman filter)
% St = S0; sqrtPt = sqrtP0; sqrtQt = chol(Qt)';
% for t=1:T
%     % prediction using transition equation (prior)
%     St_1 = St; sqrtH = diag(sqrt(Ht(t,:)));
%     
%     % use observation equation to compute forecast error and Kalman gain
%     vt = y(:,t) - St_1;              % conditional forecast error
%     
%     % setup qr for Gain and Posterior variance
%     M                = zeros(N,N+N+N);
%     M(1:N,:)         = [sqrtPt, sqrtQt, sqrtH];
%     M(N+(1:N),1:N+N) = [sqrtPt, sqrtQt];
%     [~,R] = qr(M',0);
%     R = R';
%     
%     sqrtVarvt = R(1:N, 1:N);
%     Kt        = R(N+(1:N), 1:N) / sqrtVarvt;
%     sqrtPt    = R(N+(1:N), N+(1:N));
%     
%     % update
%     St = St_1 + Kt*vt;
%     % store
%     St_collect(:,t) = St; sqrtPt_collect(:,:,t) = sqrtPt;
% end
% 
% % Backward simulation
% 
% % for better performance: draw all random variables at once
% zdraws = randn(rndStream,N,T+1); % T+1 to include also time zero
% 
% St_draw(:,T) = St_collect(:,T) + sqrtPt_collect(:,:,T) * zdraws(:,T);
% for t=T-1:-1:1
%     % compute moments
%     
%     M                = zeros(N,N+N);
%     M(1:N,:)         = [sqrtPt_collect(:,:,t), sqrtQt];
%     M(N+(1:N),1:N)   = sqrtPt_collect(:,:,t);
%     [~,R] = qr(M',0);
%     R = R';
%     
%     Kt        = R(N+(1:N), 1:N) / R(1:N, 1:N);
%     sqrtPttp1 = R(N+(1:N), N+(1:N));
%     
%     Smean = St_collect(:,t)   + Kt*(St_draw(:,t+1) - St_collect(:,t));
%     % draw and store
%     St_draw(:,t) = Smean + sqrtPttp1 * zdraws(:,t);
% end
% 
% % add time 0 smoothed state
% M                = zeros(N,N+N);
% M(1:N,:)         = [sqrtP0, sqrtQt];
% M(N+(1:N),1:N)   = sqrtP0;
% [~,R] = qr(M',0);
% R = R';
% 
% Kt        = R(N+(1:N), 1:N) / R(1:N, 1:N);
% sqrtPttp1 = R(N+(1:N), N+(1:N));
% 
% Smean = St_collect(:,t)   + Kt*(St_draw(:,t+1) - St_collect(:,t));
% 
% S0_draw = (Smean + sqrtPttp1 * zdraws(:,T+1))';
% 
% % transpose output
% St_draw = St_draw';