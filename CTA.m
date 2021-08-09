function PAI=CTA(Y,X,N,K,~,A_,sqrt_ht,iV,iVb_prior,PAI,rndStream)

% =========================================================================
% Performs a draw from the conditional posterior of the VAR conditional
% mean coefficients by using the triangular algorithm. The
% triangularization achieves computation gains of order N^2 where N is the
% number of variables in the VAR. Carriero, Clark and Marcellino (2015),
% Large Vector Autoregressions with stochastic volatility and flexible
% priors. 
%
% The model is:
%
%     Y(t) = Pai(L)Y(t-1) + v(t); Y(t) is Nx1; t=1,...,T, L=1,...,p.
%     v(t) = inv(A)*(LAMBDA(t)^0.5)*e(t); e(t) ~ N(0,I);
%                _                                         _
%               |    1          0       0       ...      0  |
%               |  A(2,1)       1       0       ...      0  |
%      A =      |  A(3,1)     A(3,2)    1       ...      0  |
%               |   ...        ...     ...      ...     ... |
%               |_ A(N,1)      ...     ...   A(N,N-1)    1 _|
%
%     Lambda(t)^0.5 = diag[sqrt_h(1,t)  , .... , sqrt_h(N,t)];
%
% INPUTS
% Data and pointers:
% Y     = (TxN) matrix of data appearing on the LHS of the VAR
% X     = (TxK) matrix of data appearing on the RHS of the VAR
% N     = scalar, #of variables in VAR 
% K     = scalar, #of regressors (=N*p+1)  
% T     = scalar, #of observations
% The matrix X needs to be ordered as: [1, y(t-1), y(t-2),..., y(t-p)]
% 
% Error variance stuff:
% invA_   = (NxN) inverse of lower triangular covariance matrix A
% sqrt_ht = (TxN) time series of diagonal elements of volatility matrix 
% For a homosckedastic system, with Sigma the error variance, one can
% perform the LDL decomposition (command [L,D]=LDL(Sigma)) and set inv_A=L
% and sqrt_ht=repmat(sqrt(diag(D)'),T,1). 
%
% Priors:
% iV          = (NKxNK) precision matrix for VAR coefficients 
% iVB_prior   = (NKx1) (prior precision)*(prior mean)
% Note 1:iV is the inverse of the prior matrix and iVB_prior is the product
% of iV and the prior mean vector, which both need to be computed only once,
% before the start of the main MCMC loop.
% Note 2:in this code, iV is assumed block-diagonal. This implies that the
% prior is independent across equations. This includes most of the priors 
% usually considered, including the Minnesota one.  To use a non-block
% diagonal iV one needs to modify the code using the recursions illustrated 
% in equations (37) and (38).  
%
% OUTPUT
% One draw from (PAI|A,Lambda,data)
% PAI=[Pai(0), Pai(1), ..., Pai(p)].
% =========================================================================

% y_til=Y*A_';
zdraws = randn(rndStream,K,N);
for j=1:N
    
    % select coefficients of equation j to remove from the LHS
    PAI(:,j)=zeros(K,1); 
    
    % build model
    lambda=vec(sqrt_ht(:,j:N));
    % Y_j=vec(y_til(:,j:N)-X*PAI*A_(j:N,:)')./lambda;
    Y_j=vec((Y-X*PAI)*A_(j:N,:)')./lambda; % quicker than the preceding line
    
    X_j=kron(A_(j:N,j),X)./lambda;
    
    % posterior moments
    index=K*(j-1)+1:(K*(j-1)+K);
    iV_post = iV(index,index)  + X_j'*X_j;
    %     V_post = iV_post\eye(K);
    %     b_post = V_post*(iVb_prior(index) + X_j'*Y_j);
    [iVchol_post, ispd]  = chol(iV_post, 'lower');
    if ispd == 0
        Vchol_post      = (iVchol_post \ eye(K))'; % note upper right factor
        V_post          = Vchol_post * Vchol_post'; % V = inv(L') * inv(L)
    else % proceed via QR
        
        warning('switching to QR routine')
        
        iVchol = chol(iV(index,index))'; % todo: could do this once and for all, not for every call of triang
        % set up Kailath's fast array matrix
        qrM = [iVchol, X_j'];
        R = triu(qr(qrM',0))';
    
        iVchol_post = R(1:K,1:K);
        Vchol_post  = (iVchol_post \ eye(K))'; % upper triangular but OK
        V_post      = Vchol_post * Vchol_post';
    end
    
    % posterior draw
    b_post   = V_post*(iVb_prior(index) + X_j'*Y_j);
    PAI(:,j) = b_post + Vchol_post * zdraws(:,j);
    
end
                    
  
