function [h, hbar, htilde, hresid, SV, outlierlog2Draws, outlierProb, outlierScaleDraws] = ...
    StochVolOutlierAR1(logy2, h, rho, hVCVsqrt, Eh0, sqrtVh0, ...
    outlierlog2Draws, outlierProb, outlieralpha, outlierbeta, outlierStates, ...
    KSC, KSCt, Nsv, T, rndStream)
% StochVolOutlierKSCcorrsqrt combines KSC Gibbs Sampling for SV with outlier model of Stock-Watson (2016, REStat)
%
% Uses Kim, Shephard and Chib normal mixtures
%
% USAGE : [h, h0, hshock, SV, outlierlog2Draws, outlierProb, outlierScaleDraws] = StochVolOutlierKSCcorrsqrt(logy2, h, rho, hVCVsqrt, Eh0, sqrtVh0, outlierlog2Draws, ...
%         outlierProb, outlieralpha, outlierbeta, outlierStates, KSC, KSCt, ...
%         Nsv, T, rndStream)
%
%
% See also getKSC7values, getKSC10values, StochVolOutlierKSC

%   Coded by  Elmar Mertens, em@elmarmertens.com


if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(sqrtVh0)
    sqrtVh0 = sqrtVh0 * eye(Nsv);
end

%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x Nmixtures
% zdraws      = bsxfun(@minus, logy2 - h - outlierlog2Draws, KSCt.mean) ./ KSCt.vol;
zdraws      = (logy2 - h - outlierlog2Draws - KSCt.mean) ./ KSCt.vol;

% construct CDF
% factor of sqrt(2 * pi) can be ommitted for kernel
pdfKernel           = KSCt.pdf ./ KSCt.vol .* exp(-.5 * zdraws.^2); 
cdf                 = cumsum(pdfKernel, 3);                % integrate
% cdf(:,:,1:end-1)    = bsxfun(@rdivide, cdf(:,:,1:end-1), cdf(:,:, end)); 
cdf(:,:,1:end-1)    = cdf(:,:,1:end-1) ./ cdf(:,:, end); % using automatic expansion 
cdf(:,:,end)        = 1;    % normalize

% draw states
% kai2States  = sum(bsxfun(@gt, rand(rndStream, Nsv, T), cdf), 3) + 1;
kai2States  = sum(rand(rndStream, Nsv, T) > cdf, 3) + 1;


%% KSC State Space
obs       = logy2 - KSC.mean(kai2States);
zerosNsv  = zeros(Nsv);
Insv      = eye(Nsv);
A     = [diag(rho) zerosNsv; zerosNsv Insv];
B     = [hVCVsqrt; zerosNsv];
C     = [Insv Insv];
sqrtR = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

sqrtVhtilde  = zeros(Nsv); % Note: need fixed prior, not depended on estimated rhos (alt: use prior rho)
x0           = [zeros(Nsv, 1); Eh0];
sqrtVx0      = [sqrtVhtilde, zerosNsv; zerosNsv sqrtVh0];
[H, Hshock, H0] = a2b2c2DisturbanceSmoothingSampler1draw(A, B, C, obs, x0, sqrtVx0, ...
    sqrtR, rndStream); 

h      = H(1:Nsv,:) + H(Nsv+1:end,:); % C * H
hbar   = H0(Nsv+1:end);
htilde = cat(2, H0(1:Nsv,:), H(1:Nsv,:));
hresid = Hshock(1:Nsv,:);

%% outlier PDF
% outlierPdf is Nsurvey times T  times Nstates
% outlierPdf2 = cat(3, repmat(1 - outlierProb, 1, T), bsxfun(@times, outlierProb, repmat(1 / outlierNgrid, 1, T, outlierNgrid)));
outlierPdf  = cat(3, repmat(1 - outlierProb, 1, T), repmat(outlierProb / outlierStates.Ngrid, 1, T, outlierStates.Ngrid));

%% outlier states
edraws      = bsxfun(@minus, logy2 - h - KSC.mean(kai2States), permute(outlierStates.log2values, [1 3 2]));
zdraws      = bsxfun(@rdivide, edraws, KSC.vol(kai2States));

% pdfKernel   = exp(-.5 * zdraws.^2);  
% division by KSC.vol is unnecessary for this kernel, since same vol would apply across outlierStates

pdfKernel   = outlierPdf .* exp(-.5 * zdraws.^2);

cdf                 = cumsum(pdfKernel, 3);                % integrate
cdf(:,:,1:end-1)    = bsxfun(@rdivide, cdf(:,:,1:end-1), cdf(:,:,end)); 
cdf(:,:,end)        = 1;    % normalize


% draw states
ndx               = sum(bsxfun(@gt, rand(rndStream, Nsv, T), cdf), 3) + 1;
outlierlog2Draws  = outlierStates.log2values(ndx);
outlierScaleDraws = outlierStates.values(ndx);

%% update outlierProb
Noutlier    = sum(ndx > 1, 2);
alpha       = outlieralpha + Noutlier;
beta        = outlierbeta + (T - Noutlier);
for n = 1 : Nsv
    outlierProb(n) = betadraw(alpha(n), beta(n), 1, rndStream);
    % re matlab's betarnd:
    % - does not seem to support randomStreams
    % - seems to be slower
    % outlierProb(n) = betarnd(alpha(n), beta(n), 1);
end

%% construct SV
SV = exp((h + outlierlog2Draws) / 2);
