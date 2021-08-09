function [h, h0, hshock, SV, SVtScalelog2Draws] = ...
    StochVoltdof(y2, logy2, h, hVCVsqrt, Eh0, sqrtVh0, ...
    SVtdof, ...
    KSC, KSCt, Nsv, T, rndStream)
% StochVolOutlierKSC combines KSC Gibbs Sampling for SV with outlier model of Stock-Watson (2016, REStat)
%
% Uses Kim, Shephard and Chib normal mixtures
%
% USAGE :[h, h0, hshock, SV, SVtScalelog2Draws, SVtdof] = ...
%     StochVoltdof(y2, logy2, h, hVCVsqrt, Eh0, sqrtVh0, ...
%     tdof, ...
%     KSC, KSCt, Nsv, T, rndStream)
%
% See also getKSC7values, getKSC10values, StochVolOutlierKSCcorrsqrt

%   Coded by  Elmar Mertens, em@elmarmertens.com


if isscalar(Eh0)
    Eh0 = repmat(Eh0, Nsv, 1);
end
if isscalar(sqrtVh0)
    sqrtVh0 = sqrtVh0 * eye(Nsv);
end

%% t-shocks via outlier draws (following JPR04)
y2scaled        = y2 .* exp(-h);

scalePosterior     = SVtdof + 1 + y2scaled; % note the explicit vector expansion of SVtdof
% note: matlab doc says stats box handles parallel streams automatically via the global stream ....
chi2draws          = chi2rnd(repmat(SVtdof + 1, 1, T));
SVtScalelog2Draws  = log(scalePosterior) - log(chi2draws);


%% draw mixture states
% zdraws are standardized draws for each component of the normal mixture 
% zdraws is thus Nsv x T x Nmixtures
% zdraws      = bsxfun(@minus, logy2 - h - outlierlog2Draws, KSCt.mean) ./ KSCt.vol;
zdraws      = (logy2 - h - SVtScalelog2Draws - KSCt.mean) ./ KSCt.vol;

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
obs   = logy2 - KSC.mean(kai2States) - SVtScalelog2Draws;
sqrtR = zeros(Nsv,Nsv,T);
for n = 1 : Nsv
    sqrtR(n,n,:) = KSC.vol(kai2States(n,:));
end

% note: for larger systems, smoothing sampler turns out to be more
% efficient than Carter-Kohn
[h, hshock, h0] = vectorRWsmoothingsampler1draw(obs, hVCVsqrt, sqrtR, Eh0, sqrtVh0, rndStream);


%% construct SV
SV = exp((h + SVtScalelog2Draws) / 2);

