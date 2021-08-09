function [obarState, obarProb] = obarGibbsdraw(resid, outlierProb, obarAlpha, obarBeta, obarStates, rndStream)
% OBARGIBBSDRAW ... 
%  
% - assumes resid is standardized (by SV and other vol factors except for obar)
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 11-Apr-2021 21:06:42 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.10.0.1602886 (R2021a) 
% FILENAME  : obarGibbsdraw.m 



[T, N] = size(resid); 

obarPrior = [1 - outlierProb, outlierProb * obarStates.uniformkernel];
% checkdiff(obarPrior, [1 - outlierProb, repmat(outlierProb / obarStates.Ngrid, 1, obarStates.Ngrid)]);
% checkdiff(sum(obarPrior) -1 );

%% obar posterior 
ssr = sum(resid.^2,2) ./ (obarStates.squaredvalues);

loglikekernel = -0.5 * (N * obarStates.log2values + ssr);
logkernel     = loglikekernel + log(obarPrior);
mlk           = max(logkernel, [], 2);
pdfKernel     = exp(logkernel - mlk);


cdf               = cumsum(pdfKernel, 2);                % integrate
cdf(:,1:end-1)    = cdf(:,1:end-1) ./ cdf(:,end); 
cdf(:,end)        = 1;    % normalize

cdfcheck = cdf;
pdfKernel = exp(-.5 * ssr) ./ (obarStates.values.^N) .* obarPrior;
cdf               = cumsum(pdfKernel, 2);                % integrate
cdf(:,1:end-1)    = cdf(:,1:end-1) ./ cdf(:,end); 
cdf(:,end)        = 1;    % normalize
checkdiff(cdf, cdfcheck);

% draw states
ndx               = sum(bsxfun(@gt, rand(rndStream, T, 1), cdf), 2) + 1;

obarState         = obarStates.values(ndx);

%% update outlierProb
Noutlier    = sum(ndx > 1, 1);
alpha       = obarAlpha + Noutlier;
beta        = obarBeta + (T - Noutlier);
obarProb    = betarnd(alpha, beta);