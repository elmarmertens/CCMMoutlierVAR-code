function [tscaleDraw, tdofDraw] = tposteriordraw(y2scaled, tdof, rndStream)
% TPOSTERIORDRAW ...
%
%   ...

% note: rndStream is currently ignored

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 03-Jul-2021 17:09:21 $
% $Revision : 1.00 $
% DEVELOPED : 9.10.0.1684407 (R2021a) Update 3
% FILENAME  : tposteriordraw.m

[N, T]        = size(y2scaled);

% y2scaled      = y2 ./ yvar;
% y2scaledGrid  = repmat(permute(y2scaled, [1 3 2]), [1 tdof.Ndof 1]);


% dataloglike   = - .5 * (tdof.values + 1) .* sum(log(tdof.values + y2scaledGrid), 3);
dataloglike   = - .5 * (tdof.values + 1) .* sum(log(tdof.values + permute(y2scaled, [1 3 2])), 3);
loglike       = tdof.loglike0 + dataloglike;

logposteriorKernel = tdof.logprior + loglike;
% note: adding prior could be dropped from kernel as long as the prior is uniform
% subtract const to avoid overflow
logposteriorKernelstar  = logposteriorKernel  - max(logposteriorKernel, [], 2);


cdf            = cumsum(exp(logposteriorKernelstar), 2);
cdf(:,1:end-1) = cdf(:,1:end-1) ./ cdf(:,end);
cdf(:,end)     = 1;
dofStates      = sum(rand(rndStream, N, 1) > cdf, 2) + 1;
tdofDraw       = tdof.values(dofStates)'; % transpose!


% draw tscales
scalePosterior  = tdofDraw + 1 + y2scaled; % note the explicit vector expansion of SVtdof
chi2draws       = chi2rnd(repmat(tdofDraw + 1, 1, T));
tscaleDraw      = scalePosterior ./ chi2draws; % scale factor for variances