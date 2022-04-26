%% load em toolboxes
warning('off','MATLAB:handle_graphics:exceptions:SceneNode')

path(pathdef)

addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/
addpath matlabtoolbox/emgibbsbox/
addpath matlabtoolbox/emeconometrics/
addpath matlabtoolbox/emstatespace/

%% Initial operations
clear; close all; clc;

resultsdir = pwd;

jumpoffDate = datenum(2021,3,1);

jumpoff     = datestr(jumpoffDate, 'yyyymm');

this = 'SVO';
matSVO        = matfile(fullfile(resultsdir, ...
    sprintf('fredMD16-2021-04-%s-p12-%s-alldraws.mat', jumpoff, this)));
this = 'SVt';
matSVt = matfile(fullfile(resultsdir, ...
    sprintf('fredMD16-2021-04-%s-p12-%s-alldraws.mat', jumpoff, this)));
this = 'SVOt';
matSVOt = matfile(fullfile(resultsdir, ...
    sprintf('fredMD16-2021-04-%s-p12-%s-alldraws.mat', jumpoff, this)));

wrap = [];
initwrap


%#ok<*UNRCH>

%% pull some objects
ydates = matSVO.ydates;
ncode  = matSVO.ncode;
N      = length(ncode);

% Ylabels = fredMDshortlabel(ncode);
Ylabels = fredMDprettylabel(ncode);



%% collect values for table
thetaMid       = NaN(N, 4);
thetaLowerTail = NaN(N, 4);
thetaUpperTail = NaN(N, 4);

lowerTail = normcdf(-1) * 100;
upperTail = normcdf(1) * 100;

% SVO
thetaMid(:,1)       = median(matSVO.SVOprob_all * 100,2);
thetaLowerTail(:,1) = prctile(matSVO.SVOprob_all * 100, lowerTail, 2);
thetaUpperTail(:,1) = prctile(matSVO.SVOprob_all * 100, upperTail, 2);

% SVt
thetaMid(:,2)       = round(median(matSVt.SVtdof_all,2));
thetaLowerTail(:,2) = round(prctile(matSVt.SVtdof_all, lowerTail, 2));
thetaUpperTail(:,2) = round(prctile(matSVt.SVtdof_all, upperTail, 2));

% SVOt
thetaMid(:,3)       = median(matSVOt.SVOprob_all * 100,2);
thetaLowerTail(:,3) = prctile(matSVOt.SVOprob_all * 100, lowerTail, 2);
thetaUpperTail(:,3) = prctile(matSVOt.SVOprob_all * 100, upperTail, 2);
thetaMid(:,4)       = round(median(matSVOt.SVtdof_all,2));
thetaLowerTail(:,4) = round(prctile(matSVOt.SVtdof_all, lowerTail, 2));
thetaUpperTail(:,4) = round(prctile(matSVOt.SVtdof_all, upperTail, 2));

%% prior
np = 12;
Ndraws = 1e4;

% SVO
SVOpriorobs = 10 * np;
SVOalpha    = 1 / (4 * np) * SVOpriorobs; % 10 years of data with 1 outlier every 4 years
SVObeta     = SVOpriorobs - SVOalpha;
priordraws  = 100 * betadraw(SVOalpha, SVObeta, Ndraws);
  
midprobSVO   = median(priordraws);
lowerprobSVO = prctile(priordraws, lowerTail);
upperprobSVO = prctile(priordraws, upperTail);


% SVOt
SVOtpriorobs = 10 * np;
SVOtalpha    = 1 / (10 * np) * SVOtpriorobs; % 10 years of data with 1 outlier every 10 years
SVOtbeta     = SVOtpriorobs - SVOtalpha;
priordraws   = 100 * betadraw(SVOtalpha, SVOtbeta, Ndraws);

midprobSVOt   = median(priordraws);
lowerprobSVOt = prctile(priordraws, lowerTail);
upperprobSVOt = prctile(priordraws, upperTail);

% tprior
tdofGrid = 3 : 40;
ndxdraws = randi(length(tdofGrid), Ndraws, 1);
priordraws = tdofGrid(ndxdraws);

midtDof   = median(priordraws);
lowertDof = prctile(priordraws, lowerTail);
uppertDof = prctile(priordraws, upperTail);

%% produce table
tabname    = 'OutlierAdjustedSVparameters.tex';
tabcaption = sprintf('Outlier-adjusted SV Parameters');

if isempty(wrap)
    tabdir = pwd;
else
    tabdir = wrap.dir;
    latexwrapper(wrap, 'add', 'tab', tabname, tabcaption)
end

% write table
fid = fopen(fullfile(tabdir, tabname), 'wt');
% fprintf(fid, '\\doublespacing\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('.3', 1, 4));
fprintf(fid, '\\toprule\n');
fprintf(fid, '& \\multicolumn{1}{c}{SVO} ');
fprintf(fid, '& \\multicolumn{1}{c}{SV-t} ');
fprintf(fid, '& \\multicolumn{2}{c}{SVO-t} ');
fprintf(fid, '\\\\');
fprintf(fid, '\\cmidrule(r){2-2}');
fprintf(fid, '\\cmidrule(lr){3-3}');
fprintf(fid, '\\cmidrule(l){4-5}');
fprintf(fid, '\n');
fprintf(fid, 'Variable ');
fprintf(fid, '& \\multicolumn{1}{c}{$p_j$} ');
fprintf(fid, '& \\multicolumn{1}{c}{$\\nu_j$} ');
fprintf(fid, '& \\multicolumn{1}{c}{$p_j$} ');
fprintf(fid, '& \\multicolumn{1}{c}{$\\nu_j$} ');
fprintf(fid, '\\\\\n');
% fprintf(fid, '\\midrule\n');
% fprintf(fid, '\\multicolumn{%d}{c}{Prior moments}', 5);
% fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
fprintf(fid, '%s ', '{Prior}');
fprintf(fid, '& %5.2f ', midprobSVO);
fprintf(fid, '& \\multicolumn{1}{c}{%d} ', midtDof);
fprintf(fid, '& %5.2f ', midprobSVOt);
fprintf(fid, '& \\multicolumn{1}{c}{%d} ', midtDof);
fprintf(fid, '\\\\\n');
fprintf(fid, '& \\multicolumn{1}{c}{\\it %5.2f -- %5.2f} ', lowerprobSVO, upperprobSVO);
fprintf(fid, '& \\multicolumn{1}{c}{\\it %5.2f -- %5.2f} ', lowertDof, uppertDof);
fprintf(fid, '& \\multicolumn{1}{c}{\\it %5.2f -- %5.2f} ', lowerprobSVOt, upperprobSVOt);
fprintf(fid, '& \\multicolumn{1}{c}{\\it %5.2f -- %5.2f} ', lowertDof, uppertDof);

fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');

% fprintf(fid, '\\multicolumn{%d}{c}{Posterior moments}\\\\\n', 5);
% fprintf(fid, '\\midrule\n');
for n = 1 : N
    fprintf(fid, '%s ', Ylabels{n});
    
    fprintf(fid, '& %5.2f ', thetaMid(n,1));
    fprintf(fid, '& \\multicolumn{1}{c}{%d} ', thetaMid(n,2));
    fprintf(fid, '& %5.2f ', thetaMid(n,3));
    fprintf(fid, '& \\multicolumn{1}{c}{%d} ', thetaMid(n,4));
    fprintf(fid, '\\\\\n');
    fprintf(fid, '& \\multicolumn{1}{c}{\\it %5.2f -- %5.2f} ', thetaLowerTail(n,1), thetaUpperTail(n,1));
    fprintf(fid, '& \\multicolumn{1}{c}{\\it %d -- %d} ', thetaLowerTail(n,2), thetaUpperTail(n,2));
    fprintf(fid, '& \\multicolumn{1}{c}{\\it %5.2f -- %5.2f} ', thetaLowerTail(n,3), thetaUpperTail(n,3));
    fprintf(fid, '& \\multicolumn{1}{c}{\\it %d -- %d} ', thetaLowerTail(n,4), thetaUpperTail(n,4));
    
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, 'Notes: Outlier probability $p_j$ (in percent), and $t$-distribution''s degrees of freedom, $\\nu_j$, of orthogonalized residuals of each variable.\n');
fprintf(fid, 'Median, 15.87\\%% and 84.14\\%% posterior quantiles from full-sample estimation using data from %s -- %s.\n', ...
    datestryymm(ydates(1)), datestryymm(jumpoffDate));
fprintf(fid, 'The priors are, $p_j \\sim Beta(%4.1f,%4.1f)$ for SVO, and $p_j \\sim Beta(%4.1f,%4.1f)$ forSVO-t.\n', SVOalpha, SVObeta, SVOtalpha, SVOtbeta);
fprintf(fid, 'and $\\nu_j \\sim U(%d,%d)$ direcretized over an integer-valued grid.\n', tdofGrid(1), tdofGrid(end));


% fprintf(fid, 'Significance of differences relative to model %s assessed by Diebold-Mariano test using Newey-West standard errors with $2$ lags.\n', ...
%     models(end).prettylabel);
fclose(fid);
type(fullfile(tabdir, tabname))



%% finish
dockAllFigures
finishwrap
finishscript