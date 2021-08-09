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

datadir = '~/jam/lager/var2021-matfiles';

MODELTYPES  = {'SVOtmax20', 'SVOmax20', 'SV', 'SVt'};
MODELLABELS = {'SVO-t', 'SVO', 'SV', 'SV-t'};

doMedian = false;
wrap = [];
initwrap

%% grab const
matCONST = matfile(fullfile(datadir, sprintf('fredMD16-2021-04-censoredYields-%s-p12.mat', 'CONST')));
medianMaxVARrootCONST   = squeeze(median(matCONST.drawsMaxlambda,1));
tailsMaxVARrootCONST = prctile(matCONST.drawsMaxlambda, [5 95], 1);


%% loop over models
for thism =  1 : length(MODELTYPES)
    modeltype  = MODELTYPES{thism};
    modellabel = MODELLABELS{thism};
    
    mat   = matfile(fullfile(datadir, sprintf('fredMD16-2021-04-censoredYields-%s-p12.mat', modeltype)));
    
    
    %#ok<*UNRCH>
    
    %% pull some objects
    Tjumpoffs = mat.Tjumpoffs;
    ydates    = mat.ydates;
    N         = mat.N;
    p         = mat.p;
    K         = N*p + 1;
    ncode     = mat.ncode;
    Ylabels   = fredMDprettylabel(ncode);
    
    %% settings
    
    medianMaxVARroot   = squeeze(median(mat.drawsMaxlambda,1));
    tailsMaxVARroot = prctile(mat.drawsMaxlambda, [5 95], 1);
    
    ndxCOVID = ydates(Tjumpoffs) > datenum(2019,12,1);
    coviddates = ydates(Tjumpoffs(ndxCOVID));
    
    %% plot maxRoots
    thisfig = figure;
    subplot(2,1,1)
    hold on
    hsv = plot(ydates(Tjumpoffs), medianMaxVARroot, 'k-', 'linewidth', 2);
    plot(ydates(Tjumpoffs), tailsMaxVARroot', 'k--', 'linewidth', 1)
    
    hconst = plot(ydates(Tjumpoffs), medianMaxVARrootCONST, 'r:', 'linewidth', 2);
    plot(ydates(Tjumpoffs), tailsMaxVARrootCONST', 'r:', 'linewidth', 1)

    xtickdates(ydates(Tjumpoffs))
    title('full sample')
    legend([hsv hconst], modellabel, 'CONST', 'location', 'best')
    set(gca, 'fontsize', 12)
    
    subplot(2,1,2)
    hold on
    plot(coviddates, medianMaxVARroot(ndxCOVID), 'k-', 'linewidth', 2)
    plot(coviddates, tailsMaxVARroot(:,ndxCOVID)', 'k--', 'linewidth', 1)
    plot(coviddates, medianMaxVARrootCONST(ndxCOVID), 'r:', 'linewidth', 2);
    plot(coviddates, tailsMaxVARrootCONST(:,ndxCOVID)', 'r:', 'linewidth', 1)
    xlim(coviddates([1 end]))
    datetick('x', 'yyyy:mmm', 'keeplimits')
    title('since 2020')
    set(gca, 'fontsize', 12)
    
    sgtitle(modellabel)
    wrapthisfigure(thisfig, sprintf('MaxVARroot-%s', modeltype), wrap);
    
    
    
end

%% finish script
finishwrap
finishscript