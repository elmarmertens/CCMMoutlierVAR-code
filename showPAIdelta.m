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

datadir = 'matPAI2020'; % '~/jam/lager/var2021-matfiles';

MODELTYPES = {'CONST-2020', 'SVOtmax20-2020',  ...
    'SVt-2020', 'SVOmax20-2020', 'SV-2020', 'SVnanO5-2020'};

doMedian = false;
doTitle  = false;

doLags        = false;
doEquations   = false;

maxIntercept  = 7;
maxLag1       = 30;
maxLagOTHER   = 5;

for thism =  1 : length(MODELTYPES)
    
    modeltype = MODELTYPES{thism};
    
    mat   = matfile(fullfile(datadir, sprintf('fredMD16-2021-04-censoredYields-%s-p12.mat', modeltype)));
    
    wrap = [];
    titlename = sprintf('PAIsinceCOVID-%s', modeltype);
    if doMedian
        titlename = strcat(titlename, '-median');
    end
    
    initwrap
    
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
    jumpoff0 = datenum(2020,1,1);
    T0       = find(ydates == jumpoff0);
    ndxT0    = find(Tjumpoffs == T0); % find to make matfile calls work
    
    %% collect PAI
    
    if ~doMedian
        PAI0   = mat.PAImean(:,:,ndxT0);
        PAI0se = mat.PAIstdev(:,:,ndxT0);
        PAIdevs = (mat.PAImean(:,:,ndxT0+1:end) - PAI0) ./ PAI0se;
    else
        PAI0   = mat.PAImedian(:,:,ndxT0);
        PAI0se = mat.PAIstdev(:,:,ndxT0);
        PAIdevs = (mat.PAImedian(:,:,ndxT0+1:end) - PAI0) ./ PAI0se;
    end
    
    maxRED = 2.5;
    
    %% plot PAI
    
    maxZ   = ceil(max(abs(PAIdevs),[], 'all') * 10) / 10;
    PAIdevs2D = reshape(PAIdevs(1:K,:,:), N * K, size(PAIdevs,3));
    thisfig = figure;
    
    hb = surf(ydates(T0+1:end),1:N*K,abs(PAIdevs2D)); %#ok<NASGU>
    
    zlim([0 maxZ])
    caxis([0 maxRED])
    
    datetick('x')
    xlabel('end of estimation window')
    ylabel('parameters')
    
    % color by height
    shading interp
    colorbar
    colormap turbo
    set(gca, 'fontsize', 12)
    wrapthisfigure(thisfig, sprintf('PAIall-%s', modeltype), wrap);
    
    %% plot PAI per equation
    if doEquations
        for n = 1 : N
            
            maxZ   = ceil(max(abs(PAIdevs),[], 'all') * 10) / 10;
            
            thisfig = figure;
            
            hb = surf(ydates(T0+1:end),1:K,squeeze(abs(PAIdevs(1:K,n,:)))); %#ok<NASGU>
            
            zlim([0 maxZ])
            caxis([0 maxRED])
            datetick('x')
            xlabel('end of estimation window')
            ylabel('parameters')
            
            
            
            % color by height
            shading interp
            colorbar
            colormap turbo
            
            wrapthisfigure(thisfig, sprintf('PAI-%s-%s', ncode{n}, modeltype), wrap, [], [], [], [], true);
            if doTitle
                title(sprintf('%s equation', Ylabels{n}))
                wrapthisfigure(thisfig, sprintf('PAI-%s-%s-WITHTITLE', ncode{n}, modeltype), wrap);
            end
        end
    end
    %% plot intercepts per equation
    thisfig = figure;
    
    thesedevs = squeeze(abs(PAIdevs(1,:,:)));
    
    % find variables with largest devs
    maxdev = max(thesedevs, [], 2);
    [maxdev,ndx] = sort(maxdev, 'desc');
    thesendx = sort(ndx(maxdev > 1.5));
    
    hb = surf(ydates(T0+1:end),1:N,thesedevs);
    
    zlim([0 maxIntercept])
    caxis([0 maxRED])
    
    datetick('x')
    yticks(thesendx)
    yticklabels(Ylabels(thesendx));
    
    xlabel('end of estimation window')
    
    
    % color by height
    shading interp
    colorbar
    colormap turbo
    wrapthisfigure(thisfig, sprintf('PAI-intercept-%s', modeltype), wrap, [], [], [], [], true);
    if doTitle
        title(sprintf('intercepts'))
        wrapthisfigure(thisfig, sprintf('PAI-intercept-%s-WITHTITLE', modeltype), wrap);
    end
    
    %% plot lag1 per equation
    thisfig = figure;
    
    thesedevs = reshape(abs(PAIdevs(1+(1:N),:,:)), N * N, []);
    
    % find variables with largest devs
    %     maxdev = max(thesedevs, [], 2);
    %     [maxdev,ndx] = sort(maxdev, 'desc');
    %     thesendx = sort(ndx(maxdev > 1.5));
    
    hb = surf(ydates(T0+1:end),1:N*N,thesedevs);
    
    zlim([0 maxLag1])
    caxis([0 maxRED])
    
    datetick('x')
    yticks(1:N:N*N)
    ylim([1 N*N])
    yticklabels(Ylabels);
    
    xlabel('end of estimation window')
    
    
    % color by height
    shading interp
    colorbar
    colormap turbo
    wrapthisfigure(thisfig, sprintf('PAI-lag1-%s', modeltype), wrap, [], [], [], [], true);
    if doTitle
        title(sprintf('intercepts'))
        wrapthisfigure(thisfig, sprintf('PAI-lag1-%s-WITHTITLE', modeltype), wrap);
    end
    
    %% plot other lags per equation
    thisfig = figure;
    
    thesedevs = reshape(abs(PAIdevs(1+N+1:end,:,:)), N * N * (p - 1), []);
    
    % find variables with largest devs
    %     maxdev = max(thesedevs, [], 2);
    %     [maxdev,ndx] = sort(maxdev, 'desc');
    %     thesendx = sort(ndx(maxdev > 1.5));
    
    hb = surf(ydates(T0+1:end),1:N*N * (p - 1),thesedevs);
    
    zlim([0 maxLagOTHER])
    caxis([0 maxRED])
    
    datetick('x')
    yticks(1:N:N*N)
    ylim([1 N*N])
    yticklabels(Ylabels);
    
    xlabel('end of estimation window')
    
    
    % color by height
    shading interp
    colorbar
    colormap turbo
    wrapthisfigure(thisfig, sprintf('PAI-lagOTHER-%s', modeltype), wrap, [], [], [], [], true);
    if doTitle
        title(sprintf('intercepts'))
        wrapthisfigure(thisfig, sprintf('PAI-lagOTHER-%s-WITHTITLE', modeltype), wrap);
    end
    
    %% plot PAI per block of each equation
    if doLags
        for n = 1 : N
            
            maxZ   = ceil(max(abs(PAIdevs(:,n,:)),[], 'all') * 10) / 10;
            
            for lag = 1 : p
                
                ndxlag = 1 + N * (lag-1) + (1 : N);
                thisfig = figure;
                
                hb = surf(ydates(T0+1:end),1:N,squeeze(abs(PAIdevs(ndxlag,n,:)))); %#ok<NASGU>
                
                zlim([0 maxZ])
                caxis([0 maxRED])
                datetick('x')
                xlabel('end of estimation window')
                %             ylabel('parameters')
                yticks(1:N)
                yticklabels(Ylabels)
                
                
                
                % color by height
                shading interp
                colorbar
                colormap turbo
                
                wrapthisfigure(thisfig, sprintf('PAI-%s-lag%d-%s', ncode{n}, lag, modeltype), wrap, [], [], [], [], true);
                if doTitle
                    title(sprintf('%s equation, lag %d', Ylabels{n}, lag))
                    wrapthisfigure(thisfig, sprintf('PAI-%s-lag%d-%s-WITHTITLE', ncode{n}, lag, modeltype), wrap);
                end
            end
            close all
        end
    end
    
    %% finish loop
    close all % dockAllFigures
    finishwrap
end

%% finish script
finishscript