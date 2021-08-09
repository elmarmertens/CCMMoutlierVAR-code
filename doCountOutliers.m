%% Count outliers in QRT data based on IQR threshold

%#ok<*NOSEL>
%#ok<*DISPLAYPROG>
%#ok<*UNRCH>

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

%% set parameters for VAR and MCMC

datalabel           = 'fredMD16-2021-04';
Ofactor             = 5;

titlename = sprintf('doOutlierCount-o%d', Ofactor);
initwrap

% fontsize = 12; % for paper
fontsize = 16; % for charts

%% load data

% load CSV file
dum=importdata(sprintf('%s.csv', datalabel),',');


ydates=dum.data(3:end,1);
% Variable names
ncode=dum.textdata(1,2:end);
% Transformation codes (data are already transformed)
tcode  =dum.data(1,2:end);
cumcode=logical(dum.data(2,2:end));
% Data
data=dum.data(3:end,2:end);

Tdata = length(ydates);

Ylabels = fredMDprettylabel(ncode);

T = length(ydates);

%% process settings
N = size(data,2);

Tjumpoffs = find(ydates >= datenum(1985,1,1));

Njumpoffs = length(Tjumpoffs);

%% count outliers
outlierNdx = NaN(Tdata, N, Njumpoffs);

for ndxT = 1 : Njumpoffs
    
    thisT = Tjumpoffs(ndxT);
    
    thisdata = data(1:thisT,:);
    dev      = abs(thisdata - median(data));
    iqr      = range(prctile(thisdata, [25 75]));
    
    outlierNdx(1:thisT,:,ndxT) = dev > Ofactor * iqr;

end


%% plot input data
for n = 1 : N
    
    sam = true(size(ydates)); 
    dev = abs(data(:,n) - median(data(sam,n)));
    iqr = range(prctile(data(sam,n), [25 75]));
    outndx = dev > Ofactor * iqr;
    
    this = figure;
    hold on
    plot(ydates, data(:,n), 'k-', 'linewidth', 2)
    plot(ydates(outndx), data(outndx,n), 'rd', 'linewidth', 2)
    ax = gca;
    set(ax, 'fontsize', fontsize)
    box(ax, 'off')
    xticks(datenum(1960:10:2020,1,1))
    xtickdates(ydates, 'keepticks')
    %     wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
    title(Ylabels{n})
    wrapthisfigure(this, sprintf('data%s-WITHTITLE', ncode{n}), wrap)
    
    if tcode(n) == 5
        thislevel = cumsum(data(:,n) / 1200) + 100;
        
        this = figure;
        hold on
        plot(ydates, thislevel, 'k-', 'linewidth', 2)
        plot(ydates(outndx), thislevel(outndx), 'ro', 'linewidth', 1)
        ax = gca;
        set(ax, 'fontsize', fontsize)
        box(ax, 'off')
        xticks(datenum(1960:10:2020,1,1))
        xtickdates(ydates, 'keepticks')
        %     wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
        title(Ylabels{n})
        wrapthisfigure(this, sprintf('leveldata%s-WITHTITLE', ncode{n}), wrap)
    end
    
end


%% 3D
outlierShare = nanmean(outlierNdx,3) * 100; %#ok<NANMEAN>

outlierShare = [zeros(T,1) outlierShare];
Ylabels2 = cat(2, ' ', Ylabels);

this1 = figure;
surf(0:N,ydates,outlierShare)
set(gca, 'fontsize', fontsize)
xticks(find(any(outlierShare > 0, 1)))
xticklabels(Ylabels2(any(outlierShare,1)))
xlim([0 N])
xtickangle(45)
ylim(ydates([1 end]))
yticks(datenum(1970:10:2020,1,1))
datetick('y', 'yyyy', 'keeplimits', 'keepticks')
shading interp
view(-26, 49)
grid off
wrapthisfigure(this1, sprintf('outlierAvgCount-o%d', Ofactor), wrap)

%% 3D for last vintage
outlierShare = outlierNdx(:,:,end);

outlierShare = [zeros(T,1) outlierShare];
Ylabels2 = cat(2, ' ', Ylabels);

this2 = figure;
surf(0:N,ydates,outlierShare)
set(gca, 'fontsize', fontsize)
xticks(find(any(outlierShare > 0, 1)))
xticklabels(Ylabels2(any(outlierShare,1)))
xtickangle(45)
xlim([0 N])
ylim(ydates([1 end]))
yticks(datenum(1970:10:2020,1,1))
zticks([0 1])
zticklabels({'no', 'yes'})
datetick('y', 'yyyy', 'keeplimits', 'keepticks')
shading interp
view(-26, 49)
grid off
wrapthisfigure(this2, sprintf('outlierLastCount-o%d', Ofactor), wrap)


%% outliers per year
wrap = diary2wrap('outlierdates.log', wrap);

for n = 1  : N

    hrulefill
    fprintf('\n%s\n', upper(Ylabels{n}))
    
    for v = [1:12:Njumpoffs Njumpoffs]
        
        thisT        = Tjumpoffs(v);
        thisCount    = sum(outlierNdx(1:thisT,n,v));
        
        if thisCount > 0

            fprintf('\t vintage %s, %d outliers in %d months (%5.2f %%) \n', ...
                datestr(ydates(Tjumpoffs(v)), 'yyyy:mmm'), thisCount, thisT, thisCount / thisT * 100)
            
            for t = find(outlierNdx(1:thisT,n,v)');
                fprintf('\t\t\t%s \n', datestr(ydates(t), 'yyyy:mmm'))
            end
        end
    end
end
    

%% wrap up
dockAllFigures
% exportgraphics(this1, sprintf('outlierAvgCount-o%d.eps', Ofactor));
finishwrap
finishscript
