%% produce charts of data

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
datalabel           = 'fredMD16levels-2021-04';

doQuarterly         = false;

samStart            = []; % datenum(1988,12,1);                 % truncate start of sample if desired (leave empty if otherwise)


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
Ylabels = strrep(Ylabels, '\', '');

Ylabelslong = Ylabels;
Ylabelslong(ismember(tcode, [2 5])) = strcat(Ylabelslong(ismember(tcode, [2 5])), ...
    sprintf(' \n (cont. compounded growth p.a. in p.p.)'));

%% process settings
N = size(data,2);

% truncate start of sample (if desired)
if ~isempty(samStart)
    ndx  = ydates >= samStart;
    data   = data(ndx,:);
    ydates = ydates(ndx);
    Tdata  = length(ydates);
end

yearStart = year(ydates(1));
yearStop  = year(ydates(end)) + 9;

if doQuarterly
    fcstdates = genrQdates(yearStart, yearStop);
else
    fcstdates = genrMdates(yearStart, yearStop, 1);
end
ndx = find(fcstdates == ydates(1));
if isempty(ndx)
    error houston
end
fcstdates = fcstdates(ndx:end);



initwrap
% wrap = [];
fontsize = 18;

%% plot input data
for n = 1 : N
    
    % determne SW outliers (based on pre COVID sample)
    %     sam = ydates < datenum(2020,1,1);
    sam = true(size(ydates)); 
    dev = abs(data(:,n) - median(data(sam,n)));
    iqr = range(prctile(data(sam,n), [25 75]));
    outndx = dev > 5 * iqr;
    
    this = figure;
    hold on
    plot(ydates, data(:,n), 'k-', 'linewidth', 2)
    plot(ydates(outndx), data(outndx,n), 'rd', 'linewidth', 2)
    ax = gca;
    set(ax, 'fontsize', fontsize)
    box(ax, 'off')
    xticks(datenum(1960:10:2020,1,1))
    xtickdates(ydates, 'keepticks')
    wrapthisfigure(this, sprintf('data%s', ncode{n}), wrap)
    title(Ylabelslong{n})
    wrapthisfigure(this, sprintf('data%s-WITHTITLE', ncode{n}), wrap)
end

%% plot recent data
% for n = 1 : N
%     this = figure;
%     
%     switch tcode(n)
%         case {2,5} % rates of change
%             bar(ydates(end-23:end), data(end-23:end,n), 1, 'facecolor', .45 * [1 1 1])
%         otherwise
%             plot(ydates(end-23:end), data(end-23:end,n), 'k-', 'linewidth', 2)
%     end
%     ax = gca;
%     set(ax, 'fontsize', fontsize)
%     box(ax, 'off')
%     
%     xtickdates([ydates(end-24:end); ydates(end) + mean(diff(ydates))])
%     datetick('x', 'yyyy:mm', 'keeplimits')
%     
%     wrapthisfigure(this, sprintf('recentdata%s', ncode{n}), wrap)
%     title(Ylabelslong{n})
%     wrapthisfigure(this, sprintf('recentdata%s-WITHTITLE', ncode{n}), wrap)
% end

%% code beamer frames
fid = fopen('foo.tex', 'wt');


for n = 1 : N
    fprintf(fid, '%s\n', '\begin{frame}');
    fprintf(fid, '%s\n', '\setlength{\skipper}{.25\baselineskip}');
    fprintf(fid, '\\frametitle{%s}\n', upper(Ylabels{n}));
    switch tcode(n)
        case 1
            % % no transformation
        case 2
            fprintf(fid, '%s\n', '\framesubtitle{\fmath{\Delta x_t}}');
        case 4
            fprintf(fid, '%s\n', '\framesubtitle{\fmath{\log(x_t)}}');
        case 5
            fprintf(fid, '%s\n', '\framesubtitle{\fmath{\Delta\log(x_t) \cdot 1200}}');
        otherwise
            error houston
    end
    fprintf(fid, '%s\n', '\vspace{-.75\baselineskip}');
    
    fprintf(fid, '%s\n', '\begin{center}');
    
    fprintf(fid, '\\includegraphics[width=\\textwidth]{data%s}\n', ncode{n});
    
    
    fprintf(fid, '%s\n', '\end{center}');
    
    fprintf(fid, '%s\n\n\n', '\end{frame}');
end

fclose(fid);
type('foo.tex')



%% stats for the paper's narrative
T = length(ydates);

display('PAYEMS')
ndx = strcmp(ncode, 'PAYEMS');
[v,d] = sort(data(:,ndx), 'descend');
v = v / 12; % deannualize
for i = [1 : 3, T - 3 : T]
    fprintf('%8.2f \t in %s\n', v(i), datestr(ydates(d(i))));
end
hrulefill
display('UNRATE')
ndx = strcmp(ncode, 'UNRATE');
[v,d] = sort(data(:,ndx), 'descend');
for i = [1 : 3, T - 3 : T]
    fprintf('%8.2f \t in %s\n', v(i), datestr(ydates(d(i))));
end
hrulefill
display('INDPRO')
ndx = strcmp(ncode, 'INDPRO');
[v,d] = sort(data(:,ndx), 'descend');
v = v / 12; % deannualize
for i = [1 : 3, T - 3 : T]
    fprintf('%8.2f \t in %s\n', v(i), datestr(ydates(d(i))));
end
hrulefill
display('RPI')
ndx = strcmp(ncode, 'RPI');
[v,d] = sort(data(:,ndx), 'descend');
v = v / 12; % deannualize
for i = [1 : 6, T - 6 : T]
    fprintf('%8.2f \t in %s\n', v(i), datestr(ydates(d(i))));
end

%% wrap up
dockAllFigures
finishwrap





