% generate_freddata.m
% =========================================================================
% DESCRIPTION:
% This script loads in raw data from a monthly database CSV file,
% transforms each series based on transformation code using
% prepare_missing.m, and removes outliers from the transformed data using
% remove_outliers.m.
%
% NOTE:
% The default CSV file read by this code is 2015-04.csv, which is the April
% 2015 version of the dataset. If using a different version, make sure to
% change the variable "csv_in" on line 26 to match the name of the relevant
% CSV file.
%
% =========================================================================
% CLEAR:
clear
close all
clc

%% load em toolboxes
path(pathdef)

addpath ../matlabtoolbox/emtools/
addpath ../matlabtoolbox/emtexbox/
addpath ../matlabtoolbox/emgibbsbox/
addpath ../matlabtoolbox/emeconometrics/
addpath ../matlabtoolbox/emstatespace/
addpath ..

%% settings
datalabel = 'fredMD16levels';

vintage = '2021-04';

outputlabel = strcat(datalabel, '-', vintage);

doQuarterly = false;

%#ok<*UNRCH>

% =========================================================================
% PARAMETER TO BE CHANGED:
% Update the .csv filename to match the desired version

% CSV file name
% csv_in='2020-10.csv';
csv_in= strcat(vintage, '.csv');

% =========================================================================
% LOAD AND LABEL DATA:
% Load data from CSV file
dum=importdata(csv_in,',');

% Variable names
names=dum.textdata(1,2:end);

% Transformation numbers
tcode=dum.data(1,:);

% Raw data
rawdata=dum.data(2:end,:); % first row contains tcodes

% Month of final observation
% final_month=month(dum.textdata(end,1));

% Year of final observation
% final_year=year(dum.textdata(end,1));

% =========================================================================
% SET UP DATES:
% Dates (monthly) are of the form YEAR+MONTH/12
% e.g. March 1970 is represented as 1970+3/12
% Dates go from 1959:01 to final_year:final_month (see above)
% dates = (1959+1/12:1/12:final_year+final_month/12)';

dates = datenum(dum.textdata(3:end,1), 'mm/dd/yy', 1950); % using matlab dates

% T = number of months in sample
T      = size(dates,1);
rawdata= rawdata(1:T,:);


%% do not first difference the fedfunds rate or yields
ndxFFR = strcmp(names, 'FEDFUNDS');
tcode(ndxFFR) = 1;

tcode(strcmp(names, 'GS1'))     = 1;
tcode(strcmp(names, 'GS10'))    = 1;
tcode(strcmp(names, 'GS5'))     = 1;
tcode(strcmp(names, 'PCEPI'))   = 5;
tcode(strcmp(names, 'UNRATE'))  = 1;

% do not double difference inflation, not difference unrate
tcode(tcode == 2) = 1;
tcode(tcode == 6) = 5;

% other
if contains(datalabel, 'levels')
    tcode(tcode == 5) = 4;
end


if strcmp(datalabel, 'fredMDlenzaprimiceri6')
    tcode(:) = 1;
end

%% apply levels if desired
if contains(datalabel, 'levels')
    
    tcode(tcode == 5) = 4;
    %     tcode(strcmp(names, 'RPI'))     = 4;
    %     tcode(strcmp(names, 'INDPRO'))    = 4;
    %     tcode(strcmp(names, 'DPCERA3M086SBEA'))    = 4;
    %     tcode(strcmp(names, 'PAYEMS'))    = 4;


end


%% transform into quarterly data if needed
if doQuarterly
    
    quarters = quarterlydates(dates);
    
    
    
    qdata = grpstats(rawdata, quarters);
    qdates = unique(quarters);

    % drop incomplete quarters
    qcount = arrayfun(@(x) sum(quarters == x), qdates);
    ndx    = qcount == 3;
    
    rawdata = qdata(ndx,:);
    dates   = qdates(ndx,:);
    
    T       = size(dates,1);

    if any(~ndx)
        warning('dropped %d monthly data points that belong to incomplete quarters', sum(~ndx))
    end
    
end

%% TRANSFORM RAW DATA INTO STATIONARY FORM:
% Use function prepare_missing.m
%   Output yt: matrix containing data after transformation
%
%     case 1, % Level (i.e. no transformation): x(t)
%     case 2, % First difference: x(t)-x(t-1)
%     case 3, % Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
%     case 4, % Natural log: ln(x)
%     case 5, % First difference of natural log: ln(x)-ln(x-1)
%     case 6, % Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
%     case 7, % First difference of percent change: (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)



yt  = prepare_missing(rawdata,tcode);
ndx = tcode == 5;
if doQuarterly
    yt(:,ndx) = yt(:,ndx) * 400;
else
    yt(:,ndx) = yt(:,ndx) * 1200;
end

ndx = tcode == 4;
yt(:,ndx) = yt(:,ndx) * 100;


% =========================================================================
% REMOVE OUTLIERS: EM dropped
% Use function remove_outliers.m (see for definition of outliers)
%   Output data: matrix containing transformed series after removal of outliers
%   Output n: matrix containing number of outliers removed for each series
% [~,n]=remove_outliers(yt);
%
% disp('there are quite a few outliers:')
% display(n)

data = yt;

% =========================================================================
% SELECT 20 VARIABLES AND STORE IN CSV

switch datalabel
    case {'fredMD3'}
        % pick Unrate, PCEdeflator and FFR
        codeVariableSelection = {'UNRATE', 'PCEPI', 'FEDFUNDS'};
    case {'fredMD16', 'fredMD16levels'}
        codeVariableSelection = {'RPI', 'DPCERA3M086SBEA', 'INDPRO', ...
            'CUMFNS', 'UNRATE', 'PAYEMS', 'CES0600000007', 'CES0600000008', 'WPSFD49207', ...
            'PCEPI', 'HOUST', 'S&P 500', 'EXUSUKx', ...
            'GS5', 'GS10', 'BAAFFM'};
    case {'fredMD8'}
        codeVariableSelection = {'RPI', 'UNRATE', 'PAYEMS', 'WPSFD49207', ...
            'PCEPI', 'GS5', 'GS10', 'BAAFFM'};
    otherwise
        error('datalabel %s not recognized', datalabel)
end



% map FRED-MD into list of 20
[~,ndxVariableSelection] = ismember(codeVariableSelection, names);

% collect data and names
varnames  = names(ndxVariableSelection);
tcode     = tcode(ndxVariableSelection);

cumcode   = false(1,length(ndxVariableSelection));
cumcode(tcode == 2) = true;
% cumcode(tcode == 5) = true; per Todd's suggestion from Dec 9 2019
cumcode(tcode == 6) = true;

tabledata = data(:,ndxVariableSelection);



varnames = matlab.lang.makeValidName(varnames, 'ReplacementStyle', 'delete'); % to get rid of '&' and other unwanted signs, seems necessary at least under linux

%% REDUCE SAMPLE TO USABLE DATES:
% Remove first two months because some series have been second differenced
if any(tcode == 3 | tcode == 6)
    tabledata = yt(3:T,:);
    dates     = dates(3:T,:);
elseif any(tcode == 2 | tcode == 5 | tcode == 7)
    tabledata = yt(2:T,:);
    dates     = dates(2:T,:);
end

% normalize log levels
ndx              = tcode == 4;
tabledata(:,ndx) = tabledata(:,ndx) - tabledata(1,ndx); % + 100;


%% check for outliers
close all
tableNoOutliers = remove_outliers(tabledata);
ndx = isnan(tableNoOutliers);
for n = 1 : size(tabledata,2)
    thisdata = tableNoOutliers(:,n);
    nanny = isnan(thisdata);
    % if any(nanny)
    theseOutliers = NaN(size(tabledata,1),1);
    theseOutliers(nanny) = tabledata(nanny,n);
    
    figure
    hold on
    plot(dates, thisdata)
    if any(nanny)
        plot(dates, theseOutliers, 'rx', 'linewidth', 2)
        plot(dates, theseOutliers, 'ro', 'linewidth', 2)
    end
    xtickdates(dates)
	titletxt = sprintf('%s (%d)', varnames{n}, tcode(n));
    title(titletxt)
    set(gcf, 'name', titletxt)
    % end
end

%% clean missing values
nanny = any(isnan(tabledata), 2);

if ~iscompact(nanny)
    error('missing data inside sample')
end

if any(nanny)
    warning('\n there is some missing data: %s', datestr(dates(nanny)))
    warning('\n data is missing for %s', varnames{any(isnan(tabledata), 1)})
end
   
%% prepare table
tabledata = tabledata(~nanny,:);
dates     = dates(~nanny);
% T         = length(dates);

%% construct yield code
yieldcode = ismember(varnames, {'FEDFUNDS', 'GS5', 'GS10', 'GS1'});

%% prepend dates, tcode etc and store table
if doQuarterly
    datalabel = strcat(datalabel, '-quarterly');
end


tabledata = [tcode; cumcode; tabledata];

datatable   = array2table(tabledata,  'VariableNames', varnames);

% check
if any(any(ismissing(datatable)))
    error('there are missing observations')
end

tabledates  = [NaN;NaN;dates]; % note: recycle the variable name dates
output = cat(2, table(tabledates, 'VariableNames', {'dates'}), datatable);
writetable(output, sprintf('%s.csv', outputlabel))

%% define minnesota prior means
N = length(varnames);
minnesotaPriorMean = ones(N,1);
% for n = 1 : N
%     switch varnames{n}
%         case {'CUMFNS', 'UNRATE', ...
%                 'WPSFD49207',  'PPICMM', 'PCEPI', ...
%                 'FEDFUNDS', 'HOUST', 'GS5', 'GS10', 'BAAFFM', 'WUXIASHADOWRATE'}
%             minnesotaPriorMean(n) = 1;
%         otherwise
%             minnesotaPriorMean(n) = 0;
%     end
% end


%% prepare fortran files
% output as txt file for fortran
mat2fortran(sprintf('%s.yData.txt', outputlabel), tabledata(3:end,:))
% mat2fortran(sprintf('%s.cumcode.txt', outputlabel), tabledata(2,:))
int2fortran(sprintf('%s.tcode.txt', outputlabel), tabledata(1,:))
mat2fortran(sprintf('%s.minnesotaPriorMean.txt', outputlabel), minnesotaPriorMean)
logical2fortran(sprintf('%s.yieldcode.txt', outputlabel), yieldcode)
mat2fortran(sprintf('%s.dates.txt', outputlabel), dates)

% output varnames
cellstr2fortran(sprintf('%s.ynames.txt', outputlabel), varnames)
cellstr2fortran(sprintf('%s.ylabel.txt', outputlabel), fredMDprettylabel(varnames))

%% generates table of variables
varlabels = fredMDprettylabel(varnames);
N = length(varlabels);

filename = sprintf('datalist-%s.tex', outputlabel);
fid = fopen(filename, 'wt');

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{lll}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Variable & FRED-MD code & tcode ');
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for n = 1 : N
    fprintf(fid, '%s ', varlabels{n});
    fprintf(fid, '& %s ', varnames{n});
    fprintf(fid, '& %d ', tcode(n));
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');
fprintf(fid, 'Note: ');
fprintf(fid, 'Data obtained from the %s vintage of FRED-MD. ', strtok(csv_in, '.'));
if doQuarterly
    fprintf(fid, 'Quarterly observations (constructed from monthly averages) ');
else
    fprintf(fid, 'Monthly observations ');
end
fprintf(fid, 'from %s to %s. ', datestryymm(dates(1)), datestryymm(dates(end)));
fprintf(fid, '%s \n', fredMDtcodeNote);
fclose(fid);
type(filename)

%% generates alt table of variables
varlabels = fredMDprettylabel(varnames);
N = length(varlabels);

filename = sprintf('datalist2-%s.tex', outputlabel);
fid = fopen(filename, 'wt');

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{lll}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Variable & FRED-MD code & transformation ');
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for n = 1 : N
    fprintf(fid, '%s ', varlabels{n});
    fprintf(fid, '& %s ', varnames{n});
    
     switch tcode(n)
        case 1
            % % no transformation
        case 2
            fprintf(fid, '& %s\n', '\ensuremath{\Delta x_t}');
        case 4
            fprintf(fid, '& %s\n', '\ensuremath{\log(x_t)}');
        case 5
            fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 1200}');
        otherwise
            error houston
    end
    
    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');
fprintf(fid, 'Note: ');
fprintf(fid, 'Data obtained from the %s vintage of FRED-MD. ', strtok(csv_in, '.'));
if doQuarterly
    fprintf(fid, 'Quarterly observations (constructed from monthly averages) ');
else
    fprintf(fid, 'Monthly observations ');
end
fprintf(fid, 'from %s to %s. ', datestryymm(dates(1)), datestryymm(dates(end)));
% fprintf(fid, '%s \n', fredMDtcodeNote);
fclose(fid);
type(filename)

%% generates alt table of variables: transformation and minnesota mean
varlabels = fredMDprettylabel(varnames);
N = length(varlabels);

filename = sprintf('datalist3-%s.tex', outputlabel);
fid = fopen(filename, 'wt');

fprintf(fid, '\\begin{center}\n');
fprintf(fid, '\\begin{tabular}{lllc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Variable & FRED-MD code & transformation & Minnesota prior');
fprintf(fid, '\\\\\n');
fprintf(fid, '\\midrule\n');
for n = 1 : N
    fprintf(fid, '%s ', varlabels{n});
    fprintf(fid, '& %s ', varnames{n});
    
     switch tcode(n)
        case 1
            fprintf(fid, ' & ');
        case 2
            fprintf(fid, ' & %s', '\ensuremath{\Delta x_t}');
        case 4
            fprintf(fid, ' & %s', '\ensuremath{\log(x_t)}');
        case 5
            if doQuarterly
                fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 400}');
            else
                fprintf(fid, '& %s\n', '\ensuremath{\Delta\log(x_t) \cdot 1200}');
            end
        otherwise
            fprintf(fid, ' & ');
     end
    
     fprintf(fid, ' & %d ', minnesotaPriorMean(n));

    fprintf(fid, '\\\\\n');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{center}\n');
fprintf(fid, '\n');
fprintf(fid, 'Note: ');
fprintf(fid, 'Data obtained from the %s vintage of FRED-MD. ', strtok(csv_in, '.'));
if doQuarterly
    fprintf(fid, 'Quarterly observations (constructed from monthly averages) ');
else
    fprintf(fid, 'Monthly observations ');
end
fprintf(fid, 'from %s to %s.\n', datestryymm(dates(1)), datestryymm(dates(end)));
fprintf(fid, 'Entries in the column ``Minnesota prior'''' report the prior mean on the first own-lag coefficient of the corresponding variable in each BVAR. Prior means on all other VAR coefficients are set to zero.\n');
% fprintf(fid, '%s \n', fredMDtcodeNote);
fclose(fid);
type(filename)

%% generates some stuff to collect latex tables
varlabels = fredMDprettylabel(varnames);
N = length(varlabels);

filename = 'scratch.tex';
fid = fopen(filename, 'wt');

fprintf(fid, 'LABELS:\n');
for n = 1 : N
    fprintf(fid, '%s\n', varlabels{n});
end
fprintf(fid, '\n');

fprintf(fid, 'NAMES:\n');
for n = 1 : N
    fprintf(fid, '%s\n', varnames{n});
end

fprintf(fid, 'FLOAT:\n');
for n = 1 : N
    fprintf(fid, '\\subfloat[%s]{\\includegraphics[width=\\picwid]{pdf%s-\\jumpoff}\\label{subfig:%s-\\jumpoff}}\n', varlabels{n}, varnames{n}, varnames{n});
    if mod(n,3) == 0
        fprintf(fid, '\\\\\n');
    else
        fprintf(fid, '\\quad\n');
    end
end




fclose(fid);
type(filename)


%% finish
dockAllFigures