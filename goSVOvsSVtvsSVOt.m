clear
close all

%% SVO
np = 12; % annual
p  = 1 / (4 * np); % once every four years
maxScale = 20;


SVOsupport = [1 2 : maxScale];
No = length(SVOsupport) - 1;
SVOpdf = [(1 - p), p * repmat(1/No, 1, No)];

%% t via inverse gamma * normal
% dof = 2 + 24 / ((maxScale - 2)^2 - 12) / p; % dof of t choosen to match variance of outlier state
% 
% % check
% if abs(dof / (dof - 2) - ((1-p) + p * (maxScale - 2)^2 / 12)) > 1e-10
%     error houston
% end

dof = 5;

% sample chi2
Ndraws = 1e5;
chi2draws = chi2rnd(repmat(dof, Ndraws, 1)) / dof; 

sqrtLambdaDraws = sqrt(1 ./ chi2draws);

[SVtpdf, SVtsupport] = ksdensity(sqrtLambdaDraws + eps, 0:.001:maxScale, 'support', 'positive');
SVtsupport = SVtsupport - eps;


%% SVO-t
Ndraws = 1e8;
% Ndraws = 1e3;
np = 12; % annual
p  = 1 / (10 * np); % once every four years
maxScale = 20;
dof = 9;

udraws = rand(Ndraws,1);
SVOstates = 1 : maxScale;
SVOtpdf = [(1 - p), p * repmat(1/No, 1, No)];
SVOtcdf = cumsum(SVOtpdf);
ndx    = sum(udraws > SVOtcdf,2) + 1;
outliers = SVOstates(ndx);

chi2draws = chi2rnd(repmat(dof, Ndraws, 1)) / dof; 

sqrtLambdaDraws = sqrt(1 ./ chi2draws) .* outliers';

[SVOtpdf, SVOtsupport] = ksdensity(sqrtLambdaDraws + eps, 0:.01:maxScale, 'support', 'positive');
SVOtsupport = SVOtsupport - eps;

%% clear mem
clear sqrtLambdaDraws chi2draws udraws ndx outliers

%% plot
theseTicks = [0 1  2 : 2 : maxScale];


fontsize = 18;


this = figure;
bar(SVOsupport, SVOpdf);
hold on
plot(SVtsupport, SVtpdf, 'k-', 'linewidth', 2);
plot(SVOtsupport, SVOtpdf, 'r:', 'linewidth', 2);
xticks(theseTicks)
ax1 = gca;
set(ax1, 'fontsize', fontsize)
orient landscape
figurename = 'SVOvsSVtvsSVOt';
print(this, '-depsc', '-r300', '-loose', figurename);

this = figure;
h1 = plot(SVOsupport(2:end), SVOpdf(2:end), 'b-.', 'linewidth', 2);
hold on
h2 = plot(SVtsupport, SVtpdf, 'k-', 'linewidth', 2);
h3 = plot(SVOtsupport, SVOtpdf, 'r:', 'linewidth', 2);
xticks(theseTicks)
xlim([3.5 maxScale+.5])
ylim([0 .01])
ax2 = gca;
set(ax2, 'fontsize', fontsize)
legend([h1 h2 h3], 'SVO', 'SV-t', 'SVO-t', 'box', 'off', 'fontsize', 28)
orient landscape
figurename = 'SVOvsSVtvsSVOt-zoomed';
print(this, '-depsc', '-r300', '-loose', figurename);
