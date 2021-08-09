clear
close all

%% SVO
np = 12; % annual
p  = 1 / ( 4 * np); % once every four years
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
Ndraws = 1e4;
chi2draws = chi2rnd(repmat(dof, Ndraws, 1)) / dof; 

sqrtLambdaDraws = sqrt(1 ./ chi2draws);

[l, lx] = ksdensity(sqrtLambdaDraws + eps, 0:.001:maxScale, 'support', 'positive');
lx = lx - eps;

%% plot
theseTicks = [0 1  2 : 2 : maxScale];


fontsize = 14;

% this = figure;
% subplot(1,2,1)
% bar(SVOsupport, SVOpdf);
% hold on
% plot(lx, l, 'r-', 'linewidth', 1);
% xticks(theseTicks)
% ax1 = gca;
% set(ax1, 'fontsize', fontsize)
% 
% 
% subplot(1,2,2)
% h1 = bar(SVOsupport, SVOpdf);
% hold on
% h2 = plot(lx, l, 'r-', 'linewidth', 1);
% xticks(theseTicks)
% xlim([1.5 maxScale+.5])
% ylim([0 .05])
% ax2 = gca;
% set(ax2, 'fontsize', fontsize)
% legend([h1 h2], 'SVO', 'SV-t', 'box', 'off')
% % title('(zoomed)')
% 
% orient landscape
% figurename = 'SVOvsSVt-twopanel';
% print(this, '-depsc', '-r300', '-loose', figurename);

this = figure;
bar(SVOsupport, SVOpdf);
hold on
plot(lx, l, 'r-', 'linewidth', 2);
xticks(theseTicks)
ax1 = gca;
set(ax1, 'fontsize', fontsize)
orient landscape
figurename = 'SVOvsSVt';
print(this, '-depsc', '-r300', '-loose', figurename);

this = figure;
h1 = bar(SVOsupport, SVOpdf);
hold on
h2 = plot(lx, l, 'r-', 'linewidth', 2);
xticks(theseTicks)
xlim([1.5 maxScale+.5])
ylim([0 .05])
ax2 = gca;
set(ax2, 'fontsize', fontsize)
legend([h1 h2], 'SVO', 'SV-t', 'box', 'off')
orient landscape
figurename = 'SVOvsSVt-zoomed';
print(this, '-depsc', '-r300', '-loose', figurename);
