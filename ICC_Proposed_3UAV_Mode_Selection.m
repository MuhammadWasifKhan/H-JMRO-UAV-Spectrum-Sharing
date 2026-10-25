function ResultsProposed = ICC_Proposed_3UAV_Mode_Selection()
% Load shared and experiment-specific MATLAB functions.
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot,'shared'), fullfile(projectRoot,'flight_time'));
% Proposed heterogeneous C/L/U 3-UAV only; NO single-C-UAV baseline.
% Figures: (1) actual average sum rate, (2) observed mode selection (%).
% Mode fractions are measured from optimized per-slot decisions across
% every Monte Carlo realization and spatial drop; no fabricated percentages.
% Requirements: MATLAB with CVX installed and configured on the path.
close all;
clc;
if exist('cvx_begin','file') ~= 2
    error('CVX is required. Install and run cvx_setup before running.');
end
[proposedRates,modePercent] = run_proposed_icc_3uav();
T = (200:20:300).';
modeLabels = {'C/C','L/U','C/U','L/C'};
assert(isequal(size(modePercent),[numel(T),4]),'Unexpected mode percentage matrix size.');
assert(all(abs(sum(modePercent,2)-100)<1e-6), ...
    'Mode percentages should sum to 100 percent for every flight time.');
ResultsProposed = table(T,proposedRates(:),modePercent(:,1),modePercent(:,2), ...
    modePercent(:,3),modePercent(:,4), ...
    'VariableNames',{'Flight_Time_s','Proposed_3UAV_Mbps', ...
    'Mode_CC_pct','Mode_LU_pct','Mode_CU_pct','Mode_LC_pct'});
fprintf('\n========= PROPOSED 3-UAV RATE AND MODE SELECTION =========\n');
disp(ResultsProposed);
writetable(ResultsProposed,'ICC_Proposed_3UAV_Mode_Selection.csv');
save('ICC_Proposed_3UAV_Mode_Selection.mat', ...
    'ResultsProposed','T','proposedRates','modePercent','modeLabels');
% Figure 1: proposed 3-UAV rate only.
figure('Color','w','Position',[100 100 620 470]);
plot(T,proposedRates(:),'-o','LineWidth',2,'MarkerSize',8);
grid on; box on;
ax=gca;
set(ax,'FontName','Times New Roman','FontSize',16,'FontWeight','normal', ...
    'LineWidth',1.2,'GridLineStyle',':','GridAlpha',0.30, ...
    'Position',[0.13 0.13 0.84 0.68]);
xlabel('Flight time (s)','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal');
ylabel('Average sum rate (Mbps)','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal');
legend('Proposed heterogeneous 3-UAV','Location','northoutside', ...
    'Orientation','horizontal','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal','Box','on');
xticks(T); xlim([min(T) max(T)]);
% Figure 2: actual selected mode percentages, grouped by flight time.
figure('Color','w','Position',[100 100 620 470]);
b=bar(T,modePercent,'grouped','BarWidth',0.82);
% Match mode palette in supplied example: blue, orange, yellow, purple.
colors=[0.0000 0.4470 0.7410; ...
        0.8500 0.3250 0.0980; ...
        0.9290 0.6940 0.1250; ...
        0.4940 0.1840 0.5560];
for m=1:4
    set(b(m),'FaceColor',colors(m,:),'EdgeColor','k','LineWidth',0.5);
end
grid on; box on;
ax=gca;
set(ax,'FontName','Times New Roman','FontSize',16,'FontWeight','normal', ...
    'LineWidth',1.2,'GridLineStyle',':','GridAlpha',0.30, ...
    'Position',[0.13 0.13 0.84 0.68]);
xlabel('Flight time period, T (s)','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal');
ylabel('Selected mode fraction (%)','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal');
legend(b,modeLabels,'Location','northwest', ...
    'Orientation','horizontal','NumColumns',2, ...
    'FontName','Times New Roman','FontSize',13, ...
    'FontWeight','normal','Box','on');
xticks(T); xlim([min(T)-12 max(T)+12]);
ylim([0 max(40,5*ceil(max(modePercent(:))/5))]);
end
