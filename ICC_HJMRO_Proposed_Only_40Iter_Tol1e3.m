function ResultsConvergence = ICC_HJMRO_Proposed_Only_40Iter_Tol1e3()
% Load shared and experiment-specific MATLAB functions.
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot,'shared'), fullfile(projectRoot,'convergence'));
% PROPOSED H-JMRO ONLY: true AO/SCA convergence at fixed T = 220 s.
% 40 maximum trajectory iterations; stopping tolerance = 1e-3.
% Each accepted 3-UAV trajectory is evaluated using 300 Monte Carlo samples
% and averaged across 30 spatial drops (the original simulation settings).
% Early termination repeats the final accepted trajectory through iteration 40.
% Individual per-PU and per-WU interference protection is unchanged.
%
% Requires MATLAB and CVX. The file contains all proposed-scheme helpers.
close all; clc;
if exist('cvx_begin','file') ~= 2
    error('CVX is required. Run cvx_setup before starting.');
end

fprintf('\n===== PROPOSED H-JMRO: 40 AO/SCA ITERATIONS | tolerance = 1e-3 =====\n');
proposedRates = run_proposed_icc_3uav();
Iteration = (1:numel(proposedRates)).';
assert(numel(proposedRates) == 40, ...
    'Expected 40 proposed-scheme rate observations.');

ResultsConvergence = table(Iteration,proposedRates(:), ...
    'VariableNames',{'Iteration','Proposed_3UAV_Mbps'});
disp(ResultsConvergence);
writetable(ResultsConvergence,'ICC_HJMRO_Proposed_Only_40Iter_Convergence.csv');
save('ICC_HJMRO_Proposed_Only_40Iter_Convergence.mat','ResultsConvergence');

%% ===================== Proposed H-JMRO convergence plot =====================
figure('Color','w','Position',[100 100 620 470]);
plot(Iteration,proposedRates,'-o','LineWidth',2,'MarkerSize',4, ...
    'MarkerIndices',1:numel(Iteration));
grid on; box on;
set(gca,'FontName','Times New Roman','FontSize',16,'FontWeight','normal', ...
    'LineWidth',1.2,'GridLineStyle',':','GridAlpha',0.30, ...
    'Position',[0.13 0.13 0.84 0.68]);
xlabel('Number of iterations','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal');
ylabel('Average sum rate (Mbps)','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal');
xlim([1 40]);
xticks([1 5:5:40]);
legend('Proposed H-JMRO','Location','northoutside', ...
    'Orientation','horizontal','FontName','Times New Roman', ...
    'FontSize',15.2,'FontWeight','normal','Box','on');

fprintf('\nFinal proposed H-JMRO average sum rate: %.6f Mbps\n', proposedRates(end));
fprintf('Results saved to ICC_HJMRO_Proposed_Only_40Iter_Convergence.csv and .mat\n');
end
