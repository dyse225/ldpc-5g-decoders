%% BF, WBF, MWBF вместе с NMS и OMS (5G NR BG2, K = 520, R = 1/2, BPSK, AWGN)
%   BF   — ldpc_bitflip_decode, правило 'gallager' = bitflipping из пакета
%          bf-tests-awgn-egldpc (Lcrypto), 50 итераций;
%   WBF  — ldpc_wbf_decode, alpha = 0 (совпадает с ldpc_wbf.m из Lcrypto);
%   MWBF — ldpc_wbf_decode, alpha = 0.7 (Zhang, Fossorier; порт из AFF3CT).
%   У WBF/MWBF за итерацию инвертируется один бит, лимит 300 итераций.
% NMS и OMS берутся из matlab/results/e1_algorithms.mat (run_e1_algorithms).
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

S = load(fullfile(here, '..', '..', 'results', 'e1_algorithms.mat'));
isMS = ismember(S.labels, {'NMS \alpha=0.75', 'OMS \beta=0.5'});
results = S.results(isMS);
labels = S.labels(isMS);

cfg = {
    'bit-flip', 0,   'BF'
    'wbf',      0,   'WBF'
    'wbf',      0.7, 'MWBF \alpha=0.7'
    };
for i = 1:size(cfg, 1)
    p = ldpc_default_params();
    p.algorithm = cfg{i, 1};
    p.wbfAlpha = cfg{i, 2};
    p.EbN0dB = 3:0.5:9;
    results = [results, ldpc_simulate(p)]; %#ok<AGROW>
    labels{end+1} = cfg{i, 3}; %#ok<SAGROW>
end
save(fullfile(resDir, 'bf_family.mat'), 'results', 'labels');

% --- BER/FER/итерации; у WBF и MWBF итерации на отдельном графике
plotRes = results;
for i = 4:5, plotRes(i).avgIter(:) = NaN; end
c = get(groot, 'defaultAxesColorOrder');
o = struct('layout', 'column', 'xlim', [0 9], 'xbreak', [2.25 2.75], ...
    'berLim', [1e-5 0.3], 'ferLim', [1e-4 1], 'colors', [c(3:4, :); 0.15 0.15 0.15; c(5, :); c(7, :)]);
fig = ldpc_plot_results(plotRes, labels, '', o);
exportgraphics(fig, fullfile(resDir, 'bf_family.png'), 'Resolution', 150);

fig = figure('Color', 'w', 'Position', [100 100 800 450]);
plot(results(4).EbN0dB, results(4).avgIter, '-^', 'Color', c(5, :), 'LineWidth', 2, 'MarkerSize', 7); hold on
plot(results(5).EbN0dB, results(5).avgIter, '-v', 'Color', c(7, :), 'LineWidth', 2, 'MarkerSize', 7);
grid on; xlabel('E_b/N_0, дБ'); ylabel('Ср. число итераций (= инверсий)');
legend(labels(4:5), 'Location', 'northeast'); set(gca, 'FontSize', 13);
title('WBF и MWBF: одна инверсия за итерацию, лимит 300');
exportgraphics(fig, fullfile(resDir, 'bf_family_wbf_iters.png'), 'Resolution', 150);
