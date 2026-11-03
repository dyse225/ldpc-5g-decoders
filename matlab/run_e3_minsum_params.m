%% E3: подбор alpha (NMS) и beta (OMS) при фиксированных Eb/N0
clear; clc;
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

EbN0dB = [1.25 1.75];
alphas = 0.5:0.05:1.0;
betas  = 0:0.1:1.0;

ferAlpha = zeros(numel(EbN0dB), numel(alphas));
ferBeta  = zeros(numel(EbN0dB), numel(betas));

p = ldpc_default_params();
p.EbN0dB = EbN0dB;
p.stopAtZeroFer = false;
p.verbose = false;

p.algorithm = 'norm-min-sum';
for i = 1:numel(alphas)
    p.scaling = alphas(i);
    r = ldpc_simulate(p);
    ferAlpha(:, i) = r.fer.';
    fprintf('alpha = %.2f: FER = %s\n', alphas(i), mat2str(r.fer, 3));
end

p.algorithm = 'offset-min-sum';
for i = 1:numel(betas)
    p.offset = betas(i);
    r = ldpc_simulate(p);
    ferBeta(:, i) = r.fer.';
    fprintf('beta  = %.2f: FER = %s\n', betas(i), mat2str(r.fer, 3));
end

save(fullfile(resDir, 'e3_minsum_params.mat'), 'EbN0dB', 'alphas', 'betas', ...
    'ferAlpha', 'ferBeta');

fig = figure('Color', 'w', 'Position', [100 100 900 360]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
title(tl, sprintf('E3: параметры min-sum, %s, BPSK, maxIter = %d', r.codeName, p.maxIter));
lbl = arrayfun(@(e) sprintf('E_b/N_0 = %.1f дБ', e), EbN0dB, 'UniformOutput', false);
nexttile; semilogy(alphas, ferAlpha, '-o', 'LineWidth', 1.3); grid on;
xlabel('\alpha (NMS)'); ylabel('FER'); legend(lbl);
nexttile; semilogy(betas, ferBeta, '-s', 'LineWidth', 1.3); grid on;
xlabel('\beta (OMS)'); ylabel('FER'); legend(lbl);
exportgraphics(fig, fullfile(resDir, 'e3_minsum_params.png'), 'Resolution', 150);
