%% Подбор alpha (NMS) и beta (OMS) на 16-QAM, как run_e3_minsum_params для BPSK
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

quick = false;                     % true - пробный прогон (до 2e4 кадров на точку)

EbN0dB = [4.25 4.5];               % водопад NMS/OMS для 16-QAM
alphas = 0.5:0.05:1.0;
betas  = 0:0.1:1.0;

base = ldpc_default_params();      % BG2, Zc = 52, K = 520, R = 1/2, AWGN
base.modulation = '16qam';
base.EbN0dB = EbN0dB;
base.stopAtZeroFer = false;
base.verbose = false;
base.maxFrames = 1e5;              % при 4.5 дБ это ~40 ошибочных кадров у NMS
if quick, base.maxFrames = 2e4; end

nA = numel(alphas); nB = numel(betas);
cfg = repmat(base, 1, nA + nB);
for i = 1:nA
    cfg(i).algorithm = 'norm-min-sum'; cfg(i).scaling = alphas(i);
end
for i = 1:nB
    cfg(nA + i).algorithm = 'offset-min-sum'; cfg(nA + i).offset = betas(i);
end
res = cell(1, nA + nB);
parfor i = 1:nA + nB
    res{i} = simulate_quantized(cfg(i));
end

ferAlpha = cell2mat(cellfun(@(r) r.fer.', res(1:nA), 'UniformOutput', false));
berAlpha = cell2mat(cellfun(@(r) r.ber.', res(1:nA), 'UniformOutput', false));
ferBeta  = cell2mat(cellfun(@(r) r.fer.', res(nA+1:end), 'UniformOutput', false));
berBeta  = cell2mat(cellfun(@(r) r.ber.', res(nA+1:end), 'UniformOutput', false));
nfeAlpha = cell2mat(cellfun(@(r) r.numFrameErrors.', res(1:nA), 'UniformOutput', false));
nfeBeta  = cell2mat(cellfun(@(r) r.numFrameErrors.', res(nA+1:end), 'UniformOutput', false));
for i = 1:nA
    fprintf('alpha = %.2f: FER = %s (ошибочных кадров %s)\n', alphas(i), ...
        mat2str(ferAlpha(:, i).', 3), mat2str(nfeAlpha(:, i).'));
end
for i = 1:nB
    fprintf('beta  = %.2f: FER = %s (ошибочных кадров %s)\n', betas(i), ...
        mat2str(ferBeta(:, i).', 3), mat2str(nfeBeta(:, i).'));
end

save(fullfile(resDir, 'minsum_params_16qam.mat'), 'EbN0dB', 'alphas', 'betas', ...
    'ferAlpha', 'berAlpha', 'ferBeta', 'berBeta', 'nfeAlpha', 'nfeBeta', 'base');

fig = plot_minsum_params_16qam(fullfile(resDir, 'minsum_params_16qam.mat'));
exportgraphics(fig, fullfile(resDir, 'minsum_params_16qam.png'), 'Resolution', 150);
