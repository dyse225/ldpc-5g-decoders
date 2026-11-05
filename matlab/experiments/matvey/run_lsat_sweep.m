%% Выбор уровня насыщения L_sat для 5-битных LLR (16-QAM, QMS)
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

quick = false;                     % true - пробный прогон (до 2e4 кадров на точку)

EbN0dB = [4.5 5.25];                % водопад min-sum для 16-QAM (FER ~ 1e-1 и ~ 1e-2)
LsatList = [2 3 4 6 8 10 12 16 24];

fer = zeros(numel(EbN0dB), numel(LsatList));
ber = zeros(numel(EbN0dB), numel(LsatList));
for i = 1:numel(LsatList)
    p = ldpc_default_params();     % BG2, Zc = 52, K = 520, R = 1/2, AWGN
    p.modulation = '16qam';
    p.algorithm = 'min-sum';
    p.llrBits = 5;
    p.llrSat = LsatList(i);
    p.EbN0dB = EbN0dB;
    p.stopAtZeroFer = false;
    if quick, p.maxFrames = 2e4; else, p.maxFrames = 2e5; end
    r = simulate_quantized(p);
    fer(:, i) = r.fer.';
    ber(:, i) = r.ber.';
end

% эталон: min-sum без квантования на тех же кадрах
p.llrBits = Inf;
rRef = simulate_quantized(p);

save(fullfile(resDir, 'lsat_sweep.mat'), 'EbN0dB', 'LsatList', 'fer', 'ber', 'rRef');

fig = plot_lsat_sweep(fullfile(resDir, 'lsat_sweep.mat'));
exportgraphics(fig, fullfile(resDir, 'lsat_sweep.png'), 'Resolution', 150);
