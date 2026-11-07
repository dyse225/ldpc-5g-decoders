%% Канал Rayleigh: BP, Layered BP, NMS, OMS, BF
% 5G NR BG2, K = 520, R = 1/2, BPSK. Плоские замирания, независимые для
% каждого символа (идеальное перемежение), приёмник знает канал точно
% (p.channel = 'rayleigh' в ldpc_simulate). Для сравнения — некодированная
% BPSK в Rayleigh по теории (berfading).
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

algs   = {'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum', 'bit-flip'};
labels = {'BP (flooding)', 'Layered BP', 'NMS \alpha=0.75', 'OMS \beta=0.5', 'Bit-flip'};
results = [];
for i = 1:numel(algs)
    p = ldpc_default_params();
    p.channel = 'rayleigh';
    p.algorithm = algs{i};
    if strcmp(algs{i}, 'bit-flip')
        p.EbN0dB = 8:1:20;
        p.maxFrames = 2e5;
    else
        p.EbN0dB = 2:0.25:5;
    end
    results = [results, ldpc_simulate(p)]; %#ok<AGROW>
end
save(fullfile(resDir, 'rayleigh.mat'), 'results', 'labels');

c = get(groot, 'defaultAxesColorOrder');
o = struct('layout', 'column', 'xlim', [2 20], 'xbreak', [5.25 7.5], ...
    'berLim', [1e-5 0.3], 'ferLim', [1e-4 1], 'colors', [c(1:4, :); 0.15 0.15 0.15]);
fig = ldpc_plot_results(results, labels, '', o);
exportgraphics(fig, fullfile(resDir, 'rayleigh.png'), 'Resolution', 150);

% --- некодированная BPSK в Rayleigh (теория) для оценки выигрыша кодирования
eb = 0:0.5:40;
berUnc = berfading(eb, 'psk', 2, 1);
save(fullfile(resDir, 'rayleigh.mat'), 'eb', 'berUnc', '-append');
