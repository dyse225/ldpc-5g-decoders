%% E1b: E1 + bit-flipping декодеры на одном графике
% Код тот же, что в E1: 5G NR BG2, Zc = 52, K = 520, R = 1/2, BPSK, AWGN.
% Кривые BP / Layered BP / NMS / OMS берутся из results/e1_algorithms.mat
% (run_e1_algorithms), bit-flip декодеры считаются здесь в своём диапазоне.
%   Bit-flip            — ldpc_bitflip_decode, правило Галлагера,
%                         при неудаче выдаётся итоговое состояние;
%   Bit-flip (оригинал) — то же, но при неудаче выдаётся решение канала,
%                         как в исходной C-функции PH_bitflipping;
%   MWBF                — ldpc_wbf_decode, alpha = 0.7, до 300 инверсий,
%                         выколотые биты не участвуют в весах проверок.
clear; clc;
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');

S = load(fullfile(resDir, 'e1_algorithms.mat'));
results = S.results;
labels  = S.labels;

cfgs = {
    'bit-flip', 'final',   'Bit-flip'
    'bit-flip', 'channel', 'Bit-flip (как в оригинале)'
    'wbf',      '',        sprintf('MWBF \\alpha=%.1f', 0.7)
    };
for i = 1:size(cfgs, 1)
    p = ldpc_default_params();
    p.algorithm = cfgs{i, 1};
    if ~isempty(cfgs{i, 2}), p.bfOnFail = cfgs{i, 2}; end
    p.EbN0dB = 3:0.5:9;
    results = [results, ldpc_simulate(p)]; %#ok<AGROW>
    labels{end+1} = cfgs{i, 3}; %#ok<SAGROW>
end
save(fullfile(resDir, 'e1b_bitflip.mat'), 'results', 'labels');

% --- график: широкие панели друг под другом, ось Eb/N0 с разрывом.
% У MWBF за итерацию инвертируется один бит (сотни итераций), поэтому на
% панели итераций он не показывается.
plotRes = results;
plotRes(end).avgIter(:) = NaN;
c = get(groot, 'defaultAxesColorOrder');
c = [c(1:4, :); 0.15 0.15 0.15; 0.6 0.6 0.6; c(7, :)];
o = struct('layout', 'column', 'xlim', [0 9], 'xbreak', [2.25 2.75], ...
    'berLim', [1e-5 0.3], 'ferLim', [1e-4 1], 'colors', c);
fig = ldpc_plot_results(plotRes, labels, '', o);
exportgraphics(fig, fullfile(resDir, 'e1_all_column.png'), 'Resolution', 150);
