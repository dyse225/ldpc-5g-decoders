%% E1: сравнение алгоритмов декодирования (5G NR BG2, K=520, R=1/2, BPSK, AWGN)
clear; clc;
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

algs   = {'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum'};
labels = {'BP (flooding)', 'Layered BP', 'NMS \alpha=0.75', 'OMS \beta=0.5'};

results = [];
for i = 1:numel(algs)
    p = ldpc_default_params();
    p.algorithm = algs{i};
    results = [results, ldpc_simulate(p)]; %#ok<AGROW>
end

save(fullfile(resDir, 'e1_algorithms.mat'), 'results', 'labels');
fig = ldpc_plot_results(results, labels, ...
    sprintf('E1: алгоритмы, %s, BPSK, maxIter = 50', results(1).codeName));
exportgraphics(fig, fullfile(resDir, 'e1_algorithms.png'), 'Resolution', 150);
