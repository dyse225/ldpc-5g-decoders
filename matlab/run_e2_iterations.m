%% E2: влияние максимального числа итераций (BP flooding и NMS layered)
clear; clc;
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

algs = {'bp', 'norm-min-sum'};
iterList = [5 10 20 50];

results = []; labels = {};
for a = 1:numel(algs)
    for it = iterList
        p = ldpc_default_params();
        p.algorithm = algs{a};
        p.maxIter = it;
        results = [results, ldpc_simulate(p)]; %#ok<AGROW>
        labels{end+1} = sprintf('%s, I_{max}=%d', algs{a}, it); %#ok<SAGROW>
    end
end

save(fullfile(resDir, 'e2_iterations.mat'), 'results', 'labels');
fig = ldpc_plot_results(results, labels, ...
    sprintf('E2: число итераций, %s, BPSK', results(1).codeName));
exportgraphics(fig, fullfile(resDir, 'e2_iterations.png'), 'Resolution', 150);
