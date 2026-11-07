function run_e1_bg1_algorithms(bgn)
% Четыре мягких декодера; по умолчанию BG1, для сравнения можно задать bgn=2.
if nargin < 1, bgn = 1; end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(fileparts(here))));
resDir = fullfile(here, 'results_1e6');
if ~exist(resDir, 'dir'), mkdir(resDir); end
algs = {'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum'};
labels = {'BP (flooding)', 'Layered BP', 'NMS \alpha=0.75', 'OMS \beta=0.5'};
results = [];
for i = 1:numel(algs)
    p = ldpc_default_params();
    p.maxFrames = 1e6;
    p.bgn = bgn;
    if bgn == 1, p.Zc = 24; end
    p.algorithm = algs{i};
    results = [results, ldpc_simulate(p)]; %#ok<AGROW>
    save(fullfile(resDir, sprintf('e1_bg%d_algorithms.mat', bgn)), 'results', 'labels');
end
fig = ldpc_plot_results(results, labels, sprintf('BG%d: %s', bgn, results(1).codeName));
exportgraphics(fig, fullfile(resDir, sprintf('e1_bg%d_algorithms.png', bgn)), 'Resolution', 180);
close(fig);
end
