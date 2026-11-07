function run_e1b_bg1_bitflip(bgn)
% Добавляем кривую bit-flip для того же базового графа, что и загруженные кривые мягких декодеров.
if nargin < 1, bgn = 1; end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(fileparts(here))));
resDir = fullfile(here, 'results_1e6');
S = load(fullfile(resDir, sprintf('e1_bg%d_algorithms.mat', bgn)));
assert(numel(S.results) == 4, 'Run the four soft algorithms first.');
p = ldpc_default_params();
p.bgn = bgn;
if bgn == 1, p.Zc = 24; end
p.maxFrames = 1e6;
p.algorithm = 'bit-flip';
p.EbN0dB = 3:0.5:9;
results = [S.results, ldpc_simulate(p)];
labels = [S.labels, {'Bit-flip'}];
save(fullfile(resDir, sprintf('e1b_bg%d_bitflip.mat', bgn)), 'results', 'labels');
fig = ldpc_plot_results(results, labels, sprintf('BG%d, maxFrames=1e6', bgn));
exportgraphics(fig, fullfile(resDir, sprintf('e1_bg%d_all.png', bgn)), 'Resolution', 180);
close(fig);
end
