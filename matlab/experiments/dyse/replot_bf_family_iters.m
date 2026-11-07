%% BF, WBF, MWBF, NMS, OMS: итерации всех декодеров на одной логарифмической оси
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');

S = load(fullfile(resDir, 'bf_family.mat'));
c = get(groot, 'defaultAxesColorOrder');
o = struct('layout', 'column', 'xlim', [0 9], 'xbreak', [2.25 2.75], ...
    'berLim', [1e-5 0.3], 'ferLim', [1e-4 1], 'colors', [c(3:4, :); 0.15 0.15 0.15; c(5, :); c(7, :)]);
fig = ldpc_plot_results(S.results, S.labels, '', o);

% нижний ряд: подпись оси x есть только у него
ax = findobj(fig, 'Type', 'axes');
ax = ax(arrayfun(@(a) ~isempty(a.XLabel.String), ax));
[~, order] = sort(arrayfun(@(a) a.XLim(1), ax));
ax = ax(order);
ticks = [1 5 10 20 50 100 300];
set(ax, 'YScale', 'log', 'YLim', [1 400], 'YTick', ticks);
set(ax(1), 'YTickLabel', string(ticks));
set(ax(2), 'YTickLabel', []);
text(ax(2), 0.03, 0.10, ['WBF, MWBF: одна инверсия за итерацию, лимит 300' newline ...
    'BF, NMS, OMS: лимит 50'], 'Units', 'normalized', 'FontSize', 12);

exportgraphics(fig, fullfile(resDir, 'bf_family_log_iters.png'), 'Resolution', 150);
