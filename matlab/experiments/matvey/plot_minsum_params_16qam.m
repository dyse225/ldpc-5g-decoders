function fig = plot_minsum_params_16qam(matFile)
%PLOT_MINSUM_PARAMS_16QAM FER от alpha (NMS) и beta (OMS) по RUN_MINSUM_PARAMS_16QAM.
%   fig = PLOT_MINSUM_PARAMS_16QAM(matFile) - перерисовка без пересчёта.
%   Вертикальный пунктир - значения проекта по умолчанию (alpha = 0.75, beta = 0.5).

S = load(matFile);
fig = figure('Color', 'w', 'Position', [100 100 900 420]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
title(tl, sprintf('Параметры min-sum на 16-QAM, 5G NR BG2 K = 520 R = 1/2, I_{max} = %d', ...
    S.base.maxIter));
lbl = arrayfun(@(e) sprintf('E_b/N_0 = %.2f дБ', e), S.EbN0dB, 'UniformOutput', false);
co = get(groot, 'defaultAxesColorOrder');
x = {S.alphas, S.betas}; v = {S.ferAlpha, S.ferBeta};
xl = {'\alpha (NMS)', '\beta (OMS)'}; def = [S.base.scaling S.base.offset];
for k = 1:2
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'YScale', 'log');
    for j = 1:numel(S.EbN0dB)
        f = v{k}(j, :); f(f == 0) = NaN;
        plot(ax, x{k}, f, '-o', 'Color', co(j, :), 'LineWidth', 1.3, 'DisplayName', lbl{j});
    end
    xline(ax, def(k), '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'по умолчанию');
    xlabel(ax, xl{k}); ylabel(ax, 'FER');
end
lg = legend(ax, 'NumColumns', 3);
lg.Layout.Tile = 'south';
end
