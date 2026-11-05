function fig = plot_lsat_sweep(matFile)
%PLOT_LSAT_SWEEP BER и FER от L_sat по результатам RUN_LSAT_SWEEP.
%   fig = PLOT_LSAT_SWEEP(matFile) - перерисовка без пересчёта.
%   Пунктир - min-sum без квантования при том же Eb/N0.

S = load(matFile);
fig = figure('Color', 'w', 'Position', [100 100 900 420]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
title(tl, sprintf('L_{sat} для 5-битных LLR, %s, 16-QAM, min-sum', S.rRef.codeName));
lbl = arrayfun(@(e) sprintf('E_b/N_0 = %.2f дБ', e), S.EbN0dB, 'UniformOutput', false);
co = get(groot, 'defaultAxesColorOrder');
for k = 1:2
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'YScale', 'log');
    if k == 1, v = S.ber; vr = S.rRef.ber; yl = 'BER'; else, v = S.fer; vr = S.rRef.fer; yl = 'FER'; end
    for j = 1:numel(S.EbN0dB)
        plot(ax, S.LsatList, v(j, :), '-o', 'Color', co(j, :), 'LineWidth', 1.3, ...
            'DisplayName', [lbl{j} ', 5 бит']);
        yline(ax, vr(j), '--', 'Color', co(j, :), 'LineWidth', 1.1, ...
            'DisplayName', [lbl{j} ', без квантования']);
    end
    xlabel(ax, 'L_{sat}'); ylabel(ax, yl);
end
lg = legend(ax, 'NumColumns', 2);
lg.Layout.Tile = 'south';
end
