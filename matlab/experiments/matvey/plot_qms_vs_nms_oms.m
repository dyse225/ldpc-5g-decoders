function fig = plot_qms_vs_nms_oms(matFile)
%PLOT_QMS_VS_NMS_OMS График BER по результатам RUN_QMS_VS_NMS_OMS.
%   fig = PLOT_QMS_VS_NMS_OMS(matFile) - перерисовка без пересчёта.
%   Стрелки на целевых BER - выигрыш кодирования QMS относительно
%   некодированной 16-QAM (по теории).

S = load(matFile);
coded = {S.rNms, S.rOms, S.rMs, S.rQms};

fig = figure('Color', 'w', 'Position', [100 100 900 600]);
ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
set(ax, 'YScale', 'log', 'FontSize', 11);
co = get(groot, 'defaultAxesColorOrder');
grey = [0.4 0.4 0.4];
plot(ax, S.rUnc.EbN0dB, S.berTh, '-', 'Color', grey, 'LineWidth', 1.3, ...
    'DisplayName', 'Некодированная 16-QAM, теория');
b = S.rUnc.ber; b(b == 0) = NaN;
plot(ax, S.rUnc.EbN0dB, b, 'o', 'Color', grey, 'LineWidth', 1.3, ...
    'DisplayName', 'Некодированная 16-QAM, модель');
style = {'-s', '-d', '--^', '-o'};
for i = 1:numel(coded)
    r = coded{i};
    b = r.ber; b(b == 0) = NaN;
    plot(ax, r.EbN0dB, b, style{i}, 'Color', co(i, :), 'LineWidth', 1.6, ...
        'MarkerSize', 6, 'DisplayName', S.labels{i});
end

% выигрыш кодирования QMS на целевых BER (строка 5 в ebAt - QMS, 1 - некодированная)
for t = 1:numel(S.targetBer)
    yline(ax, S.targetBer(t), ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    x1 = S.ebAt(5, t); x2 = S.ebAt(1, t);
    if any(isnan([x1 x2])), continue, end
    plot(ax, [x1 x2], S.targetBer(t) * [1 1], '-', 'Color', co(4, :), ...
        'LineWidth', 0.8, 'HandleVisibility', 'off');
    text(ax, (x1 + x2) / 2, S.targetBer(t) * 1.25, ...
        sprintf('BER = 10^{%d}: QMS выигрывает %.1f дБ', round(log10(S.targetBer(t))), ...
        x2 - x1), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'Color', co(4, :));
end

xlabel(ax, 'E_b/N_0, дБ'); ylabel(ax, 'BER');
xlim(ax, [2 S.rUnc.EbN0dB(end)]); ylim(ax, [1e-6 0.5]);
title(ax, sprintf('Quantized Min-Sum, %s, 16-QAM, AWGN, I_{max} = %d', ...
    S.rNms.codeName, S.rNms.params.maxIter));
legend(ax, 'Location', 'northeast');
end
