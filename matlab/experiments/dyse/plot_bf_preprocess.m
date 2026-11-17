function fig = plot_bf_preprocess(matFile)
%PLOT_BF_PREPROCESS График BF с мягкой предобработкой.
%   fig = PLOT_BF_PREPROCESS(matFile) — перерисовка из bf_preprocess.mat.
%
%   Нулевые BER/FER не выбрасываются: они показываются маркерами на уровне
%   верхней оценки 0.5/N (для BER — 0.5/(K*N)) и слегка раздвинуты по
%   вертикали, чтобы было видно все конфигурации. Настоящие значения —
%   в bf_preprocess.csv.

S = load(matFile);
results = S.results;
labels = S.labels;
n = numel(results);
co = get(groot, 'defaultAxesColorOrder');
if size(co, 1) < n, co = lines(n); end
mk = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

fig = figure('Color', 'w', 'Position', [40 40 1500 1000]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axF = nexttile(tl, [1 2]);
axB = nexttile(tl);
axI = nexttile(tl);
for ax = [axF axB axI]
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12);
    xlim(ax, [3 6]); xticks(ax, 3:0.5:6);
end
set([axF axB], 'YScale', 'log');
ylim(axF, [1e-5 2]);  ylabel(axF, 'FER');
ylim(axB, [1e-8 1]);  ylabel(axB, 'BER');
ylim(axI, [0 55]);    ylabel(axI, 'Ср. число итераций');
xlabel(axF, 'E_b/N_0, дБ');
xlabel(axB, 'E_b/N_0, дБ');
xlabel(axI, 'E_b/N_0, дБ');

for i = 1:n
    r = results(i);
    ok = r.numFrames > 0;
    x = r.EbN0dB(ok);
    % FER: нижняя оценка для нуля — 0.5/N
    draw_metric(axF, x, r.fer(ok), 0.5 ./ r.numFrames(ok), i, co, mk, labels{i});
    % BER: нижняя оценка для нуля — 0.5/(K*N)
    draw_metric(axB, x, r.ber(ok), 0.5 ./ (r.K * r.numFrames(ok)), i, co, mk, labels{i});
    % итерации
    plot(axI, x, r.avgIter(ok), '-', 'Color', co(i, :), 'Marker', mk{i}, ...
        'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', labels{i});
end

text(axF, 0.01, 0.04, ['0 ошибок: маркеры на уровне 0.5/N, ' ...
    'слегка раздвинуты по вертикали'], 'Units', 'normalized', ...
    'FontSize', 10, 'Color', [0.35 0.35 0.35]);
text(axB, 0.01, 0.04, '0 ошибок: маркеры на уровне 0.5/(K·N)', ...
    'Units', 'normalized', 'FontSize', 10, 'Color', [0.35 0.35 0.35]);

title(tl, sprintf('BF после мягкой предобработки, %s, BPSK', results(1).codeName), ...
    'FontSize', 14);
subtitle(tl, sprintf(['3–6 дБ, до %d ошибочных кадров, лимит %d кадров; ' ...
    'пунктир — оценка нуля сверху'], results(1).params.minFrameErrors, ...
    results(1).params.maxFrames), 'FontSize', 11);

lg = legend(axF, 'Location', 'eastoutside');
lg.FontSize = 11;
end

% -------------------------------------------------------------------------
function draw_metric(ax, x, v, floor, i, co, mk, label)
% Линия по ненулевым точкам + открытые маркеры на уровне floor для нулей.
v = v(:).'; floor = floor(:).';
zero = (v == 0) & isfinite(floor);
% Фиктивный handle, чтобы конфигурация попала в легенду даже при всех нулях
plot(ax, NaN, NaN, '-', 'Color', co(i, :), 'Marker', mk{i}, ...
    'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', label);
if any(~zero)
    plot(ax, x(~zero), v(~zero), '-', 'Color', co(i, :), 'Marker', mk{i}, ...
        'LineWidth', 1.6, 'MarkerSize', 6, 'HandleVisibility', 'off');
end
if any(zero)
    % небольшой сдвиг по вертикали, чтобы нули разных конфигураций не слились
    yZero = floor(zero) .* (1 + 0.12 * (i - 4));
    plot(ax, x(zero), yZero, 'Marker', mk{i}, 'LineStyle', 'none', ...
        'Color', co(i, :), 'MarkerSize', 7, 'HandleVisibility', 'off');
    plot(ax, x, max(v, floor), ':', 'Color', co(i, :), ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');
end
end
