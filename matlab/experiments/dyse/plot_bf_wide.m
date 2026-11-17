function fig = plot_bf_wide(matFile)
%PLOT_BF_WIDE BF до 12 дБ и мягкие декодеры до 5 дБ.
%   fig = PLOT_BF_WIDE(matFile) — перерисовка из bf_wide.mat.
%
%   Нулевые FER показаны маркерами на уровне 0.5/N и слегка раздвинуты,
%   иначе на логарифмической оси они пропадают.

S = load(matFile);
bf = S.bf;
R = S.results;

% series: результат, FER, итерации, подпись, до какого Eb/N0 рисовать
rList = {bf, R(1), R(2), R(2), R(3), R(4), R(5)};
fer   = {bf.fer, R(1).fer, R(2).fer, R(2).ferPre, ...
         R(3).ferPre, R(4).ferPre, R(5).ferPre};
iter  = {bf.avgIter, R(1).avgIter, R(2).avgIter, R(2).avgIterPre, ...
         R(3).avgIterPre, R(4).avgIterPre, R(5).avgIterPre};
labels = {'BF', 'BP 2 + BF', 'BP 5 + BF', 'BP 5', 'BP 10', ...
          'Layered BP 5', 'NMS 5'};
maxX = [12 5 5 5 5 5 5];

n = numel(rList);
co = get(groot, 'defaultAxesColorOrder');
if size(co, 1) < n, co = lines(n); end
mk = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

fig = figure('Color', 'w', 'Position', [40 40 1500 950]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
axF = nexttile(tl);
axI = nexttile(tl);
for ax = [axF axI]
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 12);
    xlim(ax, [0 12]); xticks(ax, 0:1:12);
end
set(axF, 'YScale', 'log');
ylim(axF, [1e-6 2]); ylabel(axF, 'FER');
ylim(axI, [0 55]);   ylabel(axI, 'Ср. число итераций');
xlabel(axI, 'E_b/N_0, дБ');

for i = 1:n
    r = rList{i};
    ok = r.numFrames > 0 & r.EbN0dB <= maxX(i);
    x = r.EbN0dB(ok);
    draw_fer(axF, x, fer{i}(ok), 0.5 ./ r.numFrames(ok), i, co, mk, labels{i});
    plot(axI, x, iter{i}(ok), '-', 'Color', co(i, :), 'Marker', mk{i}, ...
        'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', labels{i});
end

xline(axF, 5, ':', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
xline(axI, 5, ':', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
text(axF, 5.15, 1.15, 'мягкие декодеры и BP + BF — до 5 дБ', ...
    'FontSize', 10, 'Color', [0.35 0.35 0.35]);
text(axF, 0.01, 0.04, ['0 ошибок: маркеры на уровне 0.5/N, ' ...
    'слегка раздвинуты'], 'Units', 'normalized', ...
    'FontSize', 10, 'Color', [0.35 0.35 0.35]);

title(tl, sprintf('BF и мягкие декодеры, %s, BPSK', bf.codeName), 'FontSize', 14);
subtitle(tl, 'BF — 0…12 дБ; BP / layered BP / NMS и BP + BF — 0…5 дБ', ...
    'FontSize', 11);

lg = legend(axF, 'Location', 'eastoutside');
lg.FontSize = 11;
end

% -------------------------------------------------------------------------
function draw_fer(ax, x, v, floor, i, co, mk, label)
v = v(:).'; floor = floor(:).';
zero = (v == 0) & isfinite(floor);
plot(ax, NaN, NaN, '-', 'Color', co(i, :), 'Marker', mk{i}, ...
    'LineWidth', 1.6, 'MarkerSize', 6, 'DisplayName', label);
if any(~zero)
    plot(ax, x(~zero), v(~zero), '-', 'Color', co(i, :), 'Marker', mk{i}, ...
        'LineWidth', 1.6, 'MarkerSize', 6, 'HandleVisibility', 'off');
end
if any(zero)
    yZero = floor(zero) .* (1 + 0.08 * (i - 4));
    plot(ax, x(zero), yZero, 'Marker', mk{i}, 'LineStyle', 'none', ...
        'Color', co(i, :), 'MarkerSize', 7, 'HandleVisibility', 'off');
    plot(ax, x, max(v, floor), ':', 'Color', co(i, :), ...
        'LineWidth', 1.0, 'HandleVisibility', 'off');
end
end
