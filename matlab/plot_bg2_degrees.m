%% Базовый граф 5G NR BG2: где биты степени 1
% Закрашенная клетка — блок Zc x Zc (циклически сдвинутая единичная
% матрица), белая — нулевой блок. Цвет — группа столбцов.
clear; clc;
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');
Z = 52;

cfg = nr5g.internal.ldpc.getParams(2, Z);
Hf = cfg.ParityCheckMatrix;
MB = size(Hf, 1) / Z;
NB = size(Hf, 2) / Z;
B = false(MB, NB);
for i = 1:MB
    for j = 1:NB
        B(i, j) = any(any(Hf((i-1)*Z + (1:Z), (j-1)*Z + (1:Z))));
    end
end
colDeg = sum(B, 1);

% группы столбцов: 1 — инфо выколотые, 2 — инфо, 3 — ядро проверочных,
% 4 — хвост проверочных (степень 1)
grp = zeros(1, NB);
grp(1:2) = 1; grp(3:10) = 2; grp(11:14) = 3; grp(15:NB) = 4;
C = [1 1 1; 0.55 0.55 0.55; 0.15 0.15 0.15; 0.12 0.47 0.71; 0.85 0.33 0.10];
img = zeros(MB, NB);
for j = 1:NB
    img(B(:, j), j) = grp(j);
end

fig = figure('Color', 'w', 'Position', [50 50 1500 1050]);
image(img + 1); colormap(C); axis image; hold on
set(gca, 'XTick', 1:NB, 'YTick', 1:MB, 'FontSize', 9, 'TickLength', [0 0], ...
    'XAxisLocation', 'top');
for x = 0.5:1:NB+0.5, plot([x x], [0.5 MB+0.5], 'Color', [0.85 0.85 0.85]); end
for y = 0.5:1:MB+0.5, plot([0.5 NB+0.5], [y y], 'Color', [0.85 0.85 0.85]); end

% часть, используемая при R = 1/2: 12 строк-блоков x 22 столбца-блока
rectangle('Position', [0.5 0.5 22 12], 'EdgeColor', 'k', 'LineWidth', 3, 'LineStyle', '--');
text(23, 6.5, {'используется при R = 1/2', '(12 x 22 блоков)'}, 'FontSize', 12);

for j = 1:NB
    text(j, MB + 1.3, num2str(colDeg(j)), 'HorizontalAlignment', 'center', 'FontSize', 9);
end
text(0.3, MB + 1.3, 'степень:', 'HorizontalAlignment', 'right', 'FontSize', 10);
ylim([0.5 MB + 2]);
xlabel('столбцы-блоки (каждый = Z_c = 52 бита)', 'FontSize', 12);
ylabel('строки-блоки (каждая = 52 проверки)', 'FontSize', 12);

names = {'1–2: инфо, выколотые', '3–10: инфо', '11–14: ядро проверочных бит', ...
    '15–52: хвост проверочных бит — степень 1'};
h = gobjects(1, 4);
for k = 1:4, h(k) = patch(NaN, NaN, C(k+1, :)); end
legend(h, names, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 12);
title('5G NR BG2: каждая клетка — блок 52x52', 'FontSize', 13);
exportgraphics(fig, fullfile(resDir, 'bg2_degree1.png'), 'Resolution', 130);
fprintf('степени столбцов-блоков:\n'); disp(colDeg);
