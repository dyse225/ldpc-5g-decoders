%% Графики и таблицы: 256-QAM, BG1 и BG2
clear; clc;
here = fileparts(mfilename('fullpath')); addpath(here);
resDir = fullfile(here, 'results');
A = load(fullfile(resDir, 'a_params.mat'));
B = load(fullfile(resDir, 'b_curves.mat'));
codeNames = {'BG1 (K=528)', 'BG2 (K=520)'};
algs = B.algs; iterList = B.iterList;
algNames = {'BP (flooding)', 'Layered BP', 'NMS (опт. \alpha)', 'OMS (опт. \beta)'};
algShort = {'BP (flooding)', 'Layered BP', 'NMS', 'OMS'};
acol = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56];
icol = turbo(numel(iterList) + 2); icol = icol(2:end-1, :);
mk = @(w, h) figure('Color', 'w', 'Position', [50 50 w h]);

%% Рис. A: подбор параметров
for metric = ["fer", "ber"]
    f = mk(1200, 800); tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for c = 1:2
        for aIdx = 1:2
            nexttile; hold on; grid on;
            algNs = {'norm-min-sum', 'offset-min-sum'}; algName = algNs{aIdx};
            for i = 1:numel(A.iterList)
                sel = find(cellfun(@(x) x == c, A.jobs(:,1)) & cellfun(@(x) x == A.iterList(i), A.jobs(:,2)) & strcmp(A.jobs(:,3), algName));
                par = cell2mat(A.jobs(sel, 4)); y = cellfun(@(r) r.(metric), A.out(sel));
                semilogy(par, y, '-o', 'Color', icol(i,:), 'LineWidth', 1.4, 'MarkerSize', 4);
            end
            set(gca, 'YScale', 'log'); ylabel(upper(metric));
            if aIdx == 1, xlabel('\alpha (NMS)'); else, xlabel('\beta (OMS)'); xlim([0 1.5]); end
            title(sprintf('%s, %s', codeNames{c}, algShort{aIdx+2}));
        end
    end
    lg = legend(arrayfun(@(n) sprintf('I_{max}=%d', n), A.iterList, 'UniformOutput', false)); lg.Layout.Tile = 'east';
    title(tl, sprintf('%s от параметра min-sum, 256-QAM, E_b/N_0 = %.1f дБ', upper(metric), A.ebno));
    set_light(f); exportgraphics(f, fullfile(resDir, sprintf('a_params_%s.png', metric)), 'Resolution', 150);
end

%% Рис. B: сравнение декодеров при одинаковом I_max (строки: FER, BER; столбцы: I_max)
for c = 1:2
    f = mk(1500, 700); tl = tiledlayout(2, numel(iterList), 'TileSpacing', 'compact', 'Padding', 'compact');
    for mi = 1:2
        mm = {'fer', 'ber'}; metric = mm{mi};
        for i = 1:numel(iterList)
            nexttile; hold on; grid on;
            for a = 1:numel(algs)
                r = B.res{c, a, i}; ok = r.(metric) > 0 & ~isnan(r.(metric));
                semilogy(r.EbN0dB(ok), r.(metric)(ok), '-o', 'Color', acol(a,:), 'LineWidth', 1.4, 'MarkerSize', 4);
            end
            set(gca, 'YScale', 'log'); xlim([8.5 11]); xlabel('E_b/N_0, дБ'); ylabel(upper(metric));
            if mi == 1, ylim([1e-3 1]); title(sprintf('I_{max} = %d', iterList(i))); else, ylim([1e-5 0.5]); end
        end
    end
    lg = legend(algNames); lg.Layout.Tile = 'east';
    title(tl, sprintf('%s: сравнение декодеров при одинаковом числе итераций (256-QAM, R=1/2)', codeNames{c}));
    set_light(f); exportgraphics(f, fullfile(resDir, sprintf('b_decoders_bg%d.png', c)), 'Resolution', 150);
end

%% Рис. C: эффект числа итераций для каждого декодера (строки: BG1, BG2; столбцы: декодеры), FER
f = mk(1500, 650); tl = tiledlayout(2, numel(algs), 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:2
    for a = 1:numel(algs)
        nexttile; hold on; grid on;
        for i = 1:numel(iterList)
            r = B.res{c, a, i}; ok = r.fer > 0 & ~isnan(r.fer);
            semilogy(r.EbN0dB(ok), r.fer(ok), '-o', 'Color', icol(i,:), 'LineWidth', 1.4, 'MarkerSize', 4);
        end
        set(gca, 'YScale', 'log'); xlim([8.5 11]); ylim([1e-3 1]); xlabel('E_b/N_0, дБ'); ylabel('FER');
        title(sprintf('%s, %s', codeNames{c}, algShort{a}));
    end
end
lg = legend(arrayfun(@(n) sprintf('I_{max}=%d', n), iterList, 'UniformOutput', false)); lg.Layout.Tile = 'east';
title(tl, 'FER при разном числе итераций (256-QAM)');
set_light(f); exportgraphics(f, fullfile(resDir, 'c_iterations_effect.png'), 'Resolution', 150);

%% Рис. D: среднее число итераций, I_max = 50
f = mk(1000, 450); tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:2
    nexttile; hold on; grid on;
    for a = 1:numel(algs)
        r = B.res{c, a, iterList == 50}; plot(r.EbN0dB, r.avgIter, '-o', 'Color', acol(a,:), 'LineWidth', 1.4, 'MarkerSize', 4);
    end
    xlabel('E_b/N_0, дБ'); ylabel('среднее число итераций'); title(codeNames{c}); xlim([8.5 11]);
end
lg = legend(algNames); lg.Layout.Tile = 'east';
title(tl, 'Среднее число итераций при ранней остановке, I_{max}=50, 256-QAM');
set_light(f); exportgraphics(f, fullfile(resDir, 'd_avg_iterations.png'), 'Resolution', 150);

%% Таблица: Eb/N0 при FER = 1e-2 и BER = 1e-3, средние итерации и время декодирования при 10 дБ
rows = {};
for c = 1:2
    for a = 1:numel(algs)
        for i = 1:numel(iterList)
            r = B.res{c, a, i};
            k10 = find(abs(r.EbN0dB - 10) < 1e-9);
            if isempty(k10), it10 = NaN; tm10 = NaN; fer10 = NaN; else, it10 = r.avgIter(k10); tm10 = r.timePerFrame(k10) * 1e3; fer10 = r.fer(k10); end
            rows(end+1, :) = {codeNames{c}, algShort{a}, iterList(i), ebno_at(r.EbN0dB, r.fer, 1e-2), ...
                ebno_at(r.EbN0dB, r.ber, 1e-3), fer10, it10, tm10}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'code', 'decoder', 'Imax', 'EbN0_FER1e2', 'EbN0_BER1e3', 'FER_at_10dB', 'avgIter_10dB', 'ms_per_frame_10dB'});
writetable(T, fullfile(resDir, 'summary_table.csv'));
disp(T);

% оптимальные параметры
[aO, bO, TP] = pick_params(fullfile(resDir, 'a_params.mat')); %#ok<ASGLU>
writetable(TP, fullfile(resDir, 'optimal_params.csv'));
save(fullfile(resDir, 'summary.mat'), 'T', 'TP', 'aO', 'bO');

function x = ebno_at(ebno, y, target)
ok = y > 0 & ~isnan(y); e = ebno(ok); l = log10(y(ok)); x = NaN;
for k = 1:numel(e)-1
    if l(k) >= log10(target) && l(k+1) <= log10(target) && l(k) ~= l(k+1)
        x = e(k) + (log10(target) - l(k)) * (e(k+1) - e(k)) / (l(k+1) - l(k)); return
    end
end
end

function set_light(f)
if isprop(f, 'Theme'), f.Theme = 'light'; end
drawnow;
for a = findall(f, 'Type', 'axes').', try, a.Toolbar.Visible = 'off'; end, end %#ok<TRYNC>
end
