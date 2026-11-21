function qli_plots(runs, cfg)
%QLI_PLOTS Все графики задания (для каждого BG отдельно), PNG в <resDir>/figures:
%   fig1_ber_decoders_BG*      BER, 4 декодера при I_max = baseIter + Uncoded QPSK
%   fig2_fer_decoders_BG*      FER, 4 декодера при I_max = baseIter
%   fig3_ber_nms_imax_BG*      BER NMS при I_max = iterList
%   fig4_fer_nms_imax_BG*      FER NMS
%   fig5_ber_lbp_imax_BG*      BER Layered BP
%   fig6_fer_lbp_imax_BG*      FER Layered BP
%   fig7_avgiter_decoders_BG*  среднее число итераций, 4 декодера при I_max = baseIter
figDir = fullfile(cfg.resDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
col = get(groot, 'defaultAxesColorOrder');
style = {'-s', '-d', '--^', '-o'};
nd = numel(cfg.decoders);

for bgn = cfg.bgList
    sfx = sprintf('BG%d', bgn);
    K = 0;
    for i = 1:numel(runs)
        if runs{i}.meta.bgn == bgn && ~strcmp(runs{i}.meta.algorithm, 'uncoded')
            K = runs{i}.res.K; break
        end
    end
    ttl = sprintf('%s, 5G NR LDPC BG%d (Zc=%d, K=%d), R=%g, AWGN', upper(cfg.modulation), bgn, ...
        cfg.Zc.(sprintf('bg%d', bgn)), K, cfg.rate);
    unc = find_run(runs, bgn, 'uncoded', 0);

    % ---- 1-2: декодеры при baseIter ----------------------------------------
    metNames = {'ber', 'fer'};
    for fgi = 1:2
        met = metNames{fgi}; yl = upper(met);
        fig = newfig(); ax = newax(fig, yl);
        for d = 1:nd
            r = find_run(runs, bgn, cfg.decoders{d}, cfg.baseIter);
            if isempty(r), continue; end
            [x, y] = curve(r.res, met);
            plot(ax, x, y, style{d}, ...
                'Color', col(d, :), ...
                'LineWidth', 1.6, ...
                'MarkerSize', 6, ...
                'DisplayName', cfg.decoderNames{d});
        end
        % Некодированная QPSK: моделируемая и теоретическая кривые
        if fgi == 1 && ~isempty(unc)

            % 1. Моделируемая BER
             [x, y] = curve(unc.res, 'ber');

            plot(ax, x, y, 'o', ...
                'Color', [0.3 0.3 0.3], ...
                'MarkerSize', 5, ...
                'LineWidth', 1.2, ...
                'DisplayName', 'Uncoded QPSK (модель)');

            % 2. Теоретическая BER
            xt = cfg.EbN0uncoded;

            berTheory = 0.5 .* erfc(sqrt(10 .^ (xt ./ 10)));

            plot(ax, xt, berTheory, '-', ...
                'Color', [0.4 0.4 0.4], ...
                'LineWidth', 1.3, ...
                'DisplayName', 'Uncoded QPSK (теория)');

        end
        title(ax, sprintf('%s, I_{max} = %d', ttl, cfg.baseIter), 'FontSize', 11, 'FontWeight', 'bold');
        lg = legend(ax, 'Location', 'best'); set(lg, 'FontSize', 9, 'Box', 'on');
        if fgi == 1, ylim(ax, [1e-7 1]); end
        save_fig(fig, fullfile(figDir, sprintf('fig%d_%s_decoders_%s.png', fgi, met, sfx)));
    end

    % ---- 3-6: влияние I_max (NMS и Layered BP) ---------------------------------
    cases = {'norm-min-sum', 'NMS', 3; 'layered-bp', 'Layered BP', 5};
    for c = 1:2
        for fgi = 1:2
            met = metNames{fgi};
            fig = newfig(); ax = newax(fig, upper(met)); nPl = 0;
            for j = 1:numel(cfg.iterList)
                r = find_run(runs, bgn, cases{c, 1}, cfg.iterList(j));
                if isempty(r), continue; end
                [x, y] = curve(r.res, met); nPl = nPl + 1;
                plot(ax, x, y, style{j}, ...
                    'Color', col(j, :), ...
                    'LineWidth', 1.6, ...
                    'MarkerSize', 6, ...
                    'DisplayName', sprintf('I_{max} = %d', cfg.iterList(j)));
            end
            title(ax, sprintf('%s: %s, %s', cases{c, 2}, upper(met), ttl), 'FontSize', 11, 'FontWeight', 'bold');
            if nPl > 0, lg = legend(ax, 'Location', 'best'); set(lg, 'FontSize', 9, 'Box', 'on'); end
            save_fig(fig, fullfile(figDir, sprintf('fig%d_%s_%s_imax_%s.png', cases{c, 3} + fgi - 1, ...
                met, strrep(strrep(lower(cases{c, 2}), ' ', ''), 'layeredbp', 'lbp'), sfx)));
        end
    end

    % ---- 7: среднее число итераций -------------------------------------------
    fig = newfig(); ax = newax(fig, 'Среднее число итераций');
    set(ax, 'YScale', 'linear');
    for d = 1:nd
        r = find_run(runs, bgn, cfg.decoders{d}, cfg.baseIter);
        if isempty(r), continue; end
        [x, y] = curve(r.res, 'avgIter');
        plot(ax, x, y, style{d}, ...
            'Color', col(d, :), ...
            'LineWidth', 1.6, ...
            'MarkerSize', 6, ...
            'DisplayName', cfg.decoderNames{d});
    end
    title(ax, sprintf('%s, I_{max} = %d', ttl, cfg.baseIter), 'FontSize', 11, 'FontWeight', 'bold');
    lg = legend(ax, 'Location', 'best'); set(lg, 'FontSize', 9, 'Box', 'on');
    save_fig(fig, fullfile(figDir, sprintf('fig7_avgiter_decoders_%s.png', sfx)));
end
fprintf('Графики: %s\n', figDir);
end

% -------------------------------------------------------------------------
function r = find_run(runs, bgn, alg, it)
r = [];
for i = 1:numel(runs)
    m = runs{i}.meta;
    if m.bgn == bgn && strcmp(m.algorithm, alg) && m.maxIter == it, r = runs{i}; return; end
end
end

function [x, y] = curve(res, field)
ok = res.numFrames > 0;
x = res.EbN0dB(ok); y = res.(field)(ok);
if any(strcmp(field, {'ber', 'fer'})), y(y == 0) = NaN; end
end


function fig = newfig()
% Единый размер и белый фон, как в matvey_qpsk.
fig = figure('Visible', 'off', ...
    'Color', 'w', ...
    'Position', [100 100 900 600]);
end

function ax = newax(fig, ylab)
% Единое оформление осей.
ax = axes('Parent', fig);

hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

set(ax, ...
    'YScale', 'log', ...
    'FontSize', 11);

xlabel(ax, 'E_b/N_0, дБ');
ylabel(ax, ylab);
end

function save_fig(fig, file)
% Сохранение графика в PNG и FIG.

% PNG для отчёта
print(fig, file, '-dpng', '-r200');

% FIG для последующего редактирования в MATLAB
figFile = strrep(file, '.png', '.fig');
savefig(fig, figFile);

close(fig);
end
