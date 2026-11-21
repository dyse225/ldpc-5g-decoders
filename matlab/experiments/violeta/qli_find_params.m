function opt = qli_find_params(cfg)
%QLI_FIND_PARAMS Поиск оптимальных alpha (NMS) и beta (OMS) для каждого BG.
%   Сканируется сетка cfg.param.alphaList / betaList при I_max = cfg.param.iter в рабочих
%   точках cfg.param.EbN0dB. Критерий: средний log10(FER) по точкам; чтобы не выбрать
%   случайный шумовой минимум, берётся минимум скользящего среднего по трём соседним значениям.
%   Все значения считаются на одних и тех же кадрах (один seed).
%   Выход: results/opt_params.mat, param_sweep.csv, param_optimum.csv, fig0_params_BG*.png
pc = cfg.param;
if ~exist(cfg.resDir, 'dir'), mkdir(cfg.resDir); end
cacheDir = fullfile(cfg.resDir, 'cache_param');
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end

opt = struct(); rowsAll = {}; rowsOpt = {};
for bgn = cfg.bgList
    f = sprintf('bg%d', bgn);
    Zc = cfg.Zc.(f); eb = pc.EbN0dB.(f);
    % список прогонов
    specs = {};
    for a = pc.alphaList, specs{end+1} = mk(bgn, Zc, 'norm-min-sum', pc, eb, a, cfg.defaultBeta); end %#ok<AGROW>
    for b = pc.betaList,  specs{end+1} = mk(bgn, Zc, 'offset-min-sum', pc, eb, cfg.defaultAlpha, b); end %#ok<AGROW>
    n = numel(specs); out = cell(1, n);
    if cfg.parallel
        parfor i = 1:n, out{i} = qli_run_one(specs{i}, cfg, cacheDir); end
    else
        for i = 1:n, out{i} = qli_run_one(specs{i}, cfg, cacheDir); end
    end
    nA = numel(pc.alphaList);
    [alphaBest, scA, ferA, plA] = best_of(out(1:nA), pc.alphaList, 'alpha');
    [betaBest,  scB, ferB, plB] = best_of(out(nA+1:end), pc.betaList, 'beta');
    opt.(f).alpha = alphaBest; opt.(f).beta = betaBest;
    opt.(f).alphaPlateau = plA; opt.(f).betaPlateau = plB;

    fprintf('\nBG%d (K = %d): оптимальный alpha = %g (плато %s), beta = %g (плато %s)\n', bgn, ...
        out{1}.res.K, alphaBest, mat2str(plA), betaBest, mat2str(plB));
    fprintf('  FER при alpha_opt: %s,  при beta_opt: %s  (Eb/N0 = %s дБ)\n', ...
        mat2str(ferA(:, pc.alphaList == alphaBest).', 3), mat2str(ferB(:, pc.betaList == betaBest).', 3), mat2str(eb));
    if max([ferA(:); ferB(:)]) < 1e-3 || min([ferA(:); ferB(:)]) > 0.5
        warning('qli:param', ['BG%d: FER в рабочих точках %s дБ вне 1e-3...0.5, оптимум ненадёжен. ' ...
            'Сдвиньте cfg.param.EbN0dB.'], bgn, mat2str(eb));
    end

    for i = 1:n
        r = out{i}.res; m = out{i}.meta;
        if strcmp(m.algorithm, 'norm-min-sum'), pn = 'alpha'; pv = m.scaling; else, pn = 'beta'; pv = m.offset; end
        for k = 1:numel(r.EbN0dB)
            if r.numFrames(k) == 0, continue; end
            rowsAll(end+1, :) = {bgn, r.K, m.algorithm, m.maxIter, pn, pv, r.EbN0dB(k), r.ber(k), ...
                r.fer(k), r.avgIter(k), r.numFrames(k), r.numFrameErrors(k)}; %#ok<AGROW>
        end
    end
    rowsOpt(end+1, :) = {bgn, out{1}.res.K, 'norm-min-sum', pc.iter, 'alpha', alphaBest, strjoin(arrayfun(@(x) sprintf('%g', x), plA, 'UniformOutput', false), ';')}; %#ok<AGROW>
    rowsOpt(end+1, :) = {bgn, out{1}.res.K, 'offset-min-sum', pc.iter, 'beta', betaBest, strjoin(arrayfun(@(x) sprintf('%g', x), plB, 'UniformOutput', false), ';')}; %#ok<AGROW>

    plot_params(cfg, bgn, out{1}.res.K, eb, pc.alphaList, ferA, alphaBest, pc.betaList, ferB, betaBest);
end
save(fullfile(cfg.resDir, 'opt_params.mat'), 'opt');
qli_csv(fullfile(cfg.resDir, 'param_sweep.csv'), {'BG','K','algorithm','maxIter','param','value', ...
    'EbN0dB','BER','FER','avgIter','numFrames','numFrameErrors'}, rowsAll);
qli_csv(fullfile(cfg.resDir, 'param_optimum.csv'), {'BG','K','algorithm','maxIter','param', ...
    'best_value','plateau_within_1sigma'}, rowsOpt);
end

% -------------------------------------------------------------------------
function s = mk(bgn, Zc, alg, pc, eb, a, b)
s.bgn = bgn; s.Zc = Zc; s.algorithm = alg; s.maxIter = pc.iter;
s.scaling = a; s.offset = b; s.EbN0dB = eb;
s.maxFrames = pc.maxFrames; s.minFrameErrors = pc.minFrameErrors;
end

function [best, score, fer, plateau] = best_of(out, vals, ~)
nE = numel(out{1}.res.EbN0dB);
fer = nan(nE, numel(out)); nf = fer; score = nan(1, numel(out));
for i = 1:numel(out)
    r = out{i}.res;
    fer(:, i) = r.fer(:); nf(:, i) = r.numFrames(:);
end
lf = log10(max(fer, 1 ./ max(nf, 1)));          % FER = 0 -> нижняя оценка 1/nFrames
raw = mean(lf, 1);
sm = raw;                                        % скользящее среднее по 3 соседям
for i = 1:numel(raw)
    j = max(1, i-1):min(numel(raw), i+1); sm(i) = mean(raw(j));
end
[~, ib] = min(sm); best = vals(ib);
score = sm;
% плато: значения, у которых средний FER в пределах 1 sigma от лучшего
mf = mean(fer, 1);
sig = sqrt(max(mf(ib) * (1 - mf(ib)), 0) / max(mean(nf(:, ib)) * size(fer, 1), 1));
plateau = vals(mf <= mf(ib) + sig);
end

function plot_params(cfg, bgn, K, eb, aL, ferA, aBest, bL, ferB, bBest)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 480]);
data = {aL, ferA, aBest, '\alpha (NMS)'; bL, ferB, bBest, '\beta (OMS)'};
for k = 1:2
    ax = subplot(1, 2, k, 'Parent', fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    set(ax, 'FontSize', 11, 'LineWidth', 1.0, 'TickDir', 'out', ...
        'XMinorGrid', 'on', 'YMinorGrid', 'on');
    F = data{k, 2}; F(F == 0) = NaN;
    for j = 1:size(F, 1)
        plot(ax, data{k, 1}, F(j, :), '-o', 'LineWidth', 1.7, 'MarkerSize', 6, ...
            'MarkerFaceColor', 'w', 'DisplayName', sprintf('E_b/N_0 = %g дБ', eb(j)));
    end
    set(ax, 'YScale', 'log'); xlabel(ax, data{k, 4}, 'FontSize', 11); ylabel(ax, 'FER', 'FontSize', 11);
    yl = ylim(ax);
    plot(ax, data{k, 3} * [1 1], yl, '--', 'Color', [0.4 0.4 0.4], 'DisplayName', 'оптимум');
    title(ax, sprintf('BG%d, K = %d, I_{max} = %d', bgn, K, cfg.param.iter), ...
        'FontSize', 11, 'FontWeight', 'bold');
    if k == 1, lg = legend(ax, 'Location', 'best'); set(lg, 'FontSize', 9, 'Box', 'on'); end
end
print(fig, fullfile(cfg.resDir, sprintf('fig0_params_BG%d.png', bgn)), '-dpng', '-r200');
close(fig);
end
