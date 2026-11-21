function runs = qli_run(cfg)
%QLI_RUN Основной свип: BG x декодер x I_max + некодированная QPSK.
%   Для NMS и OMS берутся оптимальные alpha/beta (qli_find_params / cfg.alphaOpt / cfg.betaOpt).
%   Декодеры из cfg.iterDecoders считаются при всех I_max из cfg.iterList, остальные — только
%   при cfg.baseIter. Выход: results/all_results.mat, raw_results.csv, summary.csv.
if ~exist(cfg.resDir, 'dir'), mkdir(cfg.resDir); end
cacheDir = fullfile(cfg.resDir, 'cache');
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end

specs = {};
for bgn = cfg.bgList
    f = sprintf('bg%d', bgn);
    [alpha, beta] = qli_load_params(cfg, bgn);
    base.bgn = bgn; base.Zc = cfg.Zc.(f); base.EbN0dB = cfg.EbN0dB.(f);
    base.maxFrames = cfg.maxFrames.(f); base.minFrameErrors = cfg.minFrameErrors;
    base.scaling = alpha; base.offset = beta;
    fprintf('BG%d: alpha = %g, beta = %g\n', bgn, alpha, beta);

    % Некодированная QPSK: отдельная сетка Eb/N0 от 0 до 11 дБ
    s = base;
    s.algorithm = 'uncoded';
    s.maxIter = 0;

    s.EbN0dB = cfg.EbN0uncoded;
    s.maxFrames = cfg.uncodedMaxFrames;
    s.minFrameErrors = 200;

    s.scaling = 0;
    s.offset = 0;

    specs{end+1} = s;
    for d = 1:numel(cfg.decoders)
        alg = cfg.decoders{d};
        if any(strcmp(alg, cfg.iterDecoders)), its = cfg.iterList; else, its = cfg.baseIter; end
        for it = its
            s = base; s.algorithm = alg; s.maxIter = it;
            if ~strcmp(alg, 'norm-min-sum'), s.scaling = cfg.defaultAlpha; end   % не влияет на декодер,
            if ~strcmp(alg, 'offset-min-sum'), s.offset = cfg.defaultBeta; end   % но задаёт ключ кэша
            specs{end+1} = s; %#ok<AGROW>
        end
    end
end

n = numel(specs); fprintf('Прогонов: %d\n', n);
runs = cell(1, n);
if cfg.parallel
    parfor i = 1:n, runs{i} = qli_run_one(specs{i}, cfg, cacheDir); end
else
    for i = 1:n
        fprintf('--- %d / %d ---\n', i, n);
        runs{i} = qli_run_one(specs{i}, cfg, cacheDir);
    end
end
save(fullfile(cfg.resDir, 'all_results.mat'), 'runs', 'cfg');

% raw_results.csv и summary.csv
rows = {}; sm = {};
for i = 1:n
    m = runs{i}.meta; r = runs{i}.res;
    for k = 1:numel(r.EbN0dB)
        if r.numFrames(k) == 0, continue; end
        rows(end+1, :) = {m.bgn, m.Zc, r.K, r.E, r.K / r.E, cfg.modulation, m.algorithm, m.maxIter, ...
            m.scaling, m.offset, r.EbN0dB(k), r.ber(k), r.fer(k), r.avgIter(k), r.numFrames(k), ...
            r.numBitErrors(k), r.numFrameErrors(k)}; %#ok<AGROW>
    end
    if ~strcmp(m.algorithm, 'uncoded')
        sm(end+1, :) = {m.bgn, r.K, m.algorithm, m.maxIter, qli_ebno_at(r, 'fer', 1e-2), ...
            qli_ebno_at(r, 'fer', 1e-3), qli_ebno_at(r, 'ber', 1e-4)}; %#ok<AGROW>
    end
end
qli_csv(fullfile(cfg.resDir, 'raw_results.csv'), {'BG','Zc','K','E','rate','modulation','algorithm', ...
    'maxIter','alpha','beta','EbN0dB','BER','FER','avgIter','numFrames','numBitErrors', ...
    'numFrameErrors'}, rows);
qli_csv(fullfile(cfg.resDir, 'summary.csv'), {'BG','K','algorithm','maxIter','EbN0dB_at_FER1e-2', ...
    'EbN0dB_at_FER1e-3','EbN0dB_at_BER1e-4'}, sm);
end
