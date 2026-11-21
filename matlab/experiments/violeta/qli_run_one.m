function item = qli_run_one(s, cfg, cacheDir)
%QLI_RUN_ONE Один прогон ldpc_simulate с кэшем (отдельная функция — чтобы save работал в parfor).
p = qli_spec_params(s, cfg);
key = sprintf('%s_bg%d_Zc%d_%s_it%d_a%g_b%g', cfg.modulation, s.bgn, s.Zc, s.algorithm, ...
    s.maxIter, s.scaling, s.offset);
if strcmp(s.algorithm, 'uncoded'), key = sprintf('%s_uncoded_bg%d', cfg.modulation, s.bgn); end
file = fullfile(cacheDir, [key '.mat']);
if cfg.resume && exist(file, 'file') == 2
    S = load(file);
    if isfield(S, 'item')
        pa = S.item.res.params; pa.verbose = p.verbose;
        if isequal(pa, p)
            item = S.item; fprintf('[кэш] %s\n', key); return
        end
    end
end
fprintf('[расчёт] %s\n', key);
res = ldpc_simulate(p);
item.meta = s; item.res = res;
save(file, 'item');
end
