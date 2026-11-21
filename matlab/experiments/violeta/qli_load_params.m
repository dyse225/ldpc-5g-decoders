function [alpha, beta] = qli_load_params(cfg, bgn)
%QLI_LOAD_PARAMS Оптимальные alpha/beta: из cfg.alphaOpt/betaOpt, иначе из results/opt_params.mat,
%   иначе значения по умолчанию (с предупреждением).
f = sprintf('bg%d', bgn);
alpha = []; beta = [];
if ~isempty(cfg.alphaOpt), alpha = cfg.alphaOpt.(f); end
if ~isempty(cfg.betaOpt),  beta = cfg.betaOpt.(f); end
file = fullfile(cfg.resDir, 'opt_params.mat');
if (isempty(alpha) || isempty(beta)) && exist(file, 'file') == 2
    S = load(file);
    if isempty(alpha) && isfield(S.opt, f), alpha = S.opt.(f).alpha; end
    if isempty(beta)  && isfield(S.opt, f), beta = S.opt.(f).beta; end
end
if isempty(alpha)
    alpha = cfg.defaultAlpha;
    warning('qli:alpha', 'BG%d: оптимальный alpha не найден, берём %g (сначала qli_find_params)', bgn, alpha);
end
if isempty(beta)
    beta = cfg.defaultBeta;
    warning('qli:beta', 'BG%d: оптимальный beta не найден, берём %g (сначала qli_find_params)', bgn, beta);
end
end
