function [alphaOpt, betaOpt, T] = pick_params(matFile)
%PICK_PARAMS Лучшие alpha (NMS) и beta (OMS) для каждого кода и I_max из a_params.mat.
%   Критерий: минимум FER; если при данном I_max FER >= 0.9 (кадры почти не
%   исправляются и FER не различает параметры), то минимум BER.
%   alphaOpt(c, i), betaOpt(c, i): c - номер кода (1 - BG1, 2 - BG2), i - номер I_max.
%   T - таблица со всеми значениями для отчёта.
L = load(matFile);
nC = size(L.codes, 1); nI = numel(L.iterList);
alphaOpt = nan(nC, nI); betaOpt = nan(nC, nI);
rows = {};
for c = 1:nC
    for i = 1:nI
        for alg = ["norm-min-sum", "offset-min-sum"]
            sel = find(cellfun(@(x) x == c, L.jobs(:,1)) & cellfun(@(x) x == L.iterList(i), L.jobs(:,2)) ...
                & strcmp(L.jobs(:,3), alg));
            par = cell2mat(L.jobs(sel, 4));
            fer = cellfun(@(r) r.fer, L.out(sel)); ber = cellfun(@(r) r.ber, L.out(sel));
            if min(fer) >= 0.9, crit = ber; else, crit = fer; end
            [~, k] = min(crit);
            if alg == "norm-min-sum", alphaOpt(c, i) = par(k); else, betaOpt(c, i) = par(k); end
            rows(end+1, :) = {c, L.iterList(i), char(alg), par(k), fer(k), ber(k)}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'code', 'Imax', 'alg', 'best_param', 'FER', 'BER'});
end
