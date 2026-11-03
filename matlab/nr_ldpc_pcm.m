function H = nr_ldpc_pcm(bgn, Zc, E)
%NR_LDPC_PCM Проверочная матрица 5G NR LDPC под передаваемую длину E.
%   H = NR_LDPC_PCM(bgn, Zc, E) — sparse logical матрица для слова
%   [2*Zc выколотых систематических бит; первые E бит выхода nrLDPCEncode].
%
%   Полная H по TS 38.212 берётся из 5G Toolbox (внутренняя функция
%   nr5g.internal.ldpc.getParams — её API может измениться в новых версиях).
%   Строки расширения (после 4 базовых) имеют по одному проверочному биту,
%   поэтому для непереданных проверочных бит соответствующие строки и
%   столбцы просто отбрасываются: H = Hfull(1:nRows*Zc, 1:nCols*Zc).

if bgn == 1, Kb = 22; else, Kb = 10; end
if mod(E, Zc) ~= 0
    error('nr_ldpc_pcm:E', 'E = %d должно быть кратно Zc = %d', E, Zc);
end

persistent cache
key = sprintf('bg%d_z%d_e%d', bgn, Zc, E);
if isempty(cache), cache = struct(); end
if isfield(cache, key)
    H = cache.(key);
    return
end

cfg = nr5g.internal.ldpc.getParams(bgn, Zc);
Hfull = cfg.ParityCheckMatrix;
nCols = 2 + E / Zc;
nRows = nCols - Kb;
H = sparse(Hfull(1:nRows*Zc, 1:nCols*Zc));
cache.(key) = H;
end
