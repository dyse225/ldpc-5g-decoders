function [uHat, iters] = ldpc_wbf_decode(llr, H, K, maxFlips, alpha, ignoreErasures)
%LDPC_WBF_DECODE Взвешенный bit-flipping декодер (WBF / MWBF).
%   [uHat, iters] = LDPC_WBF_DECODE(llr, H, K, maxFlips, alpha, ignoreErasures)
%     llr      — N x B мягкие значения (LLR, > 0 -> бит 0)
%     H        — M x N проверочная матрица (sparse)
%     K        — число информационных бит (первые K позиций слова)
%     maxFlips — максимум итераций; за итерацию инвертируется ОДИН бит
%     alpha    — коэффициент MWBF; 0 -> обычный WBF (по умолчанию 0)
%     ignoreErasures — true (по умолчанию): биты с LLR = 0 (выколотые)
%                не участвуют в расчёте веса проверки
%   uHat  — K x B int8 жёсткие решения (итоговое состояние декодера)
%   iters — 1 x B число выполненных итераций
%
%   Порт Decoder_LDPC_bit_flipping (BF_decode) и
%   Decoder_LDPC_bit_flipping_OMWBF (BF_process) из AFF3CT (MIT License):
%   https://github.com/aff3ct/aff3ct,
%   src/Module/Decoder/LDPC/BF/Decoder_LDPC_bit_flipping.cpp,
%   src/Module/Decoder/LDPC/BF/OMWBF/Decoder_LDPC_bit_flipping_OMWBF.cpp.
%   Алгоритм: Kou, Lin, Fossorier (2001) — WBF; Zhang, Fossorier (2004) — MWBF.
%
%   Алгоритм (как в AFF3CT):
%     1. вес проверки  Ymin(m) = min |y(n)| по битам проверки m;
%     2. жёсткое решение z = (y < 0);
%     3. итерация: синдром s = H*z mod 2,
%        энергия бита E(n) = сумма по его проверкам (2*s(m) - 1)*Ymin(m)
%                            - alpha*|y(n)|,
%        инвертируется один бит с максимальной E;
%     4. повтор до нулевого синдрома или лимита итераций, на выходе — z.
%
%   Отличия от оригинала:
%     - ignoreErasures: в AFF3CT минимум берётся по всем битам, и любой
%       выколотый бит (|y| = 0) обнуляет вес проверки; у 5G NR выколотые
%       столбцы входят во все проверки, и WBF перестаёт работать;
%     - обрабатывается сразу пачка из B кадров.

if nargin < 5, alpha = 0; end
if nargin < 6, ignoreErasures = true; end

[M, N] = size(H);
B = size(llr, 2);
Hd = double(H);
HdT = Hd.';
a = abs(llr);

% --- веса проверок Ymin (M x B)
[rowIdx, colIdx] = find(H);
vals = a(colIdx, :);
if ignoreErasures
    vals(vals == 0) = Inf;
end
Ymin = zeros(M, B);
for b = 1:B
    Ymin(:, b) = accumarray(rowIdx, vals(:, b), [M 1], @min, Inf);
end
Ymin(isinf(Ymin)) = 0;

% --- декодирование
z = double(llr < 0);
iters = zeros(1, B);
active = true(1, B);
for it = 1:maxFlips
    idx = find(active);
    s = mod(Hd * z(:, idx), 2);
    ok = ~any(s, 1);
    active(idx(ok)) = false;
    idx = idx(~ok);
    if isempty(idx)
        break
    end
    s = s(:, ~ok);

    E = HdT * ((2*s - 1) .* Ymin(:, idx)) - alpha * a(:, idx);
    [~, nMax] = max(E, [], 1);                 % один бит на кадр
    lin = sub2ind([N, B], nMax, idx);
    z(lin) = 1 - z(lin);
    iters(idx) = it;
end
iters(active) = maxFlips;

uHat = int8(z(1:K, :));
end
