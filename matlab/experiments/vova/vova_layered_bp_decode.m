function [uHat, iters, Lout] = vova_layered_bp_decode(llr, H, Zc, K, maxIter, order, termination)
%VOVA_LAYERED_BP_DECODE Layered belief propagation (sum-product) для 5G NR LDPC
%   с явно заданным порядком обхода слоёв.
%
%   [uHat, iters] = VOVA_LAYERED_BP_DECODE(llr, H, Zc, K, maxIter, order, termination)
%     llr         - N x B мягкие значения, log(P(0)/P(1)); выколотые биты LLR = 0.
%                   N = size(H, 2): [2*Zc выколотых бит; E переданных]
%     H           - M x N проверочная матрица (sparse), например из NR_LDPC_PCM
%     Zc          - размер подъёма; слой = Zc строк H (одна строка базового графа)
%     K           - число информационных бит (первые K позиций слова)
%     maxIter     - максимум итераций (одна итерация = обход всех слоёв)
%     order       - порядок слоёв: вектор номеров слоёв 1..M/Zc
%                   (см. VOVA_LAYER_ORDER)
%     termination - 'early' (останов при нулевом синдроме, по умолчанию) | 'max'
%   uHat  - K x B int8 жёсткие решения
%   iters - 1 x B число выполненных итераций
%   Lout  - N x B апостериорные LLR
%
%   Внутри слоя строки не пересекаются по столбцам (циркулянты), поэтому
%   слой обновляется целиком за один шаг. Обновление проверочного узла -
%   точное правило sum-product в форме phi(x) = -log(tanh(x/2)):
%     Q   = L - R_old                          (сообщения бит -> проверка)
%     |R| = phi( sum_j phi(|Q_j|) - phi(|Q_i|) ),  sign(R) = prod_{j~=i} sign(Q_j)
%     L   = Q + R_new
%   Синдром проверяется после каждой полной итерации, как Termination="early"
%   в nrLDPCDecode.

if nargin < 7 || isempty(termination), termination = 'early'; end
earlyStop = strcmpi(termination, 'early');

[M, N] = size(H);
B = size(llr, 2);
if size(llr, 1) ~= N
    error('vova_layered_bp_decode:size', 'size(llr,1) = %d, а у H %d столбцов', size(llr, 1), N);
end
nLayers = M / Zc;
lay = layer_structure(H, Zc);

L = double(llr);
Rmsg = cell(1, nLayers);
for l = 1:nLayers
    Rmsg{l} = zeros(lay.deg(l), Zc, B);
end

iters  = zeros(1, B);
active = true(1, B);

for it = 1:maxIter
    act = find(active);
    if isempty(act), break; end
    La = L(:, act);
    nA = numel(act);
    for l = order(:).'
        idx = lay.idx{l};                 % deg x Zc номера столбцов
        d = lay.deg(l);
        Rold = Rmsg{l}(:, :, act);
        Q = reshape(La(idx(:), :), d, Zc, nA) - Rold;
        Rnew = check_update(Q);
        La(idx(:), :) = reshape(Q + Rnew, d * Zc, nA);
        Rmsg{l}(:, :, act) = Rnew;
    end
    L(:, act) = La;
    iters(act) = it;
    if earlyStop
        hard = double(La < 0);
        ok = ~any(mod(H * hard, 2), 1);
        active(act(ok)) = false;
    end
end

uHat = int8(L(1:K, :) < 0);
Lout = L;
end

% -------------------------------------------------------------------------
function R = check_update(Q)
% Q: deg x Zc x B. Точное правило sum-product без своего ребра.
xMin = 1e-12; xMax = 50;
a = min(max(abs(Q), xMin), xMax);
P = phi(a);
S = sum(P, 1);
mag = phi(min(max(S - P, xMin), xMax));
s = sign(Q); s(s == 0) = 1;
sgn = prod(s, 1) .* s;                    % произведение знаков без своего
R = sgn .* mag;
end

function y = phi(x)
y = -log(tanh(x / 2));
end

function lay = layer_structure(H, Zc)
% Для каждого слоя: deg (число единиц в строке) и idx (deg x Zc) - номера
% столбцов, по столбцу idx на каждую строку слоя.
persistent cache
key = sprintf('m%d_n%d_z%d_nnz%d', size(H, 1), size(H, 2), Zc, nnz(H));
if ~isempty(cache) && isfield(cache, key)
    lay = cache.(key);
    return
end
nLayers = size(H, 1) / Zc;
lay.deg = zeros(1, nLayers);
lay.idx = cell(1, nLayers);
for l = 1:nLayers
    rows = (l - 1) * Zc + (1:Zc);
    [cols, rLoc] = find(H(rows, :).');    % отсортировано по строке слоя
    d = numel(cols) / Zc;
    if any(accumarray(rLoc, 1) ~= d)
        error('vova_layered_bp_decode:deg', 'Слой %d: строки разной степени', l);
    end
    lay.deg(l) = d;
    lay.idx{l} = reshape(cols, d, Zc);
end
if isempty(cache), cache = struct(); end
cache.(key) = lay;
end
