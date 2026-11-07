function [uHat, iters] = ldpc_bitflip_decode(llr, H, K, maxIter, rule, onFail, thr, kMax)
%LDPC_BITFLIP_DECODE Параллельный жёсткий bit-flipping декодер (PHBF).
%   [uHat, iters] = LDPC_BITFLIP_DECODE(llr, H, K, maxIter, rule, onFail, thr, kMax)
%     llr     — N x B мягкие значения (используется только знак;
%               LLR = 0, т.е. выколотые биты, стартуют как 0)
%     H       — M x N проверочная матрица (sparse)
%     K       — число информационных бит (первые K позиций слова)
%     maxIter — максимум итераций
%     rule    — правило выбора инвертируемых бит:
%               'phbf'     (по умолчанию) как в оригинале: минимум
%                          F = сумма(1 - 2*S) по проверкам бита;
%               'gallager' максимум числа невыполненных проверок бита.
%               Для регулярных кодов (как EG-LDPC в оригинале) правила
%               совпадают; для нерегулярного 5G NR — нет.
%     onFail  — что выдавать, если синдром не обнулился:
%               'channel' (по умолчанию) как в оригинале — исходное
%                         жёсткое решение канала;
%               'final'   итоговое состояние декодера (часть бит уже
%                         исправлена, FER по информационным битам ниже).
%     thr     — (только 'gallager') порог инвертирования относительно
%               максимума: инвертируются биты, у которых число невыполненных
%               проверок не меньше max - thr (и больше нуля); 0 (по
%               умолчанию) — только биты с максимумом, как в оригинале.
%               Фиксированный порог (U >= T) и порог по большинству
%               проверок (Sipser–Spielman) на 5G NR не сходятся: FER = 1.
%     kMax    — (только 'gallager') не больше kMax инвертируемых бит за
%               итерацию, берутся с наибольшим числом невыполненных
%               проверок; Inf (по умолчанию) — без ограничения.
%               Если за итерацию нечего инвертировать, декодирование
%               кадра прекращается.
%   uHat  — K x B int8 жёсткие решения
%   iters — 1 x B число выполненных итераций
%
%   Порт функции PH_bitflipping из пакета bf-tests-awgn-egldpc 0.2
%   (Fernando Pujaico Rivera, GPL-3), src/extras/simula.h:
%   https://github.com/Lcrypto/LDPC-Iterative-Bit-Flipping-family-decoders
%   (каталог bf-tests-awgn-egldpc_0.2.orig).
%
%   Отличия от оригинала:
%     - оригинал работает с циклическими EG-LDPC кодами (H задаётся
%       сдвигами «линий-образцов»), здесь — произвольная разреженная H;
%     - лимит итераций задаётся числом maxIter, а не долей p*N;
%     - обрабатывается сразу пачка из B кадров (векторизация по столбцам);
%     - добавлено правило 'gallager' для нерегулярных кодов;
%     - добавлен режим выдачи onFail = 'final';
%     - добавлены порог thr и ограничение kMax.
%
%   Алгоритм (как в оригинале):
%     1. r = жёсткое решение канала, синдром S = H*r mod 2;
%     2. вес бита F(n) = сумма по его проверкам (1 - 2*S(m)) —
%        «выполненные минус невыполненные»;
%     3. инвертируются все биты с минимальным F;
%     4. повтор до нулевого синдрома или лимита итераций;
%     5. если синдром так и не обнулился — возвращается исходное
%        жёсткое решение канала.

if nargin < 5, rule = 'phbf'; end
if nargin < 6, onFail = 'channel'; end
if nargin < 7, thr = 0; end
if nargin < 8, kMax = Inf; end
Hd = double(H);
HdT = Hd.';
r0 = double(llr < 0);           % жёсткое решение канала (HardDecision)
r = r0;
B = size(r, 2);
iters = zeros(1, B);
active = true(1, B);
failed = false(1, B);           % остановились, не обнулив синдром

for it = 1:maxIter
    idx = find(active);
    S = mod(Hd * r(:, idx), 2);                 % GeraSindrome
    ok = ~any(S, 1);
    active(idx(ok)) = false;
    idx = idx(~ok);
    if isempty(idx)
        break
    end
    S = S(:, ~ok);

    if strcmp(rule, 'gallager')
        U = HdT * S;                            % невыполненные проверки бита
        flip = U >= max(U, [], 1) - thr & U > 0;
        if isfinite(kMax)
            % оставляем kMax бит с наибольшим U (при равенстве — меньший индекс)
            Uc = U; Uc(~flip) = -Inf;
            [~, order] = sort(Uc, 1, 'descend');
            keep = false(size(U));
            n = min(kMax, size(U, 1));
            cols = repmat(1:size(U, 2), n, 1);
            keep(sub2ind(size(U), order(1:n, :), cols)) = true;
            flip = flip & keep;
        end
        stuck = ~any(flip, 1);                  % нечего инвертировать
        if any(stuck)
            active(idx(stuck)) = false;
            failed(idx(stuck)) = true;
            iters(idx(stuck)) = it;
        end
    else
        F = HdT * (1 - 2*S);                    % веса синдрома
        flip = F == min(F, [], 1);              % все биты с минимальным весом
    end
    r(:, idx) = xor(r(:, idx), flip);
    iters(idx) = it;
end

% проверка после последней инверсии
idx = find(active);
if ~isempty(idx)
    ok = ~any(mod(Hd * r(:, idx), 2), 1);
    active(idx(ok)) = false;
end
iters(active) = maxIter;
failed = failed | active;
if strcmp(onFail, 'channel')
    r(:, failed) = r0(:, failed);   % неудача: исходное решение канала
end

uHat = int8(r(1:K, :));
end
