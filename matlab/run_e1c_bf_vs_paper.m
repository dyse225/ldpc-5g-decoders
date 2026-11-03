%% E1c: проверка BF и BP на коде из статьи Kou, Lin, Fossorier (2001)
% Y. Kou, S. Lin, M. Fossorier, "Low-Density Parity-Check Codes Based on
% Finite Geometries: A Rediscovery and New Results", IEEE Trans. Inf.
% Theory, vol. 47, no. 7, pp. 2711-2736, Nov. 2001 — Fig. 1:
% type-I 2-D EG-LDPC (255,175), регулярный, вес строк и столбцов 16.
% По рисунку при BER = 1e-5: SPA ~3.7 дБ, BF ~5.7 дБ (разрыв ~2 дБ).
%
% Код циклический, H — циркулянт 255x255 из вектора инцидентности прямой
% EG(2, 2^4), не проходящей через начало координат. Канал и декодеры
% симметричны, поэтому передаётся нулевое кодовое слово.
clear; clc;
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');

% --- H: прямая {1 + t*alpha : t из GF(16)} в GF(256) = EG(2, 2^4)
m = 8; n = 2^m - 1;
alpha = gf(2, m);
pts = gf(1, m);                                  % t = 0
for k = 0:14
    pts(end+1) = 1 + alpha^(17*k) * alpha;       %#ok<SAGROW> % t = alpha^(17k)
end
expo = log(pts);                                 % показатели alpha^e
row0 = zeros(1, n); row0(double(expo) + 1) = 1;
H = false(n);
for i = 0:n-1
    H(i+1, :) = circshift(row0, i);
end
H = sparse(H);
K = n - gfrank(double(full(H)), 2);
fprintf('EG-LDPC (%d,%d), веса: столбцы %d, строки %d\n', n, K, ...
    full(max(sum(H,1))), full(max(sum(H,2))));
R = K / n;

% ldpcDecoderConfig требует столбцов больше, чем строк, а H квадратная
% (255 проверок, ранг 80 — избыточность важна для BP). Добавляем фиктивный
% бит в одну проверку и подаём его как известный ноль (большой LLR).
Hbp = [H, sparse(1, 1, true, n, 1)];
decBP = ldpcDecoderConfig(Hbp, 'bp');
EbN0 = 2:0.5:6;
algs = {'BP (SPA)', 'BF (PHBF)'};
ber = nan(2, numel(EbN0)); fer = ber;
for a = 1:2
    for i = 1:numel(EbN0)
        s = RandStream('mt19937ar', Seed=i);
        snr = convertSNR(EbN0(i), "ebno", "snr", BitsPerSymbol=1, CodingRate=R);
        be = 0; fe = 0; nb = 0; nf = 0;
        while be < 200 && nf < 2e5
            B = 256;
            x = ones(n*B, 1);                    % нулевое слово -> BPSK +1
            [y, nv] = awgn(complex(x), snr, 0, s);
            llr = reshape(pskdemod(y, 2, 0, 'gray', OutputType='llr', ...
                NoiseVariance=nv), n, B);
            % pskmod(0) = +1, поэтому LLR > 0 -> бит 0
            if a == 1
                c = ldpcDecode([llr; 1e3*ones(1, B)], decBP, 50, OutputFormat='whole');
                c = c(1:n, :);
            else
                c = bf_whole(llr, H, 50);
            end
            e = sum(c ~= 0, 1);
            be = be + sum(e); fe = fe + nnz(e); nb = nb + n*B; nf = nf + B;
        end
        ber(a, i) = be / nb; fer(a, i) = fe / nf;
        fprintf('%-10s Eb/N0 = %.1f dB: BER = %.2e, FER = %.2e\n', ...
            algs{a}, EbN0(i), ber(a, i), fer(a, i));
        if be == 0, break, end
    end
end
save(fullfile(resDir, 'e1c_bf_vs_paper.mat'), 'EbN0', 'ber', 'fer', 'algs', 'n', 'K');

fig = figure('Color', 'w');
semilogy(EbN0, ber(1,:), '-o', EbN0, ber(2,:), '-s', 'LineWidth', 1.3); hold on; grid on;
semilogy(EbN0, berawgn(EbN0, 'psk', 2, 'nondiff'), 'k-');
semilogy([3.7 5.7], [1e-5 1e-5], 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
xlabel('E_b/N_0, дБ'); ylabel('BER'); ylim([1e-6 1]);
legend([algs, {'BPSK без кодирования', 'Статья, BER 10^{-5}: SPA 3.7 / BF 5.7 дБ'}], ...
    'Location', 'southwest');
title(sprintf('EG-LDPC (%d,%d): сверка со статьёй Kou-Lin-Fossorier, рис. 1', n, K));
exportgraphics(fig, fullfile(resDir, 'e1c_bf_vs_paper.png'), 'Resolution', 150);

% -------------------------------------------------------------------------
function c = bf_whole(llr, H, maxIter)
% ldpc_bitflip_decode возвращает первые K бит; здесь нужно всё слово
N = size(H, 2);
c = ldpc_bitflip_decode(llr, H, N, maxIter, 'phbf');
end
