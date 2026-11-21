function ok = qli_check_llr()
%QLI_CHECK_LLR Проверка пересчёта Eb/N0 -> SNR (convertSNR) и LLR для QPSK.
%   Те же вызовы, что в ldpc_simulate: convertSNR(..., BitsPerSymbol=2, CodingRate=R),
%   awgn(x, snr, 0, stream), qamdemod(..., 'llr', UnitAveragePower, NoiseVariance=nVar).
%   Проверяется:
%     1) Es/N0 = Eb/N0 + 10*log10(bps*R): QPSK, R = 1/2 -> Es/N0 = Eb/N0; uncoded -> +3.01 дБ;
%     2) nVar из awgn = 10^(-SNR/10) = измеренной дисперсии шума;
%     3) LLR совпадает с формулой QPSK 2*sqrt(2)*Re(y)/nVar (для Re и Im) вплоть до порядка
%        битов, знак LLR > 0 <=> бит 0;
%     4) согласованность LLR: E[exp(-LLR) | бит = 0] = 1 (ловит ошибку масштаба в 2 раза);
%     5) BER некодированной QPSK = 0.5*erfc(sqrt(Eb/N0)).
ok = true;
bps = 2;

% ---- 1) convertSNR ----------------------------------------------------------
fprintf('1) Пересчёт Eb/N0 -> Es/N0 (convertSNR)\n');
for tc = [1/2 1]
    for eb = [0 3]
        s = convertSNR(eb, "ebno", "snr", BitsPerSymbol=bps, CodingRate=tc);
        ref = eb + 10*log10(bps*tc);
        good = abs(s - ref) < 1e-9; ok = ok && good;
        fprintf('   Eb/N0 = %g дБ, R = %g: SNR = %.4f дБ, ожидается %.4f  %s\n', eb, tc, s, ref, tf(good));
    end
end

% ---- 2)-4) LLR на кодированной рабочей точке (R = 1/2, Eb/N0 = 2 дБ) ---------
N = 4e5; R = 1/2; eb = 2;
snr = convertSNR(eb, "ebno", "snr", BitsPerSymbol=bps, CodingRate=R);
stream = RandStream('mt19937ar', Seed=1);
bits = randi(stream, [0 1], 2*N, 1, 'int8');
x = qammod(bits, 4, 'gray', InputType='bit', UnitAveragePower=true);
[y, nVar] = awgn(x, snr, 0, stream);

fprintf('2) Дисперсия шума\n');
nTh = 10^(-snr/10); nMeas = mean(abs(y - x).^2);
good = abs(nVar/nTh - 1) < 1e-9 && abs(nMeas/nTh - 1) < 0.02; ok = ok && good;
fprintf('   nVar (awgn) = %.5f, 10^(-SNR/10) = %.5f, измерено = %.5f  %s\n', nVar, nTh, nMeas, tf(good));

llr = qamdemod(y, 4, 'gray', OutputType='llr', UnitAveragePower=true, NoiseVariance=nVar);
llr = double(llr(:));
lI = 2*sqrt(2)*real(y)/nVar; lQ = 2*sqrt(2)*imag(y)/nVar;
cand = [lI, lQ, -lI, -lQ];                % кандидаты для каждого бита символа
cname = {'+Re(y)', '+Im(y)', '-Re(y)', '-Im(y)'};
fprintf('3) LLR против формулы 2*sqrt(2)*Re(y)/nVar (каждый бит независимо: ось и знак)\n');
axisUsed = zeros(1, 2);
for k = 1:2
    lk = llr(k:2:end);                    % k-й бит каждого символа
    err = max(abs(cand - lk), [], 1);
    [e, i] = min(err);
    axisUsed(k) = mod(i - 1, 2);          % 0 - ось I, 1 - ось Q
    good = e < 1e-6 * max(1, max(abs(lk))); ok = ok && good;
    fprintf('   бит %d символа: LLR = %s * 2*sqrt(2)/nVar, max|разность| = %.2e  %s\n', k, cname{i}, e, tf(good));
    if ~good
        r = median(abs(lk) ./ max(abs(lI), eps));
        fprintf('      отношение |LLR_matlab| / |формула| ~ %.3f (1 - норма; 0.5 - nVar трактуется по оси)\n', r);
    end
end
good = axisUsed(1) ~= axisUsed(2); ok = ok && good;
fprintf('   биты символа лежат на разных осях (I и Q)  %s\n', tf(good));
% знак относительно настоящих битов: LLR > 0 <=> бит 0
for k = 1:2
    hk = int8(llr(k:2:end) < 0); bk = bits(k:2:end);
    good = mean(hk ~= bk) < 0.25; ok = ok && good;
    fprintf('   знак бита %d: LLR > 0 <=> бит 0 (как ждёт nrLDPCDecode), сырой BER = %.3f  %s\n', k, mean(hk ~= bk), tf(good));
end

fprintf('4) Согласованность LLR\n');
b0 = bits == 0;
c = mean(exp(-llr(b0)));
good = abs(c - 1) < 0.1; ok = ok && good;
fprintf('   E[exp(-LLR) | бит=0] = %.3f (ожидается 1)  %s\n', c, tf(good));

% ---- 5) uncoded BER ---------------------------------------------------------
fprintf('5) BER некодированной QPSK (Gray)\n');
for eb = [2 4]
    s = convertSNR(eb, "ebno", "snr", BitsPerSymbol=bps, CodingRate=1);
    [yu, nu] = awgn(x, s, 0, stream);
    lu = qamdemod(yu, 4, 'gray', OutputType='llr', UnitAveragePower=true, NoiseVariance=nu);
    ber = mean(int8(lu(:) < 0) ~= bits);
    th = 0.5*erfc(sqrt(10^(eb/10)));
    good = abs(ber/th - 1) < 0.1; ok = ok && good;
    fprintf('   Eb/N0 = %g дБ: BER модели = %.3e, теория = %.3e  %s\n', eb, ber, th, tf(good));
end
if ok, fprintf('Проверка пройдена.\n'); else, fprintf('ЕСТЬ ПРОБЛЕМЫ: см. строки FAIL выше.\n'); end
end

function s = tf(b)
if b, s = 'OK'; else, s = 'FAIL'; end
end
