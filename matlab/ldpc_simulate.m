function res = ldpc_simulate(p)
%LDPC_SIMULATE Монте-Карло моделирование 5G NR LDPC линии (AWGN или Rayleigh).
%   res = LDPC_SIMULATE(p), p — структура из LDPC_DEFAULT_PARAMS.
%
%   Цепочка: источник -> nrLDPCEncode -> согласование скорости -> модулятор
%            -> awgn -> LLR-демодулятор -> восстановление скорости
%            -> nrLDPCDecode -> подсчёт BER / FER / числа итераций.
%
%   Согласование скорости упрощённое (rv = 0, без filler-бит и без битового
%   перемежения): в канал идут первые E = K/R бит выхода nrLDPCEncode
%   (первые 2*Zc систематических бит кодер уже выколол), на приёме
%   непереданные биты получают LLR = 0.
%
%   Канал p.channel:
%     'awgn'     — только шум;
%     'rayleigh' — плоские замирания, независимые для каждого символа
%                  (идеальное перемежение), h ~ CN(0,1), E|h|^2 = 1;
%                  приёмник знает h точно: z = y/h, дисперсия шума
%                  символа nVar/|h|^2 передаётся в демодулятор.
%
%   res — структура с полями EbN0dB, ber, fer, avgIter, numFrames,
%   numBitErrors, numFrameErrors, timePerFrame, codeName, K, N, E, params.

isCoded = ~strcmpi(p.code, 'none');
bps = bits_per_symbol(p.modulation);

if isCoded
    if p.bgn == 1
        K = 22 * p.Zc; N = 66 * p.Zc;       % N — длина выхода nrLDPCEncode
    else
        K = 10 * p.Zc; N = 50 * p.Zc;
    end
    E = round(K / p.rate);                  % бит в канал
    if E > N
        error('ldpc_simulate:rate', 'R = %g ниже материнской скорости кода', p.rate);
    end
    codeName = sprintf('5G NR BG%d Zc=%d K=%d R=%g', p.bgn, p.Zc, K, K/E);
    isBitFlip = any(strcmp(p.algorithm, {'bit-flip', 'wbf'}));
    if isBitFlip
        Hbf = nr_ldpc_pcm(p.bgn, p.Zc, E);   % свой декодер работает с явной H
    else
        decArgs = decoder_args(p);
    end
else
    K = 1040; N = K; E = K;                 % делится на 1, 2, 4
    codeName = 'uncoded';
end
R = K / E;
if mod(E, bps) ~= 0
    error('ldpc_simulate:frameLen', 'E = %d не кратно log2(M) = %d', E, bps);
end

nPts = numel(p.EbN0dB);
res.EbN0dB         = p.EbN0dB(:).';
res.ber            = nan(1, nPts);
res.fer            = nan(1, nPts);
res.avgIter        = nan(1, nPts);
res.numFrames      = zeros(1, nPts);
res.numBitErrors   = zeros(1, nPts);
res.numFrameErrors = zeros(1, nPts);
res.timePerFrame   = nan(1, nPts);
res.codeName       = codeName;
res.K = K;
res.N = N;
res.E = E;
res.params = p;

if p.verbose
    fprintf('%s | %s | %s, maxIter=%d\n', codeName, upper(p.modulation), ...
        p.algorithm, p.maxIter);
end

for iPt = 1:nPts
    ebno = p.EbN0dB(iPt);
    snrdB = convertSNR(ebno, "ebno", "snr", BitsPerSymbol=bps, CodingRate=R);
    stream = RandStream('mt19937ar', Seed=p.seed + iPt);

    nFrames = 0; nBitErr = 0; nFrameErr = 0; sumIter = 0; tDec = 0;
    while nFrameErr < p.minFrameErrors && nFrames < p.maxFrames
        B = min(p.batchSize, p.maxFrames - nFrames);

        u = randi(stream, [0 1], K, B, 'int8');
        if isCoded
            cw = nrLDPCEncode(u, p.bgn);
            c = cw(1:E, :);                 % согласование скорости, rv = 0
        else
            c = u;
        end

        x = modulate(c(:), p.modulation);
        % мощность сигнала до замираний 0 дБВт (единичная мощность созвездия),
        % E|h|^2 = 1, поэтому SNR задаётся относительно 0 дБВт в обоих каналах
        if strcmpi(p.channel, 'rayleigh')
            h = (randn(stream, size(x)) + 1i*randn(stream, size(x))) / sqrt(2);
            [y, nVar] = awgn(h .* x, snrdB, 0, stream);
            llr = demodulate(y ./ h, p.modulation, nVar ./ abs(h).^2);
        else
            [y, nVar] = awgn(x, snrdB, 0, stream);
            llr = demodulate(y, p.modulation, nVar);
        end
        llr = reshape(llr, E, B);

        if isCoded
            t0 = tic;
            if isBitFlip
                % слово для Hbf: [2*Zc выколотых бит; E переданных]
                llrBf = [zeros(2*p.Zc, B); llr];
                if strcmp(p.algorithm, 'wbf')
                    [uHat, iters] = ldpc_wbf_decode(llrBf, Hbf, K, p.wbfMaxFlips, p.wbfAlpha, true);
                else
                    [uHat, iters] = ldpc_bitflip_decode(llrBf, Hbf, K, p.maxIter, ...
                        p.bfRule, p.bfOnFail, p.bfThreshold, p.bfMaxFlips);
                end
            else
                llrFull = [llr; zeros(N - E, B)];  % непереданные биты: LLR = 0
                [uHat, iters] = nrLDPCDecode(llrFull, p.bgn, p.maxIter, decArgs{:});
            end
            tDec = tDec + toc(t0);
            sumIter = sumIter + sum(iters);
        else
            uHat = int8(llr < 0);
        end

        errs = sum(uHat ~= u, 1);
        nBitErr   = nBitErr + sum(errs);
        nFrameErr = nFrameErr + nnz(errs);
        nFrames   = nFrames + B;
    end

    res.numFrames(iPt)      = nFrames;
    res.numBitErrors(iPt)   = nBitErr;
    res.numFrameErrors(iPt) = nFrameErr;
    res.ber(iPt) = nBitErr / (nFrames * K);
    res.fer(iPt) = nFrameErr / nFrames;
    if isCoded
        res.avgIter(iPt)      = sumIter / nFrames;
        res.timePerFrame(iPt) = tDec / nFrames;
    end

    if p.verbose
        fprintf('  Eb/N0 = %5.2f dB: BER = %.3e, FER = %.3e, iter = %5.2f, frames = %d\n', ...
            ebno, res.ber(iPt), res.fer(iPt), res.avgIter(iPt), nFrames);
    end
    if p.stopAtZeroFer && nFrameErr == 0
        break
    end
end
end

% -------------------------------------------------------------------------
function args = decoder_args(p)
switch p.algorithm
    case 'bp',             alg = 'Belief propagation';
    case 'layered-bp',     alg = 'Layered belief propagation';
    case 'norm-min-sum',   alg = 'Normalized min-sum';
    case 'offset-min-sum', alg = 'Offset min-sum';
    otherwise, error('ldpc_simulate:alg', 'Неизвестный алгоритм: %s', p.algorithm);
end
args = {'Algorithm', alg, 'Termination', p.termination};
switch p.algorithm
    case 'norm-min-sum'
        args = [args, {'ScalingFactor', p.scaling}];
    case 'offset-min-sum'
        args = [args, {'Offset', p.offset}];
end
end

function k = bits_per_symbol(mod)
switch lower(mod)
    case 'bpsk',  k = 1;
    case 'qpsk',  k = 2;
    case '16qam', k = 4;
    case '64qam', k = 6;
    otherwise, error('ldpc_simulate:mod', 'Неизвестная модуляция: %s', mod);
end
end

function x = modulate(bits, mod)
% Gray, единичная средняя мощность
k = bits_per_symbol(mod);
if k == 1
    x = pskmod(bits, 2, 0, 'gray', InputType='bit');
else
    x = qammod(bits, 2^k, 'gray', InputType='bit', UnitAveragePower=true);
end
end

function llr = demodulate(y, mod, nVar)
% Точные LLR: log(P(b=0)/P(b=1)), положительное -> бит 0 (соглашение nrLDPCDecode).
% nVar — скаляр или вектор (своя дисперсия для каждого символа). Вектор
% демодуляторы принимают только для символов, записанных строкой, поэтому
% в этом случае подаём строку и разворачиваем результат (биты символа подряд).
k = bits_per_symbol(mod);
if ~isscalar(nVar)
    y = y(:).'; nVar = nVar(:).';
end
if k == 1
    llr = pskdemod(y, 2, 0, 'gray', OutputType='llr', NoiseVariance=nVar);
else
    llr = qamdemod(y, 2^k, 'gray', OutputType='llr', ...
        UnitAveragePower=true, NoiseVariance=nVar);
end
llr = llr(:);
end
