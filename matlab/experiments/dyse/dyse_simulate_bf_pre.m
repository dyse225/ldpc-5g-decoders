function res = dyse_simulate_bf_pre(p)
%DYSE_SIMULATE_BF_PRE BF с предобработкой входа мягким декодером.
%   res = DYSE_SIMULATE_BF_PRE(p)
%
%   Цепочка как в LDPC_SIMULATE (5G NR, rv = 0, AWGN), но перед жёстким
%   bit-flip декодером LLR прогоняются через несколько итераций мягкого
%   декодера ldpcDecode. На выходе мягкого этапа берутся LLR всего слова
%   (OutputFormat = 'whole', DecisionType = 'soft'), которые затем подаются
%   в LDPC_BITFLIP_DECODE; BF использует только их знак.
%
%   Дополнительные поля p (к LDPC_DEFAULT_PARAMS):
%     p.preprocessAlg  - 'none' | 'bp' | 'layered-bp' | 'norm-min-sum' |
%                        'offset-min-sum'. 'none' — обычный BF без
%                        предобработки (базовый вариант).
%     p.preprocessIter - число итераций мягкого декодера; 0 при 'none'.
%   p.maxIter — лимит итераций BF, как в LDPC_SIMULATE.
%
%   res — структура как у LDPC_SIMULATE, плюс:
%     avgIter      — среднее суммарное число итераций (мягкие + BF);
%     avgIterPre   — среднее число итераций мягкого этапа;
%     avgIterBF    — среднее число итераций BF;
%     berPre, ferPre — ошибки после мягкого этапа, до BF;
%     numPreBitErrors, numPreFrameErrors — их счётчики;
%     numPreSyndromeErrors — кадров с ненулевым синдромом после мягкого этапа;
%     numPreInfoOkSynBad  — из них кадров, где информация уже верна, а синдром нет;
%     numFramesPreOk     — кадров без ошибок в информации после мягкого этапа;
%     numFramesFixedByBF — кадров, где мягкий этап ошибся, а BF исправил.

if ~strcmpi(p.code, '5gnr')
    error('dyse_simulate_bf_pre:code', 'Поддерживается только p.code = ''5gnr''');
end
if ~strcmpi(p.channel, 'awgn')
    error('dyse_simulate_bf_pre:channel', 'Поддерживается только p.channel = ''awgn''');
end
if ~isfield(p, 'preprocessAlg'),  p.preprocessAlg = 'none'; end
if ~isfield(p, 'preprocessIter'), p.preprocessIter = 0; end
p.preprocessAlg = lower(p.preprocessAlg);
if ~any(strcmp(p.preprocessAlg, {'none', 'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum'}))
    error('dyse_simulate_bf_pre:alg', 'Неизвестный мягкий декодер: %s', p.preprocessAlg);
end
doPre = ~strcmpi(p.preprocessAlg, 'none') && p.preprocessIter > 0;

bps = bits_per_symbol(p.modulation);
if p.bgn == 1
    K = 22 * p.Zc; N = 66 * p.Zc;
else
    K = 10 * p.Zc; N = 50 * p.Zc;
end
E = round(K / p.rate);
if E > N
    error('dyse_simulate_bf_pre:rate', 'R = %g ниже материнской скорости кода', p.rate);
end
if mod(E, bps) ~= 0
    error('dyse_simulate_bf_pre:frameLen', 'E = %d не кратно log2(M) = %d', E, bps);
end
if mod(E, p.Zc) ~= 0
    error('dyse_simulate_bf_pre:E', 'E = %d должно быть кратно Zc = %d для явной H', E, p.Zc);
end
R = K / E;
H = nr_ldpc_pcm(p.bgn, p.Zc, E);

if doPre
    cfg = ldpcDecoderConfig(H, p.preprocessAlg);
    preArgs = {'OutputFormat', 'whole', 'DecisionType', 'soft', 'Termination', 'max'};
    switch p.preprocessAlg
        case 'norm-min-sum'
            preArgs = [preArgs, {'MinSumScalingFactor', p.scaling}];
        case 'offset-min-sum'
            preArgs = [preArgs, {'MinSumOffset', p.offset}];
    end
    preName = sprintf('%s, %d итер.', p.preprocessAlg, p.preprocessIter);
else
    preName = 'без предобработки';
end
codeName = sprintf('5G NR BG%d Zc=%d K=%d R=%g', p.bgn, p.Zc, K, R);

nPts = numel(p.EbN0dB);
res.EbN0dB            = p.EbN0dB(:).';
res.ber               = nan(1, nPts);
res.fer               = nan(1, nPts);
res.avgIter           = nan(1, nPts);
res.avgIterPre        = nan(1, nPts);
res.avgIterBF         = nan(1, nPts);
res.berPre            = nan(1, nPts);
res.ferPre            = nan(1, nPts);
res.numFrames         = zeros(1, nPts);
res.numBitErrors      = zeros(1, nPts);
res.numFrameErrors    = zeros(1, nPts);
res.numPreBitErrors   = zeros(1, nPts);
res.numPreFrameErrors = zeros(1, nPts);
res.numPreSyndromeErrors = zeros(1, nPts);
res.numPreInfoOkSynBad   = zeros(1, nPts);
res.numFramesPreOk    = zeros(1, nPts);
res.numFramesFixedByBF = zeros(1, nPts);
res.timePerFrame      = nan(1, nPts);
res.codeName          = codeName;
res.preprocessName    = preName;
res.K = K; res.N = N; res.E = E;
res.params = p;

if p.verbose
    fprintf('%s | %s | BF, %s, maxIter=%d\n', codeName, upper(p.modulation), ...
        preName, p.maxIter);
end

for iPt = 1:nPts
    ebno = p.EbN0dB(iPt);
    snrdB = convertSNR(ebno, "ebno", "snr", BitsPerSymbol=bps, CodingRate=R);
    stream = RandStream('mt19937ar', Seed=p.seed + iPt);

    nFrames = 0; nBitErr = 0; nFrameErr = 0;
    sumIterPre = 0; sumIterBF = 0;
    nPreBitErr = 0; nPreFrameErr = 0; nPreOk = 0; nFixed = 0;
    nPreSynBad = 0; nPreInfoOkSynBad = 0;
    tDec = 0;
    while nFrameErr < p.minFrameErrors && nFrames < p.maxFrames
        B = min(p.batchSize, p.maxFrames - nFrames);

        u = randi(stream, [0 1], K, B, 'int8');
        cw = nrLDPCEncode(u, p.bgn);
        c = cw(1:E, :);
        x = modulate(c(:), p.modulation);
        [y, nVar] = awgn(x, snrdB, 0, stream);
        llr = reshape(demodulate(y, p.modulation, nVar), E, B);
        llrFull = [zeros(2 * p.Zc, B); llr];

        t0 = tic;
        if doPre
            [llrPre, itPre] = ldpcDecode(llrFull, cfg, p.preprocessIter, preArgs{:});
        else
            llrPre = llrFull;
            itPre = zeros(1, B);
        end
        hardPre = llrPre < 0;
        ePre = sum(hardPre(1:K, :) ~= u, 1);
        preSynBad = any(mod(H * double(hardPre), 2), 1);

        [uHat, itBF] = ldpc_bitflip_decode(llrPre, H, K, p.maxIter, ...
            p.bfRule, p.bfOnFail, p.bfThreshold, p.bfMaxFlips);
        tDec = tDec + toc(t0);

        errs = sum(uHat ~= u, 1);
        nBitErr   = nBitErr + sum(errs);
        nFrameErr = nFrameErr + nnz(errs);
        nPreBitErr   = nPreBitErr + sum(ePre);
        nPreFrameErr = nPreFrameErr + nnz(ePre);
        nPreOk       = nPreOk + nnz(ePre == 0);
        nPreSynBad   = nPreSynBad + nnz(preSynBad);
        nPreInfoOkSynBad = nPreInfoOkSynBad + nnz(ePre == 0 & preSynBad);
        nFixed       = nFixed + nnz(ePre > 0 & errs == 0);
        sumIterPre   = sumIterPre + sum(itPre);
        sumIterBF    = sumIterBF + sum(itBF);
        nFrames      = nFrames + B;
    end

    res.numFrames(iPt)         = nFrames;
    res.numBitErrors(iPt)      = nBitErr;
    res.numFrameErrors(iPt)    = nFrameErr;
    res.ber(iPt)               = nBitErr / (nFrames * K);
    res.fer(iPt)               = nFrameErr / nFrames;
    res.numPreBitErrors(iPt)   = nPreBitErr;
    res.numPreFrameErrors(iPt) = nPreFrameErr;
    res.numPreSyndromeErrors(iPt) = nPreSynBad;
    res.numPreInfoOkSynBad(iPt)   = nPreInfoOkSynBad;
    res.berPre(iPt)            = nPreBitErr / (nFrames * K);
    res.ferPre(iPt)            = nPreFrameErr / nFrames;
    res.numFramesPreOk(iPt)    = nPreOk;
    res.numFramesFixedByBF(iPt) = nFixed;
    res.avgIterPre(iPt)        = sumIterPre / nFrames;
    res.avgIterBF(iPt)         = sumIterBF / nFrames;
    res.avgIter(iPt)           = res.avgIterPre(iPt) + res.avgIterBF(iPt);
    res.timePerFrame(iPt)      = tDec / nFrames;

    if p.verbose
        fprintf(['  Eb/N0 = %5.2f dB: FER = %.3e (до BF %.3e), BER = %.3e, ' ...
            'iter = %.2f (pre %.2f + BF %.2f), fixed %d/%d, frames = %d\n'], ...
            ebno, res.fer(iPt), res.ferPre(iPt), res.ber(iPt), res.avgIter(iPt), ...
            res.avgIterPre(iPt), res.avgIterBF(iPt), nFixed, nPreFrameErr, nFrames);
    end
    if p.stopAtZeroFer && nFrameErr == 0
        break
    end
end
end

% -------------------------------------------------------------------------
function k = bits_per_symbol(mod)
switch lower(mod)
    case 'bpsk',  k = 1;
    case 'qpsk',  k = 2;
    case '16qam', k = 4;
    case '64qam', k = 6;
    otherwise, error('dyse_simulate_bf_pre:mod', 'Неизвестная модуляция: %s', mod);
end
end

function x = modulate(bits, mod)
% Gray, единичная средняя мощность (как в LDPC_SIMULATE)
k = bits_per_symbol(mod);
if k == 1
    x = pskmod(bits, 2, 0, 'gray', InputType='bit');
else
    x = qammod(bits, 2^k, 'gray', InputType='bit', UnitAveragePower=true);
end
end

function llr = demodulate(y, mod, nVar)
% Точные LLR: положительное -> бит 0 (как в LDPC_SIMULATE)
k = bits_per_symbol(mod);
if k == 1
    llr = pskdemod(y, 2, 0, 'gray', OutputType='llr', NoiseVariance=nVar);
else
    llr = qamdemod(y, 2^k, 'gray', OutputType='llr', ...
        UnitAveragePower=true, NoiseVariance=nVar);
end
llr = llr(:);
end
