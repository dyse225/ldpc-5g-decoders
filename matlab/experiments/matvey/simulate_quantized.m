function res = simulate_quantized(p)
%SIMULATE_QUANTIZED Монте-Карло как в LDPC_SIMULATE, но с квантованием LLR.
%   res = SIMULATE_QUANTIZED(p)
%   Цепочка та же (nrLDPCEncode -> rv = 0 -> модулятор -> awgn -> LLR,
%   те же RandStream на точку), перед декодером LLR канала квантуются
%   QUANTIZE_LLR. Сообщения внутри декодера не квантуются, как и
%   в task3 исходного репозитория.
%   Дополнительные поля p (к LDPC_DEFAULT_PARAMS):
%     p.llrBits   - число бит на LLR, Inf - без квантования
%     p.llrSat    - уровень насыщения L_sat
%     p.algorithm - помимо алгоритмов LDPC_SIMULATE ещё 'min-sum'
%                   (обычный min-sum = NMS с alpha = 1, как в task3)
%     p.minBer    - не моделировать дальше, если BER на точке < minBer
%                   (по умолчанию 0 - как в LDPC_SIMULATE)
%   При llrBits = Inf результаты покадрово совпадают с LDPC_SIMULATE.
%   res - структура как у LDPC_SIMULATE; ещё поле llrMax - средний
%         максимум |LLR| в кадре до квантования (для выбора L_sat)

if ~isfield(p, 'llrBits'), p.llrBits = Inf; end
if ~isfield(p, 'llrSat'),  p.llrSat = 8; end
if ~isfield(p, 'minBer'),  p.minBer = 0; end

bps = bits_per_symbol(p.modulation);
if p.bgn == 1
    K = 22 * p.Zc; N = 66 * p.Zc;
else
    K = 10 * p.Zc; N = 50 * p.Zc;
end
E = round(K / p.rate);
if E > N
    error('simulate_quantized:rate', 'R = %g ниже материнской скорости кода', p.rate);
end
if mod(E, bps) ~= 0
    error('simulate_quantized:frameLen', 'E = %d не кратно log2(M) = %d', E, bps);
end
R = K / E;
codeName = sprintf('5G NR BG%d Zc=%d K=%d R=%g', p.bgn, p.Zc, K, R);
decArgs = decoder_args(p);

nPts = numel(p.EbN0dB);
res.EbN0dB         = p.EbN0dB(:).';
res.ber            = nan(1, nPts);
res.fer            = nan(1, nPts);
res.avgIter        = nan(1, nPts);
res.numFrames      = zeros(1, nPts);
res.numBitErrors   = zeros(1, nPts);
res.numFrameErrors = zeros(1, nPts);
res.timePerFrame   = nan(1, nPts);
res.llrMax         = nan(1, nPts);
res.codeName = codeName;
res.K = K; res.N = N; res.E = E;
res.params = p;

if p.verbose
    if isinf(p.llrBits), qStr = 'без квантования';
    else, qStr = sprintf('%d бит, L_sat = %g', p.llrBits, p.llrSat); end
    fprintf('%s | %s | %s, maxIter=%d | LLR: %s\n', codeName, ...
        upper(p.modulation), p.algorithm, p.maxIter, qStr);
end

for iPt = 1:nPts
    ebno = p.EbN0dB(iPt);
    snrdB = convertSNR(ebno, "ebno", "snr", BitsPerSymbol=bps, CodingRate=R);
    stream = RandStream('mt19937ar', Seed=p.seed + iPt);

    nFrames = 0; nBitErr = 0; nFrameErr = 0; sumIter = 0; tDec = 0; sumMax = 0;
    while nFrameErr < p.minFrameErrors && nFrames < p.maxFrames
        B = min(p.batchSize, p.maxFrames - nFrames);

        u = randi(stream, [0 1], K, B, 'int8');
        cw = nrLDPCEncode(u, p.bgn);
        c = cw(1:E, :);
        x = modulate(c(:), p.modulation);
        [y, nVar] = awgn(x, snrdB, 0, stream);
        llr = reshape(demodulate(y, p.modulation, nVar), E, B);
        sumMax = sumMax + sum(max(abs(llr), [], 1));
        llr = quantize_llr(llr, p.llrBits, p.llrSat);

        t0 = tic;
        llrFull = [llr; zeros(N - E, B)];      % непереданные биты: LLR = 0
        [uHat, iters] = nrLDPCDecode(llrFull, p.bgn, p.maxIter, decArgs{:});
        tDec = tDec + toc(t0);

        errs = sum(uHat ~= u, 1);
        nBitErr   = nBitErr + sum(errs);
        nFrameErr = nFrameErr + nnz(errs);
        sumIter   = sumIter + sum(iters);
        nFrames   = nFrames + B;
    end

    res.numFrames(iPt)      = nFrames;
    res.numBitErrors(iPt)   = nBitErr;
    res.numFrameErrors(iPt) = nFrameErr;
    res.ber(iPt)            = nBitErr / (nFrames * K);
    res.fer(iPt)            = nFrameErr / nFrames;
    res.avgIter(iPt)        = sumIter / nFrames;
    res.timePerFrame(iPt)   = tDec / nFrames;
    res.llrMax(iPt)         = sumMax / nFrames;

    if p.verbose
        fprintf('  Eb/N0 = %5.2f dB: BER = %.3e, FER = %.3e, iter = %5.2f, frames = %d\n', ...
            ebno, res.ber(iPt), res.fer(iPt), res.avgIter(iPt), nFrames);
    end
    if (p.stopAtZeroFer && nFrameErr == 0) || res.ber(iPt) < p.minBer
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
    case 'min-sum',        alg = 'Normalized min-sum';
    otherwise, error('simulate_quantized:alg', 'Неизвестный алгоритм: %s', p.algorithm);
end
args = {'Algorithm', alg, 'Termination', p.termination};
switch p.algorithm
    case 'norm-min-sum'
        args = [args, {'ScalingFactor', p.scaling}];
    case 'offset-min-sum'
        args = [args, {'Offset', p.offset}];
    case 'min-sum'
        args = [args, {'ScalingFactor', 1}];
end
end

function k = bits_per_symbol(mod)
switch lower(mod)
    case 'bpsk',  k = 1;
    case 'qpsk',  k = 2;
    case '16qam', k = 4;
    case '64qam', k = 6;
    otherwise, error('simulate_quantized:mod', 'Неизвестная модуляция: %s', mod);
end
end

function x = modulate(bits, mod)
% как в LDPC_SIMULATE: Gray, единичная средняя мощность
k = bits_per_symbol(mod);
if k == 1
    x = pskmod(bits, 2, 0, 'gray', InputType='bit');
else
    x = qammod(bits, 2^k, 'gray', InputType='bit', UnitAveragePower=true);
end
end

function llr = demodulate(y, mod, nVar)
% как в LDPC_SIMULATE: точные LLR, положительное -> бит 0
k = bits_per_symbol(mod);
if k == 1
    llr = pskdemod(y, 2, 0, 'gray', OutputType='llr', NoiseVariance=nVar);
else
    llr = qamdemod(y, 2^k, 'gray', OutputType='llr', ...
        UnitAveragePower=true, NoiseVariance=nVar);
end
end
