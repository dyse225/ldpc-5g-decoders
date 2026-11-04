function [res, resRef] = vova_simulate(p)
%VOVA_SIMULATE Монте-Карло для своего layered BP (VOVA_LAYERED_BP_DECODE).
%   [res, resRef] = VOVA_SIMULATE(p)
%   Цепочка та же, что в LDPC_SIMULATE (nrLDPCEncode -> rv = 0 -> модулятор
%   -> awgn -> LLR), те же RandStream на точку, но декодер свой.
%   Дополнительные поля p (к LDPC_DEFAULT_PARAMS):
%     p.layerOrder   - 'forward' | 'reverse' | 'degree' | 'degree-asc'
%                      или вектор номеров слоёв
%     p.compareRef   - true: на тех же кадрах дополнительно декодировать
%                      nrLDPCDecode 'Layered belief propagation' (эталон)
%   res    - структура как у LDPC_SIMULATE (свой декодер);
%            при compareRef ещё поля numMismatchFrames (кадры, где решения
%            своего и эталонного декодера расходятся) и numIterMismatch
%   resRef - то же для nrLDPCDecode (только при compareRef)

if ~isfield(p, 'layerOrder'), p.layerOrder = 'forward'; end
if ~isfield(p, 'compareRef'), p.compareRef = false; end

bps = bits_per_symbol(p.modulation);
if p.bgn == 1
    K = 22 * p.Zc; N = 66 * p.Zc;
else
    K = 10 * p.Zc; N = 50 * p.Zc;
end
E = round(K / p.rate);
if E > N
    error('vova_simulate:rate', 'R = %g ниже материнской скорости кода', p.rate);
end
if mod(E, bps) ~= 0
    error('vova_simulate:frameLen', 'E = %d не кратно log2(M) = %d', E, bps);
end
R = K / E;
H = nr_ldpc_pcm(p.bgn, p.Zc, E);
if ischar(p.layerOrder) || isstring(p.layerOrder)
    orderName = char(p.layerOrder);
    order = vova_layer_order(H, p.Zc, orderName);
else
    order = p.layerOrder(:).';
    orderName = mat2str(order);
end
codeName = sprintf('5G NR BG%d Zc=%d K=%d R=%g', p.bgn, p.Zc, K, R);

nPts = numel(p.EbN0dB);
res = init_res(p, nPts, codeName, K, N, E);
res.layerOrder = order;
res.layerOrderName = orderName;
res.numMismatchFrames = zeros(1, nPts);
res.numIterMismatch   = zeros(1, nPts);
resRef = init_res(p, nPts, codeName, K, N, E);

if p.verbose
    fprintf('%s | %s | own layered BP, order = %s %s, maxIter=%d\n', codeName, ...
        upper(p.modulation), orderName, mat2str(order), p.maxIter);
end

for iPt = 1:nPts
    ebno = p.EbN0dB(iPt);
    snrdB = convertSNR(ebno, "ebno", "snr", BitsPerSymbol=bps, CodingRate=R);
    stream = RandStream('mt19937ar', Seed=p.seed + iPt);

    nFrames = 0; nBitErr = 0; nFrameErr = 0; sumIter = 0; tDec = 0;
    rBit = 0; rFrame = 0; rIter = 0; rT = 0; nMis = 0; nIterMis = 0;
    while nFrameErr < p.minFrameErrors && nFrames < p.maxFrames
        B = min(p.batchSize, p.maxFrames - nFrames);

        u = randi(stream, [0 1], K, B, 'int8');
        cw = nrLDPCEncode(u, p.bgn);
        c = cw(1:E, :);
        x = modulate(c(:), p.modulation);
        [y, nVar] = awgn(x, snrdB, 0, stream);
        llr = reshape(demodulate(y, p.modulation, nVar), E, B);

        t0 = tic;
        llrH = [zeros(2 * p.Zc, B); llr];
        [uHat, iters] = vova_layered_bp_decode(llrH, H, p.Zc, K, p.maxIter, ...
            order, p.termination);
        tDec = tDec + toc(t0);

        errs = sum(uHat ~= u, 1);
        nBitErr   = nBitErr + sum(errs);
        nFrameErr = nFrameErr + nnz(errs);
        sumIter   = sumIter + sum(iters);

        if p.compareRef
            t1 = tic;
            llrFull = [llr; zeros(N - E, B)];
            [uRef, itRef] = nrLDPCDecode(llrFull, p.bgn, p.maxIter, ...
                'Algorithm', 'Layered belief propagation', 'Termination', p.termination);
            rT = rT + toc(t1);
            eRef = sum(uRef ~= u, 1);
            rBit = rBit + sum(eRef);
            rFrame = rFrame + nnz(eRef);
            rIter = rIter + sum(itRef);
            nMis = nMis + nnz(any(uRef ~= uHat, 1));
            nIterMis = nIterMis + nnz(itRef(:).' ~= iters);
        end
        nFrames = nFrames + B;
    end

    res = fill_pt(res, iPt, nFrames, nBitErr, nFrameErr, sumIter, tDec, K);
    res.numMismatchFrames(iPt) = nMis;
    res.numIterMismatch(iPt)   = nIterMis;
    if p.compareRef
        resRef = fill_pt(resRef, iPt, nFrames, rBit, rFrame, rIter, rT, K);
    end

    if p.verbose
        fprintf('  Eb/N0 = %5.2f dB: BER = %.3e, FER = %.3e, iter = %5.2f, frames = %d', ...
            ebno, res.ber(iPt), res.fer(iPt), res.avgIter(iPt), nFrames);
        if p.compareRef
            fprintf(' | ref FER = %.3e, iter = %5.2f | разных кадров: %d, разных итераций: %d', ...
                resRef.fer(iPt), resRef.avgIter(iPt), nMis, nIterMis);
        end
        fprintf('\n');
    end
    if p.stopAtZeroFer && nFrameErr == 0
        break
    end
end
end

% -------------------------------------------------------------------------
function res = init_res(p, nPts, codeName, K, N, E)
res.EbN0dB         = p.EbN0dB(:).';
res.ber            = nan(1, nPts);
res.fer            = nan(1, nPts);
res.avgIter        = nan(1, nPts);
res.numFrames      = zeros(1, nPts);
res.numBitErrors   = zeros(1, nPts);
res.numFrameErrors = zeros(1, nPts);
res.timePerFrame   = nan(1, nPts);
res.codeName = codeName;
res.K = K; res.N = N; res.E = E;
res.params = p;
end

function res = fill_pt(res, iPt, nFrames, nBitErr, nFrameErr, sumIter, tDec, K)
res.numFrames(iPt)      = nFrames;
res.numBitErrors(iPt)   = nBitErr;
res.numFrameErrors(iPt) = nFrameErr;
res.ber(iPt)            = nBitErr / (nFrames * K);
res.fer(iPt)            = nFrameErr / nFrames;
res.avgIter(iPt)        = sumIter / nFrames;
res.timePerFrame(iPt)   = tDec / nFrames;
end

function k = bits_per_symbol(mod)
switch lower(mod)
    case 'bpsk',  k = 1;
    case 'qpsk',  k = 2;
    case '16qam', k = 4;
    case '64qam', k = 6;
    otherwise, error('vova_simulate:mod', 'Неизвестная модуляция: %s', mod);
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
