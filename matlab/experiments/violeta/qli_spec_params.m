function p = qli_spec_params(s, cfg)
%QLI_SPEC_PARAMS Параметры ldpc_simulate по описанию прогона s.
%   s: bgn, Zc, algorithm ('uncoded' — без кодирования), maxIter, scaling, offset,
%      EbN0dB, maxFrames, minFrameErrors
p = ldpc_default_params();
p.modulation = cfg.modulation;
p.channel = cfg.channel;
p.EbN0dB = s.EbN0dB;
p.maxFrames = s.maxFrames;
p.minFrameErrors = s.minFrameErrors;
p.batchSize = cfg.batchSize;
p.seed = cfg.seed;
p.termination = cfg.termination;
p.verbose = cfg.verbose && ~cfg.parallel;
if strcmp(s.algorithm, 'uncoded')
    p.code = 'none';
    return
end
p.code = '5gnr';
p.bgn = s.bgn; p.Zc = s.Zc; p.rate = cfg.rate;
p.algorithm = s.algorithm; p.maxIter = s.maxIter;
p.scaling = s.scaling; p.offset = s.offset;
if s.bgn == 1, K = 22 * s.Zc; else, K = 10 * s.Zc; end
E = round(K / cfg.rate);
if mod(E, 2) ~= 0
    error('qli:E', 'Для QPSK E = %d должно быть чётным (K = %d, R = %g)', E, K, cfg.rate);
end
end
