function p = ldpc_default_params()
%LDPC_DEFAULT_PARAMS Параметры моделирования по умолчанию.
%   p = LDPC_DEFAULT_PARAMS() возвращает структуру, которую можно
%   модифицировать и передать в LDPC_SIMULATE.

% Код: 5G NR LDPC (TS 38.212), nrLDPCEncode / nrLDPCDecode
p.code = '5gnr';              % '5gnr' | 'none' — без кодирования
p.bgn  = 2;                   % базовый граф: 1 (K = 22*Zc) или 2 (K = 10*Zc)
p.Zc   = 52;                  % размер подъёма -> K = 520 для BG2
p.rate = 1/2;                 % скорость после согласования (E = K/R бит в канал)

% Модуляция: 'bpsk' | 'qpsk' | '16qam' | '64qam'
p.modulation = 'bpsk';

% Канал
p.EbN0dB = 0:0.25:2;          % Eb/N0 на информационный бит, дБ (FER до ~1e-4)

% Декодер
p.algorithm   = 'bp';         % 'bp' | 'layered-bp' | 'norm-min-sum' | 'offset-min-sum'
                              % | 'bit-flip' (ldpc_bitflip_decode, порт C-кода PHBF)
                              % | 'wbf' (ldpc_wbf_decode, порт WBF/MWBF из AFF3CT)
p.maxIter     = 50;
p.scaling     = 0.75;         % alpha для norm-min-sum
p.offset      = 0.5;          % beta для offset-min-sum
p.termination = 'early';      % 'early' (синдром = 0) | 'max'
p.bfRule      = 'gallager';   % bit-flip: 'gallager' | 'phbf' (как в оригинале)
p.bfOnFail    = 'final';      % bit-flip при неудаче: 'final' (итоговое состояние)
                              % | 'channel' (решение канала, как в оригинале)
p.wbfAlpha    = 0.7;          % WBF: коэффициент MWBF (0 -> обычный WBF)
p.wbfMaxFlips = 300;          % WBF: максимум итераций (1 инверсия за итерацию)

% Монте-Карло
p.minFrameErrors = 100;       % останов на точке по числу ошибочных кадров
p.maxFrames      = 1e6;       % ... или по общему числу кадров (страховка)
p.batchSize      = 64;        % кадров за один вызов кодера/декодера
p.stopAtZeroFer  = true;      % не моделировать дальше, если FER = 0
p.seed           = 1;
p.verbose        = true;
end
