function cfg = qli_config(preset)
%QLI_CONFIG Параметры эксперимента: 5G NR LDPC, QPSK, AWGN, BP / Layered BP / NMS / OMS.
%   cfg = QLI_CONFIG('standard')  — основной прогон
%   cfg = QLI_CONFIG('quick')     — быстрая проверка всей цепочки (минуты)
if nargin < 1, preset = 'standard'; end

% ---- модуляция, канал, код --------------------------------------------------
cfg.modulation = 'qpsk';              % 2 бита на символ, Gray, единичная мощность
cfg.channel    = 'awgn';
cfg.rate       = 1/2;                 % скорость кода (E = K/R)
cfg.bgList     = [2 1];               % порядок: сначала BG2; BG1 уберите, если не хватает времени
cfg.Zc         = struct('bg2', 52, 'bg1', 176);   % K = 10*52 = 520; K = 22*176 = 3872

% ---- декодеры и итерации ----------------------------------------------------
cfg.decoders     = {'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum'};
cfg.decoderNames = {'BP (flooding)', 'Layered BP', 'NMS', 'OMS'};
cfg.iterList     = [5 10 20 50];      % I_max
cfg.baseIter     = 50;                % I_max для базового сравнения декодеров
cfg.iterDecoders = {'norm-min-sum', 'layered-bp'};  % для них считаются все I_max;
                                      % остальные — только при baseIter (экономия времени)

% ---- Eb/N0 и Монте-Карло (по BG) -------------------------------------------
cfg.EbN0dB   = struct('bg2', 0:0.25:3,  'bg1', 0:0.25:2.5);
cfg.maxFrames = struct('bg2', 3e5,       'bg1', 5e4);     % BG1 (K = 3872) считается дольше
cfg.minBer = 1e-6;
cfg.minFrameErrors = 100;
cfg.batchSize      = 64;
cfg.seed           = 1;               % один seed: кривые считаются на одних и тех же кадрах
cfg.termination    = 'early';
cfg.verbose        = true;
cfg.parallel       = false;           % true: parfor по конфигурациям
cfg.resume         = true;            % кэш готовых прогонов в results/cache

% ---- оптимальные alpha (NMS) и beta (OMS) ----------------------------------
cfg.param.iter     = 50;                           % I_max, при котором ищется оптимум
cfg.param.alphaList = 0.5:0.05:1.0;
cfg.param.betaList  = 0:0.05:1.0;
cfg.param.EbN0dB   = struct('bg2', [1.5 1.75], 'bg1', [1.0 1.25]);  % рабочие точки: FER ~ 1e-1...1e-2
cfg.param.maxFrames = 5e4;
cfg.param.minFrameErrors = 100;
cfg.alphaOpt = [];                    % задать вручную struct('bg2',0.8,'bg1',0.8) — тогда поиск не нужен
cfg.betaOpt  = [];                    % то же для beta
cfg.defaultAlpha = 0.75;              % если оптимум не найден и не задан
cfg.defaultBeta  = 0.5;

% ---- uncoded ---------------------------------------------------------------
cfg.EbN0uncoded = 0:1:11;       % отдельная сетка для некодированной передачи
cfg.uncodedMaxFrames = 2e5;     % максимальное количество кадров

here = fileparts(mfilename('fullpath'));
cfg.resDir = fullfile(here, 'results');

switch lower(preset)
    case 'quick'
        cfg.iterList = [5 50];
        cfg.EbN0dB = struct('bg2', 0.5:0.5:3, 'bg1', 0.5:0.5:2.5);
        cfg.maxFrames = struct('bg2', 2000, 'bg1', 500);
        cfg.minFrameErrors = 20;
        cfg.param.alphaList = 0.6:0.1:1.0;
        cfg.param.betaList = 0:0.25:1.0;
        cfg.param.maxFrames = 20000; cfg.param.minFrameErrors = 20;
        cfg.EbN0uncoded = 0:2:10;       % сокращённая сетка для быстрой проверки
        cfg.uncodedMaxFrames = 2000;
        cfg.resDir = fullfile(here, 'results');
    case 'standard'
    otherwise
        error('qli_config:preset', 'Неизвестный пресет: %s', preset);
end
end
