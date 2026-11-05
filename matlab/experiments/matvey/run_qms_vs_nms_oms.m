%% Quantized Min-Sum (5-битные LLR) против NMS, OMS и некодированной 16-QAM
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

quick = false;                     % true - пробный прогон (до 2e4 кадров на точку)

llrBits = 5;
llrSat  = 6;                       % по run_lsat_sweep
EbN0coded   = 3:0.25:6.5;
EbN0uncoded = 0:1:14;
targetBer   = [1e-4 1e-5];

base = ldpc_default_params();      % BG2, Zc = 52, K = 520, R = 1/2, AWGN
base.modulation = '16qam';
base.EbN0dB = EbN0coded;
base.maxFrames = 3e5;              % для BER до ~1e-6 этого хватает
base.minBer = 1e-6;                % ниже целевых BER не считаем
if quick, base.maxFrames = 2e4; end

% некодированная 16-QAM
p = base;
p.code = 'none';
p.EbN0dB = EbN0uncoded;
p.minFrameErrors = 200;
p.maxFrames = 2e5;                 % 2e8 бит: BER до ~1e-6
rUnc = ldpc_simulate(p);
berTh = berawgn(EbN0uncoded, 'qam', 16);

% кодированные кривые (считаются параллельно, по кривой на воркер):
%   NMS и OMS без квантования (alpha, beta - значения проекта по умолчанию);
%   min-sum без квантования - чтобы отделить потери квантования от потерь min-sum;
%   Quantized Min-Sum: LLR канала llrBits бит, насыщение L_sat.
% При llrBits = Inf simulate_quantized совпадает с ldpc_simulate покадрово.
cfg = repmat(base, 1, 4);
cfg(1).algorithm = 'norm-min-sum';   cfg(1).llrBits = Inf; cfg(1).llrSat = llrSat;
cfg(2).algorithm = 'offset-min-sum'; cfg(2).llrBits = Inf; cfg(2).llrSat = llrSat;
cfg(3).algorithm = 'min-sum';        cfg(3).llrBits = Inf; cfg(3).llrSat = llrSat;
cfg(4).algorithm = 'min-sum';        cfg(4).llrBits = llrBits; cfg(4).llrSat = llrSat;
coded = cell(1, 4);
parfor i = 1:4
    coded{i} = simulate_quantized(cfg(i));
end
[rNms, rOms, rMs, rQms] = coded{:};
labels = {sprintf('NMS, \\alpha = %.2f', base.scaling), ...
          sprintf('OMS, \\beta = %.2f', base.offset), ...
          'MS без квантования', ...
          sprintf('QMS, %d бит, L_{sat} = %g', llrBits, llrSat)};

% Eb/N0 при заданном BER и выигрыш кодирования
nT = numel(targetBer);
ebAt = nan(numel(coded) + 1, nT);
for t = 1:nT
    ebAt(1, t) = ebno_at_ber(EbN0uncoded, berTh, targetBer(t));   % теория
    for i = 1:numel(coded)
        ebAt(i + 1, t) = ebno_at_ber(coded{i}.EbN0dB, coded{i}.ber, targetBer(t));
    end
end
names = [{'Uncoded 16-QAM'}, {'NMS', 'OMS', 'MS', 'QMS'}];
gain = ebAt(1, :) - ebAt;
tbl = table(names(:), ebAt(:, 1), gain(:, 1), ebAt(:, 2), gain(:, 2), ...
    'VariableNames', {'curve', 'EbN0_at_1e4', 'gain_1e4', 'EbN0_at_1e5', 'gain_1e5'});
disp(tbl);
fprintf('Потери квантования (QMS - MS):  %s дБ\n', mat2str(ebAt(5, :) - ebAt(4, :), 3));
fprintf('QMS - NMS:                      %s дБ\n', mat2str(ebAt(5, :) - ebAt(2, :), 3));
fprintf('QMS - OMS:                      %s дБ\n', mat2str(ebAt(5, :) - ebAt(3, :), 3));
writetable(tbl, fullfile(resDir, 'qms_vs_nms_oms_gain.csv'));

save(fullfile(resDir, 'qms_vs_nms_oms.mat'), 'rUnc', 'berTh', 'rNms', 'rOms', 'rMs', 'rQms', ...
    'labels', 'targetBer', 'ebAt', 'gain', 'llrBits', 'llrSat');

fig = plot_qms_vs_nms_oms(fullfile(resDir, 'qms_vs_nms_oms.mat'));
exportgraphics(fig, fullfile(resDir, 'qms_vs_nms_oms.png'), 'Resolution', 150);
