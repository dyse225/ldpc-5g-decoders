function out = run_bf_wide(quick)
%RUN_BF_WIDE BF до 12 дБ против мягких декодеров до 5 дБ.
%   out = RUN_BF_WIDE(quick)
%
%   Строит общий график: обычный BF считается до 12 дБ, мягкие декодеры
%   (BP 5/10, layered BP 5, NMS 5) и варианты BP + BF — только до 5 дБ,
%   потому что выше 5 дБ они уже дают FER = 0.
%
%   quick = true — меньшие лимиты кадров для быстрой проверки.

if nargin < 1, quick = false; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

% --- обычный BF: 0..12 дБ
p = ldpc_default_params();
p.algorithm = 'bit-flip';
p.EbN0dB = 0:0.5:12;
p.minFrameErrors = 50;
p.maxFrames = 1e5;
p.stopAtZeroFer = true;
p.verbose = true;
if quick
    p.minFrameErrors = 20;
    p.maxFrames = 1e4;
end
fprintf('\n=== BF, 0..12 дБ ===\n');
bf = ldpc_simulate(p);

% --- мягкие декодеры и BP + BF: 0..5 дБ
p0 = ldpc_default_params();
p0.algorithm = 'bit-flip';
p0.EbN0dB = 0:0.5:5;
p0.minFrameErrors = 50;
p0.maxFrames = 1e4;
p0.stopAtZeroFer = false;           % нужны все точки до 5 дБ, даже нулевые
if quick
    p0.minFrameErrors = 20;
    p0.maxFrames = 3e3;
end

cfg = {
    'bp',           2, 'BP 2 + BF'
    'bp',           5, 'BP 5 + BF'
    'bp',          10, 'BP 10 + BF'
    'layered-bp',   5, 'Layered BP 5 + BF'
    'norm-min-sum', 5, 'NMS 5 + BF'
    };
results = [];
labels = {};
fprintf('\n=== мягкие декодеры и BP + BF, 0..5 дБ ===\n');
for i = 1:size(cfg, 1)
    q = p0;
    q.preprocessAlg = cfg{i, 1};
    q.preprocessIter = cfg{i, 2};
    results = [results, dyse_simulate_bf_pre(q)]; %#ok<AGROW>
    labels{end+1} = cfg{i, 3}; %#ok<SAGROW>
end

save(fullfile(resDir, 'bf_wide.mat'), 'bf', 'results', 'labels', 'quick');
fig = plot_bf_wide(fullfile(resDir, 'bf_wide.mat'));
exportgraphics(fig, fullfile(resDir, 'bf_wide.png'), 'Resolution', 150);
out = struct('bf', bf, 'results', results, 'labels', {labels});
end
