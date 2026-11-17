function results = run_bf_preprocess(quick)
%RUN_BF_PREPROCESS BF после нескольких итераций мягкого декодера.
%   results = RUN_BF_PREPROCESS(quick)
%
%   Сравнивает обычный bit-flip с вариантом, где вход BF — это мягкий
%   выход BP / layered BP / NMS после 1, 2, 5 или 10 итераций. Код и канал
%   те же, что в остальных экспериментах dyse: 5G NR BG2, Zc = 52, K = 520,
%   R = 1/2, BPSK, AWGN.
%
%   quick = true — пробный прогон (1e4 кадров на точку, 50 ошибочных
%   кадров), по умолчанию полный (3e4 кадров, 100 ошибочных кадров).
%
%   Результаты: results/bf_preprocess.mat, bf_preprocess.png,
%   bf_preprocess.csv.

if nargin < 1, quick = false; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

cfg = {
    'none',         0, 'BF'
    'bp',           1, 'BP 1 + BF'
    'bp',           2, 'BP 2 + BF'
    'bp',           5, 'BP 5 + BF'
    'bp',          10, 'BP 10 + BF'
    'layered-bp',   5, 'Layered BP 5 + BF'
    'norm-min-sum', 5, 'NMS 5 + BF'
    };

p0 = ldpc_default_params();        % BG2, Zc = 52, K = 520, R = 1/2, BPSK, AWGN
p0.algorithm = 'bit-flip';
p0.EbN0dB = 3:0.5:6;               % область, где BF ещё не бесполезен
p0.minFrameErrors = 100;
p0.maxFrames = 3e4;                % хватает и для нулевого FER (верхняя оценка ~1e-4)
p0.stopAtZeroFer = false;          % нужны все точки для сравнения
if quick
    p0.minFrameErrors = 50;
    p0.maxFrames = 1e4;
end

results = [];
labels = {};
partialFile = fullfile(resDir, 'bf_preprocess_partial.mat');
for i = 1:size(cfg, 1)
    p = p0;
    p.preprocessAlg = cfg{i, 1};
    p.preprocessIter = cfg{i, 2};
    results = [results, dyse_simulate_bf_pre(p)]; %#ok<AGROW>
    labels{end+1} = cfg{i, 3}; %#ok<SAGROW>
    save(partialFile, 'results', 'labels', 'quick');   % чекпойнт после каждой кривой
end

save(fullfile(resDir, 'bf_preprocess.mat'), 'results', 'labels', 'quick');

% --- график: FER, BER и суммарное число итераций
fig = plot_bf_preprocess(fullfile(resDir, 'bf_preprocess.mat'));
exportgraphics(fig, fullfile(resDir, 'bf_preprocess.png'), 'Resolution', 150);

% --- таблица и CSV
rows = [];
for i = 1:numel(results)
    r = results(i);
    rows = [rows; [repmat(i, numel(r.EbN0dB), 1), r.EbN0dB(:), ...
        r.ber(:), r.fer(:), r.berPre(:), r.ferPre(:), ...
        r.avgIter(:), r.avgIterPre(:), r.avgIterBF(:), ...
        r.numFrames(:), r.numFrameErrors(:), r.numPreFrameErrors(:), ...
        r.numPreSyndromeErrors(:), r.numPreInfoOkSynBad(:), ...
        r.numFramesPreOk(:), r.numFramesFixedByBF(:), r.timePerFrame(:)]]; %#ok<AGROW>
end
T = array2table(rows, 'VariableNames', {'config', 'EbN0dB', 'BER', 'FER', ...
    'BER_pre', 'FER_pre', 'iterTotal', 'iterPre', 'iterBF', 'frames', ...
    'frameErrors', 'preFrameErrors', 'preSyndromeErrors', 'preInfoOkSynBad', ...
    'preOk', 'fixedByBF', 'timePerFrame'});
writetable(T, fullfile(resDir, 'bf_preprocess.csv'));

% --- краткий вывод по точкам 3..6 дБ
probe = 3:6;
fprintf('\n%-22s', 'конфигурация');
for e = probe, fprintf('  FER %.0f дБ', e); end
fprintf('   |  FER до BF (3..6 дБ)\n');
for i = 1:numel(results)
    r = results(i);
    fprintf('%-22s', labels{i});
    for e = probe
        j = find(r.EbN0dB == e, 1);
        if isempty(j), fprintf('  %10s', '—'); else, fprintf('  %10.3e', r.fer(j)); end
    end
    fprintf('   |');
    for e = probe
        j = find(r.EbN0dB == e, 1);
        if isempty(j), fprintf('  %10s', '—'); else, fprintf('  %10.3e', r.ferPre(j)); end
    end
    fprintf('\n');
end
fprintf('\n');
eLast = probe(end);
for i = 2:numel(results)
    r = results(i);
    j = find(r.EbN0dB == eLast, 1);
    if isempty(j), continue; end
    fprintf('%-22s: BF исправил %d из %d кадров, где мягкий этап ошибся (%g дБ)\n', ...
        labels{i}, r.numFramesFixedByBF(j), r.numPreFrameErrors(j), eLast);
end
end
