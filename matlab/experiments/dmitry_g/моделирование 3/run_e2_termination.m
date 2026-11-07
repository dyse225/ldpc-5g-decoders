%% E2 (termination): критерий остановки BP: 'early' (синдром = 0) против 'max'
% BP flooding, maxIter = 50, остальное по умолчанию (BG2, K = 520).
% Для обоих прогонов seed одинаковый (seed = p.seed + номер точки), т.е. шум
% и информационные биты в кадрах те же: разница — только от критерия остановки.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(fileparts(here))), fullfile(fileparts(here), 'моделирование 1'));
maxFrames = 1e6;      % финальный лимит кадров на точку
resDir = fullfile(here, 'results_1e6');
if ~exist(resDir, 'dir'), mkdir(resDir); end
% Повторно используем только завершённый расчёт с заданным лимитом кадров (CSV записывается последним).
if exist(fullfile(resDir, 'e2_termination.csv'), 'file') && exist(fullfile(resDir, 'e2_termination.mat'), 'file')
    cached = load(fullfile(resDir, 'e2_termination.mat'), 'results');
    if numel(cached.results) == 2 && all(arrayfun(@(r) r.params.maxFrames == maxFrames, cached.results))
        fprintf('Complete e2_termination at maxFrames=%g already saved.\n', maxFrames);
        return
    end
end


terminationList = {'early', 'max'};
labels = {'early (синдром = 0)', 'max (всегда 50 итераций)'};
results = [];
valuesToRun = 1:numel(terminationList);
cachePath = fullfile(resDir, 'e2_termination.mat');
if exist(cachePath, 'file')
    cached = load(cachePath);
    if numel(cached.results) == 2 && ...
            all(arrayfun(@(r) r.params.maxFrames == maxFrames && r.params.maxIter == 50, cached.results)) && ...
            strcmp(cached.results(1).params.termination,'early') && strcmp(cached.results(2).params.termination,'max')
        results = cached.results; valuesToRun = [];
    end
end
for i = valuesToRun
    p = ldpc_default_params();
    p.maxFrames   = maxFrames;
    p.algorithm   = 'bp';
    p.maxIter     = 50;
    p.termination = terminationList{i};
    results = [results, ldpc_simulate(p)]; %#ok<AGROW>
    save(fullfile(resDir, 'e2_termination_checkpoint.mat'), 'results');
end
save(fullfile(resDir, 'e2_termination.mat'), 'results', 'labels');

fig = ldpc_plot_results(results, labels, ...
    sprintf('E2: критерий остановки BP, %s, BPSK, maxIter = 50', results(1).codeName));
exportgraphics(fig, fullfile(resDir, 'e2_termination.png'), 'Resolution', 150);

% --- числовое сравнение ---
n = min(numel(results(1).ber), numel(results(2).ber));
fprintf('\n%8s | %12s %12s | %10s %10s\n', 'Eb/N0', 'BER early', 'BER max', 'iter early', 'iter max');
for k = 1:n
    fprintf('%8.2f | %12.3e %12.3e | %10.2f %10.2f\n', results(1).EbN0dB(k), ...
        results(1).ber(k), results(2).ber(k), results(1).avgIter(k), results(2).avgIter(k));
end
fprintf('\nЧисло ошибочных кадров (early vs max); BER проверяется по numBitErrors:\n');
disp([results(1).numFrameErrors(1:n); results(2).numFrameErrors(1:n)]);

assert(all(results(2).avgIter(results(2).numFrames > 0) == 50));

% Численные результаты для экспорта, включая фактический объём выборки Монте-Карло.
rows = [];
for i=1:numel(results)
    r=results(i);
    rows=[rows; [repmat(i,numel(r.ber),1),r.EbN0dB(:),r.ber(:),r.fer(:), ...
        r.avgIter(:),r.numFrames(:),r.numBitErrors(:),r.numFrameErrors(:)]]; %#ok<AGROW>
end
T=array2table(rows,'VariableNames',{'curve','EbN0dB','BER','FER','avgIter','numFrames','numBitErrors','numFrameErrors'});
writetable(T,fullfile(resDir,'e2_termination.csv'));
close all;
