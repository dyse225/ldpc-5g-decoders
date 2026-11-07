%% E2 (BP): влияние maxIter на BP flooding, maxIter = 5, 10, 20, 50, 100
% Меняется только p.maxIter, всё остальное по умолчанию (BG2, K = 520, early).
% Образец: run_e2_iterations.m.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(fileparts(here))), fullfile(fileparts(here), 'моделирование 1'));
maxFrames = 1e6;      % финальный лимит кадров на точку
resDir = fullfile(here, 'results_1e6');
if ~exist(resDir, 'dir'), mkdir(resDir); end
% Повторно используем только завершённый расчёт с заданным лимитом кадров (CSV записывается последним).
if exist(fullfile(resDir, 'e2_bp_iterations.csv'), 'file') && exist(fullfile(resDir, 'e2_bp_iterations.mat'), 'file')
    cached = load(fullfile(resDir, 'e2_bp_iterations.mat'), 'results');
    if numel(cached.results) == 5 && all(arrayfun(@(r) r.params.maxFrames == maxFrames, cached.results))
        fprintf('Complete e2_bp_iterations at maxFrames=%g already saved.\n', maxFrames);
        return
    end
end


iterList = [5 10 20 50 100];
results = []; labels = {};
valuesToRun = iterList;
cachePath = fullfile(resDir, 'e2_bp_iterations.mat');
if exist(cachePath, 'file')
    cached = load(cachePath);
    if numel(cached.results) == numel(iterList) && ...
            all(arrayfun(@(r) r.params.maxFrames == maxFrames, cached.results)) && ...
            isequal(arrayfun(@(r) r.params.maxIter, cached.results), iterList)
        results = cached.results; labels = cached.labels; valuesToRun = [];
    end
end
for it = valuesToRun
    p = ldpc_default_params();
    p.maxFrames = maxFrames;
    p.algorithm = 'bp';
    p.maxIter   = it;
    results = [results, ldpc_simulate(p)]; %#ok<AGROW>
    save(fullfile(resDir, 'e2_bp_iterations_checkpoint.mat'), 'results');
    labels{end+1} = sprintf('I_{max}=%d', it); %#ok<SAGROW>
end
save(fullfile(resDir, 'e2_bp_iterations.mat'), 'results', 'labels', 'iterList');

% три панели: BER, FER, среднее число итераций
fig = ldpc_plot_results(results, labels, ...
    sprintf('E2: BP flooding, число итераций, %s, BPSK', results(1).codeName));
exportgraphics(fig, fullfile(resDir, 'e2_bp_iterations.png'), 'Resolution', 150);

% один график BER (для слайда)
fig2 = figure('Color', 'w', 'Position', [100 100 700 500]);
mk = {'o', 's', 'd', '^', 'v'};
for i = 1:numel(results)
    r = results(i); ok = r.numFrames > 0; b = r.ber; b(b == 0) = NaN;
    semilogy(r.EbN0dB(ok), b(ok), ['-' mk{i}], 'LineWidth', 1.8, 'MarkerSize', 7, ...
        'DisplayName', labels{i}); hold on
end
grid on; xlabel('E_b/N_0, дБ'); ylabel('BER'); legend('Location', 'southwest');
title(sprintf('BP flooding: BER, %s', results(1).codeName));
exportgraphics(fig2, fullfile(resDir, 'e2_bp_iterations_ber.png'), 'Resolution', 150);

% --- сколько итераций достаточно: Eb/N0 при FER = 1e-2 и 1e-3 ---
fprintf('\nEb/N0 (дБ) при заданном FER; потеря относительно Imax = %d\n', iterList(end));
fprintf('%6s | %10s %10s | %10s %10s\n', 'Imax', 'FER=1e-2', 'потеря', 'FER=1e-3', 'потеря');
e2 = arrayfun(@(r) ebno_at_target(r, 'fer', 1e-2), results);
e3 = arrayfun(@(r) ebno_at_target(r, 'fer', 1e-3), results);
for i = 1:numel(iterList)
    fprintf('%6d | %10.2f %+10.2f | %10.2f %+10.2f\n', iterList(i), ...
        e2(i), e2(i) - e2(end), e3(i), e3(i) - e3(end));
end

% Численные результаты для экспорта, включая фактический объём выборки Монте-Карло.
rows = [];
for i=1:numel(results)
    r=results(i);
    rows=[rows; [repmat(i,numel(r.ber),1),r.EbN0dB(:),r.ber(:),r.fer(:), ...
        r.avgIter(:),r.numFrames(:),r.numBitErrors(:),r.numFrameErrors(:)]]; %#ok<AGROW>
end
T=array2table(rows,'VariableNames',{'curve','EbN0dB','BER','FER','avgIter','numFrames','numBitErrors','numFrameErrors'});
writetable(T,fullfile(resDir,'e2_bp_iterations.csv'));
close all;
