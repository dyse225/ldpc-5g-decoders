%% Параметры BF: число инвертируемых бит и порог инвертирования
% 5G NR BG2, K = 520, R = 1/2, BPSK, AWGN, правило 'gallager'.
%   Число бит за итерацию kMax: 1 (= bitflipping_serial из bf-tests-awgn-egldpc),
%     5, 20, Inf (= bitflipping: все биты с максимумом). Лимит 300 итераций,
%     чтобы последовательному варианту хватило итераций.
%   Порог относительно максимума d: инвертируются биты, у которых невыполненных
%     проверок >= max - d; d = 0 — обычный BF (только максимум), d = 1, 2, 3.
%     Фиксированный порог (>= T, T = 1..3) и порог по большинству проверок
%     (Sipser–Spielman) проверены отдельно: на 5G NR FER = 1 при любых Eb/N0.
%     Готовой реализации порогового варианта в Lcrypto нет (Sipser–Spielman
%     вызывается в simula.c, но исходник в репозитории отсутствует).
%   Части можно запускать по отдельности: doMaxFlips / doThreshold.
clearvars -except doMaxFlips doThreshold; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
EbN0 = 3:0.5:9;
if ~exist('doMaxFlips', 'var'), doMaxFlips = true; end
if ~exist('doThreshold', 'var'), doThreshold = true; end

% --- число инвертируемых бит
if doMaxFlips
kList = [1 5 20 Inf];
resK = []; labK = {};
for k = kList
    p = ldpc_default_params();
    p.algorithm = 'bit-flip'; p.bfMaxFlips = k; p.maxIter = 300; p.EbN0dB = EbN0;
    resK = [resK, ldpc_simulate(p)]; %#ok<AGROW>
    if isinf(k), labK{end+1} = 'все с максимумом'; else, labK{end+1} = sprintf('k = %d', k); end %#ok<SAGROW>
end
save(fullfile(resDir, 'bf_maxflips.mat'), 'resK', 'labK');
fig = ldpc_plot_results(resK, labK, 'BF: число инвертируемых бит за итерацию (лимит 300 итераций)', ...
    struct('layout', 'column', 'xlim', [3 9], 'berLim', [1e-5 0.3], 'ferLim', [1e-3 1]));
exportgraphics(fig, fullfile(resDir, 'bf_maxflips.png'), 'Resolution', 150);
end

% --- порог инвертирования относительно максимума
if doThreshold
tList = [0 1 2 3];
resT = []; labT = {};
for t = tList
    p = ldpc_default_params();
    p.algorithm = 'bit-flip'; p.bfThreshold = t; p.EbN0dB = EbN0;
    resT = [resT, ldpc_simulate(p)]; %#ok<AGROW>
    if t == 0, labT{end+1} = 'max (обычный BF)'; else, labT{end+1} = sprintf('max - %d', t); end %#ok<SAGROW>
end
save(fullfile(resDir, 'bf_threshold.mat'), 'resT', 'labT');
fig = ldpc_plot_results(resT, labT, 'BF: порог инвертирования относительно максимума (50 итераций)', ...
    struct('layout', 'column', 'xlim', [3 9], 'berLim', [1e-5 0.3], 'ferLim', [1e-3 1]));
exportgraphics(fig, fullfile(resDir, 'bf_threshold.png'), 'Resolution', 150);
end
