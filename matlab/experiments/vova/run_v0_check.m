%% V0: сверка своего layered BP (прямой порядок слоёв) с nrLDPCDecode
%  'Layered belief propagation'. Декодеры работают на одних и тех же кадрах
%  и шуме, поэтому кроме кривых считается число кадров с разными решениями.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

quick = false;                     % true - пробный прогон (до 2e4 кадров на точку)

p = ldpc_default_params();         % BG2, Zc = 52, K = 520, R = 1/2, BPSK, AWGN
p.maxIter = 50;
p.layerOrder = 'forward';
p.compareRef = true;
if quick, p.maxFrames = 2e4; end

[rOwn, rRef] = vova_simulate(p);

% у rOwn есть лишние поля (порядок слоёв, счётчики расхождений) - для графика убираем
extra = setdiff(fieldnames(rOwn), fieldnames(rRef));
results = [rRef, rmfield(rOwn, extra)];
labels = {'nrLDPCDecode, layered BP', 'свой layered BP, прямой порядок'};
save(fullfile(resDir, 'v0_check.mat'), 'results', 'labels', 'rOwn', 'rRef');
fig = ldpc_plot_results(results, labels, ...
    sprintf('V0: свой layered BP против nrLDPCDecode, %s, BPSK, I_{max} = %d', ...
    rOwn.codeName, p.maxIter));
if isprop(fig, 'Theme'), fig.Theme = 'light'; end   % светлая тема независимо от настроек MATLAB
drawnow; for a = findall(fig, 'Type', 'axes').', a.Toolbar.Visible = 'off'; end   % без панели инструментов на картинке
exportgraphics(fig, fullfile(resDir, 'v0_check.png'), 'Resolution', 150);

fprintf('\nКадров с разными решениями: %d из %d\n', ...
    sum(rOwn.numMismatchFrames), sum(rOwn.numFrames));
fprintf('Кадров с разным числом итераций: %d из %d\n', ...
    sum(rOwn.numIterMismatch), sum(rOwn.numFrames));
