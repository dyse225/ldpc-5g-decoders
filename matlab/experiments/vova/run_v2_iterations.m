%% V2: влияние числа итераций на layered BP (прямой порядок слоёв)
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

quick = false;                     % true - пробный прогон (до 2e4 кадров на точку)

iterList = [5 10 20 50];

results = []; labels = {};
for it = iterList
    p = ldpc_default_params();     % BG2, Zc = 52, K = 520, R = 1/2, BPSK, AWGN
    p.maxIter = it;
    p.layerOrder = 'forward';
    if quick, p.maxFrames = 2e4; end
    results = [results, vova_simulate(p)]; %#ok<AGROW>
    labels{end+1} = sprintf('I_{max} = %d', it); %#ok<SAGROW>
end

save(fullfile(resDir, 'v2_iterations.mat'), 'results', 'labels');
fig = ldpc_plot_results(results, labels, ...
    sprintf('V2: число итераций, layered BP (прямой порядок), %s, BPSK', ...
    results(1).codeName));
if isprop(fig, 'Theme'), fig.Theme = 'light'; end   % светлая тема независимо от настроек MATLAB
drawnow; for a = findall(fig, 'Type', 'axes').', a.Toolbar.Visible = 'off'; end   % без панели инструментов на картинке
exportgraphics(fig, fullfile(resDir, 'v2_iterations.png'), 'Resolution', 150);
