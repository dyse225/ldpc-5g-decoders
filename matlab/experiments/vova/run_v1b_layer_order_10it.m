%% V1b: порядок обновления слоёв при малом числе итераций (I_max = 10)
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

quick = false;                     % true - пробный прогон (до 2e4 кадров на точку)

orders = {'forward', 'reverse', 'degree'};
names  = {'прямой', 'обратный', 'по степени строк'};

results = []; labels = {};
for k = 1:numel(orders)
    p = ldpc_default_params();     % BG2, Zc = 52, K = 520, R = 1/2, BPSK, AWGN
    p.maxIter = 10;
    p.layerOrder = orders{k};
    if quick, p.maxFrames = 2e4; end
    r = vova_simulate(p);
    results = [results, r]; %#ok<AGROW>
    labels{end+1} = sprintf('%s: %s', names{k}, mat2str(r.layerOrder)); %#ok<SAGROW>
end

save(fullfile(resDir, 'v1b_layer_order_10it.mat'), 'results', 'labels');
fig = ldpc_plot_results(results, labels, ...
    sprintf('V1b: порядок слоёв, layered BP, %s, BPSK, I_{max} = %d', ...
    results(1).codeName, results(1).params.maxIter));
if isprop(fig, 'Theme'), fig.Theme = 'light'; end   % светлая тема независимо от настроек MATLAB
drawnow; for a = findall(fig, 'Type', 'axes').', a.Toolbar.Visible = 'off'; end   % без панели инструментов на картинке
exportgraphics(fig, fullfile(resDir, 'v1b_layer_order_10it.png'), 'Resolution', 150);
