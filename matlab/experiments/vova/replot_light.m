%% Перерисовать сохранённые результаты (results/*.mat) в светлой теме без пересчёта
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));
resDir = fullfile(here, 'results');

items = {
    'v0_check',       'V0: свой layered BP против nrLDPCDecode, %s, BPSK, I_{max} = 50'
    'v1_layer_order', 'V1: порядок слоёв, layered BP, %s, BPSK, I_{max} = 50'
    'v1b_layer_order_10it', 'V1b: порядок слоёв, layered BP, %s, BPSK, I_{max} = 10'
    'v2_iterations',  'V2: число итераций, layered BP (прямой порядок), %s, BPSK'
    };
for k = 1:size(items, 1)
    f = fullfile(resDir, [items{k, 1} '.mat']);
    if ~exist(f, 'file'), continue; end
    S = load(f, 'results', 'labels');
    fig = ldpc_plot_results(S.results, S.labels, sprintf(items{k, 2}, S.results(1).codeName));
    if isprop(fig, 'Theme'), fig.Theme = 'light'; end
    drawnow; for a = findall(fig, 'Type', 'axes').', a.Toolbar.Visible = 'off'; end   % без панели инструментов на картинке
    exportgraphics(fig, fullfile(resDir, [items{k, 1} '.png']), 'Resolution', 150);
    close(fig);
end
