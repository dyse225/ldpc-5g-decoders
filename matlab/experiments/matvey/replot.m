%% Перерисовать графики из сохранённых results/*.mat без пересчёта
here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
fig = plot_lsat_sweep(fullfile(resDir, 'lsat_sweep.mat'));
exportgraphics(fig, fullfile(resDir, 'lsat_sweep.png'), 'Resolution', 150);
fig = plot_qms_vs_nms_oms(fullfile(resDir, 'qms_vs_nms_oms.mat'));
exportgraphics(fig, fullfile(resDir, 'qms_vs_nms_oms.png'), 'Resolution', 150);
fig = plot_minsum_params_16qam(fullfile(resDir, 'minsum_params_16qam.mat'));
exportgraphics(fig, fullfile(resDir, 'minsum_params_16qam.png'), 'Resolution', 150);
