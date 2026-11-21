%% QPSK, 5G NR LDPC (BG2 Zc=52/K=520 и BG1 Zc=176/K=3872): BP, Layered BP, NMS, OMS, I_max = 5,10,20,50
% Запуск:  >> run_all
% Нужны: MATLAB, Communications Toolbox, 5G Toolbox; общие ldpc_simulate.m и ldpc_default_params.m
% лежат в matlab/ (папка добавляется в путь ниже).
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..'));

cfg = qli_config('standard');      % основной прогон для итоговых графиков
% Для быстрой проверки всей цепочки используйте: cfg = qli_config('quick');
% cfg = qli_config('quick');

% Список графов берётся из cfg.bgList ([2 1] для полного прогона).
cfg.parallel = true;              % parfor (Parallel Computing Toolbox)

% Этап 0. Проверка convertSNR и LLR для QPSK (секунды)
ok = qli_check_llr();
if ~ok, error('Проверка LLR/SNR не пройдена: см. вывод выше.'); end

% Этап 1. Оптимальные alpha (NMS) и beta (OMS)
qli_find_params(cfg);

% Этап 2. Основной свип: 4 декодера x I_max, uncoded QPSK
runs = qli_run(cfg);

% Этап 3. Графики
qli_plots(runs, cfg);
% Перерисовать без пересчёта:  S = load(fullfile(cfg.resDir,'all_results.mat')); qli_plots(S.runs, S.cfg);

%%
cfg = qli_config('standard');

figDir = fullfile(cfg.resDir, 'figures');

disp(cfg.resDir);
disp(figDir);

if ~exist(figDir, 'dir')
    [status, msg] = mkdir(figDir);
    disp(status);
    disp(msg);
end

[fid, msg] = fopen(fullfile(figDir, 'test_write.txt'), 'w');

if fid == -1
    fprintf('ОШИБКА записи: %s\n', msg);
else
    fprintf(fid, 'test');
    fclose(fid);
    delete(fullfile(figDir, 'test_write.txt'));
    disp('Проверка записи пройдена.');
end
