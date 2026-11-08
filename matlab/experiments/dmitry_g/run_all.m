function run_all(maxFrames, profile, useParallel, bfMaxIter)
%RUN_ALL Последовательно запускает три независимых эксперимента.
% run_all(2e5,'long')          BG1 3872 / BG2 3840 и моделирования 2, 3
% run_all(2e5,'short')         BG1 528 / BG2 520 и моделирования 2, 3
% run_all(2e5,'long',false)    Без параллельного пула
% run_all(2e5,'long',true,50)  Проверка BF с пределом 50 вместо 300

if nargin<1, maxFrames=2e5; end
if nargin<2, profile='long'; end
if nargin<3, useParallel=true; end
if nargin<4, bfMaxIter=[]; end
here=fileparts(mfilename('fullpath'));
addpath(fullfile(here,'моделирование 1'),fullfile(here,'моделирование 2'),fullfile(here,'моделирование 3'));
run_e1_bg_compare(maxFrames,profile,useParallel,bfMaxIter);
run_e2_bp_iterations(maxFrames,useParallel);
run_e3_termination(maxFrames,useParallel);
end
