function data = run_e3_oms_iterations(mode)
% OMS [5 10 20 50 100] iterations, beta selected from matching-mode sweep.
if nargin<1, mode='trial'; end
data=maxim_run('oms_iterations',mode);
end
