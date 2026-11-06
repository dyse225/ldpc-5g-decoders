function data = run_e3_nms_iterations(mode)
% NMS [5 10 20 50 100] iterations, alpha selected from matching-mode sweep.
if nargin<1, mode='trial'; end
data=maxim_run('nms_iterations',mode);
end
