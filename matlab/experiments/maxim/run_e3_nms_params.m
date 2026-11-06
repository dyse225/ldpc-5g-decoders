function data = run_e3_nms_params(mode)
% NMS alpha sweep. Run with 'trial' first, then 'final'.
if nargin<1, mode='trial'; end
data=maxim_run('nms_alpha',mode);
end
