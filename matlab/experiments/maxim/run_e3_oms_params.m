function data = run_e3_oms_params(mode)
% OMS beta sweep, separate from NMS.
if nargin<1, mode='trial'; end
data=maxim_run('oms_beta',mode);
end
