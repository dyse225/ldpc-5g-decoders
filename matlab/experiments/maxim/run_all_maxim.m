function run_all_maxim(mode)
%RUN_ALL_MAXIM Four requested experiments, resume saved complete curves.
if nargin<1, mode='trial'; end
run_e3_nms_params(mode);
run_e3_oms_params(mode);
run_e3_nms_iterations(mode);
run_e3_oms_iterations(mode);
fprintf('ALL_FOUR_COMPLETED %s\n',mode);
end
