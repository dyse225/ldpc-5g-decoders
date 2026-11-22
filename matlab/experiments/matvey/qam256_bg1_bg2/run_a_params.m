%% A: подбор alpha (NMS) и beta (OMS) на 256-QAM для разного I_max, BG1 и BG2
% Одна точка Eb/N0 = 9.5 дБ (FER порядка 0.1 при I_max = 50).
clear; clc;
here = fileparts(mfilename('fullpath'));
sharedDir = fullfile(here, '..', '..', '..');
addpath(sharedDir, here);
resDir = fullfile(here, 'results'); if ~exist(resDir, 'dir'), mkdir(resDir); end

codes    = [1 24; 2 52];            % [bgn Zc]: BG1 K=528, BG2 K=520
iterList = [5 10 20 50];
alphaList = 0.5:0.05:1.0;           % alpha = 1 - обычный min-sum
betaList  = [0 0.2 0.3:0.05:0.7 0.8 1 1.5];   % beta = 0 - обычный min-sum; шаг 0.05 вокруг минимума
ebno = 9.5;

jobs = {};   % {код, I, алгоритм, параметр}
for c = 1:size(codes,1)
  for I = iterList
    for a = alphaList, jobs(end+1,:) = {c, I, 'norm-min-sum', a}; end %#ok<SAGROW>
    for b = betaList,  jobs(end+1,:) = {c, I, 'offset-min-sum', b}; end %#ok<SAGROW>
  end
end
nJobs = size(jobs,1);
ps = cell(nJobs,1);
for k = 1:nJobs
    p = ldpc_default_params();
    p.bgn = codes(jobs{k,1},1); p.Zc = codes(jobs{k,1},2);
    p.modulation = '256qam'; p.rate = 1/2;
    p.maxIter = jobs{k,2}; p.algorithm = jobs{k,3};
    if strcmp(p.algorithm,'norm-min-sum'), p.scaling = jobs{k,4}; else, p.offset = jobs{k,4}; end
    p.EbN0dB = ebno; p.minFrameErrors = 200; p.maxFrames = 8000;
    p.stopAtZeroFer = false; p.verbose = false;
    ps{k} = p;
end
if isempty(gcp('nocreate')), parpool('Processes', 6); end
pctRunOnAll(sprintf('addpath(''%s'', ''%s'')', strrep(sharedDir,'\','/'), strrep(here,'\','/')));
out = cell(nJobs,1);
parfor k = 1:nJobs
    out{k} = ldpc_simulate_ext(ps{k});
end
save(fullfile(resDir, 'a_params.mat'), 'jobs', 'out', 'codes', 'iterList', 'alphaList', 'betaList', 'ebno');
disp('A done');
