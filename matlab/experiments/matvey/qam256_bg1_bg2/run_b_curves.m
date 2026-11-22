%% B: BER/FER от Eb/N0 для BP, layered BP, NMS, OMS на 256-QAM, BG1 и BG2, I_max = 5, 10, 20, 50
% alpha и beta для NMS/OMS берутся оптимальные для каждого I_max из эксперимента A.
clear; clc;
here = fileparts(mfilename('fullpath'));
sharedDir = fullfile(here, '..', '..', '..');
addpath(sharedDir, here);
resDir = fullfile(here, 'results');

[alphaOpt, betaOpt] = pick_params(fullfile(resDir, 'a_params.mat'));
A = load(fullfile(resDir, 'a_params.mat'));
codes = A.codes; iterList = A.iterList;
algs = {'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum'};

quick = false;
[C, G, I] = ndgrid(1:size(codes,1), 1:numel(algs), 1:numel(iterList));
nJobs = numel(C);
ps = cell(nJobs,1);
for k = 1:nJobs
    p = ldpc_default_params();
    p.bgn = codes(C(k),1); p.Zc = codes(C(k),2);
    p.modulation = '256qam'; p.rate = 1/2;
    p.algorithm = algs{G(k)}; p.maxIter = iterList(I(k));
    p.scaling = alphaOpt(C(k), I(k)); p.offset = betaOpt(C(k), I(k));
    p.EbN0dB = 8.5:0.25:11;
    p.minFrameErrors = 100; p.maxFrames = 3e4; if quick, p.maxFrames = 500; end
    p.verbose = false;
    ps{k} = p;
end
if isempty(gcp('nocreate')), parpool('Processes', 6); end
pctRunOnAll(sprintf('addpath(''%s'', ''%s'')', strrep(sharedDir,'\','/'), strrep(here,'\','/')));
out = cell(nJobs,1);
parfor k = 1:nJobs
    out{k} = ldpc_simulate_ext(ps{k});
end
res = reshape(out, size(C));          % res{код, алгоритм, номер I_max}
save(fullfile(resDir, 'b_curves.mat'), 'res', 'algs', 'iterList', 'codes', 'alphaOpt', 'betaOpt');
disp('B done');
