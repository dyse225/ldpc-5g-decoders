function data = maxim_run(kind, mode)
%MAXIM_RUN Reproducible parameter/iteration sweep using the shared simulator.
% mode: 'trial' (2e4 frames), 'final' (1e6), 'smoke' (128).
% Shared ldpc_*.m files are deliberately left unchanged.
if nargin < 2, mode = 'trial'; end
mode = validatestring(mode, {'smoke','trial','final'});
kind = validatestring(kind, {'nms_alpha','oms_beta','nms_iterations','oms_iterations'});
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));
resDir = fullfile(here, 'results', mode);
if ~exist(resDir,'dir'), mkdir(resDir); end

p = ldpc_default_params();
p.code = '5gnr'; p.bgn = 2; p.Zc = 52; p.rate = 1/2;
p.modulation = 'bpsk'; p.EbN0dB = 0:0.5:2;
p.maxIter = 50; p.termination = 'early';
p.stopAtZeroFer = false; p.seed = 1; p.verbose = true;
p.minFrameErrors = 100;
switch mode
    case 'smoke', p.maxFrames = 128; p.minFrameErrors = 10;
    case 'trial', p.maxFrames = 2e4;
    case 'final', p.maxFrames = 1e6;
end

selection = struct();
switch kind
    case 'nms_alpha'
        p.algorithm = 'norm-min-sum'; field = 'scaling';
        values = [0.50 0.60 0.70 0.75 0.80 0.85 0.90 0.95];
        labelPrefix = '\alpha';
    case 'oms_beta'
        p.algorithm = 'offset-min-sum'; field = 'offset';
        values = [0.1 0.2 0.3 0.4 0.5]; labelPrefix = '\beta';
    otherwise
        values = [5 10 20 50 100]; field = 'maxIter'; labelPrefix = 'I_{max}';
        if startsWith(kind,'nms')
            sourceKind = 'nms_alpha'; p.algorithm = 'norm-min-sum';
            paramField = 'scaling';
        else
            sourceKind = 'oms_beta'; p.algorithm = 'offset-min-sum';
            paramField = 'offset';
        end
        source = fullfile(resDir,[sourceKind '.mat']);
        assert(isfile(source), 'Run the matching parameter sweep first: %s',source);
        loaded = load(source,'data'); ref = loaded.data;
        assert(all(ref.completed), 'The parameter sweep is incomplete.');
        selection = maxim_select(ref);
        p.(paramField) = selection.value;
        fprintf('Selected %s = %g (minimum mean log10 BER over the five points).\n', ...
            paramField,selection.value);
end

data.kind = kind; data.mode = mode; data.values = values;
data.baseParams = p; data.EbN0dB = p.EbN0dB;
data.selection = selection; data.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
data.matlabVersion = version; data.repositoryCommit = '9c2f1ca6c735b8c8bbc3ccc2eab96bdcb8afd258';
data.results = cell(1,numel(values)); data.labels = cell(1,numel(values));
data.completed = false(1,numel(values));
out = fullfile(resDir,[kind '.mat']);
if isfile(out)
    old = load(out,'data');
    assert(isequaln(old.data.baseParams,p) && isequal(old.data.values,values), ...
        'Settings differ from saved results. Move the old file before rerunning.');
    data = old.data;
end
for i = 1:numel(values)
    if data.completed(i), continue; end
    q = p; q.(field) = values(i);
    fprintf('\n=== %s / %s: %s = %g (%d/%d) ===\n', ...
        mode,kind,field,values(i),i,numel(values));
    t = tic;
    reused = false; reusedParameter = false;
    if endsWith(kind,'iterations') && values(i)==50
        candidate = ref.results{selection.index};
        same = isequaln(rmfield(candidate.params,'maxFrames'),rmfield(q,'maxFrames'));
        if same && (candidate.params.maxFrames==q.maxFrames || ...
                all(candidate.numFrameErrors>=q.minFrameErrors))
            r=candidate; reusedParameter=true;
            fprintf('Reusing the identical 50-iteration curve from parameter selection.\n');
        end
    end
    trialFile = fullfile(here,'results','trial',[kind '.mat']);
    if ~reusedParameter && strcmp(mode,'final') && isfile(trialFile)
        trial = load(trialFile,'data');
        if trial.data.completed(i)
            candidate = trial.data.results{i};
            % A larger frame cap cannot affect a point that already reached
            % the error-count stopping criterion. Preserve original params.
            same = isequaln(rmfield(candidate.params,'maxFrames'),rmfield(q,'maxFrames'));
            if same && all(candidate.numFrameErrors >= q.minFrameErrors)
                r = candidate; reused = true;
                fprintf('Reusing trial curve: every point already reached %d frame errors.\n',q.minFrameErrors);
            end
        end
    end
    if ~reused && ~reusedParameter, r = ldpc_simulate(q); end
    assert(r.K==520 && r.E==1040 && all(r.numFrames>0));
    assert(all(isfinite(r.ber)) && all(r.numBitErrors <= r.numFrameErrors*r.K));
    assert(all(r.numBitErrors >= r.numFrameErrors));
    assert(all(abs(r.ber-r.numBitErrors./(r.numFrames*r.K))<1e-14));
    assert(all(abs(r.fer-r.numFrameErrors./r.numFrames)<1e-14));
    assert(all(r.avgIter>=1 & r.avgIter<=q.maxIter));
    data.results{i} = r;
    data.labels{i} = sprintf('%s = %g',labelPrefix,values(i));
    data.completed(i) = true;
    data.reusedTrial(i) = reused;
    data.reusedParameter(i) = reusedParameter;
    data.elapsedSeconds(i) = toc(t);
    save(out,'data','-v7'); % Checkpoint after every complete curve.
end
data.ber = cell2mat(cellfun(@(r) r.ber(:),data.results,'UniformOutput',false));
data.fer = cell2mat(cellfun(@(r) r.fer(:),data.results,'UniformOutput',false));
data.avgIter = cell2mat(cellfun(@(r) r.avgIter(:),data.results,'UniformOutput',false));
if endsWith(kind,'alpha') || endsWith(kind,'beta'), data.selection = maxim_select(data); end
save(out,'data','-v7');
maxim_plot(data,resDir);
fprintf('COMPLETED: %s / %s\n',mode,kind);
end
