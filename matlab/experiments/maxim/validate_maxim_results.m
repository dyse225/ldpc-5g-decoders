function validate_maxim_results(mode)
%VALIDATE_MAXIM_RESULTS Check saved counters/settings, not curve appearance.
if nargin<1, mode='final'; end
here=fileparts(mfilename('fullpath'));
kinds={'nms_alpha','oms_beta','nms_iterations','oms_iterations'};
expected={[.50 .60 .70 .75 .80 .85 .90 .95],[.1 .2 .3 .4 .5], ...
    [5 10 20 50 100],[5 10 20 50 100]};
total=0;
for k=1:numel(kinds)
    f=fullfile(here,'results',mode,[kinds{k} '.mat']);
    loaded=load(f,'data'); d=loaded.data;
    assert(all(d.completed) && isequal(d.values,expected{k}));
    assert(isequal(d.EbN0dB,0:.5:2));
    for i=1:numel(d.values)
        r=d.results{i}; q=r.params;
        assert(r.K==520 && r.N==2600 && r.E==1040);
        assert(q.bgn==2 && q.Zc==52 && q.rate==.5 && strcmp(q.modulation,'bpsk'));
        assert(strcmp(q.termination,'early') && ~q.stopAtZeroFer && q.seed==1);
        assert(all(isfinite(r.ber)) && all(isfinite(r.fer)));
        assert(all(r.numFrames>0 & r.numFrames<=d.baseParams.maxFrames));
        assert(all(r.numFrameErrors>=q.minFrameErrors | ...
            r.numFrames==d.baseParams.maxFrames));
        assert(all(r.numBitErrors>=r.numFrameErrors));
        assert(all(abs(r.ber-r.numBitErrors./(520*r.numFrames))<1e-14));
        assert(all(abs(r.fer-r.numFrameErrors./r.numFrames)<1e-14));
        assert(all(abs(d.ber(:,i)-r.ber(:))<1e-14));
        assert(all(r.avgIter>=1 & r.avgIter<=q.maxIter));
        if endsWith(kinds{k},'iterations'), assert(q.maxIter==d.values(i));
        else, assert(q.maxIter==50); end
        if startsWith(kinds{k},'nms')
            assert(strcmp(q.algorithm,'norm-min-sum'));
            if k==1, assert(q.scaling==d.values(i));
            else, assert(q.scaling==d.selection.value); end
        else
            assert(strcmp(q.algorithm,'offset-min-sum'));
            if k==2, assert(q.offset==d.values(i));
            else, assert(q.offset==d.selection.value); end
        end
        total=total+numel(r.ber);
    end
    for name={'ber','fer','avgIter'}
        assert(isfile(fullfile(here,'results',mode,[kinds{k} '_' name{1} '.png'])));
    end
    fprintf('PASS %s: %d complete curves\n',kinds{k},numel(d.values));
end
assert(total==115);
fprintf('RESULT_VALIDATION_PASSED: %d points / 23 curves / 4 experiments (%s)\n',total,mode);
end
