function run_all()
    addpath('../../')
    resDir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~exist(resDir, 'dir'), mkdir(resDir); end
    
    if isempty(gcp('nocreate'))
        parpool(4);
    end

    p = ldpc_default_params();
    p.modulation = '64qam';
    p.EbN0dB = 2:0.25:8;   
    % BG1
    p.bgn = 1;
    p.Zc = 24;
    run_iter_alg(p, resDir)
    %BG2
    p.bgn = 2; 
    p.Zc = 48; 
    run_iter_alg(p, resDir)
end

function run_iter_alg(p, resDir)
    algs   = {'bp', 'layered-bp', 'norm-min-sum', 'offset-min-sum'};
    labels = {'BP (flooding)', 'Layered BP', 'NMS \alpha=0.75', 'OMS \beta=0.5'};
    
    resCell = cell(1, numel(algs));
    parfor i = 1:numel(algs)
        p_copy = p; 
        p_copy.algorithm = algs{i};
        resCell{i} = ldpc_simulate(p_copy); % тут лажа
    end
    results = [resCell{:}];
    titleStr = sprintf('%s, 5G NR LDPC BG%d (Zc=%d, K=%d), R=%.2g, AWGN, I_{max} = %d', ...
    upper(p.modulation), p.bgn, p.Zc, results(1).K, results(1).params.rate, p.maxIter);

    fig = ldpc_plot_results(results, labels, titleStr);
    tag = sprintf('bg%d_%s', p.bgn, p.modulation);
    save(fullfile(resDir, sprintf('e1_algorithms_%s.mat', tag)), 'results', 'labels');

    exportgraphics(fig, fullfile(resDir, sprintf('e1_algorithms_%s.png', tag)), 'Resolution', 150);
    close(fig);
end