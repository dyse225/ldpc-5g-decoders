%% Моделирование 1: новый расчёт BG1 и BG2 с maxFrames=1e6.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(here, fileparts(fileparts(fileparts(here))));
resDir = fullfile(here, 'results_1e6');
for bgn = [1 2]
    run_e1_bg1_algorithms(bgn);
    run_e1b_bg1_bitflip(bgn);
end
S1 = load(fullfile(resDir, 'e1b_bg1_bitflip.mat'));
S2 = load(fullfile(resDir, 'e1b_bg2_bitflip.mat'));
labels = S1.labels;
markers = {'o','s','d','^','v'};
for metric = {'ber', 'fer'}
    fig = figure('Color','w','Position',[60 60 1400 560]);
    tl = tiledlayout(fig,1,2,'TileSpacing','compact');
    for bg = 1:2
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        set(ax,'YScale','log','FontSize',12);
        if bg == 1, rr = S1.results; else, rr = S2.results; end
        for i = 1:numel(rr)
            v = rr(i).(metric{1}); v(v == 0) = NaN;
            plot(ax,rr(i).EbN0dB,v,['-' markers{i}], ...
                'LineWidth',1.6,'DisplayName',labels{i});
        end
        xlim(ax,[0 9]); ylim(ax,[1e-6 1]);
        xlabel(ax,'E_b/N_0, dB'); ylabel(ax,upper(metric{1}));
        title(ax,sprintf('BG%d, K=%d, R=1/2',bg,rr(1).K));
        legend(ax,'Location','southwest');
    end
    title(tl,'BPSK / AWGN; maxIter=50; maxFrames=1e6');
    exportgraphics(fig,fullfile(resDir,['bg1_vs_bg2_' metric{1} '.png']),'Resolution',180);
    close(fig);
end
% Логарифмическая интерполяция между соседними точками без экстраполяции и подстановки нулевых BER.
targets = {'ber',1e-2; 'ber',1e-3; 'ber',1e-4; 'fer',1e-2};
fid = fopen(fullfile(resDir,'bg_comparison.csv'),'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'metric,target,algorithm,BG1_dB,BG2_dB,BG1_minus_BG2_dB\n');
for t = 1:size(targets,1)
    for i = 1:5
        e1 = ebno_at_target(S1.results(i),targets{t,1},targets{t,2});
        e2 = ebno_at_target(S2.results(i),targets{t,1},targets{t,2});
        fprintf(fid,'%s,%g,%s,%.6g,%.6g,%.6g\n',targets{t,1},targets{t,2}, ...
            labels{i},e1,e2,e1-e2);
    end
end

% Исходные численные результаты для всех пяти алгоритмов и обоих базовых графов.
rows = [];
for bg=1:2
    if bg==1, rr=S1.results; else, rr=S2.results; end
    for i=1:numel(rr)
        r=rr(i);
        rows=[rows; [repmat(bg,numel(r.ber),1),repmat(i,numel(r.ber),1), ...
            r.EbN0dB(:),r.ber(:),r.fer(:),r.avgIter(:),r.numFrames(:), ...
            r.numBitErrors(:),r.numFrameErrors(:)]]; %#ok<AGROW>
    end
end
T=array2table(rows,'VariableNames',{'BG','algorithm','EbN0dB','BER','FER','avgIter','numFrames','numBitErrors','numFrameErrors'});
writetable(T,fullfile(resDir,'bg_raw_measurements.csv'));
