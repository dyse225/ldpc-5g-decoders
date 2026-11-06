function maxim_plot(data,resDir)
%MAXIM_PLOT Separate BER figure for every experiment; retain FER/iterations.
for metric = {'ber','fer','avgIter'}
    name = metric{1};
    fig = figure('Visible','off','Color','w','Position',[100 100 1100 720]);
    ax = axes(fig,'Color','w','XColor','k','YColor','k', ...
        'GridColor',[0.65 0.65 0.65],'MinorGridColor',[0.8 0.8 0.8]); hold(ax,'on');
    colors = [0 .35 .7; .85 .3 0; .15 .55 .2; .55 .15 .65; ...
        .65 .45 0; 0 .55 .6; .8 .15 .4; .25 .25 .25];
    markers = {'o','s','^','d','v','>','<','p'};
    anyZero = false;
    for i=1:numel(data.values)
        r=data.results{i}; y=r.(name); zero=y==0;
        anyZero=anyZero||any(zero);
        y(zero)=NaN; % Never draw observed zero as an exact nonzero BER.
        if strcmp(name,'avgIter')
            plot(ax,r.EbN0dB,y,['-' markers{i}], 'Color',colors(i,:), ...
                'LineWidth',1.5,'MarkerSize',6,'DisplayName',data.labels{i});
        else
            semilogy(ax,r.EbN0dB,y,['-' markers{i}], 'Color',colors(i,:), ...
                'LineWidth',1.5,'MarkerSize',6,'DisplayName',data.labels{i});
        end
    end
    % hold-on axes need an explicit logarithmic scale.
    if ~strcmp(name,'avgIter'), set(ax,'YScale','log'); end
    grid(ax,'on'); ax.XMinorGrid='off'; ax.YMinorGrid='on'; ax.FontSize=12;
    xlabel(ax,'E_b/N_0 (dB)');
    if strcmp(name,'avgIter'), ylabel(ax,'Mean actual decoder iterations');
    else, ylabel(ax,upper(name)); end
    switch data.kind
        case 'nms_alpha', topic='NMS: normalization coefficient';
        case 'oms_beta', topic='OMS: offset';
        case 'nms_iterations'
            topic=sprintf('NMS: iteration limit, alpha = %g',data.selection.value);
        otherwise
            topic=sprintf('OMS: iteration limit, beta = %g',data.selection.value);
    end
    title(ax,{topic,'5G NR BG2 | K = 520 | R = 1/2 | BPSK | AWGN'},'Interpreter','none','Color','k');
    subtitle(ax,sprintf('%s | early stop | frame cap %.0f | target %d frame errors/point', ...
        data.mode,data.baseParams.maxFrames,data.baseParams.minFrameErrors),'Interpreter','none','Color','k');
    legend(ax,'Location','eastoutside','Color','w','TextColor','k','EdgeColor',[.6 .6 .6]);
    xlim(ax,[0 2]); xticks(ax,0:0.5:2);
    if anyZero
        text(ax,0.02,0.04,'Zero-error estimates omitted; see .mat counts.', ...
            'Units','normalized','FontSize',10,'Interpreter','none','Color','k');
    end
    exportgraphics(fig,fullfile(resDir,[data.kind '_' name '.png']),'Resolution',180);
    close(fig);
end
rows=[];
for i=1:numel(data.values)
    r=data.results{i};
    rows=[rows; [repmat(data.values(i),numel(r.ber),1),r.EbN0dB(:), ...
        r.ber(:),r.fer(:),r.avgIter(:),r.numFrames(:),r.numBitErrors(:),r.numFrameErrors(:)]]; %#ok<AGROW>
end
T=array2table(rows,'VariableNames',{'parameter','EbN0dB','BER','FER','meanIterations','frames','bitErrors','frameErrors'});
writetable(T,fullfile(resDir,[data.kind '_counts.csv']));
end
