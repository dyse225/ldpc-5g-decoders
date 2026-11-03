function fig = ldpc_plot_results(results, labels, titleStr, opts)
%LDPC_PLOT_RESULTS BER, FER и среднее число итераций для набора прогонов.
%   fig = LDPC_PLOT_RESULTS(results, labels, titleStr, opts)
%     results — массив структур из LDPC_SIMULATE
%     labels  — cell-массив подписей кривых
%     opts    — (необязательно) структура:
%       .layout   'row' (по умолчанию, три графика в ряд) |
%                 'column' (друг под другом, широкие)
%       .xlim     пределы по Eb/N0, например [0 9]
%       .berLim   пределы по BER, например [1e-5 1]
%       .ferLim   пределы по FER, например [1e-4 1]
%       .xbreak   (только 'column') разрыв оси Eb/N0 [x1 x2]: слева
%                 xlim(1)..x1 на ~60% ширины, справа x2..xlim(2)
%       .colors   (необязательно) N x 3 цвета кривых

if nargin < 3, titleStr = ''; end
if nargin < 4, opts = struct(); end
if ~isfield(opts, 'layout'), opts.layout = 'row'; end
if strcmp(opts.layout, 'column') && isfield(opts, 'xbreak')
    fig = plot_column_broken(results, labels, titleStr, opts);
    return
end
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};

if strcmp(opts.layout, 'column')
    fig = figure('Color', 'w', 'Position', [50 50 1600 1000]);
    tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    lw = 2; ms = 7; fs = 14;
else
    fig = figure('Color', 'w', 'Position', [100 100 1200 380]);
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact');
    lw = 1.3; ms = 6; fs = 10;
end
if ~isempty(titleStr), title(tl, titleStr, 'FontSize', fs + 2); end

ax = gobjects(1, 3);
ylab = {'BER', 'FER', 'Ср. число итераций'};
for k = 1:3
    ax(k) = nexttile(tl);
    hold(ax(k), 'on'); grid(ax(k), 'on'); box(ax(k), 'on');
    ylabel(ax(k), ylab{k});
    set(ax(k), 'FontSize', fs);
end
xlabel(ax(3), 'E_b/N_0, дБ');
if strcmp(opts.layout, 'row')
    xlabel(ax(1), 'E_b/N_0, дБ'); xlabel(ax(2), 'E_b/N_0, дБ');
end
set(ax(1:2), 'YScale', 'log');
if isfield(opts, 'xlim'), set(ax, 'XLim', opts.xlim); end
if isfield(opts, 'berLim'), ylim(ax(1), opts.berLim); end
if isfield(opts, 'ferLim'), ylim(ax(2), opts.ferLim); end

for i = 1:numel(results)
    r = results(i);
    m = markers{mod(i-1, numel(markers)) + 1};
    ok = r.numFrames > 0;
    berPlot = r.ber; berPlot(berPlot == 0) = NaN;
    ferPlot = r.fer; ferPlot(ferPlot == 0) = NaN;
    args = {['-' m], 'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', labels{i}};
    plot(ax(1), r.EbN0dB(ok), berPlot(ok), args{:});
    plot(ax(2), r.EbN0dB(ok), ferPlot(ok), args{:});
    plot(ax(3), r.EbN0dB(ok), r.avgIter(ok), args{:});
end
if strcmp(opts.layout, 'column')
    lg = legend(ax(1), 'Location', 'eastoutside');
else
    lg = legend(ax(2), 'Location', 'southeast');
end
lg.FontSize = fs;
end

% -------------------------------------------------------------------------
function fig = plot_column_broken(results, labels, titleStr, opts)
% Три широких графика друг под другом, ось Eb/N0 с разрывом.
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
colors = get(groot, 'defaultAxesColorOrder');
if isfield(opts, 'colors'), colors = opts.colors; end
lw = 2; ms = 7; fs = 14;
nL = 3; nR = 2;                                   % доли ширины: 3/5 и 2/5

fig = figure('Color', 'w', 'Position', [50 50 1600 1000]);
tl = tiledlayout(fig, 3, nL + nR, 'TileSpacing', 'compact', 'Padding', 'compact');
if ~isempty(titleStr), title(tl, titleStr, 'FontSize', fs + 2); end

ylab = {'BER', 'FER', 'Ср. число итераций'};
yl = {[], [], [0 52]};
if isfield(opts, 'berLim'), yl{1} = opts.berLim; end
if isfield(opts, 'ferLim'), yl{2} = opts.ferLim; end
xl = [opts.xlim(1) opts.xbreak(1); opts.xbreak(2) opts.xlim(2)];

for k = 1:3
    axL = nexttile(tl, [1 nL]);
    axR = nexttile(tl, [1 nR]);
    for ax = [axL axR]
        hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
        set(ax, 'FontSize', fs);
        if k < 3, set(ax, 'YScale', 'log'); end
        if ~isempty(yl{k}), ylim(ax, yl{k}); end
    end
    xlim(axL, xl(1, :)); xlim(axR, xl(2, :));
    ylabel(axL, ylab{k});
    set(axR, 'YTickLabel', []);
    if k == 3
        xlabel(axL, 'E_b/N_0, дБ'); xlabel(axR, 'E_b/N_0, дБ');
    end

    for i = 1:numel(results)
        r = results(i);
        ok = r.numFrames > 0;
        switch k
            case 1, v = r.ber;
            case 2, v = r.fer;
            otherwise, v = r.avgIter;
        end
        v(v == 0) = NaN;
        args = {['-' markers{mod(i-1, numel(markers)) + 1}], ...
            'Color', colors(mod(i-1, size(colors, 1)) + 1, :), ...
            'LineWidth', lw, 'MarkerSize', ms, 'DisplayName', labels{i}};
        plot(axL, r.EbN0dB(ok), v(ok), args{:});
        plot(axR, r.EbN0dB(ok), v(ok), args{:}, 'HandleVisibility', 'off');
    end
    if k == 3                                     % на панели итераций справа пусто
        lg = legend(axL, 'Location', 'northeast');
        lg.FontSize = fs;
    end
end
end
