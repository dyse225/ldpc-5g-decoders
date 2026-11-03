%% E0: проверка канала и демодулятора — некодированная BPSK/QPSK/16QAM vs теория
clear; clc;
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

mods = {'bpsk', 'qpsk', '16qam'};
theory = {{'psk', 2}, {'psk', 4}, {'qam', 16}};
EbN0dB = 0:0.5:5;       % рабочая область LDPC; выше — лишняя трата времени

fig = figure('Color', 'w');
semilogy(NaN, NaN); hold on; grid on;
for i = 1:numel(mods)
    p = ldpc_default_params();
    p.code = 'none';
    p.modulation = mods{i};
    p.EbN0dB = EbN0dB;
    p.minFrameErrors = 100;
    p.maxFrames = 2e3;
    r = ldpc_simulate(p);

    berTh = berawgn(EbN0dB, theory{i}{1}, theory{i}{2}, 'nondiff');
    h = semilogy(EbN0dB, berTh, '-', 'LineWidth', 1.2, ...
        'DisplayName', [upper(mods{i}) ' теория']);
    semilogy(r.EbN0dB, r.ber, 'o', 'Color', h.Color, ...
        'DisplayName', [upper(mods{i}) ' модель']);
end
xlabel('E_b/N_0, дБ'); ylabel('BER'); ylim([1e-3 1]);
legend('Location', 'southwest'); title('E0: некодированная передача в AWGN');
exportgraphics(fig, fullfile(resDir, 'e0_uncoded.png'), 'Resolution', 150);
