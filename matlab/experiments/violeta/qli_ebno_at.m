function eb = qli_ebno_at(res, metric, target)
%QLI_EBNO_AT Eb/N0 первого пересечения уровня target кривой metric ('ber'|'fer')
%   (линейная интерполяция по log10, без экстраполяции; NaN — уровень не достигнут).
v = res.(metric); eb = NaN;
for k = 1:numel(v) - 1
    if all(res.numFrames(k:k+1) > 0) && all(isfinite(v(k:k+1))) && all(v(k:k+1) > 0) ...
            && v(k) >= target && v(k+1) <= target
        if v(k) == v(k+1), eb = res.EbN0dB(k); return; end
        f = (log10(target) - log10(v(k))) / (log10(v(k+1)) - log10(v(k)));
        eb = res.EbN0dB(k) + f * (res.EbN0dB(k+1) - res.EbN0dB(k));
        return
    end
end
end
