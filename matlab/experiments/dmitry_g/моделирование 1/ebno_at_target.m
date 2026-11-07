function eb = ebno_at_target(r, metric, target)
% Первое пересечение заданного уровня сверху вниз между соседними рассчитанными точками; логарифмическая интерполяция.
v = r.(metric); eb = NaN;
for k = 1:numel(v)-1
    if all(r.numFrames(k:k+1)>0) && all(isfinite(v(k:k+1))) && ...
            all(v(k:k+1)>0) && v(k)>=target && v(k+1)<=target
        if v(k)==v(k+1), eb=r.EbN0dB(k); return; end
        f=(log10(target)-log10(v(k)))/(log10(v(k+1))-log10(v(k)));
        eb=r.EbN0dB(k)+f*(r.EbN0dB(k+1)-r.EbN0dB(k)); return
    end
end
end
