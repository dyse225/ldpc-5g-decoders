function x = ebno_at_ber(ebno, ber, target)
%EBNO_AT_BER Eb/N0, при котором BER пересекает target.
%   x = EBNO_AT_BER(ebno, ber, target) - линейная интерполяция
%   log10(BER) между соседними точками вокруг первого пересечения.
%   NaN, если кривая не опустилась до target.

ok = ber > 0 & ~isnan(ber);
ebno = ebno(ok); lb = log10(ber(ok)); lt = log10(target);
i = find(lb(1:end-1) >= lt & lb(2:end) < lt, 1);
if isempty(i)
    x = NaN;
    return
end
x = ebno(i) + (lt - lb(i)) * (ebno(i+1) - ebno(i)) / (lb(i+1) - lb(i));
end
