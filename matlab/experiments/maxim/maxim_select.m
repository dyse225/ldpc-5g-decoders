function selection = maxim_select(data)
%MAXIM_SELECT Explicit exploratory ranking, not a significance test.
% Rank by geometric mean BER across the prescribed, equally weighted grid.
% A zero-error estimate is replaced ONLY for ranking by half an error/nbits.
ber = cell2mat(cellfun(@(r) r.ber(:),data.results,'UniformOutput',false));
used = ber;
for i=1:numel(data.values)
    r = data.results{i};
    zero = used(:,i)==0;
    floorValue = 0.5./(r.K*r.numFrames(:));
    used(zero,i) = floorValue(zero);
end
scores = mean(log10(used),1);
[~,idx] = min(scores);
selection.value = data.values(idx); selection.index = idx;
selection.meanLog10BER = scores;
selection.geometricMeanBER = 10.^scores;
selection.pointwiseBest = min(ber,[],2);
selection.hasZeroEstimates = any(ber(:)==0);
selection.rule = 'Minimum equally weighted mean log10 BER on 0:0.5:2 dB; zero estimates use 0.5/nbits for ranking only. Inspect crossings and sampling uncertainty.';
end
