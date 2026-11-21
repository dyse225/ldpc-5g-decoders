function qli_csv(file, header, rows)
%QLI_CSV Запись CSV (числа — %.10g, текст как есть).
fid = fopen(file, 'w');
if fid < 0, error('qli_csv:open', 'Не удалось открыть %s', file); end
c = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(header, ','));
for r = 1:size(rows, 1)
    parts = cell(1, size(rows, 2));
    for j = 1:size(rows, 2)
        v = rows{r, j};
        if ischar(v), parts{j} = v;
        elseif isnan(v), parts{j} = 'NaN';
        else, parts{j} = sprintf('%.10g', v); end
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end
end
