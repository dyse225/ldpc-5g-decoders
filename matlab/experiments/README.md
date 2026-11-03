# experiments

Здесь лежат эксперименты, каждый в своей папке: `maxim/`, `dima/`, `matvey/`, `vova/` и т. д.

- Общие файлы `matlab/ldpc_*.m` без согласования не меняем.
- Коммитим скрипты, `.mat` с результатами и графики (`.png`). Сохраняйте их в своей папке, например в `<имя>/results/`, а не в общую `matlab/results/` (она в git не попадает).
- Чтобы общие функции были видны из своей папки, в начале скрипта добавьте:

```matlab
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));
```

Пример сохранения в конце скрипта:

```matlab
resDir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
save(fullfile(resDir, 'nms_alpha.mat'), 'results', 'labels');
exportgraphics(fig, fullfile(resDir, 'nms_alpha.png'), 'Resolution', 150);
```
