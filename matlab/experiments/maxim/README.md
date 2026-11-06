# maxim: параметры NMS и OMS

Подбор α для NMS, β для OMS и проверка влияния числа итераций.
Используется общая `ldpc_simulate` и декодеры MATLAB `nrLDPCDecode`.

5G NR BG2, Zc=52, K=520, R=1/2, BPSK, AWGN.
Eb/N0 = 0:0.5:2 дБ. Подбор параметров при maxIter=50, досрочная остановка,
seed=1.

| Файл | Что делает |
|---|---|
| [run_e3_nms_params.m](run_e3_nms_params.m) | подбор α: 0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 |
| [run_e3_oms_params.m](run_e3_oms_params.m) | подбор β: 0.1, 0.2, 0.3, 0.4, 0.5 |
| [run_e3_nms_iterations.m](run_e3_nms_iterations.m) | итерации NMS: 5, 10, 20, 50, 100 при α=0.70 |
| [run_e3_oms_iterations.m](run_e3_oms_iterations.m) | итерации OMS: 5, 10, 20, 50, 100 при β=0.5 |
| [run_all_maxim.m](run_all_maxim.m) | запуск всех четырёх исследований по порядку |
| [maxim_run.m](maxim_run.m) | настройки, перебор значений, вызов модели и сохранение |
| [maxim_select.m](maxim_select.m) | выбор параметра по минимуму среднего log10(BER) по пяти точкам |
| [maxim_plot.m](maxim_plot.m) | графики и таблицы из результатов |
| [validate_maxim_results.m](validate_maxim_results.m) | проверка сохранённых данных без пересчёта |

Четыре файла отдельных исследований вызывают общий цикл в `maxim_run.m`.
В итерационных сериях параметр берётся из подбора того же режима.

## Запуск

MATLAB R2025b, Communications Toolbox и 5G Toolbox.
Из папки `matlab/experiments/maxim/`:

~~~matlab
run_all_maxim('trial')   % лимит 2e4 кадров на точку
run_all_maxim('final')   % лимит 1e6 кадров на точку
~~~

Отдельное исследование запускается аналогично, например
`run_e3_nms_params('final')`. Без аргумента выбран `trial`.

Точка считается до 100 ошибочных кадров либо лимита кадров.
Законченные кривые сохраняются и при повторном запуске пропускаются.
Trial-кривая используется в final только при совпадающих настройках и
достижении 100 ошибочных кадров во всех точках; исходные параметры
сохранены в MAT. Если настройки изменены, старые результаты нужно
предварительно переместить в другую папку.

## Результаты

| Исследование | График BER | Данные MATLAB | Таблица |
|---|---|---|---|
| α NMS | [PNG](results/final/nms_alpha_ber.png) | [MAT](results/final/nms_alpha.mat) | [CSV](results/final/nms_alpha_counts.csv) |
| β OMS | [PNG](results/final/oms_beta_ber.png) | [MAT](results/final/oms_beta.mat) | [CSV](results/final/oms_beta_counts.csv) |
| Итерации NMS | [PNG](results/final/nms_iterations_ber.png) | [MAT](results/final/nms_iterations.mat) | [CSV](results/final/nms_iterations_counts.csv) |
| Итерации OMS | [PNG](results/final/oms_iterations_ber.png) | [MAT](results/final/oms_iterations.mat) | [CSV](results/final/oms_iterations_counts.csv) |

В `results/final/` также есть `*_fer.png` и `*_avgIter.png`.
MAT содержит структуру `data` с настройками и результатами каждой кривой;
CSV — BER, FER, итерации и счётчики кадров и ошибок.

Модель использует исходное упрощённое согласование скорости без
filler-бит и битового перемежения.
