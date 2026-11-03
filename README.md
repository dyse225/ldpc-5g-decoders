# LDPC

Моделирование декодеров LDPC в MATLAB. Код 5G NR (BG2, Zc = 52, K = 520, R = 1/2), BPSK/QPSK/16-QAM, канал AWGN.

Нужно: MATLAB R2025b, Communications Toolbox, 5G Toolbox.

Описание системы и план экспериментов — `docs/system_description.md`.

## Как запускать

Все скрипты лежат в `matlab/`. Результаты (графики и .mat) сохраняются в `matlab/results/`.

```matlab
cd matlab
run_e1_algorithms      % BP, layered BP, NMS, OMS
```

Одна точка считается до 100 ошибочных кадров (лимит `p.maxFrames = 1e6`), поэтому полный прогон идёт долго, около часа. Для пробы ставьте `p.maxFrames = 2e4`.

## Как сделать свой эксперимент

Берёте параметры по умолчанию, меняете нужные, вызываете `ldpc_simulate`:

```matlab
p = ldpc_default_params();
p.algorithm = 'norm-min-sum';
p.scaling = 0.8;
p.EbN0dB = 0:0.25:2;
r = ldpc_simulate(p);

ldpc_plot_results(r, {'NMS 0.8'}, 'NMS');
```

Что можно менять (`ldpc_default_params.m`):

| Параметр | Значения |
|---|---|
| `p.algorithm` | `'bp'`, `'layered-bp'`, `'norm-min-sum'`, `'offset-min-sum'`, `'bit-flip'`, `'wbf'` |
| `p.bgn`, `p.Zc` | базовый граф 1 или 2 и размер подъёма; K = 22·Zc для BG1, 10·Zc для BG2 |
| `p.rate` | скорость кода |
| `p.modulation` | `'bpsk'`, `'qpsk'`, `'16qam'`, `'64qam'` |
| `p.maxIter` | максимум итераций |
| `p.scaling` | α для NMS |
| `p.offset` | β для OMS |
| `p.termination` | `'early'` — останов по нулевому синдрому, `'max'` — всегда до maxIter |
| `p.EbN0dB` | точки по Eb/N0 |
| `p.minFrameErrors`, `p.maxFrames` | когда заканчивать точку |

Несколько кривых на одном графике: соберите результаты в массив и передайте подписи.

```matlab
res = [r1, r2, r3];
ldpc_plot_results(res, {'5 итер.', '10 итер.', '50 итер.'}, 'NMS');
```

## Скрипты

- `run_e0_uncoded.m` — проверка канала: BER без кодирования против теории
- `run_e1_algorithms.m` — сравнение BP, layered BP, NMS, OMS
- `run_e1b_bitflip.m` — то же плюс bit-flip и MWBF (сначала нужен `run_e1_algorithms`)
- `run_e1c_bf_vs_paper.m` — сверка BF и BP с Kou, Lin, Fossorier (2001) на EG-коде (255,175)
- `run_e2_iterations.m` — влияние числа итераций
- `run_e3_minsum_params.m` — подбор α и β
- `plot_bg2_degrees.m` — базовый граф BG2 со степенями столбцов

## Источники кода

- `ldpc_bitflip_decode.m` — порт функции `PH_bitflipping` из пакета bf-tests-awgn-egldpc 0.2 (F. Pujaico Rivera, GPL-3), https://github.com/Lcrypto/LDPC-Iterative-Bit-Flipping-family-decoders
- `ldpc_wbf_decode.m` — порт WBF/MWBF из AFF3CT (MIT), https://github.com/aff3ct/aff3ct
- `nr_ldpc_pcm.m` берёт матрицу H из внутренней функции 5G Toolbox (`nr5g.internal.ldpc.getParams`), в других версиях MATLAB она может отличаться.
