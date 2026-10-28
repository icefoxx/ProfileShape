# Consistency against tolerance margins — evaluation results

Protocol sha256 `54442bd14167ead4`; commit `960677c694`, uncommitted src files 4 (diff `5204d9da9eeb04f5`); consistency.py `5d56ff9477fbec76`; script `786b64ca73ca733f`; seeds simulation 9990000+, semi-real 9995000+; 2026-10-27T23:34:15 → 2026-10-28T00:22:44.

## A · Simulation (100 runs per cell; margins: shift 10 bp, stretch 5 %, offset 5 %)

| n | depth | alteration | estimate (mean) | 90 % CI covers truth | within | inconclusive | beyond | offset % truth / est / upper | offset within |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 12 | shift 0 | +0.12 | 0.92 | 1.00 | 0.00 | 0.00 | nan / 0.16 / 0.84 | 1.00 |
| 2 | 12 | shift 5 | +5.12 | 0.89 | 1.00 | 0.00 | 0.00 | 2.22 / 2.22 / 2.73 | 1.00 |
| 2 | 12 | shift 10 | +9.99 | 0.89 | 0.04 | 0.89 | 0.07 | 4.43 / 4.40 / 4.88 | 0.62 |
| 2 | 12 | shift 15 | +15.05 | 0.95 | 0.00 | 0.00 | 1.00 | 6.65 / 6.62 / 7.10 | 0.00 |
| 2 | 12 | shift 20 | +20.09 | 0.96 | 0.00 | 0.00 | 1.00 | 8.87 / 8.85 / 9.33 | 0.00 |
| 2 | 12 | stretch 2.5 | +2.62 | 0.87 | 1.00 | 0.00 | 0.00 | 1.60 / 1.56 / 2.10 | 1.00 |
| 2 | 12 | stretch 5 | +5.09 | 0.89 | 0.01 | 0.89 | 0.10 | 3.08 / 3.11 / 3.61 | 1.00 |
| 2 | 12 | stretch 7.5 | +7.34 | 0.91 | 0.00 | 0.00 | 1.00 | 4.52 / 4.50 / 4.99 | 0.53 |
| 2 | 12 | stretch 10 | +9.64 | 0.89 | 0.00 | 0.00 | 1.00 | 5.95 / 5.94 / 6.43 | 0.00 |
| 2 | 40 | shift 0 | +0.09 | 0.92 | 1.00 | 0.00 | 0.00 | nan / 0.18 / 0.78 | 1.00 |
| 2 | 40 | shift 5 | +4.97 | 0.87 | 1.00 | 0.00 | 0.00 | 2.22 / 2.18 / 2.62 | 1.00 |
| 2 | 40 | shift 10 | +10.07 | 0.86 | 0.08 | 0.86 | 0.06 | 4.44 / 4.43 / 4.86 | 0.65 |
| 2 | 40 | shift 15 | +15.06 | 0.87 | 0.00 | 0.00 | 1.00 | 6.65 / 6.63 / 7.06 | 0.00 |
| 2 | 40 | shift 20 | +20.07 | 0.94 | 0.00 | 0.00 | 1.00 | 8.88 / 8.84 / 9.28 | 0.00 |
| 2 | 40 | stretch 2.5 | +2.61 | 0.92 | 1.00 | 0.00 | 0.00 | 1.60 / 1.58 / 2.04 | 1.00 |
| 2 | 40 | stretch 5 | +4.95 | 0.89 | 0.07 | 0.89 | 0.04 | 3.08 / 3.03 / 3.47 | 1.00 |
| 2 | 40 | stretch 7.5 | +7.35 | 0.85 | 0.00 | 0.00 | 1.00 | 4.52 / 4.51 / 4.95 | 0.58 |
| 2 | 40 | stretch 10 | +9.64 | 0.81 | 0.00 | 0.00 | 1.00 | 5.95 / 5.92 / 6.36 | 0.00 |
| 2 | 160 | shift 0 | -0.02 | 0.91 | 1.00 | 0.00 | 0.00 | nan / 0.14 / 0.73 | 1.00 |
| 2 | 160 | shift 5 | +5.09 | 0.87 | 1.00 | 0.00 | 0.00 | 2.22 / 2.23 / 2.66 | 1.00 |
| 2 | 160 | shift 10 | +10.08 | 0.93 | 0.03 | 0.93 | 0.04 | 4.43 / 4.43 / 4.85 | 0.71 |
| 2 | 160 | shift 15 | +15.09 | 0.93 | 0.00 | 0.00 | 1.00 | 6.65 / 6.64 / 7.05 | 0.00 |
| 2 | 160 | shift 20 | +20.07 | 0.94 | 0.00 | 0.00 | 1.00 | 8.88 / 8.84 / 9.26 | 0.00 |
| 2 | 160 | stretch 2.5 | +2.61 | 0.91 | 1.00 | 0.00 | 0.00 | 1.60 / 1.57 / 2.02 | 1.00 |
| 2 | 160 | stretch 5 | +5.07 | 0.86 | 0.07 | 0.86 | 0.07 | 3.08 / 3.11 / 3.54 | 1.00 |
| 2 | 160 | stretch 7.5 | +7.47 | 0.94 | 0.00 | 0.00 | 1.00 | 4.52 / 4.58 / 5.00 | 0.48 |
| 2 | 160 | stretch 10 | +9.73 | 0.92 | 0.00 | 0.00 | 1.00 | 5.94 / 5.98 / 6.40 | 0.00 |
| 3 | 12 | shift 0 | +0.09 | 0.88 | 1.00 | 0.00 | 0.00 | 0.00 / 0.26 / 0.84 | 1.00 |
| 3 | 12 | shift 5 | +5.01 | 0.90 | 1.00 | 0.00 | 0.00 | 2.22 / 2.21 / 2.64 | 1.00 |
| 3 | 12 | shift 10 | +10.10 | 0.92 | 0.01 | 0.92 | 0.07 | 4.43 / 4.44 / 4.87 | 0.71 |
| 3 | 12 | shift 15 | +15.30 | 0.95 | 0.00 | 0.00 | 1.00 | 6.65 / 6.74 / 7.16 | 0.00 |
| 3 | 12 | shift 20 | +20.12 | 0.92 | 0.00 | 0.00 | 1.00 | 8.87 / 8.86 / 9.28 | 0.00 |
| 3 | 12 | stretch 2.5 | +2.58 | 0.90 | 1.00 | 0.00 | 0.00 | 1.59 / 1.55 / 2.01 | 1.00 |
| 3 | 12 | stretch 5 | +5.02 | 0.88 | 0.06 | 0.88 | 0.06 | 3.07 / 3.07 / 3.51 | 1.00 |
| 3 | 12 | stretch 7.5 | +7.40 | 0.93 | 0.00 | 0.00 | 1.00 | 4.50 / 4.55 / 4.98 | 0.58 |
| 3 | 12 | stretch 10 | +9.68 | 0.84 | 0.00 | 0.00 | 1.00 | 5.92 / 5.95 / 6.38 | 0.00 |
| 3 | 40 | shift 0 | +0.12 | 0.91 | 1.00 | 0.00 | 0.00 | 0.00 / 0.20 / 0.72 | 1.00 |
| 3 | 40 | shift 5 | +4.99 | 0.93 | 1.00 | 0.00 | 0.00 | 2.22 / 2.19 / 2.57 | 1.00 |
| 3 | 40 | shift 10 | +10.17 | 0.93 | 0.02 | 0.93 | 0.05 | 4.44 / 4.49 / 4.86 | 0.68 |
| 3 | 40 | shift 15 | +15.16 | 0.91 | 0.00 | 0.00 | 1.00 | 6.66 / 6.68 / 7.06 | 0.00 |
| 3 | 40 | shift 20 | +20.11 | 0.91 | 0.00 | 0.00 | 1.00 | 8.87 / 8.86 / 9.24 | 0.00 |
| 3 | 40 | stretch 2.5 | +2.68 | 0.88 | 1.00 | 0.00 | 0.00 | 1.60 / 1.63 / 2.02 | 1.00 |
| 3 | 40 | stretch 5 | +5.07 | 0.90 | 0.04 | 0.90 | 0.06 | 3.07 / 3.11 / 3.49 | 1.00 |
| 3 | 40 | stretch 7.5 | +7.47 | 0.89 | 0.00 | 0.00 | 1.00 | 4.50 / 4.58 / 4.96 | 0.57 |
| 3 | 40 | stretch 10 | +9.69 | 0.77 | 0.00 | 0.00 | 1.00 | 5.92 / 5.94 / 6.33 | 0.00 |
| 3 | 160 | shift 0 | -0.03 | 0.94 | 1.00 | 0.00 | 0.00 | 0.00 / 0.11 / 0.63 | 1.00 |
| 3 | 160 | shift 5 | +5.05 | 0.88 | 1.00 | 0.00 | 0.00 | 2.22 / 2.21 / 2.58 | 1.00 |
| 3 | 160 | shift 10 | +10.17 | 0.85 | 0.03 | 0.85 | 0.12 | 4.44 / 4.49 / 4.85 | 0.72 |
| 3 | 160 | shift 15 | +15.18 | 0.88 | 0.00 | 0.00 | 1.00 | 6.65 / 6.68 / 7.04 | 0.00 |
| 3 | 160 | shift 20 | +20.07 | 0.88 | 0.00 | 0.00 | 1.00 | 8.88 / 8.85 / 9.22 | 0.00 |
| 3 | 160 | stretch 2.5 | +2.76 | 0.81 | 1.00 | 0.00 | 0.00 | 1.59 / 1.68 / 2.06 | 1.00 |
| 3 | 160 | stretch 5 | +5.05 | 0.84 | 0.08 | 0.84 | 0.08 | 3.07 / 3.09 / 3.46 | 1.00 |
| 3 | 160 | stretch 7.5 | +7.41 | 0.92 | 0.00 | 0.01 | 0.99 | 4.50 / 4.54 / 4.91 | 0.69 |
| 3 | 160 | stretch 10 | +9.73 | 0.79 | 0.00 | 0.00 | 1.00 | 5.92 / 5.96 / 6.33 | 0.00 |
| 4 | 12 | shift 0 | +0.06 | 0.94 | 1.00 | 0.00 | 0.00 | 0.00 / 0.14 / 0.71 | 1.00 |
| 4 | 12 | shift 5 | +5.00 | 0.91 | 1.00 | 0.00 | 0.00 | 2.22 / 2.19 / 2.60 | 1.00 |
| 4 | 12 | shift 10 | +10.26 | 0.89 | 0.01 | 0.89 | 0.10 | 4.43 / 4.52 / 4.91 | 0.63 |
| 4 | 12 | shift 15 | +15.30 | 0.90 | 0.00 | 0.00 | 1.00 | 6.65 / 6.75 / 7.14 | 0.00 |
| 4 | 12 | shift 20 | +19.98 | 0.90 | 0.00 | 0.00 | 1.00 | 8.88 / 8.82 / 9.22 | 0.00 |
| 4 | 12 | stretch 2.5 | +2.62 | 0.91 | 1.00 | 0.00 | 0.00 | 1.59 / 1.58 / 2.02 | 1.00 |
| 4 | 12 | stretch 5 | +5.02 | 0.89 | 0.02 | 0.89 | 0.09 | 3.07 / 3.07 / 3.50 | 1.00 |
| 4 | 12 | stretch 7.5 | +7.41 | 0.90 | 0.00 | 0.01 | 0.99 | 4.49 / 4.54 / 4.96 | 0.56 |
| 4 | 12 | stretch 10 | +9.71 | 0.86 | 0.00 | 0.00 | 1.00 | 5.91 / 5.96 / 6.39 | 0.00 |
| 4 | 40 | shift 0 | -0.06 | 0.82 | 1.00 | 0.00 | 0.00 | 0.00 / 0.18 / 0.68 | 1.00 |
| 4 | 40 | shift 5 | +5.02 | 0.90 | 1.00 | 0.00 | 0.00 | 2.21 / 2.20 / 2.56 | 1.00 |
| 4 | 40 | shift 10 | +10.05 | 0.90 | 0.04 | 0.90 | 0.06 | 4.44 / 4.44 / 4.79 | 0.76 |
| 4 | 40 | shift 15 | +15.15 | 0.92 | 0.00 | 0.00 | 1.00 | 6.65 / 6.68 / 7.04 | 0.00 |
| 4 | 40 | shift 20 | +20.19 | 0.89 | 0.00 | 0.00 | 1.00 | 8.87 / 8.90 / 9.26 | 0.00 |
| 4 | 40 | stretch 2.5 | +2.63 | 0.86 | 1.00 | 0.00 | 0.00 | 1.59 / 1.59 / 1.97 | 1.00 |
| 4 | 40 | stretch 5 | +5.13 | 0.91 | 0.02 | 0.91 | 0.07 | 3.07 / 3.14 / 3.50 | 1.00 |
| 4 | 40 | stretch 7.5 | +7.41 | 0.89 | 0.00 | 0.00 | 1.00 | 4.49 / 4.54 / 4.91 | 0.63 |
| 4 | 40 | stretch 10 | +9.76 | 0.85 | 0.00 | 0.00 | 1.00 | 5.90 / 5.99 / 6.35 | 0.00 |
| 4 | 160 | shift 0 | +0.05 | 0.85 | 1.00 | 0.00 | 0.00 | 0.00 / 0.17 / 0.63 | 1.00 |
| 4 | 160 | shift 5 | +5.06 | 0.92 | 1.00 | 0.00 | 0.00 | 2.21 / 2.22 / 2.57 | 1.00 |
| 4 | 160 | shift 10 | +10.11 | 0.90 | 0.05 | 0.90 | 0.05 | 4.44 / 4.47 / 4.81 | 0.79 |
| 4 | 160 | shift 15 | +15.10 | 0.87 | 0.00 | 0.00 | 1.00 | 6.66 / 6.66 / 7.01 | 0.00 |
| 4 | 160 | shift 20 | +20.14 | 0.88 | 0.00 | 0.00 | 1.00 | 8.87 / 8.88 / 9.23 | 0.00 |
| 4 | 160 | stretch 2.5 | +2.59 | 0.91 | 1.00 | 0.00 | 0.00 | 1.59 / 1.57 / 1.93 | 1.00 |
| 4 | 160 | stretch 5 | +5.08 | 0.91 | 0.01 | 0.91 | 0.08 | 3.07 / 3.11 / 3.45 | 1.00 |
| 4 | 160 | stretch 7.5 | +7.42 | 0.84 | 0.00 | 0.00 | 1.00 | 4.49 / 4.55 / 4.89 | 0.70 |
| 4 | 160 | stretch 10 | +9.71 | 0.82 | 0.00 | 0.00 | 1.00 | 5.90 / 5.95 / 6.30 | 0.00 |

Pooled coverage: shift 0.903, stretch 0.877; offset upper bound ≥ truth 0.901.

## B · Semi-real (20 shuffles; margins: shift M = max(10 bp, bin), stretch 5 %)

| dataset | group | n | alteration | estimate (mean) | 90 % CI covers truth | within | inconclusive | beyond |
|---|---|---|---|---|---|---|---|---|
| atac | snf2_0h | 3 | shift 0 bp | -0.11 | 0.80 | 20 | 0 | 0 |
| atac | snf2_0h | 3 | shift 5 bp | +6.08 | 0.55 | 20 | 0 | 0 |
| atac | snf2_0h | 3 | shift 10 bp | +10.12 | 0.90 | 0 | 18 | 2 |
| atac | snf2_0h | 3 | shift 15 bp | +17.03 | 0.15 | 0 | 0 | 20 |
| atac | snf2_0h | 3 | shift 20 bp | +20.44 | 0.85 | 0 | 0 | 20 |
| atac | snf2_0h | 3 | stretch 2.5 % | +3.13 | 0.95 | 8 | 12 | 0 |
| atac | snf2_0h | 3 | stretch 5 % | +5.61 | 1.00 | 0 | 20 | 0 |
| atac | snf2_0h | 3 | stretch 7.5 % | +8.09 | 1.00 | 0 | 10 | 10 |
| atac | snf2_0h | 3 | stretch 10 % | +8.49 | 1.00 | 0 | 5 | 15 |
| atac | snf2_1h | 3 | shift 0 bp | -0.09 | 0.85 | 20 | 0 | 0 |
| atac | snf2_1h | 3 | shift 5 bp | +5.49 | 0.80 | 20 | 0 | 0 |
| atac | snf2_1h | 3 | shift 10 bp | +10.15 | 0.75 | 2 | 15 | 3 |
| atac | snf2_1h | 3 | shift 15 bp | +16.64 | 0.05 | 0 | 0 | 20 |
| atac | snf2_1h | 3 | shift 20 bp | +21.09 | 0.20 | 0 | 0 | 20 |
| atac | snf2_1h | 3 | stretch 2.5 % | +3.27 | 0.85 | 10 | 10 | 0 |
| atac | snf2_1h | 3 | stretch 5 % | +6.24 | 0.75 | 0 | 15 | 5 |
| atac | snf2_1h | 3 | stretch 7.5 % | +9.20 | 0.50 | 0 | 0 | 20 |
| atac | snf2_1h | 3 | stretch 10 % | +12.49 | 0.25 | 0 | 0 | 20 |
| atac | sth1_0h | 4 | shift 0 bp | -0.09 | 0.85 | 20 | 0 | 0 |
| atac | sth1_0h | 4 | shift 5 bp | +5.69 | 0.60 | 20 | 0 | 0 |
| atac | sth1_0h | 4 | shift 10 bp | +10.34 | 0.85 | 0 | 17 | 3 |
| atac | sth1_0h | 4 | shift 15 bp | +17.37 | 0.00 | 0 | 0 | 20 |
| atac | sth1_0h | 4 | shift 20 bp | +21.56 | 0.10 | 0 | 0 | 20 |
| atac | sth1_0h | 4 | stretch 2.5 % | +2.91 | 1.00 | 11 | 9 | 0 |
| atac | sth1_0h | 4 | stretch 5 % | +6.55 | 0.90 | 0 | 18 | 2 |
| atac | sth1_0h | 4 | stretch 7.5 % | +10.60 | 0.40 | 0 | 0 | 20 |
| atac | sth1_0h | 4 | stretch 10 % | +13.83 | 0.15 | 0 | 0 | 20 |
| atac | sth1_1h | 4 | shift 0 bp | -0.09 | 0.80 | 20 | 0 | 0 |
| atac | sth1_1h | 4 | shift 5 bp | +5.65 | 0.65 | 20 | 0 | 0 |
| atac | sth1_1h | 4 | shift 10 bp | +10.24 | 0.85 | 0 | 17 | 3 |
| atac | sth1_1h | 4 | shift 15 bp | +17.14 | 0.10 | 0 | 0 | 20 |
| atac | sth1_1h | 4 | shift 20 bp | +21.34 | 0.30 | 0 | 0 | 20 |
| atac | sth1_1h | 4 | stretch 2.5 % | +3.26 | 0.85 | 9 | 11 | 0 |
| atac | sth1_1h | 4 | stretch 5 % | +6.57 | 0.75 | 0 | 15 | 5 |
| atac | sth1_1h | 4 | stretch 7.5 % | +9.27 | 0.65 | 0 | 1 | 19 |
| atac | sth1_1h | 4 | stretch 10 % | +13.02 | 0.45 | 0 | 0 | 20 |
| k27ac_chip | CLP | 4 | shift 0 bp | +0.00 | 0.95 | 20 | 0 | 0 |
| k27ac_chip | CLP | 4 | shift 25 bp | +26.03 | 0.25 | 20 | 0 | 0 |
| k27ac_chip | CLP | 4 | shift 50 bp | +53.43 | 0.00 | 0 | 0 | 20 |
| k27ac_chip | CLP | 4 | shift 75 bp | +79.04 | 0.00 | 0 | 0 | 20 |
| k27ac_chip | CLP | 4 | shift 100 bp | +108.05 | 0.00 | 0 | 0 | 20 |
| k27ac_chip | CLP | 4 | stretch 2.5 % | +1.84 | 0.50 | 20 | 0 | 0 |
| k27ac_chip | CLP | 4 | stretch 5 % | +3.69 | 0.05 | 19 | 1 | 0 |
| k27ac_chip | CLP | 4 | stretch 7.5 % | +5.51 | 0.00 | 0 | 16 | 4 |
| k27ac_chip | CLP | 4 | stretch 10 % | +7.26 | 0.00 | 0 | 4 | 16 |
| k27ac_chip | Sham | 4 | shift 0 bp | +0.07 | 1.00 | 20 | 0 | 0 |
| k27ac_chip | Sham | 4 | shift 25 bp | +26.14 | 0.15 | 20 | 0 | 0 |
| k27ac_chip | Sham | 4 | shift 50 bp | +53.61 | 0.00 | 0 | 0 | 20 |
| k27ac_chip | Sham | 4 | shift 75 bp | +79.38 | 0.00 | 0 | 0 | 20 |
| k27ac_chip | Sham | 4 | shift 100 bp | +108.72 | 0.00 | 0 | 0 | 20 |
| k27ac_chip | Sham | 4 | stretch 2.5 % | +1.86 | 0.70 | 20 | 0 | 0 |
| k27ac_chip | Sham | 4 | stretch 5 % | +3.76 | 0.05 | 19 | 1 | 0 |
| k27ac_chip | Sham | 4 | stretch 7.5 % | +5.81 | 0.00 | 0 | 9 | 11 |
| k27ac_chip | Sham | 4 | stretch 10 % | +7.84 | 0.00 | 0 | 0 | 20 |
| k27ac_cutrun | Arm | 4 | shift 0 bp | +0.07 | 0.95 | 20 | 0 | 0 |
| k27ac_cutrun | Arm | 4 | shift 25 bp | +25.32 | 0.95 | 20 | 0 | 0 |
| k27ac_cutrun | Arm | 4 | shift 50 bp | +50.96 | 0.50 | 0 | 10 | 10 |
| k27ac_cutrun | Arm | 4 | shift 75 bp | +76.15 | 0.45 | 0 | 0 | 20 |
| k27ac_cutrun | Arm | 4 | shift 100 bp | +102.35 | 0.00 | 0 | 0 | 20 |
| k27ac_cutrun | Arm | 4 | stretch 2.5 % | +2.22 | 0.70 | 20 | 0 | 0 |
| k27ac_cutrun | Arm | 4 | stretch 5 % | +4.30 | 0.15 | 17 | 3 | 0 |
| k27ac_cutrun | Arm | 4 | stretch 7.5 % | +6.36 | 0.00 | 0 | 0 | 20 |
| k27ac_cutrun | Arm | 4 | stretch 10 % | +8.40 | 0.00 | 0 | 0 | 20 |
| k27ac_cutrun | Naive | 4 | shift 0 bp | +0.01 | 0.90 | 20 | 0 | 0 |
| k27ac_cutrun | Naive | 4 | shift 25 bp | +25.30 | 0.95 | 20 | 0 | 0 |
| k27ac_cutrun | Naive | 4 | shift 50 bp | +50.96 | 0.75 | 0 | 15 | 5 |
| k27ac_cutrun | Naive | 4 | shift 75 bp | +76.14 | 0.65 | 0 | 0 | 20 |
| k27ac_cutrun | Naive | 4 | shift 100 bp | +102.33 | 0.05 | 0 | 0 | 20 |
| k27ac_cutrun | Naive | 4 | stretch 2.5 % | +2.12 | 0.60 | 20 | 0 | 0 |
| k27ac_cutrun | Naive | 4 | stretch 5 % | +4.24 | 0.25 | 15 | 5 | 0 |
| k27ac_cutrun | Naive | 4 | stretch 7.5 % | +6.30 | 0.00 | 0 | 0 | 20 |
| k27ac_cutrun | Naive | 4 | stretch 10 % | +8.31 | 0.00 | 0 | 0 | 20 |
| k4me3_chip_summit | CLP | 4 | shift 0 bp | +0.01 | 1.00 | 20 | 0 | 0 |
| k4me3_chip_summit | CLP | 4 | shift 12.5 bp | +12.58 | 1.00 | 20 | 0 | 0 |
| k4me3_chip_summit | CLP | 4 | shift 25 bp | +25.26 | 0.95 | 0 | 19 | 1 |
| k4me3_chip_summit | CLP | 4 | shift 37.5 bp | +37.87 | 0.95 | 0 | 0 | 20 |
| k4me3_chip_summit | CLP | 4 | shift 50 bp | +50.74 | 0.70 | 0 | 0 | 20 |
| k4me3_chip_summit | CLP | 4 | stretch 2.5 % | +1.97 | 0.45 | 20 | 0 | 0 |
| k4me3_chip_summit | CLP | 4 | stretch 5 % | +3.91 | 0.00 | 20 | 0 | 0 |
| k4me3_chip_summit | CLP | 4 | stretch 7.5 % | +5.81 | 0.00 | 0 | 3 | 17 |
| k4me3_chip_summit | CLP | 4 | stretch 10 % | +7.61 | 0.00 | 0 | 0 | 20 |
| k4me3_chip_summit | Sham | 4 | shift 0 bp | -0.03 | 0.85 | 20 | 0 | 0 |
| k4me3_chip_summit | Sham | 4 | shift 12.5 bp | +12.57 | 0.80 | 20 | 0 | 0 |
| k4me3_chip_summit | Sham | 4 | shift 25 bp | +25.26 | 0.75 | 1 | 15 | 4 |
| k4me3_chip_summit | Sham | 4 | shift 37.5 bp | +37.91 | 0.80 | 0 | 0 | 20 |
| k4me3_chip_summit | Sham | 4 | shift 50 bp | +50.80 | 0.65 | 0 | 0 | 20 |
| k4me3_chip_summit | Sham | 4 | stretch 2.5 % | +1.83 | 0.30 | 20 | 0 | 0 |
| k4me3_chip_summit | Sham | 4 | stretch 5 % | +3.81 | 0.00 | 20 | 0 | 0 |
| k4me3_chip_summit | Sham | 4 | stretch 7.5 % | +5.78 | 0.00 | 0 | 2 | 18 |
| k4me3_chip_summit | Sham | 4 | stretch 10 % | +7.68 | 0.00 | 0 | 0 | 20 |
| k4me3_chip_tss | CLP | 4 | shift 0 bp | +0.09 | 0.95 | 20 | 0 | 0 |
| k4me3_chip_tss | CLP | 4 | shift 12.5 bp | +12.71 | 0.95 | 20 | 0 | 0 |
| k4me3_chip_tss | CLP | 4 | shift 25 bp | +25.49 | 0.85 | 0 | 17 | 3 |
| k4me3_chip_tss | CLP | 4 | shift 37.5 bp | +38.15 | 0.80 | 0 | 0 | 20 |
| k4me3_chip_tss | CLP | 4 | shift 50 bp | +51.20 | 0.35 | 0 | 0 | 20 |
| k4me3_chip_tss | CLP | 4 | stretch 2.5 % | +2.17 | 0.45 | 20 | 0 | 0 |
| k4me3_chip_tss | CLP | 4 | stretch 5 % | +4.23 | 0.00 | 20 | 0 | 0 |
| k4me3_chip_tss | CLP | 4 | stretch 7.5 % | +6.23 | 0.00 | 0 | 0 | 20 |
| k4me3_chip_tss | CLP | 4 | stretch 10 % | +8.17 | 0.00 | 0 | 0 | 20 |
| k4me3_chip_tss | Sham | 4 | shift 0 bp | +0.10 | 0.90 | 20 | 0 | 0 |
| k4me3_chip_tss | Sham | 4 | shift 12.5 bp | +12.72 | 0.90 | 20 | 0 | 0 |
| k4me3_chip_tss | Sham | 4 | shift 25 bp | +25.49 | 0.75 | 0 | 15 | 5 |
| k4me3_chip_tss | Sham | 4 | shift 37.5 bp | +38.15 | 0.70 | 0 | 0 | 20 |
| k4me3_chip_tss | Sham | 4 | shift 50 bp | +51.21 | 0.45 | 0 | 0 | 20 |
| k4me3_chip_tss | Sham | 4 | stretch 2.5 % | +2.17 | 0.45 | 20 | 0 | 0 |
| k4me3_chip_tss | Sham | 4 | stretch 5 % | +4.27 | 0.00 | 20 | 0 | 0 |
| k4me3_chip_tss | Sham | 4 | stretch 7.5 % | +6.32 | 0.00 | 0 | 0 | 20 |
| k4me3_chip_tss | Sham | 4 | stretch 10 % | +8.37 | 0.00 | 0 | 0 | 20 |
| k4me3_cutrun_summit | Arm | 4 | shift 0 bp | -0.00 | 0.90 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Arm | 4 | shift 12.5 bp | +12.55 | 0.95 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Arm | 4 | shift 25 bp | +25.18 | 0.95 | 0 | 19 | 1 |
| k4me3_cutrun_summit | Arm | 4 | shift 37.5 bp | +37.75 | 0.95 | 0 | 0 | 20 |
| k4me3_cutrun_summit | Arm | 4 | shift 50 bp | +50.54 | 0.95 | 0 | 0 | 20 |
| k4me3_cutrun_summit | Arm | 4 | stretch 2.5 % | +0.45 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Arm | 4 | stretch 5 % | +1.14 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Arm | 4 | stretch 7.5 % | +2.05 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Arm | 4 | stretch 10 % | +3.04 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Naive | 4 | shift 0 bp | +0.00 | 0.85 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Naive | 4 | shift 12.5 bp | +12.54 | 0.85 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Naive | 4 | shift 25 bp | +25.16 | 0.85 | 1 | 17 | 2 |
| k4me3_cutrun_summit | Naive | 4 | shift 37.5 bp | +37.72 | 0.85 | 0 | 0 | 20 |
| k4me3_cutrun_summit | Naive | 4 | shift 50 bp | +50.49 | 0.80 | 0 | 0 | 20 |
| k4me3_cutrun_summit | Naive | 4 | stretch 2.5 % | +0.56 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Naive | 4 | stretch 5 % | +0.84 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Naive | 4 | stretch 7.5 % | +1.50 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_summit | Naive | 4 | stretch 10 % | +2.34 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Arm | 4 | shift 0 bp | +0.02 | 0.80 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Arm | 4 | shift 12.5 bp | +12.61 | 0.95 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Arm | 4 | shift 25 bp | +25.30 | 0.95 | 0 | 19 | 1 |
| k4me3_cutrun_tss | Arm | 4 | shift 37.5 bp | +37.95 | 0.95 | 0 | 0 | 20 |
| k4me3_cutrun_tss | Arm | 4 | shift 50 bp | +50.83 | 0.95 | 0 | 0 | 20 |
| k4me3_cutrun_tss | Arm | 4 | stretch 2.5 % | +1.48 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Arm | 4 | stretch 5 % | +2.88 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Arm | 4 | stretch 7.5 % | +4.36 | 0.00 | 10 | 10 | 0 |
| k4me3_cutrun_tss | Arm | 4 | stretch 10 % | +5.81 | 0.00 | 0 | 4 | 16 |
| k4me3_cutrun_tss | Naive | 4 | shift 0 bp | +0.09 | 0.80 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Naive | 4 | shift 12.5 bp | +12.69 | 0.80 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Naive | 4 | shift 25 bp | +25.38 | 0.85 | 0 | 17 | 3 |
| k4me3_cutrun_tss | Naive | 4 | shift 37.5 bp | +38.03 | 0.85 | 0 | 0 | 20 |
| k4me3_cutrun_tss | Naive | 4 | shift 50 bp | +50.89 | 0.85 | 0 | 0 | 20 |
| k4me3_cutrun_tss | Naive | 4 | stretch 2.5 % | +1.46 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Naive | 4 | stretch 5 % | +2.94 | 0.00 | 20 | 0 | 0 |
| k4me3_cutrun_tss | Naive | 4 | stretch 7.5 % | +4.39 | 0.00 | 18 | 2 | 0 |
| k4me3_cutrun_tss | Naive | 4 | stretch 10 % | +5.79 | 0.00 | 0 | 0 | 20 |
| mnase | WT | 4 | shift 0 bp | +0.02 | 0.90 | 20 | 0 | 0 |
| mnase | WT | 4 | shift 5 bp | +5.62 | 0.45 | 20 | 0 | 0 |
| mnase | WT | 4 | shift 10 bp | +11.57 | 0.00 | 0 | 0 | 20 |
| mnase | WT | 4 | shift 15 bp | +19.09 | 0.00 | 0 | 0 | 20 |
| mnase | WT | 4 | shift 20 bp | +25.58 | 0.00 | 0 | 0 | 20 |
| mnase | WT | 4 | stretch 2.5 % | +2.96 | 0.50 | 20 | 0 | 0 |
| mnase | WT | 4 | stretch 5 % | +6.89 | 0.00 | 0 | 0 | 20 |
| mnase | WT | 4 | stretch 7.5 % | +11.78 | 0.00 | 0 | 0 | 20 |
| mnase | WT | 4 | stretch 10 % | +14.27 | 0.00 | 0 | 0 | 20 |
| mnase | dKO | 4 | shift 0 bp | +0.04 | 0.95 | 20 | 0 | 0 |
| mnase | dKO | 4 | shift 5 bp | +5.81 | 0.80 | 20 | 0 | 0 |
| mnase | dKO | 4 | shift 10 bp | +11.80 | 0.15 | 0 | 3 | 17 |
| mnase | dKO | 4 | shift 15 bp | +19.27 | 0.00 | 0 | 0 | 20 |
| mnase | dKO | 4 | shift 20 bp | +25.13 | 0.00 | 0 | 0 | 20 |
| mnase | dKO | 4 | stretch 2.5 % | +0.51 | 0.75 | 16 | 4 | 0 |
| mnase | dKO | 4 | stretch 5 % | +1.35 | 0.85 | 0 | 20 | 0 |
| mnase | dKO | 4 | stretch 7.5 % | +2.43 | 0.65 | 0 | 20 | 0 |
| mnase | dKO | 4 | stretch 10 % | +2.60 | 0.15 | 0 | 20 | 0 |
| mnase | isw1rtf1 | 4 | shift 0 bp | -0.06 | 0.90 | 20 | 0 | 0 |
| mnase | isw1rtf1 | 4 | shift 5 bp | +5.71 | 0.70 | 20 | 0 | 0 |
| mnase | isw1rtf1 | 4 | shift 10 bp | +11.87 | 0.10 | 0 | 2 | 18 |
| mnase | isw1rtf1 | 4 | shift 15 bp | +19.56 | 0.00 | 0 | 0 | 20 |
| mnase | isw1rtf1 | 4 | shift 20 bp | +25.69 | 0.00 | 0 | 0 | 20 |
| mnase | isw1rtf1 | 4 | stretch 2.5 % | +1.21 | 0.85 | 19 | 1 | 0 |
| mnase | isw1rtf1 | 4 | stretch 5 % | +1.92 | 0.40 | 12 | 8 | 0 |
| mnase | isw1rtf1 | 4 | stretch 7.5 % | +3.57 | 0.45 | 0 | 20 | 0 |
| mnase | isw1rtf1 | 4 | stretch 10 % | +5.22 | 0.35 | 0 | 20 | 0 |
| mnase | rtf1 | 4 | shift 0 bp | -0.02 | 0.95 | 20 | 0 | 0 |
| mnase | rtf1 | 4 | shift 5 bp | +5.64 | 0.70 | 20 | 0 | 0 |
| mnase | rtf1 | 4 | shift 10 bp | +11.56 | 0.05 | 0 | 1 | 19 |
| mnase | rtf1 | 4 | shift 15 bp | +19.14 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1 | 4 | shift 20 bp | +25.31 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1 | 4 | stretch 2.5 % | +2.82 | 0.85 | 20 | 0 | 0 |
| mnase | rtf1 | 4 | stretch 5 % | +5.70 | 0.85 | 0 | 17 | 3 |
| mnase | rtf1 | 4 | stretch 7.5 % | +9.80 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1 | 4 | stretch 10 % | +12.90 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1-11A | 4 | shift 0 bp | +0.01 | 0.85 | 20 | 0 | 0 |
| mnase | rtf1-11A | 4 | shift 5 bp | +5.68 | 0.55 | 20 | 0 | 0 |
| mnase | rtf1-11A | 4 | shift 10 bp | +11.71 | 0.05 | 0 | 1 | 19 |
| mnase | rtf1-11A | 4 | shift 15 bp | +19.32 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1-11A | 4 | shift 20 bp | +25.59 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1-11A | 4 | stretch 2.5 % | +2.57 | 0.95 | 20 | 0 | 0 |
| mnase | rtf1-11A | 4 | stretch 5 % | +5.52 | 0.95 | 0 | 19 | 1 |
| mnase | rtf1-11A | 4 | stretch 7.5 % | +9.26 | 0.15 | 0 | 0 | 20 |
| mnase | rtf1-11A | 4 | stretch 10 % | +12.18 | 0.35 | 0 | 0 | 20 |
| mnase | rtf1-8-9A | 4 | shift 0 bp | -0.08 | 0.85 | 20 | 0 | 0 |
| mnase | rtf1-8-9A | 4 | shift 5 bp | +5.62 | 0.70 | 20 | 0 | 0 |
| mnase | rtf1-8-9A | 4 | shift 10 bp | +11.66 | 0.10 | 0 | 2 | 18 |
| mnase | rtf1-8-9A | 4 | shift 15 bp | +19.24 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1-8-9A | 4 | shift 20 bp | +25.52 | 0.00 | 0 | 0 | 20 |
| mnase | rtf1-8-9A | 4 | stretch 2.5 % | +2.57 | 1.00 | 20 | 0 | 0 |
| mnase | rtf1-8-9A | 4 | stretch 5 % | +5.16 | 1.00 | 0 | 20 | 0 |
| mnase | rtf1-8-9A | 4 | stretch 7.5 % | +8.79 | 0.85 | 0 | 0 | 20 |
| mnase | rtf1-8-9A | 4 | stretch 10 % | +11.34 | 0.85 | 0 | 0 | 20 |
| proseq | NELFB | 2 | shift 0 bp | -0.03 | 0.95 | 20 | 0 | 0 |
| proseq | NELFB | 2 | shift 5 bp | +5.02 | 0.95 | 20 | 0 | 0 |
| proseq | NELFB | 2 | shift 10 bp | +10.10 | 0.90 | 0 | 18 | 2 |
| proseq | NELFB | 2 | shift 15 bp | +15.25 | 0.50 | 0 | 0 | 20 |
| proseq | NELFB | 2 | shift 20 bp | +20.50 | 0.00 | 0 | 0 | 20 |
| proseq | NELFB | 2 | stretch 2.5 % | +2.53 | 0.95 | 20 | 0 | 0 |
| proseq | NELFB | 2 | stretch 5 % | +4.70 | 0.85 | 3 | 17 | 0 |
| proseq | NELFB | 2 | stretch 7.5 % | +6.62 | 0.45 | 0 | 0 | 20 |
| proseq | NELFB | 2 | stretch 10 % | +8.55 | 0.00 | 0 | 0 | 20 |
| proseq | PAF1 | 2 | shift 0 bp | +0.04 | 0.85 | 20 | 0 | 0 |
| proseq | PAF1 | 2 | shift 5 bp | +5.10 | 0.85 | 20 | 0 | 0 |
| proseq | PAF1 | 2 | shift 10 bp | +10.18 | 0.80 | 0 | 16 | 4 |
| proseq | PAF1 | 2 | shift 15 bp | +15.34 | 0.60 | 0 | 0 | 20 |
| proseq | PAF1 | 2 | shift 20 bp | +20.59 | 0.20 | 0 | 0 | 20 |
| proseq | PAF1 | 2 | stretch 2.5 % | +2.15 | 0.95 | 20 | 0 | 0 |
| proseq | PAF1 | 2 | stretch 5 % | +4.13 | 0.55 | 9 | 11 | 0 |
| proseq | PAF1 | 2 | stretch 7.5 % | +5.80 | 0.05 | 0 | 12 | 8 |
| proseq | PAF1 | 2 | stretch 10 % | +7.63 | 0.00 | 0 | 0 | 20 |
| proseq | RTF1 | 2 | shift 0 bp | -0.07 | 0.95 | 20 | 0 | 0 |
| proseq | RTF1 | 2 | shift 5 bp | +4.98 | 1.00 | 20 | 0 | 0 |
| proseq | RTF1 | 2 | shift 10 bp | +10.06 | 1.00 | 0 | 20 | 0 |
| proseq | RTF1 | 2 | shift 15 bp | +15.21 | 0.75 | 0 | 0 | 20 |
| proseq | RTF1 | 2 | shift 20 bp | +20.45 | 0.20 | 0 | 0 | 20 |
| proseq | RTF1 | 2 | stretch 2.5 % | +2.20 | 1.00 | 20 | 0 | 0 |
| proseq | RTF1 | 2 | stretch 5 % | +4.14 | 0.55 | 9 | 11 | 0 |
| proseq | RTF1 | 2 | stretch 7.5 % | +5.87 | 0.00 | 0 | 9 | 11 |
| proseq | RTF1 | 2 | stretch 10 % | +7.61 | 0.00 | 0 | 0 | 20 |
| tf_chip | TAF15 | 2 | shift 0 bp | -0.02 | 0.85 | 20 | 0 | 0 |
| tf_chip | TAF15 | 2 | shift 10 bp | +10.02 | 0.85 | 20 | 0 | 0 |
| tf_chip | TAF15 | 2 | shift 20 bp | +20.10 | 0.85 | 1 | 17 | 2 |
| tf_chip | TAF15 | 2 | shift 30 bp | +30.23 | 0.85 | 0 | 0 | 20 |
| tf_chip | TAF15 | 2 | shift 40 bp | +40.47 | 0.90 | 0 | 0 | 20 |
| tf_chip | TAF15 | 2 | stretch 2.5 % | +2.34 | 0.95 | 20 | 0 | 0 |
| tf_chip | TAF15 | 2 | stretch 5 % | +4.70 | 0.85 | 3 | 17 | 0 |
| tf_chip | TAF15 | 2 | stretch 7.5 % | +6.94 | 0.60 | 0 | 0 | 20 |
| tf_chip | TAF15 | 2 | stretch 10 % | +9.23 | 0.25 | 0 | 0 | 20 |
| tf_chip | TCF3 | 2 | shift 0 bp | +0.01 | 1.00 | 20 | 0 | 0 |
| tf_chip | TCF3 | 2 | shift 10 bp | +10.05 | 0.95 | 20 | 0 | 0 |
| tf_chip | TCF3 | 2 | shift 20 bp | +20.13 | 0.90 | 0 | 18 | 2 |
| tf_chip | TCF3 | 2 | shift 30 bp | +30.30 | 0.85 | 0 | 0 | 20 |
| tf_chip | TCF3 | 2 | shift 40 bp | +40.56 | 0.85 | 0 | 0 | 20 |
| tf_chip | TCF3 | 2 | stretch 2.5 % | +2.62 | 0.75 | 20 | 0 | 0 |
| tf_chip | TCF3 | 2 | stretch 5 % | +4.99 | 0.80 | 2 | 16 | 2 |
| tf_chip | TCF3 | 2 | stretch 7.5 % | +7.23 | 0.75 | 0 | 0 | 20 |
| tf_chip | TCF3 | 2 | stretch 10 % | +9.46 | 0.65 | 0 | 0 | 20 |
| tf_chip | ZNF384WT | 3 | shift 0 bp | -0.13 | 0.95 | 20 | 0 | 0 |
| tf_chip | ZNF384WT | 3 | shift 10 bp | +9.90 | 0.95 | 20 | 0 | 0 |
| tf_chip | ZNF384WT | 3 | shift 20 bp | +19.98 | 1.00 | 0 | 20 | 0 |
| tf_chip | ZNF384WT | 3 | shift 30 bp | +30.09 | 0.95 | 0 | 0 | 20 |
| tf_chip | ZNF384WT | 3 | shift 40 bp | +40.33 | 0.90 | 0 | 0 | 20 |
| tf_chip | ZNF384WT | 3 | stretch 2.5 % | +2.17 | 0.75 | 20 | 0 | 0 |
| tf_chip | ZNF384WT | 3 | stretch 5 % | +4.56 | 0.65 | 7 | 13 | 0 |
| tf_chip | ZNF384WT | 3 | stretch 7.5 % | +6.87 | 0.40 | 0 | 0 | 20 |
| tf_chip | ZNF384WT | 3 | stretch 10 % | +9.16 | 0.25 | 0 | 0 | 20 |
| tf_cutrun | A2.2 | 6 | shift 0 bp | +0.32 | 0.90 | 20 | 0 | 0 |
| tf_cutrun | A2.2 | 6 | shift 5 bp | +5.36 | 0.90 | 12 | 8 | 0 |
| tf_cutrun | A2.2 | 6 | shift 10 bp | +10.27 | 0.90 | 0 | 18 | 2 |
| tf_cutrun | A2.2 | 6 | shift 15 bp | +15.09 | 0.95 | 0 | 7 | 13 |
| tf_cutrun | A2.2 | 6 | shift 20 bp | +19.91 | 1.00 | 0 | 0 | 20 |
| tf_cutrun | A2.2 | 6 | stretch 2.5 % | +2.21 | 0.85 | 5 | 15 | 0 |
| tf_cutrun | A2.2 | 6 | stretch 5 % | +3.64 | 0.95 | 1 | 19 | 0 |
| tf_cutrun | A2.2 | 6 | stretch 7.5 % | +5.09 | 0.80 | 0 | 18 | 2 |
| tf_cutrun | A2.2 | 6 | stretch 10 % | +6.53 | 0.60 | 0 | 16 | 4 |
| tf_cutrun | A673 | 6 | shift 0 bp | +0.28 | 0.85 | 20 | 0 | 0 |
| tf_cutrun | A673 | 6 | shift 5 bp | +5.24 | 0.85 | 17 | 3 | 0 |
| tf_cutrun | A673 | 6 | shift 10 bp | +10.21 | 0.95 | 0 | 19 | 1 |
| tf_cutrun | A673 | 6 | shift 15 bp | +15.24 | 0.95 | 0 | 1 | 19 |
| tf_cutrun | A673 | 6 | shift 20 bp | +20.24 | 0.90 | 0 | 0 | 20 |
| tf_cutrun | A673 | 6 | stretch 2.5 % | +1.67 | 0.85 | 13 | 7 | 0 |
| tf_cutrun | A673 | 6 | stretch 5 % | +2.50 | 0.60 | 8 | 12 | 0 |
| tf_cutrun | A673 | 6 | stretch 7.5 % | +4.03 | 0.45 | 2 | 18 | 0 |
| tf_cutrun | A673 | 6 | stretch 10 % | +4.92 | 0.25 | 0 | 20 | 0 |

Pooled semi-real coverage 0.531 (m = 5400).

## C · Real replicates as they are (reference ranges, no truth)

| dataset | group | n | offset % (upper bound) range | shift bp range | stretch % range |
|---|---|---|---|---|---|
| mnase | WT | 4 | 2.06 – 9.50 (3.03 – 10.43) | -1.24 – 1.58 | -0.52 – 0.34 |
| mnase | dKO | 4 | 5.09 – 37.20 (5.90 – 38.49) | -1.20 – 0.85 | -10.59 – 9.30 |
| mnase | rtf1-8-9A | 4 | 2.27 – 5.33 (3.24 – 6.26) | -1.20 – 0.94 | -0.73 – 0.89 |
| mnase | rtf1-11A | 4 | 2.76 – 11.50 (3.60 – 12.39) | -1.34 – 1.88 | -1.36 – 0.99 |
| mnase | rtf1 | 4 | 5.13 – 20.42 (6.06 – 21.69) | -3.75 – 2.66 | -1.19 – 1.66 |
| mnase | isw1rtf1 | 4 | 8.64 – 16.36 (9.50 – 17.46) | -1.96 – 2.20 | -3.42 – 2.58 |
| atac | sth1_0h | 4 | 12.94 – 25.46 (15.05 – 27.55) | -4.24 – 2.00 | -15.91 – 10.99 |
| atac | sth1_1h | 4 | 10.24 – 38.37 (13.60 – 41.09) | -4.47 – 2.10 | -19.29 – 12.54 |
| atac | snf2_0h | 3 | 24.33 – 48.14 (28.16 – 53.45) | -1.79 – 2.48 | -7.66 – 13.93 |
| atac | snf2_1h | 3 | 26.30 – 54.53 (29.04 – 59.06) | -1.27 – 2.48 | -8.18 – 15.13 |
| tf_chip | ZNF384WT | 3 | 1.15 – 3.01 (1.43 – 3.20) | -0.25 – 0.40 | -4.21 – 5.50 |
| tf_chip | TAF15 | 2 | 1.39 – 1.39 (1.83 – 1.83) | -0.18 – 0.18 | -1.63 – 1.63 |
| tf_chip | TCF3 | 2 | 1.96 – 1.96 (2.39 – 2.39) | -0.07 – 0.07 | -1.56 – 1.56 |
| tf_cutrun | A673 | 6 | 1.52 – 3.39 (2.31 – 4.26) | -3.60 – 2.17 | -6.67 – 4.57 |
| tf_cutrun | A2.2 | 6 | 1.30 – 7.86 (2.56 – 9.14) | -2.97 – 3.23 | -10.27 – 16.24 |
| proseq | NELFB | 2 | 4.86 – 4.86 (5.00 – 5.00) | -2.67 – 2.67 | -9.03 – 9.03 |
| proseq | PAF1 | 2 | 11.28 – 11.28 (11.56 – 11.56) | -5.45 – 5.45 | -9.15 – 9.15 |
| proseq | RTF1 | 2 | 8.91 – 8.91 (9.08 – 9.08) | -5.04 – 5.04 | -13.09 – 13.09 |
| k27ac_chip | CLP | 4 | 0.51 – 2.12 (0.66 – 2.27) | -0.89 – 0.38 | -3.00 – 1.10 |
| k27ac_chip | Sham | 4 | 0.28 – 1.55 (0.45 – 1.70) | -0.56 – 0.45 | -2.63 – 1.80 |
| k4me3_chip_summit | CLP | 4 | 0.66 – 3.04 (0.73 – 3.11) | -0.18 – 0.31 | -5.81 – 5.01 |
| k4me3_chip_summit | Sham | 4 | 1.41 – 1.72 (1.48 – 1.80) | -0.65 – 0.57 | -2.93 – 3.08 |
| k4me3_chip_tss | CLP | 4 | 1.38 – 3.70 (1.45 – 3.80) | -3.35 – 8.17 | -2.56 – 3.96 |
| k4me3_chip_tss | Sham | 4 | 0.75 – 1.69 (0.82 – 1.77) | -2.51 – 3.52 | -1.32 – 1.68 |
| k27ac_cutrun | Naive | 4 | 0.00 – 1.51 (0.27 – 1.69) | -0.35 – 0.17 | -1.97 – 1.29 |
| k27ac_cutrun | Arm | 4 | 0.21 – 0.44 (0.37 – 0.59) | -1.19 – 0.59 | -0.78 – 0.51 |
| k4me3_cutrun_summit | Naive | 4 | 0.23 – 1.44 (0.46 – 1.69) | -0.41 – 0.42 | -2.77 – 2.34 |
| k4me3_cutrun_summit | Arm | 4 | 0.00 – 0.76 (0.37 – 1.08) | -0.41 – 0.28 | -1.48 – 0.92 |
| k4me3_cutrun_tss | Naive | 4 | 0.66 – 2.10 (0.81 – 2.25) | -5.69 – 7.16 | -1.12 – 1.49 |
| k4me3_cutrun_tss | Arm | 4 | 0.40 – 0.84 (0.53 – 0.98) | -1.54 – 2.28 | -0.45 – 0.34 |

## Criteria

| id | result |
|---|---|
| M1 | pass |
| M2 | **FAIL** |
| M3 | **FAIL** |
| M4 | **FAIL** |
