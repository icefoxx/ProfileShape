# Within-group consistency — evaluation results

Protocol `docs/PROFILESHAPE_CONSISTENCY_EVALUATION.md` (sha256 `7e0a418d9a83108e`); commit `5f7b302cb4`, uncommitted src files 4 (diff `ce4afdd42a3aff64`); consistency.py `4d7146c664a1e425`; regional.py `21ffd82395ee207f`; script `3ba1da1541180030`; 2026-10-27T19:28:36 → 2026-10-27T20:19:12.

## A · Simulation, null cells (100 runs each, 2000 regions)

| n | depth | replicate 1 at 30 % depth | any flag | sample flag | split flag | offset p < 0.01 | split p < 0.01 | bp CI covers 0 | excess F20 | split π̂₁ |
|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 12 | no | 0.00 | 0.00 | 0.00 | 0.007 | nan | 0.95 | -0.001 | nan |
| 3 | 12 | yes | 0.01 | 0.01 | 0.00 | 0.013 | nan | 0.95 | -0.001 | nan |
| 3 | 40 | no | 0.01 | 0.01 | 0.00 | 0.027 | nan | 0.92 | -0.002 | nan |
| 3 | 40 | yes | 0.01 | 0.01 | 0.00 | 0.017 | nan | 0.95 | -0.001 | nan |
| 3 | 160 | no | 0.00 | 0.00 | 0.00 | 0.020 | nan | 0.92 | -0.001 | nan |
| 3 | 160 | yes | 0.00 | 0.00 | 0.00 | 0.013 | nan | 0.96 | -0.001 | nan |
| 4 | 12 | no | 0.01 | 0.01 | 0.00 | 0.018 | 0.010 | 0.95 | -0.001 | 0.004 |
| 4 | 12 | yes | 0.02 | 0.01 | 0.01 | 0.005 | 0.017 | 0.95 | +0.010 | 0.005 |
| 4 | 40 | no | 0.04 | 0.02 | 0.02 | 0.025 | 0.020 | 0.93 | -0.001 | 0.004 |
| 4 | 40 | yes | 0.03 | 0.02 | 0.01 | 0.020 | 0.010 | 0.94 | +0.003 | 0.004 |
| 4 | 160 | no | 0.01 | 0.01 | 0.00 | 0.010 | 0.010 | 0.96 | -0.001 | 0.003 |
| 4 | 160 | yes | 0.02 | 0.02 | 0.00 | 0.013 | 0.010 | 0.95 | -0.000 | 0.004 |
| 6 | 12 | no | 0.03 | 0.01 | 0.02 | 0.012 | 0.012 | 0.95 | -0.001 | 0.003 |
| 6 | 12 | yes | 0.00 | 0.00 | 0.00 | 0.013 | 0.008 | 0.95 | +0.013 | 0.044 |
| 6 | 40 | no | 0.04 | 0.01 | 0.03 | 0.010 | 0.012 | 0.94 | -0.001 | 0.002 |
| 6 | 40 | yes | 0.00 | 0.00 | 0.00 | 0.018 | 0.010 | 0.96 | +0.006 | 0.020 |
| 6 | 160 | no | 0.03 | 0.01 | 0.02 | 0.008 | 0.010 | 0.95 | -0.001 | 0.001 |
| 6 | 160 | yes | 0.01 | 0.01 | 0.00 | 0.015 | 0.004 | 0.95 | +0.000 | 0.004 |

Pooled equal-depth offset P(p < 0.01) = 0.0141 (m = 3900)

Pooled equal-depth split contrast P(p < 0.01) = 0.0114 (m = 8400)

Pooled unequal-depth offset P(p < 0.01) = 0.0144 (m = 3900)

Pooled unequal-depth split contrast P(p < 0.01) = 0.0080 (m = 8400)

## B · Semi-real null (labels permuted within blocks of 10 regions; 20 shuffles)

| dataset | group | n | regions | runs with any flag | offset p < 0.01 | split p < 0.01 | bp CI covers 0 | excess F20 | split π̂₁ |
|---|---|---|---|---|---|---|---|---|---|
| atac | snf2_0h | 3 | 4447 | 0 / 20 | 0.033 | nan | 1.00 | +0.004 | nan |
| atac | snf2_1h | 3 | 5476 | 1 / 20 | 0.017 | nan | 0.93 | -0.005 | nan |
| atac | sth1_0h | 4 | 4747 | 1 / 20 | 0.025 | 0.000 | 0.93 | -0.005 | 0.007 |
| atac (thinned) | sth1_0h | 4 | 4279 | 0 / 20 | 0.013 | 0.000 | 0.94 | +0.006 | 0.009 |
| atac | sth1_1h | 4 | 4389 | 0 / 20 | 0.000 | 0.000 | 0.91 | -0.001 | 0.007 |
| atac (thinned) | sth1_1h | 4 | 3895 | 0 / 20 | 0.000 | 0.000 | 0.94 | +0.004 | 0.010 |
| k27ac_chip | CLP | 4 | 5000 | 0 / 20 | 0.013 | 0.000 | 0.93 | -0.008 | 0.005 |
| k27ac_chip | Sham | 4 | 5000 | 0 / 20 | 0.000 | 0.000 | 0.99 | -0.003 | 0.004 |
| k27ac_cutrun | Arm | 4 | 4999 | 0 / 20 | 0.000 | 0.000 | 0.97 | -0.000 | 0.002 |
| k27ac_cutrun | Naive | 4 | 4989 | 0 / 20 | 0.000 | 0.000 | 0.99 | -0.003 | 0.004 |
| k4me3_chip_summit | CLP | 4 | 5000 | 0 / 20 | 0.013 | 0.000 | 0.97 | -0.025 | 0.004 |
| k4me3_chip_summit | Sham | 4 | 5000 | 0 / 20 | 0.000 | 0.000 | 0.95 | -0.006 | 0.004 |
| k4me3_chip_tss | CLP | 4 | 5000 | 0 / 20 | 0.000 | 0.000 | 0.97 | -0.014 | 0.003 |
| k4me3_chip_tss | Sham | 4 | 5000 | 0 / 20 | 0.013 | 0.000 | 0.99 | -0.010 | 0.004 |
| k4me3_cutrun_summit | Arm | 4 | 4999 | 0 / 20 | 0.000 | 0.000 | 0.95 | -0.002 | 0.002 |
| k4me3_cutrun_summit | Naive | 4 | 4998 | 0 / 20 | 0.013 | 0.000 | 0.91 | -0.003 | 0.004 |
| k4me3_cutrun_tss | Arm | 4 | 4994 | 1 / 20 | 0.013 | 0.000 | 0.93 | -0.003 | 0.002 |
| k4me3_cutrun_tss | Naive | 4 | 4994 | 0 / 20 | 0.000 | 0.000 | 0.96 | +0.002 | 0.002 |
| mnase | WT | 4 | 3592 | 0 / 20 | 0.000 | 0.017 | 0.85 | -0.001 | 0.006 |
| mnase | dKO | 4 | 3592 | 1 / 20 | 0.037 | 0.017 | 0.94 | -0.001 | 0.008 |
| mnase | isw1rtf1 | 4 | 3592 | 0 / 20 | 0.000 | 0.000 | 0.97 | -0.001 | 0.014 |
| mnase | rtf1 | 4 | 3592 | 0 / 20 | 0.025 | 0.000 | 0.95 | -0.000 | 0.012 |
| mnase | rtf1-11A | 4 | 3592 | 0 / 20 | 0.000 | 0.000 | 0.94 | +0.001 | 0.005 |
| mnase | rtf1-8-9A | 4 | 3592 | 0 / 20 | 0.025 | 0.000 | 0.99 | -0.002 | 0.006 |
| proseq | NELFB | 2 | 9024 | 0 / 20 | 0.000 | nan | 0.90 | +nan | nan |
| proseq | PAF1 | 2 | 9021 | 0 / 20 | 0.050 | nan | 0.95 | +nan | nan |
| proseq | RTF1 | 2 | 9023 | 0 / 20 | 0.000 | nan | 0.90 | +nan | nan |
| tf_chip | TAF15 | 2 | 4993 | 0 / 20 | 0.050 | nan | 0.85 | +nan | nan |
| tf_chip | TCF3 | 2 | 4993 | 0 / 20 | 0.000 | nan | 0.90 | +nan | nan |
| tf_chip | ZNF384WT | 3 | 4993 | 0 / 20 | 0.000 | nan | 0.98 | +0.001 | nan |
| tf_cutrun | A2.2 | 6 | 3250 | 0 / 20 | 0.008 | 0.000 | 0.97 | -0.006 | 0.002 |
| tf_cutrun | A673 | 6 | 3786 | 0 / 20 | 0.000 | 0.000 | 0.95 | -0.002 | 0.002 |

## A · Simulation, insertions (depth 40; replicate 1 or the second half altered)

| n | insertion | sample flag (replicate 1 named) | split flag | shift bp (replicate 1) | F20 − mates | F50 − mates | target split π̂₁ | target split called |
|---|---|---|---|---|---|---|---|---|
| 3 | one_broaden1.3_all | 1.00 | 0.00 | +0.0 | +0.059 | +0.043 | nan | nan |
| 3 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.106 | +0.074 | nan | nan |
| 3 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.1 | +0.036 | +0.038 | nan | nan |
| 3 | one_shift50_30pc_rand | 0.43 | 0.00 | +0.1 | +0.114 | +0.115 | nan | nan |
| 3 | one_shift50_30pc_same | 1.00 | 0.00 | +14.9 | +0.114 | +0.117 | nan | nan |
| 3 | one_shift50_all | 1.00 | 0.00 | +51.2 | +0.384 | +0.404 | nan | nan |
| 4 | half_shift20_all | 0.00 | 1.00 | -13.5 | +0.021 | +0.013 | 0.272 | 0 |
| 4 | half_shift50_20pc_rand | 0.00 | 1.00 | +0.0 | +0.015 | +0.012 | 0.211 | 4 |
| 4 | one_broaden1.3_all | 1.00 | 0.00 | +0.1 | +0.057 | +0.038 | 0.005 | 0 |
| 4 | one_shift20_all | 1.00 | 0.01 | +20.0 | +0.116 | +0.074 | 0.003 | 0 |
| 4 | one_shift50_10pc_rand | 0.03 | 0.01 | +0.0 | +0.037 | +0.038 | 0.004 | 0 |
| 4 | one_shift50_30pc_rand | 0.55 | 0.01 | +0.1 | +0.122 | +0.119 | 0.003 | 0 |
| 4 | one_shift50_30pc_same | 1.00 | 0.03 | +14.8 | +0.126 | +0.121 | 0.004 | 0 |
| 4 | one_shift50_all | 1.00 | 0.00 | +51.0 | +0.416 | +0.415 | 0.003 | 0 |
| 6 | half_shift20_all | 0.00 | 1.00 | -12.2 | +0.007 | +0.003 | 0.336 | 0 |
| 6 | half_shift50_20pc_rand | 0.00 | 1.00 | +0.0 | +0.004 | +0.004 | 0.184 | 34 |
| 6 | one_broaden1.3_all | 1.00 | 0.03 | -0.1 | +0.056 | +0.030 | 0.003 | 0 |
| 6 | one_shift20_all | 1.00 | 0.01 | +20.0 | +0.130 | +0.074 | 0.002 | 0 |
| 6 | one_shift50_10pc_rand | 0.00 | 0.01 | +0.1 | +0.045 | +0.040 | 0.001 | 0 |
| 6 | one_shift50_30pc_rand | 0.54 | 0.01 | -0.1 | +0.132 | +0.121 | 0.002 | 0 |
| 6 | one_shift50_30pc_same | 1.00 | 0.00 | +14.9 | +0.135 | +0.122 | 0.002 | 0 |
| 6 | one_shift50_all | 1.00 | 0.01 | +50.6 | +0.445 | +0.418 | 0.001 | 0 |

## C · Semi-real insertions (20 shuffles each; no criterion)

| dataset | group | n | insertion | sample flag (pseudo-replicate 1 named) | split flag | shift bp | F20 − mates | target split π̂₁ |
|---|---|---|---|---|---|---|---|---|
| atac | snf2_0h | 3 | one_broaden1.3_all | 1.00 | 0.00 | +4.2 | +0.076 | nan |
| atac | snf2_0h | 3 | one_shift20_all | 1.00 | 0.00 | +20.6 | +0.121 | nan |
| atac | snf2_0h | 3 | one_shift50_10pc_rand | 0.05 | 0.00 | +0.1 | +0.036 | nan |
| atac | snf2_0h | 3 | one_shift50_30pc_rand | 0.15 | 0.00 | +0.2 | +0.096 | nan |
| atac | snf2_0h | 3 | one_shift50_30pc_same | 1.00 | 0.00 | +11.6 | +0.094 | nan |
| atac | snf2_0h | 3 | one_shift50_all | 1.00 | 0.00 | +56.5 | +0.274 | nan |
| atac | snf2_1h | 3 | one_broaden1.3_all | 1.00 | 0.00 | +4.0 | +0.140 | nan |
| atac | snf2_1h | 3 | one_shift20_all | 1.00 | 0.00 | +21.2 | +0.218 | nan |
| atac | snf2_1h | 3 | one_shift50_10pc_rand | 0.10 | 0.00 | +0.0 | +0.044 | nan |
| atac | snf2_1h | 3 | one_shift50_30pc_rand | 0.95 | 0.00 | -0.3 | +0.138 | nan |
| atac | snf2_1h | 3 | one_shift50_30pc_same | 1.00 | 0.00 | +11.8 | +0.133 | nan |
| atac | snf2_1h | 3 | one_shift50_all | 1.00 | 0.00 | +59.7 | +0.438 | nan |
| atac | sth1_0h | 4 | half_shift20_all | 0.00 | 0.00 | -13.8 | +0.038 | 0.819 |
| atac | sth1_0h | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | +0.0 | -0.018 | 0.217 |
| atac | sth1_0h | 4 | one_broaden1.3_all | 0.80 | 0.00 | +5.0 | +0.136 | 0.010 |
| atac | sth1_0h | 4 | one_shift20_all | 1.00 | 0.00 | +21.7 | +0.214 | 0.007 |
| atac | sth1_0h | 4 | one_shift50_10pc_rand | 0.15 | 0.00 | +0.1 | +0.041 | 0.006 |
| atac | sth1_0h | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | +0.4 | +0.128 | 0.009 |
| atac | sth1_0h | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +12.0 | +0.122 | 0.008 |
| atac | sth1_0h | 4 | one_shift50_all | 1.00 | 0.00 | +59.6 | +0.425 | 0.006 |
| atac | sth1_1h | 4 | half_shift20_all | 0.00 | 0.00 | -14.0 | +0.038 | 0.806 |
| atac | sth1_1h | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.1 | -0.017 | 0.211 |
| atac | sth1_1h | 4 | one_broaden1.3_all | 0.50 | 0.00 | +4.7 | +0.128 | 0.010 |
| atac | sth1_1h | 4 | one_shift20_all | 1.00 | 0.00 | +21.5 | +0.199 | 0.009 |
| atac | sth1_1h | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | +0.2 | +0.042 | 0.009 |
| atac | sth1_1h | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | +0.4 | +0.125 | 0.007 |
| atac | sth1_1h | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +11.9 | +0.120 | 0.005 |
| atac | sth1_1h | 4 | one_shift50_all | 1.00 | 0.00 | +60.0 | +0.399 | 0.011 |
| k27ac_chip | CLP | 4 | half_shift20_all | 0.00 | 0.00 | -14.0 | +0.033 | 0.910 |
| k27ac_chip | CLP | 4 | half_shift50_20pc_rand | 0.00 | 0.15 | -0.1 | +0.001 | 0.289 |
| k27ac_chip | CLP | 4 | one_broaden1.3_all | 1.00 | 0.00 | +5.2 | +0.079 | 0.004 |
| k27ac_chip | CLP | 4 | one_shift20_all | 1.00 | 0.00 | +20.7 | +0.144 | 0.005 |
| k27ac_chip | CLP | 4 | one_shift50_10pc_rand | 0.10 | 0.00 | +0.0 | +0.033 | 0.006 |
| k27ac_chip | CLP | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | +0.0 | +0.112 | 0.009 |
| k27ac_chip | CLP | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +15.4 | +0.112 | 0.007 |
| k27ac_chip | CLP | 4 | one_shift50_all | 1.00 | 0.00 | +53.4 | +0.380 | 0.006 |
| k27ac_chip | Sham | 4 | half_shift20_all | 0.00 | 0.00 | -14.0 | +0.024 | 0.914 |
| k27ac_chip | Sham | 4 | half_shift50_20pc_rand | 0.00 | 0.15 | -0.1 | +0.010 | 0.291 |
| k27ac_chip | Sham | 4 | one_broaden1.3_all | 1.00 | 0.00 | +5.4 | +0.088 | 0.008 |
| k27ac_chip | Sham | 4 | one_shift20_all | 1.00 | 0.00 | +20.6 | +0.138 | 0.005 |
| k27ac_chip | Sham | 4 | one_shift50_10pc_rand | 0.05 | 0.00 | +0.0 | +0.033 | 0.006 |
| k27ac_chip | Sham | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | +0.2 | +0.109 | 0.006 |
| k27ac_chip | Sham | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +15.3 | +0.110 | 0.009 |
| k27ac_chip | Sham | 4 | one_shift50_all | 1.00 | 0.00 | +53.3 | +0.379 | 0.006 |
| k27ac_cutrun | Arm | 4 | half_shift20_all | 0.00 | 0.00 | -13.8 | +0.020 | 0.393 |
| k27ac_cutrun | Arm | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.2 | +0.009 | 0.238 |
| k27ac_cutrun | Arm | 4 | one_broaden1.3_all | 1.00 | 0.00 | +0.8 | +0.054 | 0.005 |
| k27ac_cutrun | Arm | 4 | one_shift20_all | 1.00 | 0.00 | +19.9 | +0.073 | 0.003 |
| k27ac_cutrun | Arm | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.2 | +0.026 | 0.002 |
| k27ac_cutrun | Arm | 4 | one_shift50_30pc_rand | 0.55 | 0.00 | -0.1 | +0.079 | 0.005 |
| k27ac_cutrun | Arm | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.8 | +0.079 | 0.004 |
| k27ac_cutrun | Arm | 4 | one_shift50_all | 1.00 | 0.00 | +50.7 | +0.267 | 0.004 |
| k27ac_cutrun | Naive | 4 | half_shift20_all | 0.00 | 0.00 | -13.5 | +0.007 | 0.333 |
| k27ac_cutrun | Naive | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | +0.0 | +0.004 | 0.217 |
| k27ac_cutrun | Naive | 4 | one_broaden1.3_all | 1.00 | 0.00 | +1.7 | +0.046 | 0.002 |
| k27ac_cutrun | Naive | 4 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.053 | 0.000 |
| k27ac_cutrun | Naive | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.0 | +0.019 | 0.002 |
| k27ac_cutrun | Naive | 4 | one_shift50_30pc_rand | 0.20 | 0.00 | -0.0 | +0.059 | 0.004 |
| k27ac_cutrun | Naive | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +15.1 | +0.058 | 0.003 |
| k27ac_cutrun | Naive | 4 | one_shift50_all | 1.00 | 0.00 | +50.9 | +0.199 | 0.000 |
| k4me3_chip_summit | CLP | 4 | half_shift20_all | 0.00 | 1.00 | -13.5 | +0.024 | 0.983 |
| k4me3_chip_summit | CLP | 4 | half_shift50_20pc_rand | 0.00 | 1.00 | -0.0 | -0.009 | 0.251 |
| k4me3_chip_summit | CLP | 4 | one_broaden1.3_all | 1.00 | 0.00 | +5.5 | +0.074 | 0.011 |
| k4me3_chip_summit | CLP | 4 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.149 | 0.003 |
| k4me3_chip_summit | CLP | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.0 | +0.017 | 0.004 |
| k4me3_chip_summit | CLP | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | -0.0 | +0.109 | 0.003 |
| k4me3_chip_summit | CLP | 4 | one_shift50_30pc_same | 1.00 | 0.05 | +14.9 | +0.110 | 0.004 |
| k4me3_chip_summit | CLP | 4 | one_shift50_all | 1.00 | 0.00 | +50.8 | +0.430 | 0.004 |
| k4me3_chip_summit | Sham | 4 | half_shift20_all | 0.00 | 1.00 | -13.3 | +0.033 | 0.993 |
| k4me3_chip_summit | Sham | 4 | half_shift50_20pc_rand | 0.00 | 1.00 | +0.1 | +0.007 | 0.227 |
| k4me3_chip_summit | Sham | 4 | one_broaden1.3_all | 1.00 | 0.00 | +5.6 | +0.055 | 0.006 |
| k4me3_chip_summit | Sham | 4 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.139 | 0.002 |
| k4me3_chip_summit | Sham | 4 | one_shift50_10pc_rand | 0.10 | 0.00 | -0.0 | +0.034 | 0.002 |
| k4me3_chip_summit | Sham | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | +0.3 | +0.119 | 0.002 |
| k4me3_chip_summit | Sham | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.9 | +0.119 | 0.003 |
| k4me3_chip_summit | Sham | 4 | one_shift50_all | 1.00 | 0.00 | +50.6 | +0.413 | 0.003 |
| k4me3_chip_tss | CLP | 4 | half_shift20_all | 0.00 | 1.00 | -13.6 | +0.093 | 0.985 |
| k4me3_chip_tss | CLP | 4 | half_shift50_20pc_rand | 0.00 | 1.00 | -0.1 | +0.009 | 0.249 |
| k4me3_chip_tss | CLP | 4 | one_broaden1.3_all | 1.00 | 0.00 | +40.7 | +0.475 | 0.008 |
| k4me3_chip_tss | CLP | 4 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.231 | 0.006 |
| k4me3_chip_tss | CLP | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.0 | +0.044 | 0.004 |
| k4me3_chip_tss | CLP | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | +0.1 | +0.153 | 0.007 |
| k4me3_chip_tss | CLP | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.8 | +0.154 | 0.006 |
| k4me3_chip_tss | CLP | 4 | one_shift50_all | 1.00 | 0.00 | +51.0 | +0.534 | 0.014 |
| k4me3_chip_tss | Sham | 4 | half_shift20_all | 0.00 | 0.00 | -13.6 | +0.063 | 0.995 |
| k4me3_chip_tss | Sham | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.0 | +0.011 | 0.238 |
| k4me3_chip_tss | Sham | 4 | one_broaden1.3_all | 1.00 | 0.00 | +39.1 | +0.440 | 0.011 |
| k4me3_chip_tss | Sham | 4 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.247 | 0.004 |
| k4me3_chip_tss | Sham | 4 | one_shift50_10pc_rand | 0.05 | 0.00 | -0.0 | +0.045 | 0.007 |
| k4me3_chip_tss | Sham | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | -0.1 | +0.151 | 0.008 |
| k4me3_chip_tss | Sham | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.8 | +0.152 | 0.007 |
| k4me3_chip_tss | Sham | 4 | one_shift50_all | 1.00 | 0.00 | +51.1 | +0.520 | 0.007 |
| k4me3_cutrun_summit | Arm | 4 | half_shift20_all | 0.00 | 0.00 | -13.5 | +0.035 | 0.742 |
| k4me3_cutrun_summit | Arm | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.0 | +0.012 | 0.281 |
| k4me3_cutrun_summit | Arm | 4 | one_broaden1.3_all | 1.00 | 0.00 | +1.0 | +0.081 | 0.001 |
| k4me3_cutrun_summit | Arm | 4 | one_shift20_all | 1.00 | 0.00 | +20.1 | +0.165 | 0.002 |
| k4me3_cutrun_summit | Arm | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | +0.0 | +0.041 | 0.003 |
| k4me3_cutrun_summit | Arm | 4 | one_shift50_30pc_rand | 0.90 | 0.00 | +0.0 | +0.126 | 0.005 |
| k4me3_cutrun_summit | Arm | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.8 | +0.126 | 0.005 |
| k4me3_cutrun_summit | Arm | 4 | one_shift50_all | 1.00 | 0.00 | +50.5 | +0.424 | 0.006 |
| k4me3_cutrun_summit | Naive | 4 | half_shift20_all | 0.00 | 0.00 | -13.6 | +0.025 | 0.600 |
| k4me3_cutrun_summit | Naive | 4 | half_shift50_20pc_rand | 0.05 | 0.00 | -0.2 | +0.005 | 0.265 |
| k4me3_cutrun_summit | Naive | 4 | one_broaden1.3_all | 1.00 | 0.00 | +0.1 | +0.077 | 0.007 |
| k4me3_cutrun_summit | Naive | 4 | one_shift20_all | 1.00 | 0.00 | +19.9 | +0.107 | 0.003 |
| k4me3_cutrun_summit | Naive | 4 | one_shift50_10pc_rand | 0.05 | 0.00 | -0.1 | +0.030 | 0.003 |
| k4me3_cutrun_summit | Naive | 4 | one_shift50_30pc_rand | 0.70 | 0.00 | -0.1 | +0.096 | 0.002 |
| k4me3_cutrun_summit | Naive | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.7 | +0.097 | 0.004 |
| k4me3_cutrun_summit | Naive | 4 | one_shift50_all | 1.00 | 0.00 | +50.3 | +0.327 | 0.003 |
| k4me3_cutrun_tss | Arm | 4 | half_shift20_all | 0.00 | 0.00 | -13.6 | +0.023 | 0.766 |
| k4me3_cutrun_tss | Arm | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.0 | +0.006 | 0.284 |
| k4me3_cutrun_tss | Arm | 4 | one_broaden1.3_all | 1.00 | 0.00 | -19.4 | +0.219 | 0.004 |
| k4me3_cutrun_tss | Arm | 4 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.117 | 0.004 |
| k4me3_cutrun_tss | Arm | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.2 | +0.027 | 0.005 |
| k4me3_cutrun_tss | Arm | 4 | one_shift50_30pc_rand | 0.75 | 0.00 | -0.0 | +0.087 | 0.006 |
| k4me3_cutrun_tss | Arm | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +15.0 | +0.088 | 0.004 |
| k4me3_cutrun_tss | Arm | 4 | one_shift50_all | 1.00 | 0.00 | +50.8 | +0.293 | 0.003 |
| k4me3_cutrun_tss | Naive | 4 | half_shift20_all | 0.00 | 0.00 | -13.6 | +0.020 | 0.623 |
| k4me3_cutrun_tss | Naive | 4 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.1 | +0.007 | 0.265 |
| k4me3_cutrun_tss | Naive | 4 | one_broaden1.3_all | 1.00 | 0.00 | -23.7 | +0.169 | 0.006 |
| k4me3_cutrun_tss | Naive | 4 | one_shift20_all | 1.00 | 0.00 | +20.1 | +0.070 | 0.002 |
| k4me3_cutrun_tss | Naive | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.2 | +0.023 | 0.002 |
| k4me3_cutrun_tss | Naive | 4 | one_shift50_30pc_rand | 0.30 | 0.00 | +0.0 | +0.067 | 0.001 |
| k4me3_cutrun_tss | Naive | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.9 | +0.066 | 0.002 |
| k4me3_cutrun_tss | Naive | 4 | one_shift50_all | 1.00 | 0.00 | +50.7 | +0.209 | 0.002 |
| mnase | WT | 4 | half_shift20_all | 0.00 | 0.90 | -13.8 | +0.127 | 1.000 |
| mnase | WT | 4 | half_shift50_20pc_rand | 0.00 | 0.90 | -0.0 | -0.005 | 0.333 |
| mnase | WT | 4 | one_broaden1.3_all | 1.00 | 0.00 | -11.4 | +0.806 | 0.011 |
| mnase | WT | 4 | one_shift20_all | 1.00 | 0.00 | +25.6 | +0.674 | 0.012 |
| mnase | WT | 4 | one_shift50_10pc_rand | 0.25 | 0.00 | +0.1 | +0.062 | 0.012 |
| mnase | WT | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | -0.4 | +0.201 | 0.007 |
| mnase | WT | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +9.5 | +0.206 | 0.009 |
| mnase | WT | 4 | one_shift50_all | 1.00 | 0.00 | +69.1 | +0.774 | 0.012 |
| mnase | dKO | 4 | half_shift20_all | 0.00 | 0.60 | -13.4 | +0.031 | 1.000 |
| mnase | dKO | 4 | half_shift50_20pc_rand | 0.00 | 0.90 | +0.0 | +0.003 | 0.286 |
| mnase | dKO | 4 | one_broaden1.3_all | 1.00 | 0.00 | -29.9 | +0.384 | 0.010 |
| mnase | dKO | 4 | one_shift20_all | 1.00 | 0.00 | +25.3 | +0.270 | 0.004 |
| mnase | dKO | 4 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.1 | +0.042 | 0.008 |
| mnase | dKO | 4 | one_shift50_30pc_rand | 0.30 | 0.00 | -0.2 | +0.128 | 0.004 |
| mnase | dKO | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +14.7 | +0.131 | 0.007 |
| mnase | dKO | 4 | one_shift50_all | 1.00 | 0.00 | +69.3 | +0.397 | 0.003 |
| mnase | isw1rtf1 | 4 | half_shift20_all | 0.00 | 1.00 | -13.7 | +0.058 | 1.000 |
| mnase | isw1rtf1 | 4 | half_shift50_20pc_rand | 0.00 | 0.75 | +0.0 | -0.002 | 0.307 |
| mnase | isw1rtf1 | 4 | one_broaden1.3_all | 1.00 | 0.00 | -24.6 | +0.556 | 0.007 |
| mnase | isw1rtf1 | 4 | one_shift20_all | 1.00 | 0.00 | +25.8 | +0.438 | 0.012 |
| mnase | isw1rtf1 | 4 | one_shift50_10pc_rand | 0.05 | 0.00 | -0.1 | +0.050 | 0.012 |
| mnase | isw1rtf1 | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | -0.4 | +0.167 | 0.012 |
| mnase | isw1rtf1 | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +13.4 | +0.171 | 0.011 |
| mnase | isw1rtf1 | 4 | one_shift50_all | 1.00 | 0.00 | +71.6 | +0.553 | 0.011 |
| mnase | rtf1 | 4 | half_shift20_all | 0.00 | 1.00 | -13.7 | +0.089 | 1.000 |
| mnase | rtf1 | 4 | half_shift50_20pc_rand | 0.00 | 0.85 | +0.0 | -0.009 | 0.321 |
| mnase | rtf1 | 4 | one_broaden1.3_all | 1.00 | 0.00 | -16.4 | +0.672 | 0.008 |
| mnase | rtf1 | 4 | one_shift20_all | 1.00 | 0.05 | +25.5 | +0.524 | 0.009 |
| mnase | rtf1 | 4 | one_shift50_10pc_rand | 0.05 | 0.00 | -0.1 | +0.059 | 0.018 |
| mnase | rtf1 | 4 | one_shift50_30pc_rand | 0.95 | 0.00 | -0.3 | +0.186 | 0.012 |
| mnase | rtf1 | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +11.1 | +0.189 | 0.009 |
| mnase | rtf1 | 4 | one_shift50_all | 1.00 | 0.00 | +71.7 | +0.637 | 0.009 |
| mnase | rtf1-11A | 4 | half_shift20_all | 0.00 | 0.75 | -13.7 | +0.091 | 1.000 |
| mnase | rtf1-11A | 4 | half_shift50_20pc_rand | 0.00 | 0.80 | +0.1 | -0.005 | 0.324 |
| mnase | rtf1-11A | 4 | one_broaden1.3_all | 1.00 | 0.00 | -20.4 | +0.687 | 0.010 |
| mnase | rtf1-11A | 4 | one_shift20_all | 1.00 | 0.00 | +25.7 | +0.535 | 0.008 |
| mnase | rtf1-11A | 4 | one_shift50_10pc_rand | 0.15 | 0.00 | -0.1 | +0.055 | 0.010 |
| mnase | rtf1-11A | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | -0.4 | +0.173 | 0.006 |
| mnase | rtf1-11A | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +11.2 | +0.179 | 0.007 |
| mnase | rtf1-11A | 4 | one_shift50_all | 1.00 | 0.00 | +70.6 | +0.644 | 0.012 |
| mnase | rtf1-8-9A | 4 | half_shift20_all | 0.00 | 0.40 | -13.9 | +0.081 | 1.000 |
| mnase | rtf1-8-9A | 4 | half_shift50_20pc_rand | 0.00 | 0.80 | -0.3 | -0.007 | 0.339 |
| mnase | rtf1-8-9A | 4 | one_broaden1.3_all | 1.00 | 0.00 | -22.3 | +0.638 | 0.007 |
| mnase | rtf1-8-9A | 4 | one_shift20_all | 1.00 | 0.00 | +25.4 | +0.485 | 0.004 |
| mnase | rtf1-8-9A | 4 | one_shift50_10pc_rand | 0.15 | 0.00 | -0.3 | +0.049 | 0.014 |
| mnase | rtf1-8-9A | 4 | one_shift50_30pc_rand | 1.00 | 0.00 | -0.4 | +0.168 | 0.014 |
| mnase | rtf1-8-9A | 4 | one_shift50_30pc_same | 1.00 | 0.00 | +11.5 | +0.171 | 0.010 |
| mnase | rtf1-8-9A | 4 | one_shift50_all | 1.00 | 0.00 | +71.2 | +0.598 | 0.008 |
| tf_chip | ZNF384WT | 3 | one_broaden1.3_all | 1.00 | 0.00 | +0.0 | +0.048 | nan |
| tf_chip | ZNF384WT | 3 | one_shift20_all | 1.00 | 0.00 | +20.2 | +0.107 | nan |
| tf_chip | ZNF384WT | 3 | one_shift50_10pc_rand | 0.05 | 0.00 | +0.0 | +0.036 | nan |
| tf_chip | ZNF384WT | 3 | one_shift50_30pc_rand | 1.00 | 0.00 | +0.1 | +0.109 | nan |
| tf_chip | ZNF384WT | 3 | one_shift50_30pc_same | 1.00 | 0.00 | +14.9 | +0.109 | nan |
| tf_chip | ZNF384WT | 3 | one_shift50_all | 1.00 | 0.00 | +50.8 | +0.356 | nan |
| tf_cutrun | A2.2 | 6 | half_shift20_all | 0.00 | 0.00 | -12.3 | -0.004 | 0.082 |
| tf_cutrun | A2.2 | 6 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.3 | -0.006 | 0.032 |
| tf_cutrun | A2.2 | 6 | one_broaden1.3_all | 1.00 | 0.00 | -0.5 | -0.000 | 0.001 |
| tf_cutrun | A2.2 | 6 | one_shift20_all | 1.00 | 0.00 | +19.8 | -0.002 | 0.001 |
| tf_cutrun | A2.2 | 6 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.1 | -0.003 | 0.001 |
| tf_cutrun | A2.2 | 6 | one_shift50_30pc_rand | 0.05 | 0.00 | -0.5 | -0.000 | 0.001 |
| tf_cutrun | A2.2 | 6 | one_shift50_30pc_same | 1.00 | 0.00 | +14.1 | +0.000 | 0.000 |
| tf_cutrun | A2.2 | 6 | one_shift50_all | 1.00 | 0.00 | +49.7 | +0.013 | 0.000 |
| tf_cutrun | A673 | 6 | half_shift20_all | 0.00 | 0.00 | -12.4 | -0.007 | 0.177 |
| tf_cutrun | A673 | 6 | half_shift50_20pc_rand | 0.00 | 0.00 | -0.1 | -0.007 | 0.075 |
| tf_cutrun | A673 | 6 | one_broaden1.3_all | 0.50 | 0.00 | -2.7 | +0.002 | 0.001 |
| tf_cutrun | A673 | 6 | one_shift20_all | 1.00 | 0.00 | +19.9 | -0.002 | 0.001 |
| tf_cutrun | A673 | 6 | one_shift50_10pc_rand | 0.00 | 0.00 | -0.0 | -0.001 | 0.001 |
| tf_cutrun | A673 | 6 | one_shift50_30pc_rand | 0.05 | 0.00 | -0.4 | +0.003 | 0.001 |
| tf_cutrun | A673 | 6 | one_shift50_30pc_same | 1.00 | 0.00 | +14.3 | +0.002 | 0.000 |
| tf_cutrun | A673 | 6 | one_shift50_all | 1.00 | 0.00 | +50.2 | +0.022 | 0.003 |

## D · Real data as is (descriptive; no truth)

**mnase · WT** (n = 4, 3592 regions): g: WT_3 stands out from its group-mates: its sample-wide shape offset is 2.8× theirs (9.6% of the aggregate profile; p = 1.3e-60, Bonferroni over 4 samples). Check fragment length, alignment and IP quality
- WT_1: offset 6.1 % of profile (1.15× mates), shift -1.2 bp [-1.7, -0.8], F20 0.11 (mates 0.10)
- WT_2: offset 3.5 % of profile (0.48× mates), shift +0.1 bp [-0.3, +0.5], F20 0.11 (mates 0.10)
- WT_3: offset 9.6 % of profile (2.81× mates), shift +1.6 bp [+1.2, +2.0], F20 0.10 (mates 0.10)
- WT_4: offset 2.6 % of profile (0.38× mates), shift -0.4 bp [-0.8, -0.0], F20 0.09 (mates 0.11)
- split WT_1+WT_2 | WT_3+WT_4: z +2.16, p 0.015, π̂₁ 0.101, called 0, shift +0.9 bp
- split WT_1+WT_4 | WT_2+WT_3: z +0.79, p 0.21, π̂₁ 0.026, called 0, shift +1.2 bp
- split WT_1+WT_3 | WT_2+WT_4: z -2.89, p 1, π̂₁ 0.000, called 0, shift -0.3 bp

**mnase · dKO** (n = 4, 3592 regions): g: replicates form two sub-groups (dKO_1+dKO_2 | dKO_3+dKO_4; per-region distance contrast p = 3.8e-28, Bonferroni over 3 splits; meta-profile across/within ratio 1.49). Contrasts or null splits along this partition are not between exchangeable replicates
- dKO_1: offset 33.0 % of profile (1.56× mates), shift -1.2 bp [-2.3, -0.1], F20 0.52 (mates 0.47)
- dKO_2: offset 9.6 % of profile (0.27× mates), shift -0.4 bp [-1.2, +0.5], F20 0.41 (mates 0.52)
- dKO_3: offset 5.3 % of profile (0.16× mates), shift +0.7 bp [-0.2, +1.6], F20 0.42 (mates 0.50)
- dKO_4: offset 37.3 % of profile (2.00× mates), shift +0.9 bp [-0.2, +1.9], F20 0.52 (mates 0.46)
- split dKO_1+dKO_2 | dKO_3+dKO_4: z +10.94, p 3.8e-28, π̂₁ 0.652, called 0, shift +1.2 bp
- split dKO_1+dKO_3 | dKO_2+dKO_4: z -2.37, p 0.99, π̂₁ 0.000, called 0, shift +0.4 bp
- split dKO_1+dKO_4 | dKO_2+dKO_3: z -13.91, p 1, π̂₁ 0.000, called 0, shift +0.3 bp

**mnase · rtf1-8-9A** (n = 4, 3592 regions): g: rtf1-8-9A_1 stands out from its group-mates: its sample-wide shape offset is 2.1× theirs (5.6% of the aggregate profile; p = 2.1e-18, Bonferroni over 4 samples). Check fragment length, alignment and IP quality
- rtf1-8-9A_1: offset 5.6 % of profile (2.10× mates), shift -1.2 bp [-1.9, -0.5], F20 0.29 (mates 0.25)
- rtf1-8-9A_2: offset 3.7 % of profile (0.83× mates), shift +0.7 bp [+0.1, +1.3], F20 0.24 (mates 0.29)
- rtf1-8-9A_3: offset 2.8 % of profile (0.73× mates), shift -0.5 bp [-1.1, +0.2], F20 0.30 (mates 0.24)
- rtf1-8-9A_4: offset 3.3 % of profile (0.72× mates), shift +0.9 bp [+0.3, +1.5], F20 0.26 (mates 0.29)
- split rtf1-8-9A_1+rtf1-8-9A_3 | rtf1-8-9A_2+rtf1-8-9A_4: z +1.03, p 0.15, π̂₁ 0.005, called 0, shift +1.2 bp
- split rtf1-8-9A_1+rtf1-8-9A_4 | rtf1-8-9A_2+rtf1-8-9A_3: z +0.03, p 0.49, π̂₁ 0.000, called 0, shift +0.2 bp
- split rtf1-8-9A_1+rtf1-8-9A_2 | rtf1-8-9A_3+rtf1-8-9A_4: z -1.77, p 0.96, π̂₁ 0.000, called 0, shift +0.4 bp

**mnase · rtf1-11A** (n = 4, 3592 regions): g: rtf1-11A_3 stands out from its group-mates: its sample-wide shape offset is 2.3× theirs (11.6% of the aggregate profile; p = 2.8e-96, Bonferroni over 4 samples). Check fragment length, alignment and IP quality; g: replicates form two sub-groups (rtf1-11A_1+rtf1-11A_2 | rtf1-11A_3+rtf1-11A_4; per-region distance contrast p = 2.7e-14, Bonferroni over 3 splits; meta-profile across/within ratio 2.24). Contrasts or null splits along this partition are not between exchangeable replicates
- rtf1-11A_1: offset 6.6 % of profile (0.65× mates), shift -1.3 bp [-1.9, -0.8], F20 0.24 (mates 0.22)
- rtf1-11A_2: offset 8.3 % of profile (0.90× mates), shift -1.1 bp [-1.6, -0.5], F20 0.22 (mates 0.24)
- rtf1-11A_3: offset 11.6 % of profile (2.34× mates), shift +1.9 bp [+1.3, +2.4], F20 0.24 (mates 0.21)
- rtf1-11A_4: offset 3.1 % of profile (0.39× mates), shift +0.5 bp [-0.0, +1.1], F20 0.20 (mates 0.24)
- split rtf1-11A_1+rtf1-11A_2 | rtf1-11A_3+rtf1-11A_4: z +7.52, p 2.7e-14, π̂₁ 0.333, called 0, shift +1.8 bp
- split rtf1-11A_1+rtf1-11A_4 | rtf1-11A_2+rtf1-11A_3: z -2.23, p 0.99, π̂₁ 0.000, called 0, shift +0.6 bp
- split rtf1-11A_1+rtf1-11A_3 | rtf1-11A_2+rtf1-11A_4: z -4.36, p 1, π̂₁ 0.000, called 0, shift -0.4 bp

**mnase · rtf1** (n = 4, 3592 regions): g: rtf1_1 stands out from its group-mates: its sample-wide shape offset is 2.0× theirs (20.5% of the aggregate profile; p = 5.5e-151, Bonferroni over 4 samples). Check fragment length, alignment and IP quality; g: replicates form two sub-groups (rtf1_1+rtf1_2 | rtf1_3+rtf1_4; per-region distance contrast p = 1.9e-56, Bonferroni over 3 splits; meta-profile across/within ratio 2.77). Contrasts or null splits along this partition are not between exchangeable replicates
- rtf1_1: offset 20.5 % of profile (2.04× mates), shift -3.7 bp [-4.4, -3.1], F20 0.31 (mates 0.21)
- rtf1_2: offset 5.4 % of profile (0.33× mates), shift -0.6 bp [-1.1, -0.0], F20 0.21 (mates 0.25)
- rtf1_3: offset 8.8 % of profile (0.45× mates), shift +1.6 bp [+1.1, +2.1], F20 0.19 (mates 0.27)
- rtf1_4: offset 16.9 % of profile (1.06× mates), shift +2.7 bp [+2.1, +3.2], F20 0.24 (mates 0.22)
- split rtf1_1+rtf1_2 | rtf1_3+rtf1_4: z +15.79, p 1.9e-56, π̂₁ 0.669, called 0, shift +3.2 bp
- split rtf1_1+rtf1_3 | rtf1_2+rtf1_4: z -5.81, p 1, π̂₁ 0.000, called 0, shift +1.6 bp
- split rtf1_1+rtf1_4 | rtf1_2+rtf1_3: z -12.72, p 1, π̂₁ 0.000, called 0, shift +0.8 bp

**mnase · isw1rtf1** (n = 4, 3592 regions): g: replicates form two sub-groups (isw1rtf1_1+isw1rtf1_2 | isw1rtf1_3+isw1rtf1_4; per-region distance contrast p = 6e-91, Bonferroni over 3 splits; meta-profile across/within ratio 3.77). Contrasts or null splits along this partition are not between exchangeable replicates
- isw1rtf1_1: offset 13.2 % of profile (1.09× mates), shift -2.0 bp [-2.7, -1.2], F20 0.36 (mates 0.34)
- isw1rtf1_2: offset 11.8 % of profile (0.93× mates), shift -2.0 bp [-2.7, -1.2], F20 0.37 (mates 0.33)
- isw1rtf1_3: offset 8.8 % of profile (0.77× mates), shift +2.2 bp [+1.6, +2.8], F20 0.29 (mates 0.37)
- isw1rtf1_4: offset 16.5 % of profile (1.89× mates), shift +1.7 bp [+1.0, +2.4], F20 0.36 (mates 0.34)
- split isw1rtf1_1+isw1rtf1_2 | isw1rtf1_3+isw1rtf1_4: z +20.19, p 6e-91, π̂₁ 0.721, called 187, shift +2.9 bp
- split isw1rtf1_1+isw1rtf1_3 | isw1rtf1_2+isw1rtf1_4: z -15.49, p 1, π̂₁ 0.000, called 0, shift -0.2 bp
- split isw1rtf1_1+isw1rtf1_4 | isw1rtf1_2+isw1rtf1_3: z -14.08, p 1, π̂₁ 0.000, called 0, shift +0.2 bp

**atac · sth1_0h** (n = 4, 4747 regions): no flag
- sth1_0h_rep1: offset 21.7 % of profile (0.89× mates), shift +2.3 bp [+1.5, +3.1], F20 0.39 (mates 0.41)
- sth1_0h_rep2: offset 26.1 % of profile (1.13× mates), shift +1.5 bp [+0.7, +2.3], F20 0.38 (mates 0.40)
- sth1_0h_rep3: offset 17.5 % of profile (0.62× mates), shift +0.3 bp [-0.5, +1.2], F20 0.41 (mates 0.37)
- sth1_0h_rep4: offset 30.9 % of profile (1.60× mates), shift -4.1 bp [-4.9, -3.3], F20 0.37 (mates 0.41)
- split sth1_0h_rep1+sth1_0h_rep4 | sth1_0h_rep2+sth1_0h_rep3: z +1.00, p 0.16, π̂₁ 0.145, called 0, shift +1.3 bp
- split sth1_0h_rep1+sth1_0h_rep3 | sth1_0h_rep2+sth1_0h_rep4: z +1.00, p 0.16, π̂₁ 0.000, called 0, shift -2.0 bp
- split sth1_0h_rep1+sth1_0h_rep2 | sth1_0h_rep3+sth1_0h_rep4: z -1.00, p 0.84, π̂₁ 0.021, called 0, shift -2.8 bp

**atac · sth1_1h** (n = 4, 4389 regions): no flag
- sth1_1h_rep1: offset 16.7 % of profile (0.40× mates), shift +1.0 bp [+0.2, +1.8], F20 0.35 (mates 0.45)
- sth1_1h_rep2: offset 40.0 % of profile (1.36× mates), shift +1.6 bp [+0.8, +2.5], F20 0.43 (mates 0.39)
- sth1_1h_rep3: offset 17.7 % of profile (0.40× mates), shift +1.7 bp [+0.8, +2.6], F20 0.42 (mates 0.41)
- sth1_1h_rep4: offset 48.6 % of profile (1.60× mates), shift -4.3 bp [-5.2, -3.4], F20 0.43 (mates 0.40)
- split sth1_1h_rep1+sth1_1h_rep4 | sth1_1h_rep2+sth1_1h_rep3: z +1.00, p 0.16, π̂₁ 0.208, called 0, shift +2.5 bp
- split sth1_1h_rep1+sth1_1h_rep3 | sth1_1h_rep2+sth1_1h_rep4: z -1.00, p 0.84, π̂₁ 0.000, called 0, shift -2.0 bp
- split sth1_1h_rep1+sth1_1h_rep2 | sth1_1h_rep3+sth1_1h_rep4: z -1.00, p 0.84, π̂₁ 0.064, called 0, shift -1.9 bp

**atac · snf2_0h** (n = 3, 4447 regions): g: snf2_0h_rep1 stands out from its group-mates: its sample-wide shape offset is 8.5× theirs (68.0% of the aggregate profile; p = 2e-54, Bonferroni over 3 samples). Check fragment length, alignment and IP quality
- snf2_0h_rep1: offset 68.0 % of profile (8.51× mates), shift +2.1 bp [+0.7, +3.4], F20 0.56 (mates 0.46)
- snf2_0h_rep2: offset 35.2 % of profile (0.55× mates), shift -1.4 bp [-2.5, -0.3], F20 0.51 (mates 0.59)
- snf2_0h_rep3: offset 33.7 % of profile (0.48× mates), shift -0.6 bp [-1.7, +0.4], F20 0.52 (mates 0.53)

**atac · snf2_1h** (n = 3, 5476 regions): g: snf2_1h_rep1 stands out from its group-mates: its sample-wide shape offset is 26.2× theirs (77.0% of the aggregate profile; p = 7.2e-109, Bonferroni over 3 samples). Check fragment length, alignment and IP quality
- snf2_1h_rep1: offset 77.0 % of profile (26.19× mates), shift +2.4 bp [+1.5, +3.4], F20 0.48 (mates 0.28)
- snf2_1h_rep2: offset 36.2 % of profile (0.51× mates), shift -1.3 bp [-1.9, -0.6], F20 0.37 (mates 0.43)
- snf2_1h_rep3: offset 40.9 % of profile (0.61× mates), shift -1.2 bp [-1.8, -0.5], F20 0.33 (mates 0.48)

**tf_chip · ZNF384WT** (n = 3, 4993 regions): g: ZNF384WT_2 stands out from its group-mates: its sample-wide shape offset is 2.3× theirs (3.0% of the aggregate profile; p = 4.8e-143, Bonferroni over 3 samples). Check fragment length, alignment and IP quality
- ZNF384WT_1: offset 2.4 % of profile (1.05× mates), shift +0.4 bp [-0.3, +1.1], F20 0.43 (mates 0.56)
- ZNF384WT_2: offset 3.0 % of profile (2.27× mates), shift -0.2 bp [-0.9, +0.6], F20 0.44 (mates 0.54)
- ZNF384WT_3: offset 1.2 % of profile (0.40× mates), shift -0.2 bp [-1.3, +0.8], F20 0.59 (mates 0.33)

**tf_chip · TAF15** (n = 2, 4993 regions): no flag
- TAF15_1: offset 1.5 % of profile (nan× mates), shift +0.2 bp [-1.1, +1.4], F20 0.63 (mates nan)
- TAF15_2: offset 1.5 % of profile (nan× mates), shift -0.2 bp [-1.4, +1.1], F20 0.63 (mates nan)

**tf_chip · TCF3** (n = 2, 4993 regions): no flag
- TCF3_1: offset 2.1 % of profile (nan× mates), shift -0.1 bp [-1.3, +1.0], F20 0.62 (mates nan)
- TCF3_2: offset 2.1 % of profile (nan× mates), shift +0.1 bp [-1.0, +1.3], F20 0.62 (mates nan)

**tf_cutrun · A673** (n = 6, 3786 regions): no flag
- A673_1: offset 2.3 % of profile (0.95× mates), shift +2.3 bp [+0.1, +4.4], F20 0.73 (mates 0.74)
- A673_2: offset 1.5 % of profile (0.68× mates), shift -0.6 bp [-2.4, +1.3], F20 0.71 (mates 0.75)
- A673_3: offset 1.8 % of profile (0.73× mates), shift -0.1 bp [-2.1, +1.9], F20 0.73 (mates 0.75)
- A673_4: offset 3.7 % of profile (1.62× mates), shift +0.3 bp [-1.8, +2.4], F20 0.73 (mates 0.76)
- A673_5: offset 2.3 % of profile (0.94× mates), shift -2.3 bp [-6.1, +1.4], F20 0.84 (mates 0.72)
- A673_6: offset 2.9 % of profile (1.09× mates), shift +0.4 bp [-3.4, +4.3], F20 0.84 (mates 0.72)
- split A673_1+A673_2+A673_5 | A673_3+A673_4+A673_6: z +1.37, p 0.085, π̂₁ 0.000, called 0, shift +0.4 bp
- split A673_3+A673_4 | A673_1+A673_2+A673_5+A673_6: z +1.32, p 0.093, π̂₁ 0.000, called 0, shift -0.1 bp
- split A673_4+A673_6 | A673_1+A673_2+A673_3+A673_5: z +1.15, p 0.12, π̂₁ 0.153, called 0, shift -0.5 bp

**tf_cutrun · A2.2** (n = 6, 3250 regions): g: A2.2_2 stands out from its group-mates: its sample-wide shape offset is 2.0× theirs (8.1% of the aggregate profile; p = 3.8e-26, Bonferroni over 6 samples). Check fragment length, alignment and IP quality; it is also much shallower than its group-mates; g: A2.2_5 stands out from its group-mates: its sample-wide shape offset is 2.2× theirs (8.1% of the aggregate profile; p = 8.5e-15, Bonferroni over 6 samples). Check fragment length, alignment and IP quality
- A2.2_1: offset 3.1 % of profile (0.67× mates), shift -0.5 bp [-2.7, +1.7], F20 0.74 (mates 0.86)
- A2.2_2: offset 8.1 % of profile (2.05× mates), shift -0.2 bp [-4.4, +4.0], F20 0.85 (mates 0.82)
- A2.2_3: offset 4.1 % of profile (0.97× mates), shift +1.3 bp [-1.4, +3.9], F20 0.79 (mates 0.85)
- A2.2_4: offset 2.6 % of profile (0.55× mates), shift -0.5 bp [-3.9, +3.0], F20 0.82 (mates 0.85)
- A2.2_5: offset 8.1 % of profile (2.23× mates), shift +2.9 bp [-1.9, +7.7], F20 0.87 (mates 0.82)
- A2.2_6: offset 4.3 % of profile (1.21× mates), shift -2.9 bp [-7.7, +1.9], F20 0.88 (mates 0.82)
- split A2.2_4+A2.2_6 | A2.2_1+A2.2_2+A2.2_3+A2.2_5: z +2.37, p 0.009, π̂₁ 0.130, called 0, shift +2.1 bp
- split A2.2_1+A2.2_3+A2.2_6 | A2.2_2+A2.2_4+A2.2_5: z +1.69, p 0.045, π̂₁ 0.000, called 0, shift +1.2 bp
- split A2.2_2+A2.2_5 | A2.2_1+A2.2_3+A2.2_4+A2.2_6: z +1.63, p 0.052, π̂₁ 0.217, called 0, shift -1.7 bp

**proseq · NELFB** (n = 2, 9024 regions): no flag
- NELFB_1: offset 4.9 % of profile (nan× mates), shift +2.7 bp [+2.5, +2.8], F20 0.02 (mates nan)
- NELFB_2: offset 4.9 % of profile (nan× mates), shift -2.7 bp [-2.8, -2.5], F20 0.02 (mates nan)

**proseq · PAF1** (n = 2, 9021 regions): no flag
- PAF1_1: offset 11.3 % of profile (nan× mates), shift -5.5 bp [-5.6, -5.3], F20 0.05 (mates nan)
- PAF1_2: offset 11.3 % of profile (nan× mates), shift +5.5 bp [+5.3, +5.6], F20 0.05 (mates nan)

**proseq · RTF1** (n = 2, 9023 regions): no flag
- RTF1_1: offset 8.9 % of profile (nan× mates), shift +5.0 bp [+4.9, +5.2], F20 0.02 (mates nan)
- RTF1_2: offset 8.9 % of profile (nan× mates), shift -5.0 bp [-5.2, -4.9], F20 0.02 (mates nan)

**k27ac_chip · CLP** (n = 4, 5000 regions): g: CLP_4 stands out from its group-mates: its sample-wide shape offset is 4.3× theirs (2.1% of the aggregate profile; p = 1.5e-101, Bonferroni over 4 samples). Check fragment length, alignment and IP quality; it is also much shallower than its group-mates
- CLP_1: offset 1.1 % of profile (0.99× mates), shift +0.2 bp [-0.3, +0.7], F20 0.27 (mates 0.32)
- CLP_2: offset 0.6 % of profile (0.47× mates), shift -0.9 bp [-1.5, -0.3], F20 0.30 (mates 0.30)
- CLP_3: offset 0.7 % of profile (0.54× mates), shift +0.4 bp [-0.2, +1.0], F20 0.29 (mates 0.30)
- CLP_4: offset 2.1 % of profile (4.28× mates), shift +0.3 bp [-0.3, +1.0], F20 0.34 (mates 0.29)
- split CLP_1+CLP_2 | CLP_3+CLP_4: z +1.04, p 0.15, π̂₁ 0.146, called 0, shift +0.5 bp
- split CLP_1+CLP_4 | CLP_2+CLP_3: z +0.17, p 0.43, π̂₁ 0.000, called 0, shift -0.4 bp
- split CLP_1+CLP_3 | CLP_2+CLP_4: z -1.02, p 0.85, π̂₁ 0.000, called 0, shift -0.4 bp

**k27ac_chip · Sham** (n = 4, 5000 regions): no flag
- Sham_1: offset 1.6 % of profile (1.47× mates), shift +0.5 bp [-0.1, +1.1], F20 0.31 (mates 0.31)
- Sham_2: offset 1.2 % of profile (0.87× mates), shift +0.3 bp [-0.2, +0.8], F20 0.30 (mates 0.32)
- Sham_3: offset 1.4 % of profile (1.22× mates), shift -0.6 bp [-1.2, +0.1], F20 0.34 (mates 0.32)
- Sham_4: offset 0.4 % of profile (0.30× mates), shift -0.2 bp [-0.9, +0.4], F20 0.32 (mates 0.31)
- split Sham_1+Sham_2 | Sham_3+Sham_4: z +1.01, p 0.16, π̂₁ 0.044, called 0, shift -0.6 bp
- split Sham_1+Sham_3 | Sham_2+Sham_4: z -1.00, p 0.84, π̂₁ 0.107, called 0, shift +0.1 bp
- split Sham_1+Sham_4 | Sham_2+Sham_3: z -1.03, p 0.85, π̂₁ 0.000, called 0, shift -0.2 bp

**k4me3_chip_summit · CLP** (n = 4, 5000 regions): g: CLP_4 stands out from its group-mates: its sample-wide shape offset is 2.1× theirs (3.0% of the aggregate profile; p = 0, Bonferroni over 4 samples). Check fragment length, alignment and IP quality; g: replicates form two sub-groups (CLP_1+CLP_2 | CLP_3+CLP_4; per-region distance contrast p = 7.9e-30, Bonferroni over 3 splits; meta-profile across/within ratio 1.39). Contrasts or null splits along this partition are not between exchangeable replicates
- CLP_1: offset 0.7 % of profile (0.22× mates), shift -0.2 bp [-0.5, +0.2], F20 0.12 (mates 0.19)
- CLP_2: offset 2.6 % of profile (1.62× mates), shift +0.3 bp [-0.2, +0.9], F20 0.29 (mates 0.11)
- CLP_3: offset 0.7 % of profile (0.25× mates), shift +0.0 bp [-0.4, +0.4], F20 0.15 (mates 0.20)
- CLP_4: offset 3.0 % of profile (2.10× mates), shift -0.2 bp [-0.6, +0.2], F20 0.17 (mates 0.18)
- split CLP_1+CLP_2 | CLP_3+CLP_4: z +11.28, p 7.9e-30, π̂₁ 0.229, called 0, shift -0.1 bp
- split CLP_1+CLP_3 | CLP_2+CLP_4: z -0.82, p 0.79, π̂₁ 0.000, called 0, shift +0.1 bp
- split CLP_1+CLP_4 | CLP_2+CLP_3: z -11.32, p 1, π̂₁ 0.000, called 0, shift +0.3 bp

**k4me3_chip_summit · Sham** (n = 4, 5000 regions): g: replicates form two sub-groups (Sham_1+Sham_2 | Sham_3+Sham_4; per-region distance contrast p = 1.3e-35, Bonferroni over 3 splits; meta-profile across/within ratio 3.55). Contrasts or null splits along this partition are not between exchangeable replicates
- Sham_1: offset 1.7 % of profile (1.25× mates), shift +0.6 bp [+0.1, +1.0], F20 0.16 (mates 0.17)
- Sham_2: offset 1.7 % of profile (1.30× mates), shift +0.1 bp [-0.3, +0.6], F20 0.15 (mates 0.16)
- Sham_3: offset 1.7 % of profile (1.41× mates), shift -0.0 bp [-0.5, +0.4], F20 0.20 (mates 0.15)
- Sham_4: offset 1.4 % of profile (1.06× mates), shift -0.7 bp [-1.1, -0.2], F20 0.16 (mates 0.16)
- split Sham_1+Sham_2 | Sham_3+Sham_4: z +12.40, p 1.3e-35, π̂₁ 0.367, called 0, shift -0.5 bp
- split Sham_1+Sham_3 | Sham_2+Sham_4: z -1.64, p 0.95, π̂₁ 0.000, called 0, shift -0.4 bp
- split Sham_1+Sham_4 | Sham_2+Sham_3: z -12.80, p 1, π̂₁ 0.000, called 0, shift +0.1 bp

**k4me3_chip_tss · CLP** (n = 4, 5000 regions): g: CLP_2 stands out from its group-mates: its sample-wide shape offset is 2.1× theirs (3.7% of the aggregate profile; p = 0, Bonferroni over 4 samples). Check fragment length, alignment and IP quality; g: replicates form two sub-groups (CLP_1+CLP_2 | CLP_3+CLP_4; per-region distance contrast p = 8.4e-13, Bonferroni over 3 splits; meta-profile across/within ratio 1.23). Contrasts or null splits along this partition are not between exchangeable replicates
- CLP_1: offset 1.4 % of profile (0.52× mates), shift -2.0 bp [-2.3, -1.7], F20 0.07 (mates 0.15)
- CLP_2: offset 3.7 % of profile (2.15× mates), shift +8.2 bp [+7.7, +8.6], F20 0.24 (mates 0.07)
- CLP_3: offset 1.4 % of profile (0.44× mates), shift -3.3 bp [-3.7, -3.0], F20 0.10 (mates 0.16)
- CLP_4: offset 2.7 % of profile (1.27× mates), shift -2.8 bp [-3.2, -2.4], F20 0.14 (mates 0.12)
- split CLP_1+CLP_2 | CLP_3+CLP_4: z +7.06, p 8.4e-13, π̂₁ 0.221, called 0, shift -4.6 bp
- split CLP_1+CLP_3 | CLP_2+CLP_4: z -0.54, p 0.71, π̂₁ 0.000, called 0, shift +4.0 bp
- split CLP_1+CLP_4 | CLP_2+CLP_3: z -4.93, p 1, π̂₁ 0.000, called 0, shift +3.6 bp

**k4me3_chip_tss · Sham** (n = 4, 5000 regions): no flag
- Sham_1: offset 1.7 % of profile (1.82× mates), shift +3.5 bp [+3.2, +3.9], F20 0.10 (mates 0.11)
- Sham_2: offset 0.8 % of profile (0.64× mates), shift +0.2 bp [-0.1, +0.6], F20 0.09 (mates 0.10)
- Sham_3: offset 1.3 % of profile (1.14× mates), shift -2.5 bp [-2.9, -2.1], F20 0.13 (mates 0.09)
- Sham_4: offset 0.8 % of profile (0.49× mates), shift -1.2 bp [-1.6, -0.9], F20 0.11 (mates 0.08)
- split Sham_1+Sham_4 | Sham_2+Sham_3: z +1.04, p 0.15, π̂₁ 0.000, called 0, shift -1.7 bp
- split Sham_1+Sham_2 | Sham_3+Sham_4: z -0.94, p 0.83, π̂₁ 0.354, called 0, shift -2.8 bp
- split Sham_1+Sham_3 | Sham_2+Sham_4: z -1.18, p 0.88, π̂₁ 0.000, called 0, shift -0.8 bp

**k27ac_cutrun · Naive** (n = 4, 4989 regions): g: Naive_2 stands out from its group-mates: its sample-wide shape offset is 3.7× theirs (1.6% of the aggregate profile; p = 5.3e-42, Bonferroni over 4 samples). Check fragment length, alignment and IP quality
- Naive_1: offset 0.8 % of profile (0.84× mates), shift +0.1 bp [-1.0, +1.2], F20 0.59 (mates 0.59)
- Naive_2: offset 1.6 % of profile (3.74× mates), shift +0.1 bp [-1.1, +1.4], F20 0.64 (mates 0.57)
- Naive_3: offset 0.7 % of profile (0.70× mates), shift -0.3 bp [-1.4, +0.7], F20 0.58 (mates 0.60)
- Naive_4: offset 0.3 % of profile (0.29× mates), shift +0.1 bp [-1.0, +1.2], F20 0.58 (mates 0.60)
- split Naive_1+Naive_4 | Naive_2+Naive_3: z +1.03, p 0.15, π̂₁ 0.000, called 0, shift -0.1 bp
- split Naive_1+Naive_2 | Naive_3+Naive_4: z -0.14, p 0.56, π̂₁ 0.000, called 0, shift -0.2 bp
- split Naive_1+Naive_3 | Naive_2+Naive_4: z -1.56, p 0.94, π̂₁ 0.119, called 0, shift +0.1 bp

**k27ac_cutrun · Arm** (n = 4, 4999 regions): no flag
- Arm_1: offset 0.3 % of profile (0.66× mates), shift +0.6 bp [-0.3, +1.5], F20 0.53 (mates 0.55)
- Arm_2: offset 0.5 % of profile (1.37× mates), shift -1.1 bp [-2.1, -0.1], F20 0.55 (mates 0.53)
- Arm_3: offset 0.4 % of profile (0.98× mates), shift +0.1 bp [-0.9, +1.0], F20 0.53 (mates 0.54)
- Arm_4: offset 0.5 % of profile (1.08× mates), shift +0.5 bp [-0.5, +1.5], F20 0.54 (mates 0.53)
- split Arm_1+Arm_4 | Arm_2+Arm_3: z +0.92, p 0.18, π̂₁ 0.000, called 0, shift -0.8 bp
- split Arm_1+Arm_2 | Arm_3+Arm_4: z +0.43, p 0.33, π̂₁ 0.000, called 0, shift +0.4 bp
- split Arm_1+Arm_3 | Arm_2+Arm_4: z -1.10, p 0.86, π̂₁ 0.000, called 0, shift -0.5 bp

**k4me3_cutrun_summit · Naive** (n = 4, 4998 regions): g: Naive_3 stands out from its group-mates: its sample-wide shape offset is 2.4× theirs (1.5% of the aggregate profile; p = 7.4e-21, Bonferroni over 4 samples). Check fragment length, alignment and IP quality
- Naive_1: offset 0.6 % of profile (0.45× mates), shift -0.4 bp [-1.2, +0.4], F20 0.46 (mates 0.48)
- Naive_2: offset 0.4 % of profile (0.30× mates), shift -0.1 bp [-0.9, +0.7], F20 0.45 (mates 0.48)
- Naive_3: offset 1.5 % of profile (2.39× mates), shift +0.4 bp [-0.6, +1.4], F20 0.55 (mates 0.45)
- Naive_4: offset 1.2 % of profile (1.35× mates), shift +0.1 bp [-0.7, +0.9], F20 0.47 (mates 0.47)
- split Naive_1+Naive_3 | Naive_2+Naive_4: z +1.74, p 0.04, π̂₁ 0.000, called 0, shift -0.0 bp
- split Naive_1+Naive_2 | Naive_3+Naive_4: z -1.66, p 0.95, π̂₁ 0.000, called 0, shift +0.4 bp
- split Naive_1+Naive_4 | Naive_2+Naive_3: z -1.74, p 0.96, π̂₁ 0.013, called 0, shift +0.2 bp

**k4me3_cutrun_summit · Arm** (n = 4, 4999 regions): g: Arm_1 stands out from its group-mates: its sample-wide shape offset is 2.8× theirs (0.8% of the aggregate profile; p = 9.3e-05, Bonferroni over 4 samples). Check fragment length, alignment and IP quality
- Arm_1: offset 0.8 % of profile (2.79× mates), shift -0.1 bp [-0.8, +0.6], F20 0.38 (mates 0.38)
- Arm_2: offset 0.5 % of profile (1.17× mates), shift +0.3 bp [-0.4, +1.0], F20 0.37 (mates 0.38)
- Arm_3: offset 0.3 % of profile (0.54× mates), shift +0.2 bp [-0.4, +0.8], F20 0.37 (mates 0.38)
- Arm_4: offset 0.3 % of profile (0.54× mates), shift -0.4 bp [-1.2, +0.4], F20 0.42 (mates 0.36)
- split Arm_1+Arm_3 | Arm_2+Arm_4: z +1.07, p 0.14, π̂₁ 0.002, called 0, shift -0.1 bp
- split Arm_1+Arm_4 | Arm_2+Arm_3: z -0.68, p 0.75, π̂₁ 0.000, called 0, shift +0.4 bp
- split Arm_1+Arm_2 | Arm_3+Arm_4: z -1.17, p 0.88, π̂₁ 0.014, called 0, shift -0.1 bp

**k4me3_cutrun_tss · Naive** (n = 4, 4994 regions): g: Naive_4 stands out from its group-mates: its sample-wide shape offset is 2.1× theirs (2.1% of the aggregate profile; p = 1.1e-125, Bonferroni over 4 samples). Check fragment length, alignment and IP quality
- Naive_1: offset 1.1 % of profile (0.66× mates), shift -0.0 bp [-1.0, +0.9], F20 0.57 (mates 0.60)
- Naive_2: offset 0.7 % of profile (0.41× mates), shift -1.4 bp [-2.4, -0.4], F20 0.58 (mates 0.60)
- Naive_3: offset 1.6 % of profile (1.18× mates), shift +7.2 bp [+6.0, +8.4], F20 0.64 (mates 0.57)
- Naive_4: offset 2.1 % of profile (2.13× mates), shift -5.7 bp [-6.8, -4.7], F20 0.59 (mates 0.60)
- split Naive_1+Naive_2 | Naive_3+Naive_4: z +1.26, p 0.1, π̂₁ 0.000, called 0, shift +1.1 bp
- split Naive_1+Naive_3 | Naive_2+Naive_4: z -1.07, p 0.86, π̂₁ 0.018, called 0, shift -5.4 bp
- split Naive_1+Naive_4 | Naive_2+Naive_3: z -1.46, p 0.93, π̂₁ 0.000, called 0, shift +4.3 bp

**k4me3_cutrun_tss · Arm** (n = 4, 4994 regions): no flag
- Arm_1: offset 0.7 % of profile (1.31× mates), shift -1.2 bp [-2.0, -0.4], F20 0.47 (mates 0.48)
- Arm_2: offset 0.5 % of profile (0.63× mates), shift -1.6 bp [-2.3, -0.9], F20 0.45 (mates 0.49)
- Arm_3: offset 0.6 % of profile (0.87× mates), shift +0.4 bp [-0.3, +1.2], F20 0.48 (mates 0.49)
- Arm_4: offset 0.9 % of profile (1.66× mates), shift +2.4 bp [+1.5, +3.3], F20 0.52 (mates 0.46)
- split Arm_1+Arm_2 | Arm_3+Arm_4: z +0.71, p 0.24, π̂₁ 0.000, called 0, shift +2.1 bp
- split Arm_1+Arm_4 | Arm_2+Arm_3: z +0.05, p 0.48, π̂₁ 0.000, called 0, shift -0.9 bp
- split Arm_1+Arm_3 | Arm_2+Arm_4: z -1.47, p 0.93, π̂₁ 0.054, called 0, shift +0.6 bp

## E · Functional consistency (real-null splits)

| dataset | group | split | contrast z | split flagged | shape_any null-call share |
|---|---|---|---|---|---|
| mnase (development) | WT | WT_1+WT_2 vs WT_3+WT_4 | +2.16 |  | 0.0000 |
| mnase (development) | WT | WT_1+WT_3 vs WT_2+WT_4 | -2.89 |  | 0.0000 |
| mnase (development) | WT | WT_1+WT_4 vs WT_2+WT_3 | +0.79 |  | 0.0000 |
| mnase (development) | dKO | dKO_1+dKO_2 vs dKO_3+dKO_4 | +10.94 | yes | 0.0554 |
| mnase (development) | dKO | dKO_1+dKO_3 vs dKO_2+dKO_4 | -2.37 |  | 0.0000 |
| mnase (development) | dKO | dKO_1+dKO_4 vs dKO_2+dKO_3 | -13.91 |  | 0.0000 |
| mnase (development) | rtf1-8-9A | rtf1-8-9A_1+rtf1-8-9A_2 vs rtf1-8-9A_3+rtf1-8-9A_4 | -1.77 |  | 0.0000 |
| mnase (development) | rtf1-8-9A | rtf1-8-9A_1+rtf1-8-9A_3 vs rtf1-8-9A_2+rtf1-8-9A_4 | +1.03 |  | 0.0000 |
| mnase (development) | rtf1-8-9A | rtf1-8-9A_1+rtf1-8-9A_4 vs rtf1-8-9A_2+rtf1-8-9A_3 | +0.03 |  | 0.0000 |
| mnase (development) | rtf1-11A | rtf1-11A_1+rtf1-11A_2 vs rtf1-11A_3+rtf1-11A_4 | +7.52 | yes | 0.0000 |
| mnase (development) | rtf1-11A | rtf1-11A_1+rtf1-11A_3 vs rtf1-11A_2+rtf1-11A_4 | -4.36 |  | 0.0000 |
| mnase (development) | rtf1-11A | rtf1-11A_1+rtf1-11A_4 vs rtf1-11A_2+rtf1-11A_3 | -2.23 |  | 0.0000 |
| mnase (development) | rtf1 | rtf1_1+rtf1_2 vs rtf1_3+rtf1_4 | +15.79 | yes | 0.0000 |
| mnase (development) | rtf1 | rtf1_1+rtf1_3 vs rtf1_2+rtf1_4 | -5.81 |  | 0.0000 |
| mnase (development) | rtf1 | rtf1_1+rtf1_4 vs rtf1_2+rtf1_3 | -12.72 |  | 0.0000 |
| mnase (development) | isw1rtf1 | isw1rtf1_1+isw1rtf1_2 vs isw1rtf1_3+isw1rtf1_4 | +20.19 | yes | 0.0880 |
| mnase (development) | isw1rtf1 | isw1rtf1_1+isw1rtf1_3 vs isw1rtf1_2+isw1rtf1_4 | -15.49 |  | 0.0000 |
| mnase (development) | isw1rtf1 | isw1rtf1_1+isw1rtf1_4 vs isw1rtf1_2+isw1rtf1_3 | -14.08 |  | 0.0000 |
| atac | sth1_0h | sth1_0h_rep1+sth1_0h_rep2 vs sth1_0h_rep3+sth1_0h_rep4 | -1.00 |  | 0.0000 |
| atac | sth1_0h | sth1_0h_rep1+sth1_0h_rep3 vs sth1_0h_rep2+sth1_0h_rep4 | +1.00 |  | 0.0000 |
| atac | sth1_0h | sth1_0h_rep1+sth1_0h_rep4 vs sth1_0h_rep2+sth1_0h_rep3 | +1.00 |  | 0.0000 |
| atac | sth1_1h | sth1_1h_rep1+sth1_1h_rep2 vs sth1_1h_rep3+sth1_1h_rep4 | -1.00 |  | 0.0000 |
| atac | sth1_1h | sth1_1h_rep1+sth1_1h_rep3 vs sth1_1h_rep2+sth1_1h_rep4 | -1.00 |  | 0.0000 |
| atac | sth1_1h | sth1_1h_rep1+sth1_1h_rep4 vs sth1_1h_rep2+sth1_1h_rep3 | +1.00 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_2+A673_3 vs A673_4+A673_5+A673_6 | +0.48 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_2+A673_4 vs A673_3+A673_5+A673_6 | -0.19 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_2+A673_5 vs A673_3+A673_4+A673_6 | +1.37 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_2+A673_6 vs A673_3+A673_4+A673_5 | +0.02 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_3+A673_4 vs A673_2+A673_5+A673_6 | +0.99 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_3+A673_5 vs A673_2+A673_4+A673_6 | +0.09 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_3+A673_6 vs A673_2+A673_4+A673_5 | -1.47 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_4+A673_5 vs A673_2+A673_3+A673_6 | -0.68 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_4+A673_6 vs A673_2+A673_3+A673_5 | -1.37 |  | 0.0000 |
| tf_cutrun | A673 | A673_1+A673_5+A673_6 vs A673_2+A673_3+A673_4 | -0.38 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_2+A2.2_3 vs A2.2_4+A2.2_5+A2.2_6 | +0.91 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_2+A2.2_4 vs A2.2_3+A2.2_5+A2.2_6 | -2.30 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_2+A2.2_5 vs A2.2_3+A2.2_4+A2.2_6 | +1.36 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_2+A2.2_6 vs A2.2_3+A2.2_4+A2.2_5 | -0.66 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_3+A2.2_4 vs A2.2_2+A2.2_5+A2.2_6 | +1.15 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_3+A2.2_5 vs A2.2_2+A2.2_4+A2.2_6 | +0.28 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_3+A2.2_6 vs A2.2_2+A2.2_4+A2.2_5 | +1.69 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_4+A2.2_5 vs A2.2_2+A2.2_3+A2.2_6 | -1.76 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_4+A2.2_6 vs A2.2_2+A2.2_3+A2.2_5 | -0.04 |  | 0.0000 |
| tf_cutrun | A2.2 | A2.2_1+A2.2_5+A2.2_6 vs A2.2_2+A2.2_3+A2.2_4 | -1.99 |  | 0.0000 |
| k27ac_chip | CLP | CLP_1+CLP_2 vs CLP_3+CLP_4 | +1.04 |  | 0.0000 |
| k27ac_chip | CLP | CLP_1+CLP_3 vs CLP_2+CLP_4 | -1.02 |  | 0.0000 |
| k27ac_chip | CLP | CLP_1+CLP_4 vs CLP_2+CLP_3 | +0.17 |  | 0.0000 |
| k27ac_chip | Sham | Sham_1+Sham_2 vs Sham_3+Sham_4 | +1.01 |  | 0.0000 |
| k27ac_chip | Sham | Sham_1+Sham_3 vs Sham_2+Sham_4 | -1.00 |  | 0.0000 |
| k27ac_chip | Sham | Sham_1+Sham_4 vs Sham_2+Sham_3 | -1.03 |  | 0.0000 |
| k4me3_chip_summit | CLP | CLP_1+CLP_2 vs CLP_3+CLP_4 | +11.28 | yes | 0.0000 |
| k4me3_chip_summit | CLP | CLP_1+CLP_3 vs CLP_2+CLP_4 | -0.82 |  | 0.0000 |
| k4me3_chip_summit | CLP | CLP_1+CLP_4 vs CLP_2+CLP_3 | -11.32 |  | 0.0000 |
| k4me3_chip_summit | Sham | Sham_1+Sham_2 vs Sham_3+Sham_4 | +12.40 | yes | 0.0000 |
| k4me3_chip_summit | Sham | Sham_1+Sham_3 vs Sham_2+Sham_4 | -1.64 |  | 0.0000 |
| k4me3_chip_summit | Sham | Sham_1+Sham_4 vs Sham_2+Sham_3 | -12.80 |  | 0.0000 |
| k4me3_chip_tss | CLP | CLP_1+CLP_2 vs CLP_3+CLP_4 | +7.06 | yes | 0.0000 |
| k4me3_chip_tss | CLP | CLP_1+CLP_3 vs CLP_2+CLP_4 | -0.54 |  | 0.0000 |
| k4me3_chip_tss | CLP | CLP_1+CLP_4 vs CLP_2+CLP_3 | -4.93 |  | 0.0000 |
| k4me3_chip_tss | Sham | Sham_1+Sham_2 vs Sham_3+Sham_4 | -0.94 |  | 0.0000 |
| k4me3_chip_tss | Sham | Sham_1+Sham_3 vs Sham_2+Sham_4 | -1.18 |  | 0.0000 |
| k4me3_chip_tss | Sham | Sham_1+Sham_4 vs Sham_2+Sham_3 | +1.04 |  | 0.0000 |
| k27ac_cutrun | Naive | Naive_1+Naive_2 vs Naive_3+Naive_4 | -0.14 |  | 0.0000 |
| k27ac_cutrun | Naive | Naive_1+Naive_3 vs Naive_2+Naive_4 | -1.56 |  | 0.0000 |
| k27ac_cutrun | Naive | Naive_1+Naive_4 vs Naive_2+Naive_3 | +1.03 |  | 0.0000 |
| k27ac_cutrun | Arm | Arm_1+Arm_2 vs Arm_3+Arm_4 | +0.43 |  | 0.0000 |
| k27ac_cutrun | Arm | Arm_1+Arm_3 vs Arm_2+Arm_4 | -1.10 |  | 0.0000 |
| k27ac_cutrun | Arm | Arm_1+Arm_4 vs Arm_2+Arm_3 | +0.92 |  | 0.0000 |
| k4me3_cutrun_summit | Naive | Naive_1+Naive_2 vs Naive_3+Naive_4 | -1.66 |  | 0.0000 |
| k4me3_cutrun_summit | Naive | Naive_1+Naive_3 vs Naive_2+Naive_4 | +1.74 |  | 0.0000 |
| k4me3_cutrun_summit | Naive | Naive_1+Naive_4 vs Naive_2+Naive_3 | -1.74 |  | 0.0000 |
| k4me3_cutrun_summit | Arm | Arm_1+Arm_2 vs Arm_3+Arm_4 | -1.17 |  | 0.0000 |
| k4me3_cutrun_summit | Arm | Arm_1+Arm_3 vs Arm_2+Arm_4 | +1.07 |  | 0.0000 |
| k4me3_cutrun_summit | Arm | Arm_1+Arm_4 vs Arm_2+Arm_3 | -0.68 |  | 0.0000 |
| k4me3_cutrun_tss | Naive | Naive_1+Naive_2 vs Naive_3+Naive_4 | +1.26 |  | 0.0000 |
| k4me3_cutrun_tss | Naive | Naive_1+Naive_3 vs Naive_2+Naive_4 | -1.07 |  | 0.0000 |
| k4me3_cutrun_tss | Naive | Naive_1+Naive_4 vs Naive_2+Naive_3 | -1.46 |  | 0.0000 |
| k4me3_cutrun_tss | Arm | Arm_1+Arm_2 vs Arm_3+Arm_4 | +0.71 |  | 0.0000 |
| k4me3_cutrun_tss | Arm | Arm_1+Arm_3 vs Arm_2+Arm_4 | -1.47 |  | 0.0000 |
| k4me3_cutrun_tss | Arm | Arm_1+Arm_4 vs Arm_2+Arm_3 | +0.05 |  | 0.0000 |

New datasets: 62 splits; Spearman(contrast z, null-call share) = +nan, one-sided p = nan; flagged splits 3 (mean share 0.0000) vs unflagged 59 (mean 0.0000).

## Criteria

| id | result |
|---|---|
| K1 | pass |
| K2 | pass |
| K3 | pass |
| K4 | **FAIL** |
| K5 | pass |
| K6 | pass |
| K7 | **FAIL** |
| K8 | pass |
| K9 | pass |
| K10 | **FAIL** |
