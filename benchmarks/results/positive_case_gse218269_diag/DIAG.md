# GSE218269 · diagnostic supplement (amendment B5) — EXPLORATORY

S3 was not met; the confirmatory case is "not evaluable". Nothing here updates that verdict. Commit `271ae31d10` (uncommitted files 0); regional.py `21ffd82395ee207f`; protocol `441d659d19b0a184`; script `79d94bf2beca096f`; 2026-10-28T09:21:08 → 2026-10-28T09:23:19. B4 A1 reproduced: True.

## D1 · A2 ablation (same rows, folds, bootstrap draws)

| readout | regions | R² base | + shift | + kernel | + both | gain shift | gain kernel | gain both | both − shift | both − kernel |
|---|---|---|---|---|---|---|---|---|---|---|
| 30min | 7863 | 0.0336 | 0.0393 | 0.0336 | 0.0392 | +0.0057 [+0.0029, +0.0085] | +0.0000 [-0.0008, +0.0006] | +0.0056 [+0.0025, +0.0085] | -0.0001 [-0.0009, +0.0003] | +0.0056 [+0.0029, +0.0084] |
| 15min | 7863 | 0.0830 | 0.0833 | 0.0833 | 0.0835 | +0.0003 [-0.0004, +0.0012] | +0.0002 [-0.0006, +0.0012] | +0.0005 [-0.0007, +0.0016] | +0.0002 [-0.0006, +0.0011] | +0.0003 [-0.0005, +0.0011] |
| 60min | 7863 | 0.0514 | 0.0559 | 0.0518 | 0.0562 | +0.0046 [+0.0020, +0.0075] | +0.0005 [-0.0005, +0.0015] | +0.0048 [+0.0021, +0.0074] | +0.0003 [-0.0006, +0.0012] | +0.0044 [+0.0017, +0.0072] |

## D2 · Per-region width (MAD of the core profile, CDK7i − DMSO)

| regions | n | diff quantiles 5/25/50/75/95 % (bp) | share narrower | all 3 CDK7i < all 3 DMSO | all 3 > | expected under exchangeability |
|---|---|---|---|---|---|---|
| `shape_any` called | 3078 | -65.1 / -23.2 / +1.0 / +25.1 / +77.9 | 0.488 | 0.119 | 0.117 | 0.050 each |
| `shape_any` not called | 4445 | -65.4 / -21.9 / +0.3 / +23.7 / +73.6 | 0.497 | 0.058 | 0.056 | 0.050 each |

`mad_width` test: 7513 tested, 0 FDR calls, share p < 0.01 0.024. Aggregate IQR: DMSO 271 bp, difference -30.3 bp → equivalent f ≈ 0.89.

Detectability (500 random shape-tested regions narrowed by f in the three CDK7i samples; share called at FDR 0.05):

| f | draw | `mad_width` without → with | `shape_kernel` without → with |
|---|---|---|---|
| 0.95 | 0 | 0.000 → 0.000 | 0.492 → 0.590 |
| 0.95 | 1 | 0.000 → 0.000 | 0.458 → 0.564 |
| 0.95 | 2 | 0.000 → 0.000 | 0.498 → 0.586 |
| 0.9 | 0 | 0.000 → 0.000 | 0.466 → 0.618 |
| 0.9 | 1 | 0.000 → 0.000 | 0.470 → 0.620 |
| 0.9 | 2 | 0.000 → 0.000 | 0.482 → 0.626 |
| 0.8 | 0 | 0.000 → 0.003 | 0.452 → 0.648 |
| 0.8 | 1 | 0.000 → 0.000 | 0.508 → 0.684 |
| 0.8 | 2 | 0.000 → 0.000 | 0.500 → 0.663 |

Reported `shift` estimate: translation_consistent: cross-correlation and projection agree in sign, translation R² ≥ 0.5, ambiguity < 0.9, not at the search edge (regional.py:848–858); 173 of 7865 regions; quantiles 5/25/50/75/95 % -86.4 / -72.3 / -45.7 / -16.7 / +48.4 bp; flags {'': 173, 'direction_mismatch': 622, 'poor_translation_fit': 7062, 'search_edge': 8}.

All regions, cross-correlation shift: observed share negative 0.614; one-sided p = 0.1 (floor 0.1, 10 assignments). Rule: share of regions with a negative (upstream) cross-correlation shift, CDK7i vs DMSO; one-sided p = rank of the observed share among the 10 assignments (ties counted against).

## D3 · Processing and replicate sensitivity

| variant | depth (median) | `shape_any` share | `shift` sig. | upstream share | centroid bp [95 %] | A2 30 min gain [95 %] |
|---|---|---|---|---|---|---|
| 20-70 midpoint | 50 | 0.276 | 0 | nan | +4.3 [-3.4, +12.0] | +0.0022 [+0.0002, +0.0039] |
| 20-70 coverage | – | 0.311 | 382 | 0.79 | +4.3 [-3.5, +12.0] | +0.0045 [+0.0006, +0.0080] |
| 20-120 midpoint | 85 | 0.400 | 0 | nan | +5.6 [-1.6, +12.9] | +0.0056 [+0.0025, +0.0085] |
| 20-120 coverage | – | 0.343 | 638 | 0.79 | +6.7 [-1.2, +14.6] | +0.0085 [+0.0035, +0.0129] |
| 70-120 midpoint | 35 | 0.000 | 0 | nan | +8.7 [-2.7, +20.1] | +0.0025 [+0.0004, +0.0042] |
| 70-120 coverage | – | 0.000 | 0 | nan | +9.6 [-3.8, +22.9] | +0.0061 [+0.0018, +0.0101] |
| drop DMSO_1 | – | 0.143 | 0 | nan | +8.0 [+1.0, +15.1] | +0.0038 [+0.0007, +0.0062] |
| drop DMSO_2 | – | 0.157 | 0 | nan | +5.0 [-13.1, +23.1] | +0.0046 [+0.0017, +0.0073] |
| drop DMSO_3 | – | 0.173 | 0 | nan | +4.9 [-12.9, +22.7] | +0.0047 [+0.0017, +0.0069] |
| drop CDK7i_1 | – | 0.165 | 0 | nan | +7.2 [-1.0, +15.4] | +0.0061 [+0.0034, +0.0090] |
| drop CDK7i_2 | – | 0.177 | 0 | nan | +4.9 [-6.8, +16.5] | +0.0053 [+0.0014, +0.0083] |
| drop CDK7i_3 | – | 0.208 | 0 | nan | +4.8 [-7.6, +17.2] | +0.0049 [+0.0019, +0.0075] |
