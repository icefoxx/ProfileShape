# Independent validation on GSE346229 (results)

Commit `c3e17ec053`; uncommitted src files 0; regional.py `87f0c6c86498c6a5`; protocol sha256 `a5234324972ab7af`; anchors {'plus1': 6089, 'minus1': 6030, 'ndr': 6030}; 2026-10-27T13:31:42 → 2026-10-27T13:34:34.

## A. Real-null splits (share of regions called, BH 5 %)

| anchor | fragments | condition | split | regions | shape_any | shape_kernel | shift | diffSplice |
|---|---|---|---|---|---|---|---|---|
| plus1 | nuc | 0h | sth1_0h_rep1+sth1_0h_rep2 vs sth1_0h_rep3+sth1_0h_rep4 | 6089 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| plus1 | nuc | 0h | sth1_0h_rep1+sth1_0h_rep3 vs sth1_0h_rep2+sth1_0h_rep4 | 6089 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| plus1 | nuc | 0h | sth1_0h_rep1+sth1_0h_rep4 vs sth1_0h_rep2+sth1_0h_rep3 | 6089 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| plus1 | nuc | 1h | sth1_1h_rep1+sth1_1h_rep2 vs sth1_1h_rep3+sth1_1h_rep4 | 6089 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| plus1 | nuc | 1h | sth1_1h_rep1+sth1_1h_rep3 vs sth1_1h_rep2+sth1_1h_rep4 | 6089 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| plus1 | nuc | 1h | sth1_1h_rep1+sth1_1h_rep4 vs sth1_1h_rep2+sth1_1h_rep3 | 6089 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| minus1 | nuc | 0h | sth1_0h_rep1+sth1_0h_rep2 vs sth1_0h_rep3+sth1_0h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| minus1 | nuc | 0h | sth1_0h_rep1+sth1_0h_rep3 vs sth1_0h_rep2+sth1_0h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| minus1 | nuc | 0h | sth1_0h_rep1+sth1_0h_rep4 vs sth1_0h_rep2+sth1_0h_rep3 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| minus1 | nuc | 1h | sth1_1h_rep1+sth1_1h_rep2 vs sth1_1h_rep3+sth1_1h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| minus1 | nuc | 1h | sth1_1h_rep1+sth1_1h_rep3 vs sth1_1h_rep2+sth1_1h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| minus1 | nuc | 1h | sth1_1h_rep1+sth1_1h_rep4 vs sth1_1h_rep2+sth1_1h_rep3 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0002 |
| ndr | sub | 0h | sth1_0h_rep1+sth1_0h_rep2 vs sth1_0h_rep3+sth1_0h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0007 |
| ndr | sub | 0h | sth1_0h_rep1+sth1_0h_rep3 vs sth1_0h_rep2+sth1_0h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| ndr | sub | 0h | sth1_0h_rep1+sth1_0h_rep4 vs sth1_0h_rep2+sth1_0h_rep3 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| ndr | sub | 1h | sth1_1h_rep1+sth1_1h_rep2 vs sth1_1h_rep3+sth1_1h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| ndr | sub | 1h | sth1_1h_rep1+sth1_1h_rep3 vs sth1_1h_rep2+sth1_1h_rep4 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| ndr | sub | 1h | sth1_1h_rep1+sth1_1h_rep4 vs sth1_1h_rep2+sth1_1h_rep3 | 6030 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

| anchor | A1: mean shape_any share ≤ 0.005 | A2: max split ≤ 0.01 |
|---|---|---|
| plus1 | 0.0000 → pass | 0.0000 → pass |
| minus1 | 0.0000 → pass | 0.0000 → pass |
| ndr | 0.0000 → pass | 0.0000 → pass |

## B. Fragment-level insertion (20 % of regions; 6 splits × 2 sets)

| anchor | insertion | shape_any power / FDP | shape_kernel power / FDP | shift power / FDP | diffSplice power / FDP | translation coverage · direction · estimate ÷ true | criterion |
|---|---|---|---|---|---|---|---|
| plus1 | shift20 | 0.04 / 0.050 | 0.00 / 0.000 | 0.08 / 0.041 | 0.00 / 0.032 | 0.08 · 1.00 · 1.06 | B1 pass |
| plus1 | shift40 | 0.10 / 0.055 | 0.00 / 0.035 | 0.13 / 0.037 | 0.07 / 0.014 | 0.12 · 1.00 · 0.96 | B1 pass |
| plus1 | broaden | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | – | B1 pass |
| plus1 | amplitude | called 0.0000 | called 0.0000 | called 0.0000 | called 0.0000 | – | B2 pass |
| minus1 | shift20 | 0.04 / 0.039 | 0.00 / 0.000 | 0.15 / 0.043 | 0.01 / 0.029 | 0.14 · 1.00 · 1.04 | B1 pass |
| minus1 | shift40 | 0.12 / 0.048 | 0.00 / 0.030 | 0.15 / 0.044 | 0.09 / 0.017 | 0.15 · 1.00 · 0.97 | B1 pass |
| minus1 | broaden | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.000 | 0.00 / 0.167 | – | B1 pass |
| minus1 | amplitude | called 0.0000 | called 0.0000 | called 0.0000 | called 0.0000 | – | B2 pass |

## C. Biology

| comparison | anchor | fragments | regions | shape_any | shape_kernel | shift | diffSplice | overlap (Jaccard) | consistent upstream / downstream | median consistent shift (bp) |
|---|---|---|---|---|---|---|---|---|---|---|
| Sth1 4v4 | plus1 | nuc | 6089 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 | 0 / 0 | nan |
| Sth1 4v4 | minus1 | nuc | 6030 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 | 0 / 0 | nan |
| Sth1 4v4 | ndr | sub | 6030 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 | 0 / 0 | nan |
| Snf2 3v3 | plus1 | nuc | 6089 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 | 0 / 0 | nan |
| Snf2 3v3 | minus1 | nuc | 6030 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 | 0 / 0 | nan |
| Snf2 3v3 | ndr | sub | 6030 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 | 0 / 0 | nan |

| criterion | result |
|---|---|
| C1: +1 nucleosomes, consistent calls mostly upstream (share > 0.5, binomial p < 0.05) | share nan, p = nan → FAIL |
| C2: −1 nucleosomes, consistent calls mostly downstream | share nan, p = nan → FAIL |
| C3: Sth1 (mean over 3 vs 3 subsets) has more shape_any calls than Snf2 3 vs 3 | +1: 0.000 vs 0.000; −1: 0.000 vs 0.000 → FAIL |

Notes from the Sth1 4 vs 4 runs: plus1: shape tests: 3057 testable regions; null in 10 spread strata of ~10370 values each; low depth: median 18 fragments per region and sample (BAM (exact)); per-region shape tests have little power below ~50 — a non-significant region is not evidence of no change / minus1: shape tests: 3015 testable regions; null in 10 spread strata of ~10234 values each; low depth: median 24 fragments per region and sample (BAM (exact)); per-region shape tests have little power below ~50 — a non-significant region is not evidence of no change / ndr: shape tests: 4368 testable regions; null in 10 spread strata of ~14824 values each; low depth: median 44 fragments per region and sample (BAM (exact)); per-region shape tests have little power below ~50 — a non-significant region is not evidence of no change

