# Independent positive case · GSE218269 (CDK7 inhibition, Pol II MNase-ChIP-seq)

**EXPLORATORY (amendment B4): S3 (depth) was not met; the confirmatory case is "not evaluable". This run is not a confirmation and is not counted as one.**

Protocol `docs/PROFILESHAPE_POSITIVE_CASE_PROTOCOL.md` (stage 2) sha256 `244049c1b4a7058f`; commit `c5defbcb7b`, uncommitted src files 0 (diff `e3b0c44298fc1c14`); regional.py `21ffd82395ee207f`; script `fc5b15fe471ae7a9`; 2026-10-28T07:13:19 → 2026-10-28T07:16:22.

## A0 · Depth and consistency

```
{
 "median_fragments_per_sample": {
  "DMSO_1": 86.0,
  "DMSO_2": 95.0,
  "DMSO_3": 84.0,
  "CDK7i_1": 77.0,
  "CDK7i_2": 82.0,
  "CDK7i_3": 93.0
 },
 "median_over_samples": 85.0,
 "S3_met": false,
 "source": "BAM, exact fragments with a midpoint in the core window (amendment B2)"
}
```

Consistency flags: DMSO: DMSO_1 stands out from its group-mates: its sample-wide shape offset is 4.7× theirs (4.2% of the aggregate profile; p = 4.9e-121, Bonferroni over 3 samples). Check fragment length, alignment and IP quality; CDK7i: CDK7i_1 stands out from its group-mates: its sample-wide shape offset is 2.7× theirs (3.3% of the aggregate profile; p = 8.3e-78, Bonferroni over 3 samples). Check fragment length, alignment and IP quality

## A1 · Expected direction

| prediction | result | met |
|---|---|---|
| P1 `shape_any` share ≥ 0.20 | 0.400 of 7865 | yes |
| P2 upstream share among `shift`-significant ≥ 0.80 | nan of 0 | **no** |
| P3 aggregate centroid upstream (Welch 95 % < 0) | +5.6 bp [-1.6, +12.9] | **no** |

Consistent translations among `shift`-significant regions: 0 (median +nan bp).

## A2 · Added explanatory value for TT-seq (primary at 30 min)

| readout | regions | CV R² amount-only | CV R² + shape | gain [95 % CI] | std. coef. signed shift | std. coef. kernel |
|---|---|---|---|---|---|---|
| TT-seq 30min | 7863 | 0.034 | 0.039 | +0.0056 [+0.0025, +0.0085] | -0.027 | +0.003 |
| TT-seq 15min | 7863 | 0.083 | 0.084 | +0.0005 [-0.0007, +0.0016] | -0.010 | +0.008 |
| TT-seq 60min | 7863 | 0.051 | 0.056 | +0.0048 [+0.0021, +0.0074] | -0.022 | +0.007 |

Primary criterion (30 min, lower bound of gain > 0): **met**.
