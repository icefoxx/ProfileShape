# Independent positive case · GSE218269 (CDK7 inhibition, Pol II MNase-ChIP-seq)

Protocol `docs/PROFILESHAPE_POSITIVE_CASE_PROTOCOL.md` (stage 2) sha256 `dbdc8ff0afcac14e`; commit `22a03cfd84`, uncommitted src files 0 (diff `e3b0c44298fc1c14`); regional.py `21ffd82395ee207f`; script `9f1546decb6470bf`; 2026-10-28T07:12:28 → 2026-10-28T07:12:30.

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

**not evaluable: S3 (depth) not met**
