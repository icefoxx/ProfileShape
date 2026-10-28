# Independent positive case · GSE218269 (CDK7 inhibition, Pol II MNase-ChIP-seq)

Protocol `docs/PROFILESHAPE_POSITIVE_CASE_PROTOCOL.md` (stage 2) sha256 `4cb2702cdbcfd0ad`; commit `960677c694`, uncommitted src files 4 (diff `5204d9da9eeb04f5`); regional.py `21ffd82395ee207f`; script `6a4fc471828ec9a8`; 2026-10-28T00:13:11 → 2026-10-28T00:13:15.

## A0 · Depth and consistency

```
{
 "median_fragments_per_sample": {
  "DMSO_1": 59.9,
  "DMSO_2": 65.0,
  "DMSO_3": 55.9,
  "CDK7i_1": 51.0,
  "CDK7i_2": 54.3,
  "CDK7i_3": 60.7
 },
 "median_over_samples": 57.918614861607004,
 "S3_met": false,
 "source": "bigWig count quantum (estimate)"
}
```

Consistency flags: DMSO: DMSO_1 stands out from its group-mates: its sample-wide shape offset is 2.3× theirs (3.3% of the aggregate profile; p = 1e-59, Bonferroni over 3 samples). Check fragment length, alignment and IP quality

**not evaluable: S3 (depth) not met**
