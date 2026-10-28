# Prerequisite R · midpoint vs coverage under real nulls (same reads)

Protocol `docs/PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md` sha256 `2ad46875beab9620` (amendment R1); commit `de4059bd89` (uncommitted 0); regional.py `21ffd82395ee207f`; script `8fa2517c9af7764d`; 2026-10-28T09:50:32 → 2026-10-28T09:51:48. Regions: {'cutrun': 801, 'yeast': 3592}.

Per split: whole-list calls at FDR 0.05 / share p < 0.01, for midpoint | coverage.

| dataset | split | shape_any mid | shape_kernel mid | shift mid | shape_any cov | shape_kernel cov | shift cov |
|---|---|---|---|---|---|---|---|
| atac | minus1|0h|sth1_0h_rep1+sth1_0h_rep2 | 0 / 0.0098 | 0 / 0.0117 | 0 / 0.0106 | 0 / 0.0068 | 0 / 0.0073 | 0 / 0.0099 |
| atac | minus1|0h|sth1_0h_rep1+sth1_0h_rep3 | 0 / 0.0073 | 0 / 0.0018 | 0 / 0.0101 | 0 / 0.0048 | 0 / 0.0038 | 0 / 0.0037 |
| atac | minus1|0h|sth1_0h_rep1+sth1_0h_rep4 | 0 / 0.0162 | 0 / 0.0190 | 0 / 0.0081 | 0 / 0.0242 | 0 / 0.0191 | 0 / 0.0183 |
| atac | minus1|1h|sth1_1h_rep1+sth1_1h_rep2 | 0 / 0.0093 | 0 / 0.0047 | 0 / 0.0093 | 0 / 0.0099 | 0 / 0.0048 | 0 / 0.0102 |
| atac | minus1|1h|sth1_1h_rep1+sth1_1h_rep3 | 0 / 0.0029 | 0 / 0.0019 | 0 / 0.0093 | 0 / 0.0044 | 0 / 0.0041 | 0 / 0.0053 |
| atac | minus1|1h|sth1_1h_rep1+sth1_1h_rep4 | 0 / 0.0296 | 0 / 0.0449 | 0 / 0.0119 | 0 / 0.0193 | 0 / 0.0266 | 0 / 0.0134 |
| atac | ndr|0h|sth1_0h_rep1+sth1_0h_rep2 | 0 / 0.0269 | 0 / 0.0350 | 0 / 0.0127 | 0 / 0.0353 | 0 / 0.0485 | 0 / 0.0208 |
| atac | ndr|0h|sth1_0h_rep1+sth1_0h_rep3 | 0 / 0.0052 | 0 / 0.0024 | 0 / 0.0070 | 0 / 0.0042 | 0 / 0.0018 | 0 / 0.0059 |
| atac | ndr|0h|sth1_0h_rep1+sth1_0h_rep4 | 0 / 0.0083 | 0 / 0.0040 | 0 / 0.0103 | 0 / 0.0039 | 0 / 0.0034 | 0 / 0.0068 |
| atac | ndr|1h|sth1_1h_rep1+sth1_1h_rep2 | 0 / 0.0127 | 0 / 0.0175 | 0 / 0.0075 | 0 / 0.0148 | 0 / 0.0172 | 0 / 0.0114 |
| atac | ndr|1h|sth1_1h_rep1+sth1_1h_rep3 | 0 / 0.0034 | 0 / 0.0011 | 0 / 0.0068 | 0 / 0.0046 | 0 / 0.0021 | 0 / 0.0062 |
| atac | ndr|1h|sth1_1h_rep1+sth1_1h_rep4 | 0 / 0.0159 | 0 / 0.0161 | 0 / 0.0143 | 0 / 0.0128 | 0 / 0.0144 | 0 / 0.0119 |
| atac | plus1|0h|sth1_0h_rep1+sth1_0h_rep2 | 0 / 0.0069 | 0 / 0.0085 | 0 / 0.0114 | 0 / 0.0101 | 0 / 0.0075 | 0 / 0.0112 |
| atac | plus1|0h|sth1_0h_rep1+sth1_0h_rep3 | 0 / 0.0074 | 0 / 0.0031 | 0 / 0.0096 | 0 / 0.0035 | 0 / 0.0040 | 0 / 0.0058 |
| atac | plus1|0h|sth1_0h_rep1+sth1_0h_rep4 | 0 / 0.0168 | 0 / 0.0227 | 0 / 0.0090 | 0 / 0.0178 | 0 / 0.0218 | 0 / 0.0191 |
| atac | plus1|1h|sth1_1h_rep1+sth1_1h_rep2 | 0 / 0.0064 | 0 / 0.0074 | 0 / 0.0101 | 0 / 0.0076 | 0 / 0.0072 | 0 / 0.0073 |
| atac | plus1|1h|sth1_1h_rep1+sth1_1h_rep3 | 0 / 0.0034 | 0 / 0.0006 | 0 / 0.0048 | 0 / 0.0025 | 0 / 0.0014 | 0 / 0.0054 |
| atac | plus1|1h|sth1_1h_rep1+sth1_1h_rep4 | 0 / 0.0227 | 0 / 0.0302 | 0 / 0.0154 | 0 / 0.0259 | 0 / 0.0271 | 0 / 0.0183 |
| cutrun | A673_1+A673_2+A673_3 | 0 / 0.0064 | 0 / 0.0037 | 0 / 0.0127 | 0 / 0.0101 | 0 / 0.0110 | 0 / 0.0127 |
| cutrun | A673_1+A673_2+A673_4 | 0 / 0.0076 | 0 / 0.0149 | 0 / 0.0064 | 0 / 0.0051 | 0 / 0.0110 | 0 / 0.0038 |
| cutrun | A673_1+A673_2+A673_5 | 0 / 0.0102 | 0 / 0.0056 | 0 / 0.0115 | 0 / 0.0190 | 0 / 0.0155 | 0 / 0.0165 |
| cutrun | A673_1+A673_2+A673_6 | 0 / 0.0255 | 0 / 0.0279 | 0 / 0.0178 | 0 / 0.0025 | 0 / 0.0022 | 0 / 0.0076 |
| cutrun | A673_1+A673_3+A673_4 | 0 / 0.0102 | 0 / 0.0112 | 0 / 0.0127 | 0 / 0.0165 | 0 / 0.0132 | 0 / 0.0190 |
| cutrun | A673_1+A673_3+A673_5 | 0 / 0.0102 | 0 / 0.0130 | 0 / 0.0076 | 0 / 0.0063 | 0 / 0.0044 | 0 / 0.0063 |
| cutrun | A673_1+A673_3+A673_6 | 0 / 0.0038 | 0 / 0.0056 | 0 / 0.0076 | 0 / 0.0089 | 0 / 0.0110 | 0 / 0.0063 |
| cutrun | A673_1+A673_4+A673_5 | 0 / 0.0025 | 0 / 0.0056 | 0 / 0.0038 | 0 / 0.0051 | 0 / 0.0088 | 0 / 0.0101 |
| cutrun | A673_1+A673_4+A673_6 | 0 / 0.0038 | 0 / 0.0019 | 0 / 0.0064 | 0 / 0.0013 | 0 / 0.0000 | 0 / 0.0025 |
| cutrun | A673_1+A673_5+A673_6 | 0 / 0.0089 | 0 / 0.0093 | 0 / 0.0076 | 0 / 0.0139 | 0 / 0.0199 | 0 / 0.0089 |
| yeast | WT_rep1+WT_rep2 | 0 / 0.0103 | 0 / 0.0077 | 0 / 0.0077 | 0 / 0.0092 | 0 / 0.0046 | 0 / 0.0105 |
| yeast | WT_rep1+WT_rep3 | 0 / 0.0135 | 0 / 0.0148 | 0 / 0.0115 | 0 / 0.0131 | 0 / 0.0184 | 0 / 0.0098 |
| yeast | WT_rep1+WT_rep4 | 0 / 0.0058 | 0 / 0.0051 | 0 / 0.0090 | 0 / 0.0098 | 0 / 0.0072 | 0 / 0.0085 |

## Verdict

| dataset | test | coverage more calls (splits) | Δ p<0.01 share (cov − mid) | other | coverage higher |
|---|---|---|---|---|---|
| atac | shape_any | 0 / 18 | +0.0001 |  | no |
| atac | shape_kernel | 0 / 18 | -0.0006 |  | no |
| atac | shift | 0 / 18 | +0.0007 |  | no |
| cutrun | shape_any | 0 / 10 | -0.0000 | zero-call splits 10/10 | no |
| cutrun | shape_kernel | 0 / 10 | -0.0001 | zero-call splits 10/10 | no |
| cutrun | shift | 0 / 10 | -0.0000 | zero-call splits 10/10 | no |
| yeast | shape_any | 0 / 3 | +0.0009 |  | no |
| yeast | shape_kernel | 0 / 3 | +0.0009 |  | no |
| yeast | shift | 0 / 3 | +0.0002 |  | no |

**Coverage eligible for stage 2: yes** (not higher than midpoints on any dataset and test). Yeast midpoints are known to exceed nominal at 2 vs 2; this comparison is relative. CUT&RUN uses a 4 M-pair, chr1 + chr2 subsample (R1).
