# Scenario C · changes that window counts and aggregate profiles miss (GSE218269) — EXPLORATORY

S3 (depth) not met; the confirmatory case is not evaluable and stays so. Plan sha256 `3e077e717e8304af`; commit `cddb56ca98` (uncommitted code files 0); regional.py `21ffd82395ee207f`; 2026-10-28T13:04:57 → 2026-10-28T13:05:17.

(i) Aggregate core profile per replicate: DMSO_1 centroid -11.5 bp, MAD 200.6, IQR 347; DMSO_2 centroid -15.1 bp, MAD 198.0, IQR 338; DMSO_3 centroid -15.5 bp, MAD 198.6, IQR 340; CDK7i_1 centroid -15.2 bp, MAD 187.2, IQR 313; CDK7i_2 centroid -11.7 bp, MAD 190.2, IQR 322; CDK7i_3 centroid -11.9 bp, MAD 189.8, IQR 322

(ii) 3 vs 3: `shape_any` 3148 / 7865; edgeR amount 2201; shape calls without an amount call 2177 (0.69); amount without shape 1230; diffSpliceDGE 1644 (1294 shared with `shape_any`; 1106 without an amount call).

(iii) Held-out sign agreement (called on replicates 1–2: 2336 / 7861; checked on replicate 3), share [95 % exact]:

- width, called: negative 1138/1617 = 0.704 [0.681, 0.726]; positive 386/719 = 0.537 [0.500, 0.574]
- centroid, called: negative 672/1166 = 0.576 [0.547, 0.605]; positive 788/1170 = 0.674 [0.646, 0.700]
- width, not called (reference): negative 2177/3370 = 0.646 [0.630, 0.662]; positive 1179/2155 = 0.547 [0.526, 0.568]
- centroid, not called (reference): negative 1445/2740 = 0.527 [0.508, 0.546]; positive 1666/2785 = 0.598 [0.580, 0.616]

**Grade: partly consistent** (criteria in the plan, clarification for C). Exploratory; below the planned depth; no mechanism claimed.
