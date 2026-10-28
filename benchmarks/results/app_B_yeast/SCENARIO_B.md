# Scenario B · complex shape change at nucleosome arrays (yeast MNase, development data)

Plan sha256 `610bdd4e5922cc64`; commit `490b14a6ee` (uncommitted 1); regional.py `21ffd82395ee207f`; seed 9910000; 2026-10-28T13:03:19 → 2026-10-28T13:03:34.

(i) WT (4) vs isw1Δ chd1Δ (4): 1433 / 1555 genes `shape_any`; `shift` 553 called, 0 with a direction; reasons for the rest {'direction_mismatch': 2, 'poor_translation_fit': 551}.

Per-gene phasing change (dKO − WT; autocorrelation at 150–190 bp), difference of medians called − not called [95 % bootstrap]:

- (ii) main, 4 vs 4: -0.122 [-0.142, -0.094] (called -0.242, n = 1433; not -0.120, n = 122)
- (iii) held-out: called on replicates 1–2 (1321 / 1476), phasing on 3–4: -0.103 [-0.143, -0.081] (called -0.271, n = 1320; not -0.169, n = 155)
- diffSpliceDGE (2549 calls, 1086 shared with `shape_any`): -0.041 [-0.059, -0.025] (called -0.247, n = 1138; not -0.206, n = 417)
- edgeR window counts, amount (82 calls; `shape_any` calls without an amount call: 1400): +0.045 [+0.010, +0.099] (called -0.191, n = 34; not -0.236, n = 1521)

Median per-sample phasing: WT_rep1 0.34, WT_rep2 0.33, WT_rep3 0.34, WT_rep4 0.36, dKO_rep1 0.14, dKO_rep2 0.12, dKO_rep3 0.10, dKO_rep4 0.09.

**Grade: consistent** (criteria in the plan). Development data; counts for the comparators are fragment equivalents from CPM coverage tracks (approximation); held-out calls use 2 vs 2, where the yeast null is known to exceed nominal.
