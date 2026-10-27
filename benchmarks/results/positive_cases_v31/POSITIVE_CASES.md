# Real positive cases with ProfileSeg v3.1 (development datasets; biological consistency, not validation)

Commit `c3e17ec053`; uncommitted src files 0; regional.py `87f0c6c86498c6a5`; 2026-10-27T15:45:06 → 2026-10-27T15:45:31.

## Yeast MNase-seq (GSE280671), ATG-anchored ORFs

| comparison | shape_any | shape_kernel | shift (called) | consistent: downstream / upstream | Changed (no direction) | median consistent shift, bp [IQR] |
|---|---|---|---|---|---|---|
| WT(4) vs isw1Δ chd1Δ(4) | 1433 / 1555 | 1453 | 553 | 0 / 0 | 553 | nan [nan, nan] |
| WT(4) vs chd1Δ(2) | 905 / 1383 | 798 | 492 | 0 / 0 | 492 | nan [nan, nan] |
| WT(4) vs isw1Δ(2) | 1043 / 1506 | 971 | 442 | 0 / 0 | 442 | nan [nan, nan] |
| null WT WT_rep1+WT_rep2 vs rest | 0 / 1490 | 0 | 0 | 0 / 0 | 0 | nan [nan, nan] |
| null WT WT_rep1+WT_rep3 vs rest | 0 / 1490 | 0 | 0 | 0 / 0 | 0 | nan [nan, nan] |
| null WT WT_rep1+WT_rep4 vs rest | 0 / 1490 | 0 | 0 | 0 / 0 | 0 | nan [nan, nan] |

Phasing (meta-profile autocorrelation at the nucleosome repeat; not used by ProfileSeg): WT [0.81, 0.82, 0.82, 0.82]; dKO [0.36, 0.17, 0.14, 0.37]; chd1 [0.67, 0.68]; isw1 [0.68, 0.67]

## PRO-seq (GSE233930), NELFB degradation vs DMSO, 2 vs 2, sense strand at TSSs

| comparison | shape_any | shape_kernel | shift (called) | consistent: downstream / upstream (share downstream) | Changed (no direction) | median consistent shift, bp [IQR] |
|---|---|---|---|---|---|---|
| DMSO vs VHL2h | 10867 / 13797 | 10605 | 9117 | 3455 / 22 (0.99) | 5640 | 16.0 [11.2, 22.2] |
| DMSO vs VHL48h | 8667 / 13923 | 8652 | 4816 | 1885 / 7 (1.00) | 2924 | 13.9 [10.9, 19.2] |

Meta-profile over the top 2,000 promoters (descriptive; not used by ProfileSeg): DMSO_rep1 peak 35 bp, centroid 63.6 bp; DMSO_rep2 peak 35 bp, centroid 57.9 bp; VHL2h_rep1 peak 45 bp, centroid 81.0 bp; VHL2h_rep2 peak 45 bp, centroid 80.4 bp; VHL48h_rep1 peak 45 bp, centroid 80.7 bp; VHL48h_rep2 peak 45 bp, centroid 74.4 bp

Notes: WT(4) vs isw1Δ chd1Δ(4): shape tests: 1555 testable regions; null in 7 spread strata of ~7548 values each; shift: 553 of 553 significant regions reported as Changed (see translation_flag: direction_mismatch, poor_translation_fit, ambiguous_match, search_edge); translation_consistent is a fit-consistency flag, not proof of a translation | WT(4) vs chd1Δ(2): shape tests: 1383 testable regions; null in 6 spread strata of ~3220 values each; shift: 492 of 492 significant regions reported as Changed (see translation_flag: direction_mismatch, poor_translation_fit, ambiguous_match, search_edge); translation_consistent is a fit-consistency flag, not proof of a translation | WT(4) vs isw1Δ(2): shape tests: 1505 testable regions; null in 7 spread strata of ~3010 values each; shift: 442 of 442 significant regions reported as Changed (see translation_flag: direction_mismatch, poor_translation_fit, ambiguous_match, search_edge); translation_consistent is a fit-consistency flag, not proof of a translation | null WT WT_rep1+WT_rep2 vs rest: shape tests: 1490 testable regions; null in 1 spread strata of ~2980 values each (too few regions to stratify: per-class calibration as in v2.4; see docs/PROFILESEG_CALIBRATION.md) | null WT WT_rep1+WT_rep3 vs rest: shape tests: 1490 testable regions; null in 1 spread strata of ~2980 values each (too few regions to stratify: per-class calibration as in v2.4; see docs/PROFILESEG_CALIBRATION.md) | null WT WT_rep1+WT_rep4 vs rest: shape tests: 1490 testable regions; null in 1 spread strata of ~2980 values each (too few regions to stratify: per-class calibration as in v2.4; see docs/PROFILESEG_CALIBRATION.md) | VHL2h: shape tests: 12289 testable regions; null in 10 spread strata of ~2456 values each; shift: 5640 of 9117 significant regions reported as Changed (see translation_flag: direction_mismatch, poor_translation_fit, ambiguous_match, search_edge); translation_consistent is a fit-consistency flag, not proof of a translation | VHL48h: shape tests: 12320 testable regions; null in 10 spread strata of ~2464 values each; shift: 2924 of 4816 significant regions reported as Changed (see translation_flag: direction_mismatch, poor_translation_fit, ambiguous_match, search_edge); translation_consistent is a fit-consistency flag, not proof of a translation

