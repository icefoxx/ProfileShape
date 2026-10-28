# Scenario A · directional redistribution at the pause (PRO-seq NELFB, development data)

Plan sha256 `facbcd6cfc90247d`; commit `43f3644134` (uncommitted 0); regional.py `21ffd82395ee207f`; named cache identical to the original: True; 2026-10-28T13:01:46 → 2026-10-28T13:03:28.

| contrast | tested | shape_any | shift called | consistent down / up | downstream share | Changed (no direction) | consistent shift bp 5/25/50/75/95 % | diffSplice calls | diffSplice among shape_any | diffSplice among downstream | shape_any among diffSplice |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DMSO vs VHL2h | 13797 | 10867 | 9117 | 3455 / 22 | 0.994 | 5640 {'direction_mismatch': 21, 'poor_translation_fit': 5612, 'search_edge': 7} | +8 / +11 / +16 / +22 / +40 | 5726 | 0.52 | 0.69 | 0.99 |
| DMSO vs VHL48h | 13923 | 8667 | 4816 | 1885 / 7 | 0.996 | 2924 {'direction_mismatch': 6, 'poor_translation_fit': 2910, 'search_edge': 8} | +8 / +11 / +14 / +19 / +31 | 2643 | 0.30 | 0.43 | 0.98 |

Per-replicate meta-profile (2000 promoters with the highest pooled signal): DMSO_rep1 peak 35 bp, centroid 63.6 bp; DMSO_rep2 peak 35 bp, centroid 57.9 bp; VHL2h_rep1 peak 45 bp, centroid 81.0 bp; VHL2h_rep2 peak 45 bp, centroid 80.4 bp; VHL48h_rep1 peak 45 bp, centroid 80.7 bp; VHL48h_rep2 peak 45 bp, centroid 74.4 bp

Representative genes (25 / 50 / 75 % of the consistent 2 h shift): RETREG3, TSKU, RAD54L (`representatives_2h.png`).

**Grade: consistent** (criteria in the plan). Development data, seen before the plan was written: the grade was largely foreseeable; this is a demonstration, not a test.
