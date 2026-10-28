# Within-group consistency — evaluation protocol (v1)

*Written 2026-10-27, before the new datasets (GSE260577, GSE285245) were downloaded and before any run with the
seeds below. The code under evaluation is `src/profileshape/consistency.py` as of the run; its sha256, the
regional.py hash, the commit and the uncommitted-diff hash are recorded in the results. This file's sha256 is
recorded too; changes after the first run are appended under "Amendments" with the reason, never edited in place.*

## Questions

1. **Calibration.** Without within-group differences, are the sample flag (offset test + ≥ 2× group-mates)
   and the split flag (per-region distance contrast) at or below their nominal level (0.01 family-wise per
   group and diagnostic), including when one replicate is shallower and when neighbouring regions are dependent?
2. **Size.** Does the translation-equivalent offset (bp) recover a known sample-wide shift, and does its
   block interval cover the truth?
3. **Extent.** How do the descriptive shares F₁₀ / F₂₀ / F₅₀ (with the group-mates' reference) and the
   experimental split π̂₁ respond to a known share of shifted regions, and how far are they from it?
4. **Sensitivity.** Which kinds of inconsistency are detected (sample-wide shift, broadening, a sub-group)?
5. **Functional consistency.** On new real data, do splits of a group that the split test marks as
   structured carry more calls when used as a real-data null for the differential tests?

Not addressed: equivalence (no margin is tested), genome-wide statements (all shares are shares of the analysed
regions), the biological origin of any real inconsistency.

## Data and regions

Regions are condition-blind in every dataset (annotation, or the pooled signal of all samples of the mark).

| id | assay | data | groups (replicates) | regions | window / bin / background / core (bp) |
|---|---|---|---|---|---|
| mnase | MNase-seq | GSE280671 bigWigs (development data) | 6 genotypes × 4 | ORF ATG, oriented (`ins_yeast.npz`, coordinates from `orf_atg.bed`) | 1500 / 10 / 300 / 1000 |
| atac | ATAC-seq | GSE346229 BAM fragments (`psv_cache.npz`) | Sth1 0 h (4), Sth1 1 h (4), Snf2 0 h (3), Snf2 1 h (3) | +1 nucleosome dyad (condition-blind, as in the independent validation), nucleosomal fragments | as in `independent_validation.py` |
| tf_chip | TF ChIP-seq | GSE112560 bigWigs | ZNF384 WT (3), TAF15 (2), TCF3 (2) | `consensus_top5000.bed` | 2000 / 20 / 500 / 1000 |
| tf_cutrun | TF CUT&RUN | GSE291134 bigWigs (`cutrun_profiles.npz`) | A673 (6), A2.2 (6) | pooled-signal summits (as in `cutrun_null.py`) | 1000 / 10 / 300 / 600 |
| proseq | PRO-seq | GSE233930 (`ins_proseq.npz`) | NELFB, PAF1, RTF1 DMSO (2 each) | expressed TSS, oriented | 1000 / 10 / 300 / 500 |
| k27ac_chip | H3K27ac ChIP-seq | GSE260577 CPM bigWigs, mouse liver | Sham (4), CLP (4) | pooled-signal summits of the 8 samples: 5000 strongest local maxima at 200 bp, ≥ 5 kb apart, refined at 10 bp | 3000 / 50 / 500 / 2500 |
| k4me3_chip | H3K4me3 ChIP-seq | GSE260577 | Sham (4), CLP (4) | (a) pooled-signal summits, ≥ 3 kb apart; (b) RefSeq TSS (assembly from the bigWig header), oriented, the 5000 with the highest pooled signal within ±1 kb | 2000 / 25 / 400 / 1500 |
| k27ac_cutrun | H3K27ac CUT&RUN | GSE285245 bigWigs, mouse CD8 T cells | Naive day 30 (4), Arm day 30 (4) | as k27ac_chip | as k27ac_chip |
| k4me3_cutrun | H3K4me3 CUT&RUN | GSE285245 | Naive day 30 (4), Arm day 30 (4) | as k4me3_chip (a) and (b) | as k4me3_chip |

Summit anchors are oriented along the genome (+). If a new dataset turns out unusable (e.g. a file missing, an
assembly that cannot be matched), it is reported as such and dropped; no replacement is chosen after seeing data.

## Arms

**A. Simulation** (`profileshape.simulate`, parametric peaks). Seeds 9 950 000 + 1000 × cell + run; 100 runs per
cell; 2000 regions. Cells: replicates n ∈ {3, 4, 6} × depth ∈ {12, 40, 160}, each as (i) exchangeable null and
(ii) null with replicate 1 simulated at 30 % depth (same shapes). Insertion cells at depth 40, n ∈ {3, 4, 6}:
replicate 1 shifted +20 and +50 bp in every region; replicate 1 shifted 50 bp in 10 % and in 30 % of regions
(random sign per region); replicate 1 shifted +50 bp in 30 % of regions (same sign); replicate 1 broadened ×1.3
in every region; for n ≥ 4, the second half of the replicates shifted +20 bp in every region and, separately,
50 bp in 20 % of regions (random sign).

**B. Semi-real null** (real replicate noise, no within-group structure by construction). Per dataset and group,
regions are cut into blocks of 10 consecutive regions along the genome and, in every block independently, the
replicate labels are permuted at random (one permutation per block, so local dependence is kept). Any
sample-wide difference between real replicates is spread over all pseudo-replicates; the expected sample
offset and split contrast are zero. 20 shuffles per group, seeds 9 960 000 + 1000 × dataset + 100 × group +
shuffle. For atac, additionally: after shuffling, pseudo-replicate 1 thinned to 30 % of its fragments.

**C. Semi-real insertion.** On the shuffled data of B (same seeds), the insertions of arm A are applied to the
signal by linear interpolation (shift) or rescaling about the anchor (broadening, background kept, mass
conserved; as `insertion.insert`), 20 shuffles each, for groups with ≥ 3 replicates (splits: ≥ 4).

**D. Real data as is** (descriptive). The consistency report of every group; no truth.

**E. Functional consistency.** For every group with 4 replicates (all 2 vs 2 splits) or 6 (all 3 vs 3 splits)
in the new datasets atac, tf_cutrun, k27ac_chip, k4me3_chip (a, b), k27ac_cutrun, k4me3_cutrun: the split's
contrast z from the consistency module and the share of regions called by `shape_any` (BH 5 %) when the split is
analysed as a real-data null with `per_region_tests`. mnase (development data) is reported separately and does
not enter the criteria.

## Criteria (fixed before running)

| id | arm | criterion |
|---|---|---|
| K1 | A null cells (i) and (ii) | share of runs with any sample or split flag in the group ≤ 0.04 in every cell |
| K2 | B | runs with any flag ≤ 2 of 20 for every group and dataset (atac thinned included) |
| K3 | A (i), (ii) pooled per kind | offset test P(p < 0.01) in [0.005, 0.02]; split contrast P(p < 0.01) in [0.005, 0.02] |
| K4 | A null, B | 95 % interval of the sample shift (bp) covers 0 in ≥ 0.90 of samples, every cell / dataset |
| K5 | A, sample-wide +20 / +50 bp | mean estimate ÷ truth in [0.8, 1.1] for each n |
| K6 | A (i), B | mean excess F₂₀ (sample − group-mates) within ±0.03 |
| K7 | A (i), (ii), B | split π̂₁ mean ≤ 0.02 |
| K8 | A, half 50 bp in 20 % | mean split π̂₁ ≤ 0.25 (no overestimation; its shortfall from 0.20 is reported) |
| K9 | A, depth ≥ 40 | sample-wide +20 bp flagged in ≥ 0.8 of runs (n ≥ 3); half +20 bp split-flagged in ≥ 0.8 (n ≥ 4) |
| K10 | E | Spearman(contrast z, `shape_any` null-call share) > 0 with one-sided p < 0.05 over all new splits; if any split is flagged, flagged splits have a higher mean null-call share than unflagged ones |

Reported without a criterion: sensitivity in C per dataset; π̂₁ and F-excess against the inserted share; the
attenuation of the bp estimate on real shapes; everything in D.

Consequences: a metric that fails its calibration criterion (K1–K4, K6, K7) is not reported as a calibrated
quantity — it is removed or labelled as failing in the report text, and the failure is stated in the results.
K8–K10 failing limits what is claimed. K10 is an association on the same samples (the contrast and the null calls
are computed from the same data), not proof that the split test predicts calibration on unseen data.

## Amendments

**A1 (2026-10-27, after run v1).** Run v1 (`benchmarks/results/consistency_v1`, seeds 9 950 000+ / 9 960 000+,
protocol sha256 `7e0a418d9a83108e`) is kept as written. Its outcome: K1, K2, K3, K5, K6, K8, K9 pass; K4, K7 and
K10 fail. Inspection found a defect in the code under evaluation: regions were usable whenever every sample's
support mass was positive, so near-empty windows (ATAC: 10 % of regions below 0.1 × the sample's median mass)
were divided by a tiny mass, reaching values of ~5 × 10¹⁶ that dominated every mean over regions. On the
affected datasets the split contrast was degenerate (z = ±1.00 in every split of atac and k27ac_chip), so v1's
B and E results for those datasets do not test the method. K7 failed in simulation (6 replicates, one at 30 %
depth: split π̂₁ 0.044), traced to splits into halves of unequal size, where a shallow replicate's noise does not
enter every label assignment alike. K10 failed with no variation in the outcome: none of the 62 new real-null
splits had a `shape_any` call. K4 failed in two cells (mnase WT 0.85 of 80 sample intervals; tf_chip TAF15 0.85
of 20).

Changes, made before any v2 run:
1. A region is used only if every sample's support mass is ≥ 0.1 × that sample's median (label-free).
2. Split π̂₁ is computed only for splits into equal halves.
3. Scoring of K10: if the null-call share is identical in every split, K10 is reported as *not evaluable*
   instead of failing (a correlation with a constant is undefined). The criterion itself is unchanged.

Run v2 repeats the protocol unchanged otherwise (same arms, cells, numbers of runs and criteria) on new seeds:
simulation 9 970 000 + 1000 × cell + run, semi-real 9 980 000 + 1000 × dataset + 100 × group + shuffle. v2 is the
run for which the criteria are judged; v1 remains on record with its failures.
