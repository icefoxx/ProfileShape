# Independent positive case, second attempt — selection protocol (draft for review)

*Written 2026-10-28, before any search result for this attempt is looked at and before any download. Same
staging as the first attempt (`PROFILESHAPE_POSITIVE_CASE_PROTOCOL.md`): stage 1 (this file) fixes the selection
procedure and the analysis framework; stage 2 (appended after a dataset is chosen, still before the full download)
fixes the dataset-specific predictions, the representation and the power check. Changes afterwards go under
"Amendments". The sha256 of this file is recorded in the results.*

## What the first attempt taught (GSE218269; B1–B5)

1. **The depth estimate before download was about four times too high.** Stage 1 estimated ≈ 250–500 fragments
   per active promoter from run sizes × an assay-typical fraction; the exact count after alignment was a median
   of 85 (20–120 bp, TSS ± 500 bp). Run size × assumed fraction is not a usable estimate.
2. **≥ 100 fragments is not a power guarantee.** Whether the expected change could be seen depended on the
   effect size, the replicate variability and the representation, none of which S3 measures.
3. **The representation decided what the `shift` test saw.** Midpoints: 0 significant regions in all three
   fragment ranges; fragment coverage: 382–638, 0.79 upstream (B5 D3, exploratory). Earlier evidence: in
   simulation both representations are calibrated and coverage is slightly more powerful (v3 plan, bigWig vs
   BAM on the same fragments); real-data calibration exists for each representation on different data
   (coverage: GSE112560 semi-real, CUT&RUN null; midpoints: yeast MNase, ATAC GSE346229) but never on the same
   reads.
4. **`mad_width` has almost no power at this depth** (B5 D2), and the A2 gain came from the signed shift alone
   with a sign opposite to the expected one (B5 D1, D3). An A2 gain without the expected sign does not support
   the stated mechanism.

## Prerequisite R · representation calibration on the same reads (before stage 2)

The two representations have not been compared under a real null on the same reads; since stage 2 chooses
one, both must be shown calibrated on the same data first.

- Data (BAMs already aligned): TF CUT&RUN GSE291134, line A673, six replicates — the ten 3 vs 3 splits (the
  design of the positive case; same regions as `cutrun_null.py`, which used the bigWigs); yeast MNase GSE280671
  WT, four replicates — the three 2 vs 2 splits (the dataset with the 2 vs 2 excess); ATAC GSE346229 — its
  real-null splits as in `independent_validation.py`.
- Representations: midpoint and fragment coverage from the same BAMs, all other options as in each dataset's
  earlier analysis.
- Outcome: number of `shape_any`, `shape_kernel` and `shift` calls per split (whole list at FDR 0.05), and the
  share of p < 0.01, per representation.
- Criterion: coverage is eligible for stage 2 only if, on each dataset, its calls per split and its p < 0.01
  share are not higher than those of midpoints by more than the Monte Carlo noise of the comparison (fixed in
  the script before running: a difference in p < 0.01 share above 0.005 on CUT&RUN or ATAC, or more calls in
  2 of the 3 yeast splits, is "higher"; on CUT&RUN, in addition, coverage's whole-list calls must be 0 in at
  least 9 of the 10 splits). Yeast midpoints are known to exceed nominal at 2 vs 2; the comparison is relative,
  not a claim that either is calibrated there. If coverage fails, stage 2 uses midpoints.
- Results go to `benchmarks/results/representation_null_v1`, reported whatever they show.

## Selection procedure

### Step 1 · Search and screening (metadata only)

Families (S1, as before): Pol II / nascent RNA redistribution within the promoter-proximal window
(pause-release, initiation or escape perturbations); H3K4me3 breadth perturbations; nucleosome positioning /
spacing perturbations measured by MNase. Human, mouse, yeast. Excluded (S5): the seven development datasets
of the first protocol, GSE218269 (first attempt) and its companion series.

Required as before: S1 assay–hypothesis match **with a quantitative expected effect from an independent
source** (another study or an orthogonal assay in the same study, not the dataset's own processed tracks),
S2 ≥ 3 replicates per condition processed alike, S4 an independent readout with ≥ 2 replicates, S6 reads public.

### Step 2 · Pilot depth (replaces run size × fraction)

For every candidate that passes step 1, before any full download: the first 2 million read pairs of one
replicate per condition (ENA ranged download), aligned and filtered exactly as the full analysis would; the
count of fragments per region in the core window, scaled to the run's total read pairs by the pilot's
alignment and filter rates. Recorded per candidate before ranking. (Processed tracks, if the study provides
them, are not used for depth: B1 showed that their fragment selection need not match the analysis.)

### Step 3 · Power by simulation (before ranking)

For each candidate: `profileshape.simulate`, 3 vs 3 (or the candidate's design), 2000 regions, depth set to the
pilot estimate, replicate noise taken from the closest assay of `consistency_v2` (mapping fixed here: Pol II /
nascent RNA → PRO-seq group; MNase → yeast MNase; H3K4me3 → H3K4me3 ChIP or CUT&RUN; the table of parameters
is written by the script before candidates are evaluated), the expected effect inserted in the share of
regions the independent source implies, both representations (coverage only if prerequisite R passed).
Outcome: power of the stage-2 primary A1 test at FDR 0.05, 100 runs, fresh seeds 9 800 000 +.
**Requirement: power ≥ 0.80 for at least one eligible representation.** That representation is fixed for the
candidate (if both reach 0.80, the one with higher power; within 0.05 of each other → coverage if eligible,
else midpoints).

### Step 4 · Ranking and order

Among candidates meeting steps 1–3: (a) simulated power, (b) more replicates, (c) a second condition or dose for
replication, (d) the independent readout closer in time to the profiled assay. The ranked list is recorded
before the first full download. Candidates are taken in order; a candidate failing a post-download check (step
5) is reported as not evaluable and the next is taken — this is fixed now, so taking the next one is not a
choice made after seeing data.

### Step 5 · After the full download, before any test

Exact checks, each failing → not evaluable (no change of thresholds or windows): (i) exact core-window depth ≥
0.8 × the pilot estimate and ≥ 100 fragments (median over samples); (ii) the consistency report: no sample whose
offset exceeds 3× its group-mates' (one such sample → the protocol's pre-specified leave-that-sample-out design
is the primary analysis, if ≥ 3 replicates remain; otherwise not evaluable); (iii) power re-run with the exact
depth ≥ 0.70.

## Analysis framework (predictions in stage 2)

- **A1 · Expected direction.** One primary prediction with a numeric criterion, chosen so that the fixed
  representation and bin size can resolve it (expected effect ≥ 2 bins and inside the core window); secondary
  predictions as needed. A prediction whose denominator can be empty (a share among significant regions) is
  paired with a condition for "not evaluable" (fewer than 50 significant regions).
- **A2 · Added explanatory value (primary endpoint).** As in the first attempt (folds by chromosome,
  chromosome bootstrap with folds kept per chromosome, B3), with two changes fixed now: (1) the expected sign of
  the shape statistic's coefficient is part of the criterion (gain lower bound > 0 **and** coefficient sign as
  predicted, with its own bootstrap interval excluding the opposite sign); (2) the four-model ablation (base,
  + shift, + kernel, + both; same rows, folds and draws) is always reported.
- **A3 · Sensitivity (reported, no criterion).** The other representation and the leave-one-replicate-out
  designs, as in B5 D3.

Every result is reported; exploratory analyses are labelled as such and never replace a confirmatory verdict.

## Amendments

**R1 (2026-10-28, before prerequisite R is run).** (a) The CUT&RUN BAMs on disk are the first 4 M read pairs per
replicate aligned to chr1 + chr2 only (made for the MMDiff2 comparison, `realbam_cutrun.py`); R uses them as
they are, with the condition-blind anchors of `cutrun_null.py` on those two chromosomes; the number of regions
and the depth are reported, and the limitation is stated with the result. The yeast BAMs are the first 3 M pairs
(genome-wide); the ATAC fragments are those of the GSE346229 validation cache (all reads), from which coverage is
rebuilt exactly (each fragment's span, from its midpoint and length). (b) The criterion applies to each of
`shape_any`, `shape_kernel` and `shift`: coverage is "higher" on a dataset if it is higher for any of the three
tests (stricter than one test). (c) Paths in each dataset: the same fragment filters for both representations
(yeast and CUT&RUN: MAPQ ≥ 10, proper pairs, `extract` defaults otherwise; ATAC: the validation's nucleosomal /
sub-nucleosomal ranges and anchors). Script `benchmarks/representation_null.py`.
