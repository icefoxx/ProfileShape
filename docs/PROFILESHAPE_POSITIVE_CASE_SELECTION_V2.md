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

## Outcome of prerequisite R (run representation_null_v1, 2026-10-28; recorded after the run)

Coverage **eligible**: on CUT&RUN (10 × 3 vs 3), yeast (3 × 2 vs 2) and ATAC (18 analyses) no test had more calls
with coverage than with midpoints and the difference in p < 0.01 share was −0.0006 to +0.0009 (margin 0.005).
Limits of the evidence: every split had 0 whole-list calls in both representations, so the comparison of calls does
not discriminate and the p < 0.01 shares carry the information; the CUT&RUN subsample is shallow (median 10
fragments per region, 801 regions); the p < 0.01 share varies between splits (0.000–0.045) alike in both
representations. "Eligible" means not worse than midpoints here, not shown calibrated in general. The criterion compares shares averaged over splits, so differences of opposite sign can cancel; single splits differ by more (e.g. ATAC, an NDR split, `shape_kernel` 0.0350 → 0.0485, +0.0135). Coverage met the pre-specified average criterion; per-split non-inferiority is not shown (added after an external review; the verdict is not changed after the fact).

## Step 1 result (2026-10-28; metadata and literature only, nothing downloaded)

GEO E-utilities, four queries (Pol II / nascent RNA with perturbation; MNase with remodelers; H3K4me3 breadth
regulators; pause-release perturbations: flavopiridol, NVP-2, CDK9, DRB, triptolide, NELF, SPT5, heat shock),
human / mouse (yeast for MNase), 2016–2026: 2725 series; automated screen of sample titles for ≥ 2 groups of the
profiled assay with ≥ 3 replicates: 113 series; read by hand for S1–S6. The title-based screen misses series whose
sample titles do not encode replicates; this is a limitation of the search, not a selection after data.

| series | perturbation · profiled assay (design) | independent readout | S1 (independent source) | status |
|---|---|---|---|---|
| GSE112379 | heat shock vs none, HEK293 · RPB1 ChIP-seq (scramble 3 vs 3) | PRO-seq 4 vs 4 (same study) | Mahat et al. 2016 (PRO-seq): heat shock represses thousands of genes by reduced pause release → Pol II concentrated at the pause relative to downstream | **long list** |
| GSE205504 | NELFB-AID +IAA, empty vector vs NELFB rescue, MEF · Pol II ChIP-seq (3 vs 3, spike-in) | RNA-seq GSE205388, but in the floxed/Cre system, not the AID lines | Aoi et al. 2020: after NELF depletion Pol II accumulates at a second, downstream pause associated with the +1 nucleosome | long list; S4 weak (readout from another system) |
| GSE199086 | SEC disruptor KL-2 vs DMSO, SHEP · TT-seq / FP-TT-seq (3 vs 3) | RNA-seq 3 vs 3 | Liang et al. 2018: impaired release from pause, slower elongation | long list; assay–hypothesis match weaker (TT-seq near the TSS) |
| GSE120201 | dominant-negative TFIIS (Dox) vs H2O · Pol II ChIP (3 vs 3) | GRO-seq / NET-seq 2 vs 2 | only the study itself | S1 not met (no independent source) |
| GSE97527 | shPAF1 vs control · Pol II ChIP (3 vs 3) | — | conflicting literature on PAF1 and pause release | S1 not met |
| GSE260852 | cohesin degradation · RPB1 ChIP | mNET-seq | — | S2 not met (2 replicates) |
| GSE288059 | CFP1 deletion, GC B · H3K4me3 ChIP (3 vs 3) | RNA-seq 3 vs 3 | the paper reports loss at a subset of genes (amount), not breadth | S1 not met |
| GSE289627 / GSE289600 | CFP1 deletion, pro-B · H3K4me3 CUT&Tag | RNA-seq | CFP1 broad domains | S2 not met (2 knockout replicates) |
| GSE295011 · GSE216057 · GSE234251 | KDM5-C70 · SMARCA4/5 KD MNase · ISWI subunit KO MNase | none found / — | — | S4 not met (no readout found) or S2 not met |

Next: step 2 (pilot depth) for the three long-list candidates, then step 3 (power simulation); ranking after both.

## Step 3 specification (2026-10-28, written before any pilot count is looked at)

Script `benchmarks/power_sim.py`; seeds 9 800 000 +; 100 runs per candidate and representation; 3 vs 3; 3000
regions drawn from the candidate's condition-blind selection.

- **Shapes:** each region's template is its pooled pilot profile (both pilot samples summed, condition-blind,
  smoothed over 3 bins), used through `simulate(..., templates=...)`; depth scaled so that the median core-window
  count equals the pilot estimate of the full run (step 2).
- **Noise (fixed for all candidates, from `consistency_v2`, PRO-seq group, the closest assay to Pol II ChIP /
  nascent RNA in that set):** replicate amplitude sd = median over the top-half regions of the within-group sd of
  log region totals; bin dispersion = median CV² of the top-decile bins across replicates (an upper bound: it
  includes the counting noise, so power is underestimated rather than overestimated); AR parameters at their
  defaults. The values are written to the result file before the candidates run.
- **Coverage representation:** the simulated midpoint counts convolved with the candidate's pilot fragment-length
  distribution (paired-end) or its estimated fragment length (single-end).
- **Expected effect and primary A1 criterion per candidate** (from the independent sources of step 1):

| candidate | inserted in the alternative group | regions changed | primary A1 criterion (met in a run) |
|---|---|---|---|
| GSE112379 heat shock | signal downstream of +300 bp halved (reduced pause release; Mahat et al. 2016) | 30 % (repression of "several thousand" genes; conservative) | `shape_any` share ≥ 0.10 **and** ≥ 50 `shift`-significant regions with ≥ 0.80 upstream |
| GSE205504 NELF depletion | profile translated +60 bp (accumulation at a second, +1-nucleosome-associated pause; Aoi et al. 2020) | 50 % | ≥ 50 `shift`-significant regions with ≥ 0.80 downstream |
| GSE199086 KL-2 | signal downstream of +300 bp reduced by 30 % (impaired pause release; Liang et al. 2018) | 50 % | as for heat shock |

Power = share of runs meeting the criterion. Requirement (step 3): ≥ 0.80 in at least one eligible
representation. These effect sizes are assumptions taken from the sources' qualitative descriptions, not
measured; the power is conditional on them, and the stage-2 predictions will use the same values.

## Steps 2–3 result (2026-10-28)

| candidate | step 2: depth, complexity-adjusted (linear) | step 3 power, midpoint / coverage (dispersion 0.15) | note |
|---|---|---|---|
| GSE112379 heat shock | 199 (1195) | 0.00 / 0.00 (0.00 / 0.00) | only ≈ 8 % of the core-window signal lies in +300…+500 bp, so halving the signal downstream of +300 bp moves ≈ 4 % of the mass: the specified effect is mostly outside the core window (TSS ± 500 bp) |
| GSE205504 NELF depletion | 86 (166) | 1.00 / 1.00 (1.00 / 1.00) | +60 bp translation in 50 % of regions; predicted-direction share 0.90 (midpoint) / 0.99 (coverage) |
| GSE199086 KL-2 | 0 (0) | not simulated | library ≈ 0.2 M distinct fragments; alignment 30 % |

No candidate meets all conditions (see R5). NELF met the simulated A1 criterion under the assumed translation model only. Open points before any full download (not decided here): (1) NELF's complexity-adjusted
depth (86, conservative for single-end) is below step 5's ≥ 100, so the case may become not evaluable after download
under the rule as written; (2) NELF's readout (RNA-seq) comes from the floxed/Cre system, not the AID lines, so A2 is
weak or not definable; (3) heat shock fails step 3 only because of the window: a design with the core window covering
the gene-body part (e.g. TSS ± 1 kb) is a new specification and would have to be written and simulated before any
full download.

## Amendments

**R5 (2026-10-28, after an external review of steps 2–3; before the region-level estimate is computed).** R3 scaled
promoter counts by the growth of distinct fragments in the whole library. That assumes the selected core windows
gain distinct fragments at the same rate as the library; at dense promoters duplicates saturate earlier, so R3 may
**overstate** depth there, not understate it as R3 claimed. Added: the same Poisson fit restricted to fragments with
a midpoint in the selected core windows (25 / 50 / 100 % points), and the region-level scale reported beside the
library-wide one. Both are model-based extrapolations from three shallow points; they are reported as estimates
with model uncertainty, not as expected depths. Also corrected in the record: (a) step 3 has **not** found a
candidate meeting all conditions — NELF met the simulated A1 criterion only under an assumed rigid +60 bp translation
in 50 % of regions, which is not a quantitative effect from an independent source as S1 requires, and its readout
comes from another system; (b) heat shock's failure is shown only for the fixed TSS ± 500 bp window — that the
window is the only cause is not shown.

**R4 (2026-10-28, after the first step-3 run).** The first run (`positive_case_v2_selection/power`, kept) gave power
0.00 for heat shock with not a single `shape_any` call in any run — implausible for halving the signal downstream of
+300 bp in 30 % of regions. Cause: an implementation error, not the specification — the script applied the effect to
the template's peak component only, while the specification says the *signal* is halved (and the profile
translated). With real Pol II templates much of the signal downstream of +300 bp sits in the template's background
component (estimated from the window ends, floored at 0.05 of the apex): only ≈ 32 % of it was in the peak component
for heat shock (≈ 56 % for NELF), so the inserted effect was far smaller than specified. Corrected to act on the whole
expected signal of the changed regions; the specification (effects, shares, criteria, seeds) is unchanged and the
full run is repeated into `positive_case_v2_selection/power_r4`. Both runs are reported.

**R3 (2026-10-28, after the first step-2 estimate; before step 3 is run).** The step-2 estimate scaled
deduplicated pilot counts linearly to the run size. Duplicate shares in the pilots grow with depth (heat shock: 3.1 /
6.0 / 11.2 % at 25 / 50 / 100 % of the pilot, as expected for a library of ≈ 7 M distinct fragments), so linear
scaling overstates the depth after duplicate removal (≈ 6× for heat shock) — the kind of error that misled the first
attempt. Changed: duplicates are removed in the full analysis of these ChIP data (standard; `extract` default), and
the scale factor is the expected number of distinct fragments at the run's size over that in the pilot, from a
Poisson library-complexity curve fitted by least squares to the 25 / 50 / 100 % subsample points of each pilot
(fit residuals reported). Where duplicates also arise without PCR (single-end reads at dense sites), the Poisson fit
overstates duplication and so understates depth — the conservative direction. The linear estimate is kept in the
result for comparison. Both estimates were seen before this amendment: heat shock 1195 (linear), NELF 166, KL-2 0.

**R2 (2026-10-28, before step 3 is run; no pilot count looked at).** The protocol's noise from the PRO-seq group
gives replicate sd 0.10 and bin dispersion 0.036 (deep PRO-seq tracks). Pol II ChIP is likely noisier, so this may
overstate power. Added, reported only (not part of the requirement): the same runs with dispersion 0.15 (the
simulator default), fresh seeds (9 800 000 + 10 000 × candidate + 1000 + run).

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
