# ProfileShape — manuscript plan (working draft, rewritten 2026-10-27)

Replaces the v2.4-era plan (its numbers — e.g. "98 % downstream, +17.5 bp", "V1–V3 pass in all cells" — are
superseded; current evidence: `docs/PROFILESHAPE_V3_REVIEW.md`). Work packages referred to as P0-x / P1-x are
in `docs/PROFILESHAPE_NEXT_PLAN.md`.

## Positioning

**Core claim.** At predefined genomic regions, ProfileShape screens for *replicate-supported spatial
redistribution of signal* between two conditions with 2–3 biological replicates, separately from changes in
amount, directly on coverage tracks or BAM files, and says what this adds beyond amount-based analysis.

**How (the contribution is the combination and its evaluation, not a new theory):**
1. a label-invariant stabilising constant in the per-region statistics;
2. a permutation null pooled over regions within strata of a label-invariant covariate;
3. one list from an omnibus (kernel) and a directional (projection) test, combined per region;
4. a translation description (bp) reported only where it fits, as a descriptive effect size;
5. replicate-consistency diagnostics (sample offsets, sub-group structure) and depth / batch diagnostics.

**What it rests on (stated in the paper).** With 2 replicates per group there are 3 label assignments per
region; inference is possible only through the cross-region assumption of the pooled null. Many regions do
not replace biological replicates. The assumption held in simulation and several real-noise settings and
failed on yeast MNase 2 vs 2 (FDP 0.10–0.13, associated with replicate structure).

**Not claimed.** A proven FDR bound; exact permutation validity across regions; per-region absolute amplitude
inference (miscalibrated, reported only for exploration); mechanism ("translation" is a fit-consistency
description of a population profile, not observed movement of a nucleosome or polymerase); localisation of the
changed sub-region (RegCFinder does this); genome-wide discovery without anchors; "amount unchanged" without an
equivalence margin or an amount-adjusted model.

## Prior work to position against

MMDiff (kernel shape test with a replicate-derived null), DANPOS (position / fuzziness / occupancy), csaw
(window counts; can detect local redistribution), RegCFinder (sub-regions of relative density change in
predefined windows, DEXSeq), IWTomics (curves as units — different target; a misuse example when replicates
are few), DiffBind / edgeR (amount). Details: `docs/PROFILESHAPE_LITERATURE.md`, `APPLICATIONS.md`.

## Evidence map

| claim | evidence in hand | still needed |
|---|---|---|
| whole-list FDR of the shape tests near nominal | simulation and semi-real grid incl. a pre-specified confirmatory run on new seeds (`calibration/confirm_v3`); PRO-seq 3 vs 3 real nulls; ATAC (GSE346229) real nulls and insertions | stress tests beyond the simulated range: 70–100 % of regions changed, same-direction global shifts, shared perturbations, multimodal / asymmetric shapes, batch, spatial dependence (P1-F) |
| per-class calibration | ≥ 3000 regions (2 vs 2), ≥ 1000 (3 vs 3) in the tested mix | a null for few regions (open) |
| real-data limits | yeast MNase 2 vs 2 above nominal; associated with replicate structure; the consistency module flags that structure (`consistency_v2`) but flagged structure did not produce null calls in 62 new splits | cause; whether a pre-analysis rule can be validated on independent data |
| sensitivity vs comparators | edgeR diffSpliceDGE, pre-specified (`round6/diffsplice`); MMDiff2 in simulation and on real BAMs (v2.x code, `PROFILESHAPE_COMPARISON.md` §2, §4, §9, §11.6) | RegCFinder; MMDiff2 with the current code (P0-C, deferred) |
| positive cases | yeast remodeler mutants (ordering matches an independent phasing measure); NELFB PRO-seq (downstream redistribution) — development data, consistency with known biology | **an independent case** with matched assay, sufficient depth, an independent readout, and added explanatory value beyond amount, depth, expression and width (P0-D) |
| translation description | direction correct in every evaluated cell for the selected subset; estimate ÷ true 0.94–1.08 for ≥ 50 bp with power ≥ 0.3; coverage of the description can be small | selection-adjusted intervals (open) |
| replicate consistency | pre-specified evaluation on 11 datasets / 7 assays (`consistency_v2`): calibration passes; interval coverage fails in 2 of 31 cells; sub-region share (π̂₁) rough | user-set tolerance margins with calibrated upper bounds (P1-E) |
| practicality | bigWig / BAM / matrix input; runtime and memory vs MMDiff2 (`COMPARISON.md` §10) | pinned environment, input-fetch script, tagged release (P1-G) |

## Ablations (one result file each)

| design element | evidence so far | to do |
|---|---|---|
| label-invariant s0 vs per-assignment s0 | development finding: 50 % of regions shifted, per-assignment s0 gave far more false calls than a common s0 (`COMPARISON.md` §4) | re-run as a clean ablation with the current code |
| stratified vs unstratified pooled null | v2.4 null vs v3 null on the same seeds (`calibration/dev_v24` vs `dev_v3`; region-count boundary `round3/mixed_n_*`) | — (exists) |
| Cauchy combination vs union of two BH lists | design argument; union not FDR-controlled | simulation result file |
| translation consistency rule vs none; stricter R² | R² 0.8 rule chosen on development seeds and rejected on new seeds (`round5/shift_r2_0.8_rejected`) | report the rule's coverage vs direction errors in one table |
| signed support normalisation vs clipped / whole-window | `COMPARISON.md` §11.3 (amplitude-only calls 1–3 % with clipping; power loss with whole-window mass) | result file with the current code |

## Figures (draft)

1. Problem: amount vs shape; why regions and reads are not replicates; what RegCFinder / MMDiff / count models target.
2. Method: per-region statistics, pooled stratified null, label-invariant s0 (ablation panel).
3. Calibration: confirmatory grid, real nulls, region-count boundary, stress tests; the yeast failure shown, not hidden.
4. Power and effect size: vs replicates and depth; translation description and its coverage.
5. Independent positive case with the added-value analysis (P0-D).
6. Replicate consistency and practical use (inputs, diagnostics, runtime).

## Threats to validity (to state)

- One simulator family developed alongside the method; mitigated by fresh-seed confirmatory runs, real noise
  and insertion, not eliminated.
- Development history (amended rules, fixed defects, a defect found by the consistency evaluation v1) is
  documented; criteria were pre-specified for confirmatory runs.
- Real-data null splits reuse the same replicates; they are candidate nulls, not gold standards.
- Shared-noise and replicate structure can make splits non-exchangeable (yeast).

## Target venues

A focused methods / software paper (e.g. NAR Genomics and Bioinformatics, Bioinformatics, BMC Bioinformatics)
once P0-A, P0-B and P0-D are complete; the independent positive case weighs more than further simulation.
