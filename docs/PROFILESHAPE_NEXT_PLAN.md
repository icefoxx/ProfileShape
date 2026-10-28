# ProfileShape — integrated plan after v3.1 (draft for review, 2026-10-27)

Sources: the within-group consistency work and its evaluation (`docs/PROFILESHAPE_CONSISTENCY_EVALUATION.md`,
`benchmarks/results/consistency_v1`, `consistency_v2`), an external-style review ("major revision"), and a
reading of RegCFinder (Weiß & Friedel 2023, *Bioinformatics Advances*; paper and `amss` module source).

Nothing below has been run yet. Items marked **run** need a pre-specified protocol (written first, reviewed,
then run on fresh seeds) as in `CLAUDE.md`.

## 0. Where things stand

| item | state |
|---|---|
| within-group consistency module (`consistency.py`, `profileshape consistency`, report 07, figure H) | implemented, 41 tests pass, **not committed** |
| consistency evaluation v1 | defect found (near-empty windows); kept on record; K4, K7, K10 failed |
| consistency evaluation v2 (new seeds, after amendment A1) | K1–K3, K5–K9 pass; **K4 fails** (interval coverage 0.85 in 2 of 31 real-noise cells; TF ChIP TAF15 pair also 0.85 in v1); **K10 not evaluable** (no null call in 62 new splits) |
| new data | GSE260577 (H3K27ac/H3K4me3 ChIP, mouse liver, mm39), GSE285245 (H3K27ac/H3K4me3 CUT&RUN, mm10), mouse RefSeq; caches `inputs/caches/cons_*.npz` |
| open from v3.1 | independent positive case; yeast 2 vs 2 excess; translation vs other redistribution; per-class calibration with few regions; per-region amplitude; graphical abstract |

Review points accepted (checked against code and docs): RegCFinder missing from the positioning; the pooled
null is the central untested-in-general assumption; yeast excess unresolved and not fixed by the consistency
module; real applications (1433/1555, 10867/13797 regions called) exceed the simulated share changed (≤ 50 %);
default outputs (per-region absolute amplitude, consistency "95 % CI") are stronger than the evidence; ATAC
failure may be assay–hypothesis mismatch, not only depth; biological value shown as consistency with known
effects, not added insight. Partly overstated: MMDiff / DANPOS / csaw are already positioned as prior work in
`PROFILESHAPE_LITERATURE.md` and `APPLICATIONS.md` (but the manuscript plan is stale).

## 1. Work packages

### P0-A · Make default outputs match the evidence (no run)

Decision taken: per-region amplitude tests **stay in the default output**.

- Report 05: state that per-region absolute amplitude (`log2_auc`, `log2_signal`) was miscalibrated under
  complete nulls (FWER = FDR 0.20 and 0.275 in two designs) and is not a claimed inference; keep the
  normalisation caveats. Same sentence in the guide and in `region_tests_summary` descriptions.
- Report 07 and TSVs: label the bp interval "approximate (coverage below nominal in 2 of 31 real-noise cells)";
  for two replicates, state that the pair's offset p-value can be too small; keep π̂₁ as experimental with the
  observed range (0.05–0.32 for a true 0.20 on real noise).
- Wording: "spatial redistribution" instead of "movement"; `translation_consistent` is a fit-consistency flag;
  a non-significant region at low depth is not evidence of no change; per-class calibration needs ≥ 3000
  (2 vs 2) / ≥ 1000 (3 vs 3) regions in the tested mix — warn below.
- Tests: assert the new caveat strings are present.

Deliverable: one commit; no change to `regional.py` statistics.

### P0-B · Positioning and manuscript plan (no run)

- Add RegCFinder to `PROFILESHAPE_LITERATURE.md` and the comparison table in `APPLICATIONS.md`, with what it
  does (MSS segmentation of replicate-averaged within-window density differences; position-permutation filter;
  DEXSeq on data-derived features) and how the question differs (localisation vs replicate-supported screening).
- Rewrite `PROFILESHAPE_MANUSCRIPT_PLAN.md` (still quotes v2.4 numbers): core claim — *replicate-supported
  screening of spatial signal redistribution at predefined regions with 2–3 replicates, via label-invariant
  stabilisation, a covariate-stratified pooled permutation null and a combined test, with coverage-track input,
  replicate-consistency diagnostics and a constrained shift description*; list what is not new; state that
  2 vs 2 inference rests on the cross-region assumption.
- Ablation list for the paper (each needs a result file; several exist): label-invariant s0 vs per-assignment
  s0; stratified vs unstratified null; Cauchy combination vs separate lists; translation consistency rule vs
  none; signed support normalisation vs clipped.

### P0-C · Comparison with RegCFinder, MMDiff2 and a count model — **deferred** (decision 2026-10-27)

Existing comparisons: edgeR diffSpliceDGE with pre-specified endpoints and the current code (`round4`–`round6`
diffsplice; also in the ATAC validation); MMDiff2 in simulation and on real yeast MNase and TF CUT&RUN BAMs
(`PROFILESHAPE_COMPARISON.md` §2, §4, §9, §11.6; v2.x code). Missing: RegCFinder, and MMDiff2 with the current
code. The design below is kept for when the package is resumed.


- Implementation: RegCFinder via its own code (`amss` R script + featureCounts + DEXSeq, defaults ρ = 1,
  1000 randomisations); the Watchdog wrapper only if the scripts cannot be run alone. MMDiff2 as in earlier
  rounds; edgeR diffSpliceDGE (existing comparator).
- Data: real-null splits (yeast MNase 2 vs 2, TF CUT&RUN 3 vs 3, histone 2 vs 2 — BAM needed for RegCFinder:
  yeast and CUT&RUN BAMs exist; histone datasets are bigWig-only, so either fetch reads or exclude); simulated
  counts with known truth (the diffSplice simulator); insertions into real-null splits.
- Endpoint at the window level (a window counts as called if any of its sub-regions is significant):
  observed FDP, power, untestable windows, runtime. Also localisation accuracy for RegCFinder (overlap of
  called sub-region with the inserted change), where ProfileShape has no counterpart.
- Pre-specified predictions: RegCFinder's window-level FDP under nulls exceeds nominal (segmentation ignores
  replicate variance; the position-permutation null does not reflect replicate variability; features are
  selected and tested on the same data). If it does not, that is reported as such.
- Fairness: comparator settings fixed before seeing results; IWTomics-type methods only as background /
  misuse examples.

### P0-D · Independent positive case with added explanatory value (**run**, selection and protocol first)

- Selection criteria written before searching: assay matched to the hypothesis (the expected change must be
  a redistribution of the analysed fragments); depth ≥ 100 fragments per region and sample in the core window
  (median); ≥ 3 replicates per condition; expected per-region effect large (e.g. a perturbation with known
  strong pausing / positioning / spreading phenotype); an **independent readout** from the same study
  (RNA-seq, PRO-seq, nascent RNA or a second assay).
- Pre-specified analyses: direction and size of shape changes; whether shape calls explain the independent
  readout **beyond** amplitude, depth, baseline expression and peak width (regression with these covariates;
  replication in a second condition or dataset); statements "amplitude unchanged" only with an equivalence
  margin or an amplitude-adjusted model.
- Deliverable: protocol → dataset choice reviewed → download → run → report including failures.

### P1-E · Consistency with user-set tolerance margins (implementation, then **run**)

- Bias-corrected sample offset (noise expectation removed) with a confidence interval, in % of the aggregate
  profile and in bp; a width-change measure (%) for broad marks.
- Config: `consistency.margin_offset_pct`, `consistency.margin_shift_bp` (and width %), **no default
  margins**; three-way verdict per sample / pair / split: *within margin* (upper bound < margin), *beyond
  margin* (lower bound > margin), *inconclusive*; the upper bound always reported ("offsets above X % are
  excluded"). Two replicates: pair-level verdict, no attribution.
- Assay presets (`assay: atac | mnase | tf_chip | tf_cutrun | histone_broad | histone_narrow | proseq`) set
  window, bin, background, mass floor, block size and the units shown — not margins. Reference ranges of real
  replicate offsets from the evaluated datasets shown as context, recomputed after bias correction.
- "Share of regions differing" margins are not offered (π̂₁ and F_δ not accurate enough); documented.
- Evaluation: coverage of the upper bound; rate of a wrong "within margin" verdict when the true offset
  exceeds the margin ≤ 0.05 (primary); inconclusive rate by depth; fresh seeds; same datasets plus a held-out one.

### P1-F · Stress tests of the pooled null (**run**)

- Share changed 70–100 %; same-direction global shifts; shared perturbations across regions; asymmetric and
  multimodal real shapes; batch structure; spatial dependence beyond duplicated windows (overlapping windows,
  arrays). Real-data analogues by insertion into real-null splits. Criteria: whole-list FDP as in the
  calibration study; report where it breaks.

### P1-G · Reproducibility (no run)

- Pinned environment file; script that fetches all inputs (GEO/SRA accessions in `docs/DATA.md`) and rebuilds
  caches; tagged release; results compressed (`*.json.gz`) before committing (v1+v2 JSON ≈ 50 MB).

### P2 · Carried over

Yeast 2 vs 2 excess (cause), separating translations from other redistributions, per-class calibration with
few regions, per-region amplitude test, graphical abstract.

## 2. Order and dependencies

1. P0-A, P0-B (no run) → commit together with the consistency module (results gzipped). **Done.**
2. Protocols for P0-C and P0-D written and reviewed (dataset search for P0-D needs only metadata).
3. P0-D download and run (P0-C deferred).
4. P1-E implementation and protocol, then run; P1-F; P1-G alongside.

## 3. Decisions needed

| # | question | proposal |
|---|---|---|
| 1 | commit the consistency module and v1/v2 results now (gzipped)? | yes, after P0-A |
| 2 | RegCFinder: own scripts vs Watchdog workflow | own scripts with its defaults; Watchdog only if needed |
| 3 | histone data for P0-C (bigWig only) | exclude from the RegCFinder arm, or fetch reads (≈ 40 GB) |
| 4 | P0-D depth threshold and readout requirement | ≥ 100 fragments / region / sample; independent nascent-RNA or RNA-seq readout |
| 5 | P1-E before or after P0-C/D | after (the review weighs the positive case higher) |
