# ZNF384 fusion versus wild type — execution protocol for HPC (readouts A, B, C)

*2026-11-12. Written before any of these analyses is run. Executed by an agent on HPC, where the data are.
Background, screening and the paused pilot: `PROFILESHAPE_ZNF384_SCENARIO.md`. This protocol does not change the
stop rule of `PROFILESHAPE_APPLICATIONS_PLAN.md`; whether the result enters the paper is the user's decision.*

> 中文概要：问题是 TCF3-ZNF384 融合蛋白与野生型 ZNF384 相比，在总量没有明显变化的区域上，信号形状是否不同。
> 分三种读出信号：A 是 TF 结合谱本身，先验最弱；B 是 H3K27ac 的扩散，先验最强；C 是 ATAC/MNase 的 NDR 和核小体位置。
> 每种读出各有一套预先写定的主分析、独立验证量、参照区域和分级标准。执行前先做数据清点并提交（第 0 步），
> 不满足门槛就停下来记录，不改标准。场景 C 的教训已经写进去：未检出区域也有一致性，所以必须和参照区域比较。

---

## Rules for the executing agent (fixed)

1. **Statistical core frozen.** Before anything else, check `sha256sum src/profileshape/regional.py | cut -c1-16`
   equals `21ffd82395ee207f`. Do not edit `src/`. New code goes in `benchmarks/znf384_*.py`, results in
   `benchmarks/results/znf384_<readout>/` (JSON + a short `.md`), each recording commit, uncommitted-diff count,
   the regional.py hash, this protocol's sha256 and start/finish times (as `benchmarks/app_scenario_c.py` does).
2. **Order:** step 0 → commit → readouts in the order B, A, C (strongest prior first) → each result recorded on this
   page under *Results* before the next readout starts. A readout whose data are absent is recorded as
   "not available" and skipped.
3. **No tuning.** Window, preset, representation, region set, thresholds and grading are fixed below. Anything
   this page leaves open is decided in a dated *Clarification* committed **before** the analysis it concerns, and
   may not use that analysis's results. Gates that fail are reported as failures; thresholds are not moved.
4. One main analysis and at most one sensitivity analysis per readout (named below). More calls are not success.
5. Fresh seeds: bootstrap and any random step use seeds 9 920 000 + k (k = 0, 1, … in the order used, recorded).
6. Records kept as written: `benchmarks/results/**` are never overwritten; a rerun goes to a new directory with a
   reason in the page.
7. Commit after step 0 and after each readout with a message ending
   `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Push only when the user asks.
8. Stop and ask the user (do not improvise) if: the design differs from anything below (e.g. no WT arm), a gate
   fails in a way that leaves the readout ungradable, or a code error would change a conclusion.

## Step 0 · Data inventory (fill in, commit, before any count is looked at)

Write a table on this page under *Step 0 result*, one row per sample: sample id, group (WT / fusion; any other
construct listed but not used), readout (A / B / C), assay and antibody or tag, cell background, biological
replicate unit (clone, culture, patient), batch / sequencing run, genome build, paired / single-end, read
pairs after filters, duplicate handling, fragment-length median and IQR (paired-end), file path on HPC. Then state
per readout which gates hold:

| gate | requirement | if not met |
|---|---|---|
| G1 design | WT and fusion profiled with the same antibody or tag, same cell background, processed alike | not gradable; report as description only |
| G2 replicates | ≥ 3 biological replicates per group | 2 vs 2: application demonstration, no grade (the 2 vs 2 null is known to exceed nominal); held-out step impossible |
| G3 input | reads (BAM or FASTQ), not only processed tracks | tracks only: midpoints are not recoverable; coverage representation from tracks, stated; depth from `bigwig_count_unit` |
| G4 depth | median over samples of core-window fragments ≥ 100 at the analysed regions (computed after the region set is fixed, condition-blind) | stop the readout; record the depth |
| G5 fragment length | per-group median fragment length differs by ≤ 10 % | not a stop: the coverage-representation sensitivity becomes mandatory and is the readout's one sensitivity |

Also run `profileshape consistency` per group and readout (report 07) and record its flags. A flagged sample is
**not removed** in the main analysis; removing it is allowed only as the readout's one sensitivity, and only if
G5 does not already use that slot (then report the flag as a limit).

## Shared definitions

- **Representation.** BAM input: fragment midpoints (`bam.mode: midpoint`), duplicates removed, MAPQ ≥ 10, proper
  pairs; size selection per readout below. Coverage only as the named sensitivity.
- **Regions (condition-blind).** ZNF384 sites = summits of peaks called by MACS2 (`-q 0.01`, paired-end `-f BAMPE`)
  on the **pooled** reads of all A-readout samples of both groups (equal reads subsampled per group before
  pooling, seed 9 920 000), main chromosomes, merged within 200 bp, blacklist removed. If readout A is not
  available, the user supplies the ZNF384 site list in a Clarification. Never use per-group peak calls to define
  regions.
- **Amount classes (after testing, not before).** edgeR QL on core-window fragment counts (TMM, design ~group,
  `glmQLFit`/`glmQLFTest`), per region log2FC and SE = |log2FC| / √F (F = QL F statistic; regions with log2FC = 0
  get SE from the median SE of their depth decile). Classes: **up / down** = BH-FDR < 0.05; **equivalent** = 95 %
  interval log2FC ± 1.96 SE within ± log2(1.25); **undetermined** = the rest. Only "equivalent" is called "no amount
  change"; non-significance alone is not.
- **Shape calls.** `shape_any` (3 vs 3 or the full design), BH-FDR 0.05, `per_region_tests(analyse(...))` as in
  `app_scenario_c.calls`, with the readout's preset. Shape tests run on all regions; amount classes are
  cross-tabulated afterwards (selection on label-dependent amount statistics before the shape test would break the
  permutation null).
- **Per-sample shape measures.** Width = mean absolute deviation, and centroid, of the mass-normalised,
  background-subtracted, clipped core profile (`app_scenario_c.shape_measures`).
- **Held-out sign agreement** (needs G2): calls on replicates 1–2 of each group (replicate order = the sample
  order of step 0, fixed before counts); direction = sign of mean(fusion 1–2) − mean(WT 1–2); held-out sign =
  mean(fusion 3…n) − mean(WT 3…n). Per direction subgroup: share agreeing, exact 95 % Clopper–Pearson interval.
- **Reference agreement** (scenario C lesson: not-called regions agreed too, a genome-wide trend). The same
  held-out agreement among **tested, not-called** regions of the same amount class and depth decile distribution
  (not-called regions resampled to the called regions' depth-decile counts, seed fixed). Statistic: difference
  called − reference, 95 % bootstrap over regions (2000 resamples).
- **Comparators** (reported, not graded): edgeR amount (above); diffSpliceDGE on the same bins
  (`app_scenario_a.diffsplice`); per-replicate aggregate profile (centroid, MAD, IQR) over all regions and over the
  equivalent class.

## Readout B · H3K27ac (strongest prior; run first)

- **Question.** At ZNF384 sites with equivalent H3K27ac amount, does H3K27ac redistribute (spread or shift) with
  the fusion, which recruits EP300/CBP through the TCF3 activation domains?
- **Regions.** ZNF384 sites (shared definition) within 1 kb of a pooled H3K27ac peak (MACS2 on pooled H3K27ac
  reads, `--broad` off, same equal-subsampling rule). Anchor: the ZNF384 summit. Orientation: + towards the
  nearest annotated TSS (RefSeq Select, ≤ 50 kb; regions without a TSS within 50 kb are unoriented and analysed
  separately, descriptively). Preset `histone_narrow` (flank 2000, bin 25, core ± 1500).
- **Size selection.** Fragments 120–600 bp (mono- and di-nucleosomal).
- **Main result.** (i) Cross-table of `shape_any` × amount class; (ii) held-out sign agreement of **width** in the
  equivalent class; (iii) called − reference agreement for width.
- **Grading.** Consistent: (a) ≥ 50 `shape_any` calls in the equivalent class **and** ≥ 25 % of all `shape_any`
  calls are in it; (b) held-out width agreement interval above 0.5 in each direction subgroup with ≥ 30 regions
  (a subgroup with < 30 regions is reported, not graded); (c) called − reference difference interval above 0 for
  width. Partly consistent: (a) and one of (b), (c). Not supported: otherwise.
- **Descriptive, not graded:** centroid towards / away from the TSS (oriented regions), share of widening versus
  narrowing; association between H3K27ac widening and ZNF384 binding change at the same site (readout A amount
  log2FC; Spearman, if A is available).
- **Sensitivity (one):** coverage representation (mandatory if G5 fails), otherwise exclusion of a sample flagged
  in report 07, otherwise none.
- **Boundaries.** Acetylation amount and breadth are linked; "equivalent" is an interval statement at ± 25 %,
  not zero change. No mechanism is claimed from the shape alone.

## Readout A · ZNF384 / fusion binding profile (weakest prior)

- **Question.** At sites with equivalent binding amount, does the TF profile shift (choice of another motif
  instance in the region) or change width (clustered or cooperative binding)?
- **Regions.** ZNF384 sites (shared definition); anchor summit, unoriented. Preset `tf_chip` (ChIP) or
  `tf_cutrun` (CUT&RUN / CUT&Tag). Size selection: ChIP 0–1000 bp; CUT&RUN ≤ 120 bp (TF footprint fragments).
- **Main result.** As readout B, with **width and centroid** as measures (each graded separately; the readout is
  consistent only if both are, partly if one is).
- **Grading.** As readout B, per measure.
- **Independent measure (descriptive, not graded; no direction predicted).** ZNF384 motif matches (JASPAR CORE
  ZNF384 matrix — ID and version recorded in the Clarification before the run; FIMO p < 10⁻⁴) in the core window:
  for regions with a held-out-reproduced centroid change, the share whose shift points towards the side with more
  motif matches, with an exact interval, and the same in the reference regions.
- **Sensitivity (one):** as readout B.
- **Boundaries.** DNA-binding specificity is expected to be shared (all zinc fingers kept), so a null result is
  plausible and is reported as such. Different expression levels of the constructs change amount; the equivalent
  class handles this only partly. With a tag, the comparison is fusion-tag versus WT-tag, not endogenous protein.

## Readout C · Accessibility / nucleosome position (ATAC or MNase)

- **Question.** At ZNF384 sites with equivalent accessibility, does the nucleosome-depleted region change width or
  do flanking nucleosomes move?
- **Regions.** ZNF384 sites; anchor summit; orientation as readout B. Preset `atac` (flank 600, bin 10, core
  ± 300) for ATAC; `mnase` with anchor `center` for MNase (stated in a Clarification).
- **Size selection (main).** ATAC: fragments ≤ 100 bp (nucleosome-free). MNase: 120–180 bp (mononucleosomes).
- **Main result.** As readout B with **width** as the measure.
- **Independent measure (not used by the main analysis).** ATAC: from fragments 180–247 bp only, per region and
  sample the distance between the strongest nucleosomal-midpoint peak upstream and downstream of the anchor
  (smoothed 30 bp, search 100–500 bp each side); per-region change = fusion − WT mean. Prediction (direction
  stated before the run): NDR widening calls have a larger nucleosome-to-nucleosome distance increase than
  reference regions — difference of medians, 95 % bootstrap. MNase: lag-150–190 bp autocorrelation over ± 1 kb as in
  scenario B. Graded as an additional criterion (d): interval excludes 0 in the predicted direction.
- **Grading.** Consistent: (a), (b), (c) as readout B and (d). Partly: (a) and two of (b)–(d). Not supported:
  otherwise.
- **Sensitivity (one):** as readout B. Tn5 insertion bias is not corrected; stated.
- **Boundaries.** ATAC NFR width depends on digestion and fragment size distribution (hence G5); the ATAC depth
  boundary (`PROFILESHAPE_INDEPENDENT_VALIDATION.md`: 18–44 fragments per region gave no per-region detection)
  applies; G4 guards it.

## Reporting

Per readout: the gates table, the cross-table, held-out and reference agreement with intervals, comparators, the
sensitivity, the grade and its limits, three representative regions chosen by a fixed rule (25th / 50th / 75th
percentile of held-out width change among called equivalent-class regions) with per-replicate profiles. Add one
row to the progress log of `PROFILESHAPE_REVISION_REPORT.md` per readout. Do not write paper text.

## Step 0 result

(to be filled on HPC, committed before any readout)

## Clarifications

(dated, each committed before the analysis it concerns)

## Results

(per readout, recorded after its run; nothing above changed)
