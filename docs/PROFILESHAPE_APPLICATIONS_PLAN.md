# Applications for the paper — scenario pages (written before any of these analyses is run)

*2026-10-28. Replaces the search for a confirmatory positive case as the route to the paper's biological section.
Goal of the paper: with few replicates, find region-level spatial redistributions that window counts or aggregate
profiles do not show, and recover interpretable patterns in a few known perturbations. Not required: a new
mechanism, the same direction in every case, superiority over every tool in every setting.*

## Rules (fixed)

- **Statistical core frozen** (`regional.py` sha256 `21ffd82395ee207f`). Only errors that change a conclusion are
  fixed; no new classifier, margin framework or general background model.
- Per scenario: one main analysis, at most one sensitivity analysis aimed at a named concern. No change of window,
  subset or representation to obtain agreement.
- Each result is graded **consistent / partly consistent / not supported** against the criteria below.
- More calls alone is not success: replicate consistency, an independent measure or locus-level support is
  required.
- Development datasets are labelled as application demonstrations, not independent confirmation.
- Stop looking for data once three complementary, interpretable cases exist; failures are kept and reported.
- A2-type prediction gains are supplementary evidence, not an entry condition.
- Kept out of the main line: per-region absolute amplitude test (miscalibrated), stretch verdict (withdrawn),
  margin framework and consistency module (supplementary tools, with their evaluations), the second positive-case
  search (`PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md`, stopped at steps 1–3 with no candidate meeting all
  conditions; record kept).
- Comparators per scenario: one strongly related method (edgeR QL diffSpliceDGE on the same bins, the round-6
  configuration) and one routine baseline (window counts with edgeR for amount; aggregate meta-profile).
  RegCFinder's local segment output is a different output; it is discussed as complementary, not raced.

---

## Scenario A · Directional redistribution (promoter-proximal pause)

- **Question.** After NELF loss, does ProfileShape recover a downstream redistribution of PRO-seq signal at the
  pause, with replicate support, and show how heterogeneous it is across genes?
- **Data.** PRO-seq GSE233930 (development dataset, labelled as such): DMSO vs VHL-mediated NELFB degradation
  2 h, 2 vs 2; the 48 h contrast as within-study replication. Existing cache and settings of
  `benchmarks/positive_cases.py`. Real-null calibration of PRO-seq exists at 3 vs 3, not 2 vs 2 — stated.
- **Main result.** (i) Among `shift`-significant genes with a consistent translation, the share downstream;
  (ii) the distribution of per-gene shift estimates (bp) and the share `Changed` with `poor_translation_fit`
  (heterogeneity: moved peak versus redistribution into the gene body); (iii) per-replicate meta-profiles and three
  representative loci chosen by a fixed rule (the genes at the 25th, 50th and 75th percentile of the shift
  estimate among consistent calls).
- **Comparators.** diffSpliceDGE on the same bins (does it detect these genes; it reports no direction); the
  aggregate meta-profile (pause peak and 0–300 bp centroid per replicate).
- **Grading.** Consistent: downstream share ≥ 0.90 at 2 h and 48 h, and the per-replicate centroid moves
  downstream in every VHL replicate relative to both DMSO replicates. Partly: one of the two. Not supported:
  neither.
- **Boundaries.** Development data; 2 vs 2; `translation_consistent` is a fit-consistency flag, not proof of a
  translation; the bp values are selected-subset estimates (selection bias near the detection limit).
- **Sensitivity (if needed).** None planned; the named concern would be the 2 vs 2 null — answered by the WT-like
  null splits only if a 3 vs 3 design exists in the study.

## Scenario B · Complex shape change not summarised by a translation (nucleosome arrays)

- **Question.** In remodeler mutants that disrupt nucleosome arrays, does ProfileShape detect per-gene changes
  that a translation does not describe, and do the detected genes agree with an independent per-gene phasing
  measure?
- **Data.** Yeast MNase-seq GSE280671 (development dataset): WT (4) vs isw1Δ chd1Δ (4), ATG ± 1.5 kb, the existing
  4 vs 4 analysis of `positive_cases.py` (4 vs 4 avoids the known 2 vs 2 excess; the three WT null splits had no
  call).
- **Independent measure (not used by ProfileShape).** Per gene and sample: autocorrelation of the dyad profile at
  the nucleosome repeat (lags 150–180 bp, the lag window of the meta-profile measure already reported), averaged
  within genotype; per-gene change = dKO − WT.
- **Main result.** (i) Calls of `shape_any` and the reasons attached to `shift` calls (expected: `poor_translation_fit`,
  no bp value forced on them); (ii) the per-gene phasing change in called versus not-called genes (difference of
  medians with a bootstrap interval over genes); (iii) replicate support without selection bias: genes called on
  replicates 1–2 of each genotype (2 vs 2), phasing change measured on replicates 3–4 only.
- **Comparators.** diffSpliceDGE on the same bins (overlap of calls and its phasing change in the same way);
  window counts with edgeR (do amount calls track phasing loss).
- **Grading.** Consistent: (ii) called genes lose more phasing (interval excludes 0) and (iii) the held-out phasing
  change agrees in sign; partly: one of them; not supported: neither.
- **Boundaries.** Development data; the phasing measure is a summary, not a gold standard; MNase digestion level
  can change phasing measures.
- **Sensitivity (if needed).** The single mutants (2 replicates each) are reported descriptively only.

## Scenario C · Changes that window counts and aggregate profiles miss

- **Question.** Can region-level shape testing show redistribution that an amount analysis and an aggregate
  profile do not express — opposite-direction subgroups hidden in the aggregate, or shape changes without an amount
  change?
- **Data.** GSE218269 (CDK7 inhibition, Pol II MNase-ChIP 3 vs 3; not a development dataset), the B2 BAM cache
  (20–120 bp, midpoints). **Exploratory**: depth criterion S3 failed (median 85 < 100), the confirmatory case is not
  evaluable and stays so; the expected upstream translation was not supported (B4, B5).
- **Main result.** (i) Aggregate: centroid and width per replicate (already: centroid +5.6 bp [−1.6, +12.9]);
  (ii) per-region: `shape_any` calls cross-tabulated against amount calls (edgeR on core-window counts, FDR 0.05);
  (iii) direction heterogeneity with replicate support free of selection: regions called on replicates 1–2 of each
  group (2 vs 2); the sign of their width change and of their centroid change measured on replicate 3 of each group
  only; agreement share against 0.5.
- **Comparators.** edgeR window counts (amount); diffSpliceDGE on the same bins; the aggregate profile.
- **Grading.** Consistent: a substantial share of shape calls without an amount call (≥ 25 % of calls) and held-out
  sign agreement > 0.5 with a binomial interval excluding 0.5, in both directions of change; partly: one of them;
  not supported: neither.
- **Boundaries.** Exploratory, below the planned depth; no mechanism claimed; B5 found the `shift` result
  representation-dependent — only midpoints (the B2 representation) are used here.
- **Sensitivity (if needed).** Coverage representation, reported alongside, not substituted.

## Optional · Failure boundary (short)

ATAC GSE346229 (independent validation): calibration held, no per-region detection at 18–44 fragments per region
with a ≈ 2 bp genome-wide effect. Reported as the depth boundary, with the GSE218269 depth story; no new analysis.

## Clarifications (before the scenario they concern is run)

- **B, phasing lag window (2026-10-28, before B is run).** The page says "lags 150–180 bp, the lag window of the
  meta-profile measure already reported"; that measure (`positive_cases.phasing`) uses lags 150–190 bp (15–19 bins).
  The stated intent — the same window as the reported measure — is followed: 150–190 bp, gene-body side (0 to
  +1200 bp), mass-normalised per gene and sample. Counts for diffSplice and edgeR come from the CPM coverage tracks
  as area / (count quantum × single-fragment footprint) per bin (fragment equivalents, rounded) — an approximation
  stated with the result. Bootstrap intervals: 2000 resamples of genes, seed 9 910 000. Held-out "agrees in sign":
  the held-out difference of medians (called − not called) has the same sign as the main one.

## Order and stop

A → B → C, each with its page completed by its result and grade. Three graded cases (whatever the grades) end
the data phase; then writing. At most one new dataset, and only if a scenario cannot be run at all on the data
above.
