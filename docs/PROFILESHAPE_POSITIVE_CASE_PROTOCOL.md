# Independent positive case — selection and analysis protocol (draft for review)

*Stage 1 (this file, 2026-10-27): selection criteria and the analysis framework, written before any candidate's
data are downloaded. Stage 2 (appended after the dataset is chosen, still before download): the
dataset-specific predictions. Changes after stage 2 go under "Amendments". The sha256 of this file is recorded
in the results.*

## Purpose

Show, on data not used in developing ProfileShape, that (1) shape calls behave as expected for a perturbation
with a known, strong redistribution of the measured signal, and (2) the shape information explains an
independent readout **beyond** amount, depth, baseline expression and peak width. (1) alone would repeat what
the development datasets already show; (2) is the point.

## Selection criteria (all required)

| id | criterion | how checked before download |
|---|---|---|
| S1 | assay–hypothesis match: the perturbation is known (from other studies) to redistribute the *measured* signal within predefined regions, at many regions, by a per-region effect large relative to the window (e.g. Pol II / nascent-RNA signal moving from the promoter-proximal pause into the gene body after pause-release perturbation; H3K4me3 domain broadening; nucleosome array disruption measured by MNase) | literature; the study's own description |
| S2 | ≥ 3 biological replicates per condition, processed alike | GEO sample list |
| S3 | depth: expected median ≥ 100 fragments per region and sample in the core window | SRA run sizes × an assay-typical fraction in regions; **verified after download and before any test**; failing → dataset dropped, next ranked candidate taken |
| S4 | an independent readout of the same conditions from the same study (RNA-seq, nascent RNA of a different kind, or a second assay), with ≥ 2 replicates | GEO series |
| S5 | not used in ProfileShape development: excludes GSE280671, GSE233930, GSE112560, GSE291134, GSE346229, GSE260577, GSE285245 | — |
| S6 | reads public (FASTQ/BAM) or bigWigs with documented processing; human, mouse or yeast | GEO / SRA |

Ranking among candidates meeting S1–S6, fixed now: (a) larger expected per-region effect, (b) more replicates,
(c) reads available (BAM path, exact depth), (d) a second condition or dose for replication within the study.
The ranked list is recorded before the first download; candidates are taken in order.

## Analyses (framework; dataset-specific predictions added in stage 2)

- **A0 · Data checks.** Depth per region (S3); within-condition consistency report (`profileshape
  consistency`); if ≥ 4 replicates, within-condition real-null splits (shape_any share called).
- **A1 · Expected direction.** The stage-2 prediction for the perturbation (e.g. share of significant regions
  with gene-body redistribution; direction of consistent shifts), with a numeric criterion.
- **A2 · Added explanatory value (primary).** Per region, the independent readout (e.g. log2 fold change of
  the RNA level of the gene) is modelled by (i) amount covariates (`log2_signal_relative` estimate, baseline
  signal), depth, baseline expression, peak / domain width, and (ii) the same plus the shape statistic (signed
  shift where oriented; −log10 p of `shape_any`). Primary endpoint: gain in out-of-sample R² (10-fold
  cross-validation over regions, folds by chromosome) of (ii) over (i), with a bootstrap interval over
  chromosomes; criterion fixed in stage 2. Replication: the same model in the study's second condition / dose
  if present.
- **A3 · Amount-unchanged statements** only with an equivalence margin for `log2_signal_relative` fixed in
  stage 2, or within the A2 model (shape effect adjusted for amount).
- **A4 · Comparators** on the same regions: amount-only (edgeR on window counts, or `log2_signal_relative`) and
  diffSpliceDGE as in round 6, for the A2 gain.

Failures are reported as failures; no analysis choice is changed after the data are seen except under
"Amendments" with the reason.

## Candidate list (stage 1 search, 2026-10-27; metadata only, nothing downloaded)

Search: GEO (E-utilities) for pause-release / elongation perturbations profiled by Pol II ChIP, PRO-seq, TT-seq,
mNET-seq or PRO-cap, and for H3K4me3 breadth perturbations (KDM5, COMPASS, CFP1); human and mouse; 481 series
screened by sample titles for ≥ 2 conditions with ≥ 3 replicates of the profiled assay plus RNA-seq. Depth from
SRA run tables; the per-region figure is an estimate (assay-typical share of fragments in the regions) to be
verified after download.

| rank | series | perturbation · assay | replicates | independent readout | reads per sample | S1–S6 | note |
|---|---|---|---|---|---|---|---|
| 1 | GSE218269 (CDK7as, HEK293) | 1-NM-PP1 30 min vs DMSO · Pol II (RPB1) MNase-ChIP-seq | 3 vs 3 | TT-seq labelled RNA, 15 / 30 / 60 min (2 each); mNET-seq; START-seq | 46–52 M pairs (≈ 250–500 fragments per active promoter, estimated) | all met | authors report Pol II retention at promoters on CDK7 inhibition (a within-window redistribution); replication within the study: SPT5 and Ser5P MNase-ChIP-seq (3 vs 3) and other TT-seq time points |
| 2 | GSE152595 (mESC) | Set1B KO + Set1A ΔSET vs WT · H3K4me3 ChIP-seq | 4 vs 4 | RNA-seq (WT 7, mutant 6) | 18–51 M single-end | all met | H3K4me3 breadth (domain narrowing) — the stretch dimension on a broad mark; replicates in pairs of submission batches (to be checked by the consistency report) |
| 3 | GSE245600 (SUDHL5) | CDK-TCIP1 2 h / 4 h vs DMSO 4 h · Pol II ChIP-seq | 3 vs 3 (two time points) | RNA-seq (3–4) | 24–29 M pairs | S3 borderline | effect concentrated on BCL6-target genes (few regions) |

Screened out: GSE112608 (the 4-replicate Pol II ChIP is in *Drosophila* S2 cells — S6); GSE180845 (ChIP-seq 2
replicates; TT-seq 3 — S2 for ChIP, weak S1 for TT-seq); GSE349745 (TT-seq 4 vs 4, but the RNA-seq readout is
at 72 h with rescue constructs — S4); GSE188417 (replicated Pol II CUT&Tag only in controls — S2 per contrast).

Proposed: rank 1 as the primary case, rank 2 as the second case (broad mark, stretch). Stage-2 predictions are
written for the chosen dataset(s) before download.

**Decision (2026-10-27).** Rank 1 (GSE218269) is the primary case. Rank 2 (GSE152595) is withdrawn as a positive
case after checking the paper (Sze et al. 2020, *Sci Adv*, doi:10.1126/sciadv.aaz4764): it reports *loss of
H3K4me3 coverage* within broad domains, not narrowing of the domains, so S1 (a known redistribution of the
measured signal) is not met. It may be analysed later as an exploratory amount-versus-shape case, labelled as such.

## Stage 2 · GSE218269 (written before download)

**Literature basis** (Velychko et al. 2024, *Mol Cell* 84:2287, doi:10.1016/j.molcel.2024.05.007, Fig. 3B,
S4B): in Pol II (RPB1) MNase-ChIP-seq of short fragments (< 70 bp), 30 min of CDK7 inhibition decreases Pol II in
the promoter-proximal pause region (TSS to +100 bp) and strongly increases it ~30 bp upstream of the TSS (the
PIC position); TT-seq shows an immediate global decrease of synthesis. The expected change is a redistribution
within the TSS window from downstream to upstream of the TSS.

**Data.** GEO processed bigWigs (coverage of < 70 bp fragments, normalised by the authors to non-transcribed
regions): Pol II MNase-ChIP-seq DMSO 30 min (3) and 1-NM-PP1 30 min (3). Readout: TT-seq labelled RNA
coverage (spike-in normalised, strand-specific), DMSO and 1-NM-PP1 at 15, 30 and 60 min (2 each). Depth is
estimated from the bigWig count quantum (no reads); S3 is checked on that estimate.

**Regions** (condition-blind): RefSeq Select (hg38, NM_, main chromosomes, one TSS per gene), genes ≥ 5 kb; of
these, the half with the highest pooled Pol II signal (all 6 samples) within TSS ± 300 bp, and pooled TT-seq
signal (all 30-min samples) > 0 in the readout window. Windows oriented by strand; flank 1000, bin 10,
background 300, core 500 (preset `proseq`).

**A1 predictions (criteria).**
- P1: `shape_any` calls in ≥ 20 % of tested regions (BH 5 %).
- P2: among regions significant in `shift`, ≥ 80 % have an upstream direction (towards lower coordinates of
  the oriented window).
- P3: the aggregate centroid moves upstream (Welch 95 % interval below 0).

**A2 (primary).** Readout y = log2 fold change (1-NM-PP1 / DMSO, 30 min) of mean sense-strand TT-seq coverage
in [TSS + 1000, TSS + 4000] (pseudocount 1 % of the median), replicate means. Model (i): baseline Pol II core
mass (log, DMSO mean), Pol II amount change (log2 ratio of core mass means), depth proxy (log pooled window
signal), baseline TT-seq (log, DMSO mean), log gene length, pooled profile width (`mad_width`). Model (ii): (i)
plus the signed `shift` statistic (direction × −log10 p) and −log10 p of `shape_kernel`. Ordinary least
squares; out-of-sample R² by 10-fold cross-validation with folds formed by chromosome; gain = R²(ii) − R²(i)
with a 95 % interval from 500 bootstrap resamples of chromosomes. **Criterion: lower bound of the gain > 0.**
Reported: the coefficient sign of the signed shift statistic (expectation: more upstream redistribution ↔
larger synthesis decrease, i.e. a positive coefficient).

**Replication (reported, no criterion):** A2 with the 15-min and 60-min TT-seq readouts.

**Checks (A0):** S3 on the estimated depth (median ≥ 100 fragments per region and sample in the core window;
failing → the case is reported as not evaluable, no replacement chosen after seeing data); consistency report
of both groups (3 replicates each; no real-null splits possible).

**Not tested:** comparisons with other methods (P0-C deferred); A3 (no "amount unchanged" claim is made).

## Amendments

**B1 (2026-10-28, after the A0 check of GSE218269; no A1 / A2 result was computed).** S3 failed: estimated
median 58 fragments per region and sample (per sample 51–65; count quantum 0.87–1.11, single-fragment footprint
58 bp, over the ±1 kb window — the core window has fewer), against the required ≥ 100. The processed tracks keep
only fragments < 70 bp, which gives the resolution the authors needed (PIC vs pause) at the cost of depth. Result
on record: `benchmarks/results/positive_case_gse218269` ("not evaluable"). Stage 1 said a failing dataset is
replaced by the next ranked candidate; stage 2 said "reported as not evaluable, no replacement chosen after
seeing data". The later, more specific rule (stage 2) is followed: no replacement is taken automatically. Rank 2
was withdrawn before download; rank 3 (GSE245600) has about half the reads per sample. Any further step (a new
search with depth estimated per region before selection, or GSE218269 re-processed from reads with a wider
fragment range) is a new decision, recorded here before it is carried out.

**B2 (2026-10-28, decision after B1; still no A1 / A2 result computed).** GSE218269 is re-processed from reads
with a wider fragment range, written here before any read is downloaded:
- Reads: ENA FASTQ of the six Pol II (RPB1) MNase-ChIP-seq runs — DMSO rep 1–3: SRR22327380, SRR22327379,
  SRR22327378; 1-NM-PP1 rep 1–3: SRR22327377, SRR22327376, SRR22327375.
- Alignment: bowtie2 2.5.5 `--very-sensitive -X 1000 --no-mixed --no-discordant` to GRCh38 (no-alt analysis set,
  prebuilt genome-idx index); Drosophila spike-in reads are not used (they do not pass MAPQ ≥ 10 on hg38 except
  by chance); samtools sort / index. No duplicate removal (MNase cuts produce identical fragments genuinely);
  the duplicate share is reported.
- ProfileShape BAM input: fragment midpoints, fragments 20–120 bp (sub-nucleosomal; excludes mononucleosome-
  protected fragments ≈ 147 bp), MAPQ ≥ 10, proper pairs, duplicates kept.
- Regions: the 8038 regions already selected condition-blind (B1 cache), unchanged.
- S3 is checked on exact fragment counts in the core window (TSS ± 500 bp). If it fails again the case is
  final as "not evaluable"; no further change.
- Predictions P1–P3, A2 and the readout are unchanged. Midpoints of fragments up to 120 bp blur the PIC / pause
  separation (~80 bp), which works against P2; this is accepted, not compensated.

**B3 (2026-10-28, from an external review of the analysis code; no A2 result had been computed).** The A2
bootstrap re-assigned folds by draw position, so a chromosome drawn more than once had its copies in different
folds: the same genes entered training and test folds together (seed 30, first draw: 6 of 22 chromosomes
affected). Corrected in `benchmarks/positive_case.py` (`bootstrap_folds`): every copy keeps its chromosome's
original fold; a resample in which all draws fall in one fold is skipped. The point estimate of the gain was not
affected. Criterion and model unchanged.

**Outcome of B2 (2026-10-28, recorded after the A0 check; no A1 / A2 result computed before this entry).** All six
samples aligned (overall alignment 86.8–94.6 %). Exact fragments (20–120 bp, midpoint in TSS ± 500 bp) per region,
median per sample: DMSO 86, 95, 84; CDK7i 77, 82, 93; median over samples 85 < 100. **S3 not met; the
confirmatory case is final as "not evaluable"** (`benchmarks/results/positive_case_gse218269_bam`). Open item 1
stays open. Duplicate removal step of the alignment script used an invalid option (`samtools markdup -n`) and
produced no statistics; the duplicate share is computed afterwards from the BAMs (same position, mate position
and fragment length) and reported in the revision report. This does not affect the BAMs, which keep duplicates
by design.

**B4 (2026-10-28, written before running; exploratory, not a confirmation).** A1 and A2 are run once on the B2
cache with the S3 gate disabled, unchanged otherwise (same predictions, model, seeds, readouts), into
`benchmarks/results/positive_case_gse218269_bam_exploratory`, whose report is headed as exploratory with S3
failed. Whatever it shows, it is not counted as independent positive confirmation, is not used to choose a
replacement dataset, and is reported in full (including predictions not met). Its purpose is to see whether the
method's output points in the expected direction at this depth, as information for planning the next search.

**B5 (2026-10-28, diagnostic supplement after an external review of B4; written before any of it is run).**
Exploratory throughout: none of it replaces or updates the confirmatory verdict ("not evaluable"), none of it is
used to choose a dataset, a fragment range or a representation for a later claim, and every variant listed below
is reported, whatever it shows. Script `benchmarks/positive_case_diag.py`; output
`benchmarks/results/positive_case_gse218269_diag`. Code under study unchanged (`regional.py` hash recorded).

- **D1 · A2 ablation.** On the B2 cache, the same rows, folds (chromosome → fold as in `added_value`) and the same
  500 chromosome bootstrap draws for all models (seed = readout minutes, as in A2): base; base + signed shift;
  base + kernel; base + both. Reported per readout (15, 30, 60 min): CV R² of each, gains over base with 95 %
  intervals, and the paired difference full − (base + shift) and full − (base + kernel). Question: which shape
  statistic carries the added prediction; no criterion.
- **D2 · Per-region width and its replicate support.** Per region and sample, the mean absolute deviation of the
  mass-normalised core profile (the `mad_width` statistic). Reported: (a) the distribution (quantiles, share < 0)
  of the group difference CDK7i − DMSO among `shape_any`-called and not-called regions; (b) the share of regions
  where all three CDK7i widths lie below all three DMSO widths, and above (each 1/20 under exchangeable labels),
  in both region sets; (c) the subset behind the reported `shift` estimate: its selection condition and size, and
  the estimate's quantiles there. (d) Detectability (real data, insertion): in 500 randomly chosen regions
  (seed 9 700 000 + k), the three CDK7i profiles are narrowed about the background by f ∈ {0.95, 0.90, 0.80}
  (`insertion.insert`, kind `broaden<f>`); reported: the share of those regions called by `mad_width` and by
  `shape_kernel` (FDR 0.05) with and without the insertion, 3 draws per f. The aggregate IQR change of B4 is
  converted to the equivalent f for orientation. Question: does the aggregate narrowing correspond to a common
  per-region change, and could the per-region test have seen one of that size?
- **D3 · Processing and replicate sensitivity.** (a) From the same BAMs and regions, six representations:
  fragment range 20–70, 20–120 (= B2), 70–120 bp × midpoint or coverage; all other BAM options as in B2. For each:
  core-window depth (median over samples), A1 quantities (P1–P3 as defined, with `shift`-significant count) and
  the A2 30-min gain with interval. (b) On the B2 cache, each of the six 2 vs 3 designs that drop one replicate
  (DMSO_1, DMSO_2, DMSO_3, CDK7i_1, CDK7i_2, CDK7i_3): the same quantities. DMSO_1 and CDK7i_1 (flagged by the
  consistency module in B4) are named but not treated differently. Question: does the failure of the expected
  direction depend on the processing or on one replicate?

The diagnostic statistics on ranks of 10 label assignments (3 vs 3) are reported with their permutation p and
the note that 0.1 is the floor.
