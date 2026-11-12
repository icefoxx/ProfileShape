# Candidate scenario · ZNF384 fusion versus wild type — data screening and pilot (written before the pilot)

*2026-10-28. Proposed by the user after the data phase of `PROFILESHAPE_APPLICATIONS_PLAN.md` ended. Question:
at ZNF384-bound regions whose signal amount does not change much, does the profile shape differ between the
TCF3-ZNF384 fusion and wild-type ZNF384? Status: screening only. Whether this scenario enters the paper (as a
replacement for scenario C, or as follow-up work) is the user's decision, taken after this record; the stop rule
of the applications plan is unchanged until then.*

## Prior expectation (before any data)

- The fusion keeps all ZNF384 zinc fingers and adds the TCF3 N-terminal activation domains (EP300/CBP
  recruitment). DNA-binding specificity is expected to be similar, so the shape of the **TF profile itself** at
  shared sites is expected to change little. Larger shape changes are more plausible in chromatin readouts
  downstream of recruitment (H3K27ac spreading, accessibility / nucleosome position).
- Possible TF-profile changes: choice of a different motif instance within a region (centroid shift),
  cooperative or clustered binding (width).
- Known traps (from earlier rounds): "no amount change" must not be defined as "edgeR not significant", and
  regions must not be selected with label-dependent amount statistics before the shape test (tests run on a
  condition-blind region set; amount classes are cross-tabulated afterwards, as in scenario C, with an equivalence
  bound for "unchanged"); fragment-length differences can mimic shape changes (B5: `shift` was representation-
  dependent); 2 vs 2 is not calibrated; core-window depth ≥ 100 fragments is necessary, not sufficient.

## Step 1 · Screening (GEO / SRA / ENA metadata, 2026-10-28; nothing downloaded)

Search: GEO `ZNF384` (37 series), SRA `ZNF384` excluding RNA-seq (183 runs, 19 studies), literature search.

| series | design | assay, reads | status |
|---|---|---|---|
| **GSE112560** (Mullighan lab; murine Arf⁻/⁻ pre-B, lentiviral HA-tagged constructs, anti-HA ChIP) | ZNF384 WT 3 (clones C10, C11, C6) · TCF3-ZNF384 2 (C16, C17) · TAF15-ZNF384 2; one input per construct | ChIP-seq, 50 bp single-end, 53–74 M reads; GEO tracks mm9, reads extended by SPP | the only fusion-versus-WT design with the same tag. **Development dataset of this project** (semi-real calibration, real templates, consistency evaluation; `consensus_top5000.bed` was derived from its peaks; a WT vs TCF3 config exists from the GenomeCanvas era). 3 vs 2 (10 label assignments); single-end; WT replicate 2 flagged by the consistency module (sample-wide shape offset 2.3× its mates, `consistency_v2`) |
| **GSE181617** (same lab, GSE181499 super-series) | TCF3-ZNF384 **+exon 8** (2 extra zinc fingers) 4 vs **−exon 8** 4, separate clones (C86–C93), input per sample | ChIP-seq anti-HA (Abcam ab9110), paired-end 2 × 51 bp, 38–47 M pairs; per-sample narrowPeak | not a development dataset. A different contrast (fusion isoforms that differ in the DNA-binding domain), 4 vs 4, paired-end. H3K27ac in the same series: WT 1, TCF3-ZNF384 2, TAF15-ZNF384 2 — no replicated WT, not testable |
| GSE181498 | H3K27ac HiChIP, WT 2 · TCF3 2 · TAF15 2 | HiChIP | not a 1-D profile assay |
| GSE197890 | human JIH-5 line and TCF3-ZNF384 PDX: ZNF384, H3K27ac, H3K4me3 CUT&Tag / ChIP, single samples; EP300-ZNF384 PDX ATAC 2 | mixed | no replicates per condition |
| GSE226400 (GSE211631) | 157 primary B-ALL ATAC-seq | ATAC | subtype labels (incl. ZNF384-rearranged) not in GEO metadata; patients are not exchangeable replicates; weakest design (application demonstration at best) |
| ENCODE ZNF384 (HEK293, HepG2, K562) | WT only | ChIP-seq | no fusion |

No dataset has a fusion-versus-WT contrast with ≥ 3 replicates in both groups, paired-end reads, and data not
used in development. The two usable designs answer neighbouring questions:

- **(i) GSE112560, fusion vs WT, 3 vs 2** — the question as asked, on development data, single-end, uncalibrated
  design size (3 vs 2 lies between the calibrated 3 vs 3 and the known 2 vs 2 excess). At best an application
  demonstration.
- **(ii) GSE181617, TCF3-ZNF384 +e8 vs −e8, 4 vs 4** — new data, paired-end, replicated; the two isoforms differ in
  two zinc fingers, which gives a stronger prior for a change in binding position or footprint than fusion vs WT.
  It does not answer fusion vs WT.

## Step 2 · Pilot depth (specification, written before any pilot count is looked at)

- **(ii) GSE181617**: the first 2 M read pairs of one +e8 HA (SRR15365234, replicate 1) and one −e8 HA run
  (SRR15365232, replicate 1), `benchmarks/pilot_depth.sh` against mm10 (bowtie2 --very-sensitive, -X 1000, proper
  pairs, MAPQ ≥ 10). Regions: summits of the two samples' GEO narrowPeak files (the same replicates), merged within
  200 bp, main chromosomes; `tf_chip` preset (flank 2 kb, bin 20 bp, core ± 1 kb); the condition-blind top half by
  pooled pilot fragments in the core. Full-run depth: region-level complexity extrapolation (amendment R5 of the
  selection protocol), linear scaling and library-wide complexity reported alongside.
- **(i) GSE112560**: depth from the local bigWigs (already on disk) at `consensus_top5000.bed`, core ± 1 kb,
  fragments = area / (count quantum × footprint) (`bigwig_count_unit`). Processed-track depth, stated as such.
- **Criterion**: median core-window fragments over pilot samples ≥ 100 (S3). Meeting it is necessary, not
  sufficient: a power simulation with a stated effect (as in step 3 of the selection protocol) would follow before
  any full download, and is not run here.

## Results

(recorded below after the pilot, unchanged above)

### Paused and moved to HPC (2026-11-12)

The pilot was stopped at the user's request on 2026-10-28 (the user holds this part of the analysis and the data
are on HPC). Nothing above was changed; no criterion has been evaluated; no grade is given.

State at the pause:

- **(i) GSE112560, processed-track depth (done, from the local bigWigs):** median core-window (± 1 kb) fragments at
  `consensus_top5000.bed`, area / (quantum × footprint):

  | sample | quantum | footprint (bp) | median | 25 % |
  |---|---|---|---|---|
  | WT C10 (GSM3073020) | 0.365 | 190 | 189 | 153 |
  | WT C11 (GSM3073021) | 0.399 | 210 | 187 | 152 |
  | WT C6 (GSM3073022) | 0.394 | 210 | 141 | 116 |
  | TAF15 C13 (GSM3073023) | 0.445 | 205 | 104 | 80 |
  | TAF15 C14 (GSM3073024) | 0.366 | 210 | 162 | 131 |
  | TCF3 C16 (GSM3073025) | 0.297 | 200 | 121 | 98 |
  | TCF3 C17 (GSM3073026) | 0.389 | 180 | 85 | 68 |

  Not a pass: the five constructs' median for WT vs TCF3 (141–189 vs 85–121) is from processed tracks (B1: their
  fragment selection need not match the analysis), at the 5000 strongest pooled regions (a favourable selection),
  and duplicate handling of the tracks is not stated (Picard marked, removal unclear). Recorded only. The number is
  printed by `gse112560()` in `benchmarks/pilot_znf384.py` (not written to a results JSON yet).
- **(ii) GSE181617 pilot (not completed):** SRR15365234 (+e8 rep 1) was aligned (2 M pairs; 97.0 % overall,
  73.6 % concordant unique; BAM in `inputs/pilot_znf384/`, its region counts **not looked at**); SRR15365232
  (−e8 rep 1) was stopped during alignment (FASTQ only). The GEO narrowPeak files have 21 685 (+e8, GSM5508181) and
  3 521 (−e8, GSM5508179) peaks: the region set is dominated by +e8 peaks, consistent with the study's report of
  increased +e8 binding; the condition-blind top-half rule stays as specified.
- GEO file names and titles disagree on replicate numbers in GSE181617 (e.g. GSM5508179 "−e8 HA Replicate 1",
  file `…-no-e8-rep2-HA-AB`); ENA maps SRR15365232 → GSM5508179 and SRR15365234 → GSM5508181. Use GSM / SRR, not
  titles.

To continue on HPC (as specified above; nothing re-tuned): run `benchmarks/pilot_depth.sh` for the two runs (mm10
bowtie2 index), put the two narrowPeak files as `PEAK_DIR/GSM5508181.narrowPeak.gz` and
`PEAK_DIR/GSM5508179.narrowPeak.gz`, then `PYTHONPATH=benchmarks python benchmarks/pilot_znf384.py PILOT_DIR PEAK_DIR
benchmarks/results/znf384_pilot/pilot.json`; record the result here under a dated heading. Part (i) needs the
GSE112560 bigWigs at `inputs/chip_gse112560/`. Local environment used: Python 3.11.5, numpy 1.26.4, scipy 1.17.1,
pysam 0.24.1, pyBigWig 0.3.25, profileshape 3.1.0 (editable). Next decisions after the pilot (the user's): power
simulation with a stated effect, and whether the scenario replaces C or stays follow-up work.

Execution protocol for the user's HPC data (readouts A TF profile, B H3K27ac, C accessibility / nucleosomes),
written 2026-11-12 before any of those analyses: `PROFILESHAPE_ZNF384_HPC_PROTOCOL.md`.
