# ProfileShape — related work and positioning

Every citation below was checked against the publisher, PubMed/PMC or arXiv record (October 2026).
For each one: what it does, and where ProfileShape agrees, borrows or differs.

## 1. Differential binding (amplitude / counts)

| work | what it does | relation to ProfileShape |
|---|---|---|
| Ross-Innes *et al.* (2012) *Nature* 481:389–393, doi:10.1038/nature10730; DiffBind package (Stark & Brown) | consensus peak sets, read counts per peak, count-model tests (edgeR/DESeq2) | answers *which peaks change in total signal*; does not look at shape. ProfileShape's per-region amplitude tests address the same question on coverage tracks; DiffBind on BAM counts remains the reference when BAMs are available |
| Lun & Smyth (2016) *Nucleic Acids Res* 44(5):e45, doi:10.1093/nar/gkv1191 (csaw); Lun & Smyth (2014) *Nucleic Acids Res* (peaks and windows, error-rate control) | sliding-window counts, edgeR tests, FDR controlled over **regions** not windows | the region-level FDR principle is the one ProfileShape follows (one test per region per metric). csaw is count-based and peak-shape-agnostic |
| Robinson & Oshlack (2010) *Genome Biol* 11:R25 (TMM) | scaling factors from trimmed log-ratios, assuming most features do not change | ProfileShape's `log2_signal_relative` uses the same assumption with a different estimator (half-sample mode + iterative exclusion). The comparison benchmark shows both control FDR similarly; window counts that include background make TMM-based tests slightly anti-conservative under uniform IP-efficiency changes |
| Tu *et al.* (2021) *Genome Res* 31(1):131–145 (MAnorm2) | hierarchical normalisation, empirical-Bayes modelling of within-group variability | closest in spirit to the per-region amplitude part (variance borrowing across regions); count/peak-based and amplitude-only |
| Allhoff *et al.* (2016) *Nucleic Acids Res* 44(20):e153 (THOR) | HMM differential peak calling with replicates, replicate normalisation (housekeeping genes) | segmentation-based; ProfileShape v2 deliberately does *not* segment (planned v3) |
| Steinhauser *et al.* (2016) *Brief Bioinform* 17(6):953–966, doi:10.1093/bib/bbv110; Eder & Grebien (2022) *Genome Biol* 23, doi:10.1186/s13059-022-02686-y | benchmarks of differential ChIP-seq tools (33 tools and approaches in Eder & Grebien) | Eder & Grebien report that tool performance depends strongly on **peak size and shape** and on the regulatory scenario, which motivates benchmarking on the user's own profile shapes (`benchmark --templates-config`, `power --templates-config`) |

## 2. Shape and position changes

| work | what it does | relation to ProfileShape |
|---|---|---|
| Schweikert *et al.* (2013) *BMC Genomics* 14:826, doi:10.1186/1471-2164-14-826 (MMDiff; Bioconductor MMDiff2) | per-peak kernel two-sample test (maximum mean discrepancy on read midpoints); null from replicate-vs-replicate comparisons binned by coverage; BH | **the closest prior art for shape**. Its authors found the plain MMD bootstrap rejected almost every replicate-vs-replicate comparison, so they built a replicate-derived null; this is the same lesson as ProfileShape's (replicates, not reads or regions, are the unit of inference). Differences: MMDiff gives one omnibus distance per peak and needs BAMs; ProfileShape reports interpretable components (centroid, spread, asymmetry, amplitude kept separate) at both aggregate and per-region level on coverage tracks. MMDiff reports poor sensitivity to count changes by design, i.e. it too separates shape from amplitude |
| Chen *et al.* (2013) *Genome Res* 23(2):341–351 (DANPOS) | nucleosome position shift, fuzziness and occupancy changes from MNase-seq | same decomposition (position / spread / occupancy ≈ centroid / width / amplitude) for nucleosomes; ProfileShape applies it to anchored windows of any assay |
| Weiß & Friedel (2023) *Bioinform Adv* 3(1):vbad085, doi:10.1093/bioadv/vbad085, PMID 37456509 (RegCFinder) | in user-defined windows (promoters, gene bodies, peaks; any short-read assay, BAM), finds sub-regions whose *relative* read density differs between two conditions: per-position densities normalised within the window and averaged over replicates, all maximal scoring subsequences of the difference (Ruzzo–Tompa), a filter against position-permuted differences, then DEXSeq differential usage of the found sub-regions versus the rest of the window; ≥ 2 replicates per condition; applied to PRO-seq and Pol II ChIP-seq | the closest prior work: same target (within-region redistribution, amplitude removed by within-window normalisation), predefined windows, replicates. It **localises** the change (sub-region coordinates), which ProfileShape does not. It differs in inference: segmentation ignores replicate variance, the permutation is over positions (not replicate labels), and sub-regions are selected and tested on the same data; no simulation or FDR calibration is reported. Not run here (comparison planned; `docs/PROFILESHAPE_NEXT_PLAN.md` P0-C) |

## 3. Functional data / local testing

| work | what it does | relation to ProfileShape |
|---|---|---|
| Pini & Vantini (2017) *J Nonparametr Stat* 29(2):407–424, doi:10.1080/10485252.2017.1306627 (interval-wise testing) | permutation tests over all intervals of the domain; adjusted p-value function with interval-wise error control | the right tool for *where along the profile* a difference sits (planned v3 "local inference") |
| Cremona *et al.* (2018) *Bioinformatics* 34:2289–2291, doi:10.1093/bioinformatics/bty090 (IWTomics) | IWT for groups of genomic regions at multiple scales | its unit is the **region (curve)**, comparing two region sets. For a WT-vs-KO design on the same regions, treating regions as units ignores replicate variability; the comparison benchmark's "regions as units" baseline (same unit choice) calls amplitude-only changes as shape changes in 100 % of runs and global rescaling in 12 % |
| Ramírez *et al.* (2016) *Nucleic Acids Res* 44:W160–W165, doi:10.1093/nar/gkw257 (deepTools2) | computeMatrix / plotProfile / plotHeatmap | ProfileShape reads computeMatrix output directly; deepTools plots are descriptive and have no replicate-level inference |

## 4. Statistical components used

| component | reference | used for |
|---|---|---|
| moderated t, prior fitted by moments of log variances | Smyth (2004) *Stat Appl Genet Mol Biol* 3(1):Art. 3, doi:10.2202/1544-6115.1027 (limma) | per-region tests; numerically identical to `limma::eBayes` (unit test) |
| robust hyperparameters | Phipson *et al.* (2016) *Ann Appl Stat* 10(2):946–963, doi:10.1214/16-AOAS920 | **not yet used**: candidate fix for the mild tail inflation of per-region tests on low-SNR real profiles |
| false discovery rate | Benjamini & Hochberg (1995) *J R Stat Soc B* 57(1):289–300 | per-region multiplicity |
| sequentially rejective Bonferroni | Holm (1979) *Scand J Stat* 6(2):65–70 | secondary aggregate metrics |
| exact permutation p-values | Phipson & Smyth (2010) *Stat Appl Genet Mol Biol* 9(1):Art. 39 | reporting the smallest attainable p (0.1 for 3 vs 3) instead of a misleading p |
| half-sample mode | Bickel & Frühwirth (2006) *Comput Stat Data Anal* 50(12):3500–3530 | per-sample offsets for relative amplitude (robust to one-sided changes) |
| hierarchical bootstrap | Saravanan, Berman & Sober (2020) *Neurons Behav Data Anal Theory* 3(5) | resampling replicates then regions (descriptive with ≤ 3 replicates) |
| pseudoreplication | Lazic (2010) *BMC Neurosci* 11:5, doi:10.1186/1471-2202-11-5 | rationale for replicates (not regions or reads) as the inferential unit |
| spike-in normalisation | Orlando *et al.* (2014) *Cell Rep* 9(3):1163–1170, doi:10.1016/j.celrep.2014.10.018 (ChIP-Rx) | the only reliable way to interpret *absolute* global amplitude changes; ProfileShape flags but cannot resolve them |

## 5. What is new, and what is not

**Not new:** empirical-Bayes per-region testing (limma, MAnorm2), region-level FDR (csaw), TMM-type relative
normalisation, the position / spread / occupancy decomposition (DANPOS), kernel-based shape testing with a
replicate-derived null (MMDiff), and testing within-window redistribution of signal at predefined windows with
replicates (RegCFinder, which also localises the change). Window-level count methods (csaw) can detect local
redistribution too; they are not "total-only" methods.

**What ProfileShape contributes** — a combination of inference details, not a new statistical theory:
replicate-supported screening of spatial signal redistribution at predefined regions with 2–3 replicates, via
(i) a label-invariant stabilising constant (a per-assignment one makes pooled-permutation shape tests
anti-conservative when many regions change), (ii) a permutation null pooled over regions within strata of a
label-invariant covariate, (iii) a per-region combination of an omnibus and a directional test into one list,
on coverage tracks (bigWig) as well as BAM, with replicate-consistency diagnostics and a constrained,
descriptive shift (bp) reported only where a translation describes the change.

What this rests on: inference with 2 replicates per group is possible only through the cross-region
assumption of the pooled null (regions comparable within strata), which is studied empirically
(`PROFILESHAPE_CALIBRATION.md`), not guaranteed, and was violated on yeast MNase 2 vs 2. Each element above
needs an ablation in the paper (`PROFILESHAPE_MANUSCRIPT_PLAN.md`).
