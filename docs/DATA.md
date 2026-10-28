# Benchmark input data

Raw data and caches live in `inputs/` (git-ignored; about 10 GB). All are public; this file records where each
comes from, what was built from it and which benchmark reads it. Paths below are relative to the repository root.

## Datasets

| directory | source | content | used by |
|---|---|---|---|
| `inputs/yeast_mnase_gse280671/` | GEO GSE280671 (Tripplehorn et al., NAR 2025) | mononucleosome CPM bigWigs (WT, isw1Δ chd1Δ, rtf1 mutants, isw1Δ rtf1Δ, chd1Δ, isw1Δ; 2–4 replicates); `sacCer3_ncbiRefSeq.txt.gz`; `workspace/` = ATG-anchored profile cache `profiles_all.npz`, manifest, `orf_atg.bed`, configs | `positive_cases.py`, `insertion.py --build-yeast` |
| `inputs/yeast_mnase_bam/` | SRA runs of GSE280671 (WT, dKO), aligned with bowtie2 | BAMs (+ first-mate-only `.R1.bam` for MMDiff2), `sacCer3.fa` and bowtie2 index | `realbam.py`, `bw_vs_bam.py` checks, alignment of new yeast data |
| `inputs/proseq_gse233930/` | GEO GSE233930 (NELFB/PAF1/RTF1-FKBP12F36V HAP1, NAR 2024) | strand-specific PRO-seq bigWigs; `hg38_ncbiRefSeqSelect.txt.gz`; `proseq_profiles.npz` (NELFB DMSO / VHL 2 h / VHL 48 h, TSS ± 1 kb) | `positive_cases.py`, `insertion.py --build-proseq` |
| `inputs/chip_gse112560/` | GEO GSE112560 (HA ChIP-seq, ZNF384 WT / TAF15 / TCF3 fusions) | bigWigs; `semireal_config/` (config + manifest that build `caches/semireal_cache.npz`) | semi-real cells of `calibration_boundary.py` |
| `inputs/cutrun_gse291134_bam/`, `inputs/cutrun_gse291134_bigwig/` | GEO GSE291134 (A673 FLI1 CUT&RUN) | BAMs (chr1/chr2) and bigWigs | `cutrun_null.py`, `realbam_cutrun.py` |
| `inputs/atac_gse346229_bam/` | GEO GSE346229 (Sth1 / Snf2 anchor-away ATAC-seq; independent validation) | 14 deduplicated BAMs (bowtie2 `--very-sensitive -X 1000`, MAPQ ≥ 10, proper pairs, markdup) with logs | `independent_validation.py extract` |
| `inputs/histone_chip_gse260577/` | GEO GSE260577 (mouse liver, Sham vs CLP) | H3K27ac and H3K4me3 ChIP-seq CPM bigWigs, 4 replicates per group (assembly mm39 by the bigWig header, although the GEO text also names mm10) | `consistency_eval.py build` |
| `inputs/histone_cutrun_gse285245/` | GEO GSE285245 (mouse P14 CD8 T cells) | H3K27ac and H3K4me3 CUT&RUN bigWigs, Naive and Arm day 30, 4 replicates each (mm10) | `consistency_eval.py build` |
| `inputs/annotation_mouse/` | UCSC `ncbiRefSeqSelect` for mm10 and mm39 | TSS anchors for H3K4me3 | `consistency_eval.py build` |
| `inputs/annotation_tss/` | derived from public annotations and GSE112560 peaks | `consensus_top5000.bed`, `tss_1bp.bed`, `promoters_tss1kb.bed` | semi-real cache build |

## Caches (`inputs/caches/`)

| file | built from | by | read by |
|---|---|---|---|
| `semireal_cache.npz` | `chip_gse112560/semireal_config/config.yaml` | `semireal.build_cache` (re-build verified identical) | `calibration_boundary.py --semireal` |
| `real_templates.npz` | GSE112560 profiles (templates for simulation) | earlier template extraction | `shift_boundary.py --templates`, region tuning |
| `templates_2cell.npz` (+ `templates_2cell_source/`) | 2-cell IP (GSM5746912/3) | earlier template extraction | `validation.py --templates-2cell` (design V4) |
| `ins_yeast.npz` | yeast MNase bigWigs | `insertion.py --build-yeast` | `insertion.py`, `replicate_structure.py`, `replicate_holdout.py`, `batch_aware.py` |
| `ins_proseq.npz` | PRO-seq DMSO bigWigs (line-balanced) | `insertion.py --build-proseq` | `insertion.py` |
| `psv_cache.npz` | `atac_gse346229_bam/` | `independent_validation.py extract` | `independent_validation.py run / aggregate` |
| `cons_*.npz` | all datasets above, with region coordinates | `consistency_eval.py build` | `consistency_eval.py run` |
| `cutrun_profiles.npz` | CUT&RUN bigWigs | `cutrun_null.py` | `cutrun_null.py` |

## Typical commands

```bash
I=inputs
python benchmarks/positive_cases.py --yeast-cache $I/yeast_mnase_gse280671/workspace/profiles_all.npz \
    --yeast-manifest $I/yeast_mnase_gse280671/workspace/metadata/sample_manifest.tsv \
    --proseq-cache $I/proseq_gse233930/proseq_profiles.npz --out OUT/positive
PYTHONPATH=benchmarks python benchmarks/calibration_boundary.py --semireal $I/caches/semireal_cache.npz --reps 100 --out OUT/cal
python benchmarks/insertion.py --yeast-cache $I/caches/ins_yeast.npz --proseq-cache $I/caches/ins_proseq.npz --out OUT/ins
python benchmarks/independent_validation.py run $I/caches/psv_cache.npz OUT/independent
python benchmarks/shift_boundary.py --templates $I/caches/real_templates.npz --reps 30 --out OUT/shift
python benchmarks/diffsplice_compare.py --reps 50 --out OUT/dsp          # needs R with edgeR and jsonlite
```

Manifests inside `inputs/` contain absolute paths to the bigWigs on this machine; if `inputs/` moves, update them
(`*/metadata/sample_manifest.tsv`).

## Tools used to produce the BAMs

bowtie2 and samtools from the conda environment `psbench` (`~/anaconda3/envs/psbench`, also holds sra-tools);
R 4.x with edgeR (diffSplice comparisons). Downloads from NCBI need the system `/usr/bin/curl` on this network
(a TLS-intercepting gateway is trusted by the macOS keychain, not by conda's CA bundle or sra-tools' prefetch).
