# ProfileShape

**Replicate-aware detection of differential signal-profile shapes at predefined genomic regions.**

ProfileShape compares signal profiles (ChIP-seq, ATAC-seq, MNase-seq, PRO-seq, CUT&RUN coverage) around predefined
anchors — TSSs, summits, nucleosome dyads — between two conditions with biological replicates. It asks whether the
*spatial distribution* of signal at a region changed (broadening, narrowing, redistribution between peaks, array
disorder, and, where it describes the change, a translation in bp), separately from how much signal there is.

Formerly **ProfileSeg**, developed inside the [GenomeCanvas](https://github.com/icefoxx/genomeCanvas) repository
(branch `profileseg-v3`, tag `profileseg-v3.1`). "Seg" never meant segmentation; the new name says what the method does.

## What it does

- **Replicates are the unit of inference.** With 2–3 replicates per group there are only 3 or 10 label assignments per
  region, so per-region p-values come from a permutation null pooled over regions, stratified by a label-invariant
  region covariate when enough regions are available.
- **Shape tests** on background-corrected, mass-normalised profiles: `shape_kernel` (kernel distance between condition
  means), `shift` (directional projection test), `shape_any` (their per-region Cauchy combination, one BH list).
- **Translation description, only where it fits.** A direction and bp value are reported when a cross-correlation
  shift is consistent with the projection test, explains most of the difference and is unambiguous; otherwise the
  region is reported as `Changed` with the reason (`translation_flag`). The bp value is a descriptive effect size.
- **Within-group consistency** (`profileshape consistency`): sample distance matrices, sample-wide offsets
  of one replicate from its group-mates, sub-group structure among replicates, divergent regions. Descriptive QC;
  not an equivalence test.
- **Diagnostics**: depth per region (exact for BAM, estimated for bigWig), replicate influence, batch confounding
  (tests refused when condition and batch coincide), fragment-length QC for BAM.
- **Inputs**: bigWig or BAM (fragment midpoints, coverage, read ends; MAPQ / duplicate / proper-pair / size filters),
  or a precomputed matrix; a BED of anchors; a sample manifest.

## Install

```bash
pip install -e ".[all]"        # numpy, scipy, matplotlib, PyYAML + pyBigWig, pysam
```

## Quick start

```bash
profileshape example --out ps_example     # synthetic data, full report
profileshape init --out config.yaml       # annotated template; edit it and the manifest
profileshape run --config config.yaml
profileshape power --replicates 2 3 4     # power guidance
```

User guide: [`docs/PROFILESHAPE.md`](docs/PROFILESHAPE.md).

## Status (v3.1) and evidence

The consolidated evidence, with what is and is not supported, is in
[`docs/PROFILESHAPE_V3_REVIEW.md`](docs/PROFILESHAPE_V3_REVIEW.md). In short:

- Whole-list FDR of the shape tests near or below 0.05 across the simulated and semi-real scenarios tested, including a
  pre-specified confirmatory run on new seeds — empirical, not a proven bound.
- Per-class calibration holds with enough regions (≥ 3000 at 2 vs 2, ≥ 1000 at 3 vs 3 in the tested mix), not with few.
- Real data: PRO-seq within nominal; yeast MNase 2 vs 2 above nominal (associated with replicate structure). An
  independent ATAC-seq dataset replicated calibration but was too shallow for positive confirmation.
- More sensitive than edgeR diffSpliceDGE for medium translations of narrow peaks in simulation; no proof of general
  superiority. Translation estimates are descriptive; `translation_consistent` is not proof of a translation.
- Per-region absolute amplitude inference is not claimed.

## Repository layout

| path | content |
|---|---|
| `src/profileshape/` | the package |
| `tests/` | unit tests |
| `benchmarks/` | evaluation scripts (calibration, insertion, comparison, validation) |
| `benchmarks/results/` | result files with provenance, produced in the GenomeCanvas repository before the move (see below) |
| `docs/` | user guide, review packet, calibration, comparison, validation, plans |

### Records kept as written

Result files under `benchmarks/results/`, `docs/PROFILESHAPE_V3_REVIEW_HISTORY.md` and
`docs/PROFILESHAPE_INDEPENDENT_VALIDATION.md` (a protocol whose hash is recorded in its results) are kept exactly as
written. They refer to the module `genomecanvas.profileseg`, to commits of the GenomeCanvas repository (e.g. `5a5be13`,
`ad1e875`, `c3e17ec`) and to script names with a `profileseg_` prefix. The code moved unchanged apart from names; the
statistical core (`regional.py`) is functionally identical to GenomeCanvas regional.py `87f0c6c86498c6a5`.

The graphical abstract in `docs/figures/` still shows the old name and needs to be regenerated.

## Licence

MIT (see `LICENSE`). `src/profileshape/bigwig_io.py` is taken from GenomeCanvas (MIT).
