# ProfileShape — notes for Claude Code sessions

## What this is

Replicate-aware detection of differential signal-profile shapes at predefined genomic regions (formerly ProfileSeg in
the GenomeCanvas repository; moved here at v3.1). Package `src/profileshape`, CLI `profileshape`, tests in `tests/`,
evaluation scripts in `benchmarks/`, results with provenance in `benchmarks/results/`.

Start with `README.md`, then `docs/PROFILESHAPE_V3_REVIEW.md` (consolidated evidence and limits) and
`docs/DATA.md` (inputs). The full history of six review rounds is `docs/PROFILESHAPE_V3_REVIEW_HISTORY.md`.

## Working conventions

- Run tests with `python -m pytest -q tests/` (installed editable: `pip install -e ".[all]"`).
- Claims are kept narrow and tied to result files; every summary number should be traceable to a JSON under
  `benchmarks/results/`. Evaluation scripts record commit, uncommitted-diff hash and regional.py hash at run start.
- New evaluations: write the protocol and criteria before running; use fresh seeds; report failures as failures;
  do not tune on the evaluation seeds.
- Records kept as written: `benchmarks/results/**`, `docs/PROFILESHAPE_V3_REVIEW_HISTORY.md`,
  `docs/PROFILESHAPE_INDEPENDENT_VALIDATION.md` (its hash is recorded in its results).
- Commit and push only when asked. The GitHub repository `icefoxx/ProfileShape` is private.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Environment notes

- Input data: `inputs/` (git-ignored, ~10 GB); see `docs/DATA.md`.
- bowtie2, samtools, sra-tools: conda env `psbench` (`~/anaconda3/envs/psbench/bin`). R with edgeR on PATH.
- Network: use `/usr/bin/curl` for downloads (TLS-intercepting gateway; conda's curl and prefetch fail).
- The shell is zsh: an unquoted `$VAR` holding several arguments is not split; write commands out or use a script.

## Open items (as of v3.1)

1. Independent positive confirmation on a deeper dataset with a larger per-region effect (depth as a selection
   criterion; the ATAC attempt replicated calibration only).
2. Yeast MNase 2 vs 2 excess (associated with replicate structure); a replicate-structure diagnostic in the output.
3. Separating translations from other redistributions (R² does not); selection-adjusted bp estimates.
4. Per-class calibration with few regions; per-region absolute amplitude test.
5. Regenerate the graphical abstract (`docs/figures/`) with the new name and corrected footer.
