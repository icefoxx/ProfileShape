# Stress tests of the pooled null (P1-F) — protocol

*Written 2026-10-28 before any run with the seeds below (9 600 000 + 10 000 × cell + 1000 × run + 100 × n; not
used before). Script `benchmarks/stress_null.py`; regional.py is not changed by this work.*

## Why

The calibration study (`benchmarks/calibration_boundary.py`) varied one factor at a time with at most 50 % of
regions changed. Real applications called far more (yeast double mutant 1433 / 1555 regions; PRO-seq 2 h
10 867 / 13 797), and the pooled null assumes regions are comparable and replicates exchangeable. These cells
probe where the whole-list FDR of `shape_kernel`, `shift` and `shape_any` breaks.

## Cells (600 regions as in the calibration baseline unless stated; 2 vs 2 and 3 vs 3; 100 runs each)

| factor | levels | what is changed |
|---|---|---|
| proportion_high | 0.7, 0.9, 1.0 | share of regions changed (half shifted +100 bp, half broadened ×1.5); at 1.0 only power is reported |
| same_direction | 0.5, 0.9 | changed regions are all shifted +100 bp (one direction, no broadening) |
| sample_effect | τ = 5, 10, 20 bp | every replicate receives its own sample-wide shift ~ N(0, τ²) in all regions (exchangeable, shared across regions); 20 % changed |
| batch_partial | undeclared, declared | 3 vs 3 with batches (A, A, B \| A, B, B), 2 vs 2 with (A, A \| A, B); batch B carries a 15 bp sample-wide shift; "declared" puts the batches in the manifest and permutes within batches |
| real_shapes | GSE112560 templates | shapes drawn from `inputs/caches/real_templates.npz` (multimodal, asymmetric), 20 % changed |
| tiling | 3 offsets | every region present three times with anchors at 0, +150 and +300 bp (overlapping windows); 20 % changed |

## Outcomes and classification (as in the calibration study)

Whole-list FDP per test (mean, 95 % Monte Carlo CI); EXCEEDS if the CI lies above 0.05, WATCH if the estimate is
above 0.075 without the CI excluding 0.05, ok otherwise; power; share of unchanged regions with p < 0.01.

## Consequences

No change to the method in this round. Cells classified EXCEEDS or WATCH are added to the documented limits
(review packet, guide) with the factor level at which they occur; "ok" is not proof of FDR control.

## Amendments

(none)
