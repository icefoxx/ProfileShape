# Consistency against user tolerance margins — evaluation protocol (v1)

*Written 2026-10-27 before any run with the seeds below. Code under evaluation: `consistency.py`
(`tolerance`, `offset_size`, `projection_part`, `verdict`) as of the run; hashes recorded in the results.
Development checks before this protocol used seeds 100–129 (simulation only). Changes after the first run go
under "Amendments".*

## Question

When a user fixes a tolerance margin before the analysis, how often does ProfileShape declare a replicate
(or a pair) **within** the margin although its true sample-wide difference is at or beyond the margin — the
error that matters most — and how often can it declare consistency when there is no difference? Do the
intervals behind the verdicts cover the truth?

Quantities: shift (bp, 90 % interval), stretch (% about the anchor, 90 % interval), offset (% of the aggregate
profile, bias-corrected, one-sided 95 % bounds). Verdict *within* = interval inside ±margin (shift, stretch) or
upper bound below margin (offset).

## Arms

**A. Simulation** (`profileshape.simulate`, 2000 regions, parametric peaks). Seeds 9 990 000 + 1000 × cell + run,
100 runs per cell. n ∈ {2, 3, 4} × depth ∈ {12, 40, 160}. Replicate 1 is altered in every region: shift
d ∈ {0, 5, 10, 15, 20} bp, or stretch (broadening as `insertion.insert`) f ∈ {2.5, 5, 7.5, 10} %. Margins:
shift 10 bp, stretch 5 %, offset 5 %. Truth: d and f as inserted; offset truth = the module's offset estimate on
the noise-free profiles (`simulate(..., noiseless=True)`, same seed, same insertion).

**B. Semi-real** (real replicate noise). The 11 datasets and groups of `consistency_v2` (caches
`inputs/caches/cons_*.npz`), replicate labels permuted within blocks of 10 neighbouring regions (as in
`consistency_eval.py`), 20 shuffles per group and insertion, seeds 9 995 000 + 1000 × dataset + 100 × group +
shuffle. Pseudo-replicate 1 altered in every region: shift d ∈ {0, 0.5, 1, 1.5, 2} × M with M = max(10 bp,
bin size), or stretch f ∈ {2.5, 5, 7.5, 10} %. Margins: shift M, stretch 5 %. Offset is not judged here (its
truth on real shapes is not known).

**C. Real data as is** (descriptive): estimates and bounds of every replicate in the 11 datasets — reference
ranges, no truth.

## Criteria (fixed before running)

| id | arm | criterion |
|---|---|---|
| M1 | A | wrong *within*: at truth = margin ≤ 0.10 in every cell; at truth ≥ 1.5 × margin ≤ 0.05 in every cell (shift, stretch, offset) |
| M2 | B | wrong *within*: at truth = margin ≤ 3 / 20; at truth ≥ 1.5 × margin ≤ 1 / 20, every group (shift, stretch) |
| M3 | A | coverage: shift and stretch 90 % intervals ≥ 0.85 per cell and ≥ 0.88 pooled; offset upper bound ≥ truth in ≥ 0.93 pooled |
| M4 | B | coverage of shift and stretch 90 % intervals ≥ 0.85 pooled over groups |

Reported without a criterion: *within* at truth 0 (how often consistency can be shown) by n and depth;
*inconclusive* rates; cross-talk (stretch read under a pure shift and vice versa); arm C ranges.

Consequences: if M1 or M2 fails for a quantity, its verdict is not offered (only the estimate and interval);
if M3 or M4 fails, the interval is labelled as failing its coverage check in the report.

## Outcome (run margin_v1, 2026-10-28; recorded after the run)

M1 pass; M2 **fail** (stretch only: 24 real-noise groups exceed the wrong-*within* limit; shift: 0 groups);
M3 **fail** (per-cell coverage down to 0.77 for stretch and 0.82 for shift; offset upper bound 0.901 < 0.93);
M4 **fail** (pooled real-noise coverage 0.53). Consequences applied as specified: the stretch verdict is
removed (estimate and interval shown, described as failing); shift and offset verdicts are kept; all three
intervals are labelled as failing their coverage checks in report 07 and the guide.

## Amendments

(none)
