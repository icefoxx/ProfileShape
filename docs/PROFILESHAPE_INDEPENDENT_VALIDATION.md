# ProfileSeg v3.1 — independent real-data validation: protocol

*Written 2026-10-27, before any read of this dataset was downloaded or inspected. Method frozen at the v3.1
code on branch `profileseg-v3` (regional.py as of commit `c3e17ec`). Results will be added in a separate
section below; changes to this protocol after data are seen will be listed as deviations.*

## Dataset

GEO GSE346229 (BioProject PRJNA1523976; released 2026; not used at any point during development; ATAC-seq,
an assay type not used during development). *S. cerevisiae*, paired-end, NextSeq 2000.

- Primary: Sth1 anchor-away (RSC catalytic subunit depletion), 0 h (rep1–4) versus 1 h (rep1–4).
- Secondary: Snf2 anchor-away (SWI/SNF), 0 h (rep1–3) versus 1 h (rep1–3).

Prior biological expectation (literature, not these data): RSC depletion narrows nucleosome-depleted regions —
the −1 and +1 nucleosomes move into the NDR — and reduces NDR accessibility genome-wide (Kubik et al. 2018 and
others). SWI/SNF depletion has more restricted effects.

## Processing (fixed)

bowtie2 `--very-sensitive -X 1000` to sacCer3; MAPQ ≥ 10, proper pairs, duplicates removed (samtools markdup);
Tn5 offset (+4 / −5). Fragment classes: sub-nucleosomal (≤ 120 bp; NDR accessibility) and nucleosomal
(140–200 bp; nucleosome positions). Signal = fragment midpoints (ProfileSeg BAM path, midpoint mode), 10 bp bins.

## Anchors (condition-blind)

Gene starts (ATG) from the sacCer3 RefSeq annotation, nuclear chromosomes, genes ≥ 500 bp apart from the
nearest other gene start on the same strand. From the nucleosomal-fragment midpoints of all eight Sth1 samples
pooled with equal weight (labels not used):

- +1 nucleosome dyad: maximum of the smoothed pooled density in [ATG − 200, ATG + 100] (oriented);
- −1 nucleosome dyad: maximum in [+1 − 450, +1 − 200];
- NDR centre: midpoint between −1 and +1.

Windows: ±300 bp around each anchor (oriented by gene strand). Analyses per anchor set: +1 and −1 with nucleosomal
fragments; NDR centre with sub-nucleosomal fragments.

## Analyses and criteria

ProfileSeg v3.1 defaults (`shape_any`, `shape_kernel`, `shift`; BH 5 %). diffSplice comparator fixed as in round 6:
edgeR QL diffSpliceDGE, support bins of the background-corrected pooled profile merged to 60 bp, gene-level F.

**A. Real-null splits.** Within-condition 2 vs 2 splits: 3 within Sth1 0 h, 3 within Sth1 1 h (6 splits), for
each anchor set. Criteria (per test, per anchor set):
- A1: mean share of regions called over the 6 splits ≤ 0.005 for `shape_any`;
- A2: no single split with more than 1 % of regions called by `shape_any`.
Reported: per-split calls for all tests and diffSplice, and the label-free-of-p-values replicate-pairing ratio.

**B. Fragment-level insertion on the null splits.** In each of the 6 splits, 20 % of regions (random, fixed seed)
are modified in the second group only, at the level of individual fragments before binning: translation of
nucleosomal-fragment midpoints by 20 or 40 bp (random sign per region); broadening (each midpoint jittered by
N(0, 30 bp)); amplitude-only (binomial thinning of fragments by 0.5, which keeps the shape). Two insertion sets
per split. Criteria for `shape_any`:
- B1: mean FDP over runs ≤ 0.075 for each translation and broadening insertion;
- B2: amplitude-only: mean share of regions called ≤ 0.01.
Reported, not judged: power of all tests; end-to-end power and FDP of diffSplice on the same counts; ProfileSeg
translation estimate ÷ true and coverage for consistent calls.

**C. Biology (Sth1 1 h vs 0 h, 4 vs 4; Snf2 secondary).** Predictions:
- C1: at +1 nucleosomes, among `shift` calls with a consistent translation, the majority move upstream (towards
  the NDR; negative in the oriented frame): share > 0.5, one-sided binomial p < 0.05;
- C2: at −1 nucleosomes, the majority of consistent calls move downstream (towards the NDR): share > 0.5,
  binomial p < 0.05;
- C3: Sth1 has more `shape_any` calls than Snf2 at +1 and −1 nucleosomes (share of regions; Snf2 is 3 vs 3,
  so Sth1 is also run on every 3 vs 3 subset of its replicates and the comparison uses the mean).
Reported descriptively: sub-nucleosomal NDR signal change (amplitude tests are not claimed), diffSplice calls and
overlap with ProfileSeg, replicate-structure diagnostics, depth.

A failure of any criterion is reported as such; the method is not changed in response to these data.

## Amendment before data (2026-10-27, before download)

The ProfileSeg BAM reader does not apply the Tn5 +4 / −5 offset. It moves fragment midpoints by 0.5 bp
(negligible at 10 bp bins) and lengths by 9 bp, so the offset is not applied and the size classes are defined on
the unshifted TLEN as written above (≤ 120 bp; 140–200 bp). Fragment-level insertion needs per-fragment
midpoints, so midpoints are extracted once per sample with the same filters as the BAM reader (MAPQ ≥ 10, proper
pairs, duplicates removed, each fragment counted once from its leftmost mate) and binned per experiment; the
unmodified counts are cross-checked against the ProfileSeg BAM path for one sample.
ProfileSeg configuration for every anchor set (fixed before data, by analogy with the earlier yeast settings):
flank 600 bp, bins 10 bp, core ±300 bp (the analysed window), background = mean of the outer 200 bp on each side,
oriented by gene strand. Fragments are extracted ±1300 bp around each gene start so that insertions of up to 40 bp
move fragments into and out of the window correctly.

---

## Results (2026-10-27; added after the run)

Run: commit `c3e17ec`, no uncommitted source changes, regional.py `87f0c6c86498c6a5`; protocol file as amended
before download (sha256 `a5234324972ab7af`; the first version, before the two amendments, was `5f2f19efe1626ff3`).
Outputs: `benchmarks/results/independent_gse346229/` (`INDEPENDENT_VALIDATION.md`, JSON, alignment summary).
Processing: 14 samples, alignment 87–96 %, 0.78–2.78 M deduplicated fragments per sample. The extracted
fragment counts were identical, bin by bin, to the ProfileSeg BAM path for the checked sample (400 regions).

| criterion | result |
|---|---|
| A1 / A2 real-null splits (`shape_any` share ≤ 0.005; no split > 1 %) | **pass** for all three anchor sets: no ProfileSeg call in any of the 18 null analyses (diffSplice: 0–0.07 % of regions) |
| B1 insertion FDP ≤ 0.075 | **pass** (`shape_any` mean FDP 0.000–0.055; the tolerance, not the nominal 0.05) |
| B2 amplitude-only calls ≤ 1 % | **pass** (0 calls) |
| C1 / C2 direction at +1 / −1 nucleosomes | **fail — not evaluable**: no region was called by any method, ProfileSeg or diffSplice |
| C3 Sth1 more calls than Snf2 | **fail** (0 vs 0) |

Power was low throughout: translations of 20 / 40 bp inserted at the fragment level were detected in 4–12 %
of regions by `shape_any` (`shift` 8–15 %; diffSplice 0–9 %), broadening in 0 %. The reason was flagged by
the depth diagnostic in every analysis: a median of 18 (+1), 24 (−1) and 44 (NDR, sub-nucleosomal) fragments
per region and sample in the ±300 bp windows, below the ~50 at which per-region shape tests have power, and
only about half of the regions testable. Translation estimates for the few consistent calls were accurate
(direction 1.00, estimate ÷ true 0.96–1.06; coverage 0.08–0.15).

What this validation shows, and what it does not:

- It supports, on an independent dataset of a new assay type, that the shape tests do not produce false calls
  on within-condition splits or under amplitude-only changes, and that FDP for inserted changes stays within the
  pre-specified tolerance (≤ 0.075; point estimates 0.000–0.055, i.e. one cell slightly above the nominal 0.05).
  With 12 runs per cell the FDP estimates are imprecise; neither a single high run nor these means establish
  or refute FDR control.
- It does not provide positive biological evidence: at this depth no method detected per-region changes. The
  protocol did not check depth before selecting the dataset, which in hindsight should have been a selection
  criterion.

### Exploratory, not pre-specified: aggregate replicate-level analysis

To ask whether the expected biology is present at all, per-sample meta-profiles over all regions were compared
between 1 h and 0 h with replicates as units (`AGGREGATE_EXPLORATORY.md`). Sth1 depletion moved the +1 dyad by
−1.9 bp (towards the NDR; Welch p = 0.07) and the −1 dyad by +1.0 bp (towards the NDR; p = 0.12); NDR
accessibility and the −1 to +1 distance did not change detectably. The direction matches the literature, but the
genome-wide effect at 1 h in these data is about 2 bp — far below what any per-region test can resolve at this
depth. The prior expectation written in the protocol (clear NDR narrowing) was too strong for this time point.

### Consequences

1. The calibration claims gain an independent, new-assay replication; the power and biology claims do not.
2. A depth requirement belongs in dataset selection for positive validation (e.g. a median of ≥ 100 fragments
   per region and sample in the analysed window, estimable from total depth before analysis).
3. A positive biological validation still needs a deeper dataset with a larger expected per-region effect.
