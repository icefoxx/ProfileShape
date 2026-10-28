"""Within-group profile consistency: are the replicates of each condition alike?

The differential tests ask whether two conditions differ; this module asks whether the samples of one
condition agree with each other, in the shape of their profiles and separately in amplitude. It uses
no condition contrast and runs on every condition in the manifest. Outputs are descriptive QC:

* sample × sample distance matrices (shape: per-region kernel distance between mass-normalised
  profiles, median over regions; meta-profile: kernel distance between the samples' mean profiles;
  amplitude: robust spread of per-region log2 signal ratios);
* per-sample systematic offset from the mean of its group-mates (shared across regions; linear in the
  signal so that depth alone does not produce it), flagged when significant and ≥ 2× the group-mates'
  own offsets (≥ 3 replicates; with two, the offset of the pair is reported but not attributed);
* replicate structure: for groups of ≥ 4, every split into two halves of ≥ 2 samples, tested by the
  per-region contrast of across-half and within-half pair distances (unequal depth cancels), with the
  across / within ratio of meta-profile distances as effect size;
* per-region within-group spread, scored against regions of similar signal mass, and for ≥ 3
  replicates the sample that drives each divergent region.

None of this is a test of equivalence: a low score means no inconsistency was seen at this depth, not
that the replicates are shown to be the same. Distances grow with sampling noise, so a shallower
sample sits further from the consensus; depth is reported next to every sample score.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

from . import metrics as M
from .data import Config, Profiles

# Same kernel as the differential shape test (regional.KERNEL_SIGMA_BP), so distances are on the scale
# the tests see. Kept as a separate constant so that regional.py stays unchanged.
KERNEL_SIGMA_BP = 60.0
# Sample offsets and sub-group splits are flagged at this family-wise level per group (Bonferroni over
# samples or splits). The region flag is a heuristic ranking cut-off, not a test.
FLAG_ALPHA = 0.01
# With thousands of regions, real replicates always show a significant offset (yeast MNase, PRO-seq: nearly
# every sample), so a sample is flagged only when its offset is also this many times the median offset of
# its group-mates: it stands out from its group, not merely differs as replicates do.
OFFSET_FLAG_RATIO = 2.0
REGION_Z_FLAG = 3.0  # robust z of log within-group spread, within signal-mass strata
MIN_REGIONS_PER_STRATUM = 200
# Neighbouring windows can overlap and share signal, so variances over regions are computed from blocks of
# this many consecutive regions along the genome (batch means) instead of treating regions as independent.
BLOCK_REGIONS = 10
# A region is used only if every sample's support mass is at least this share of that sample's median: a
# near-empty window divided by its tiny mass gives values up to ~1e16 (ATAC, evaluation v1) that dominate any
# mean over regions. Label-free; such windows carry no shape information.
MIN_RELATIVE_MASS = 0.1
# Per-region translation-equivalent shift thresholds (bp) for the descriptive shares F_δ.
SHIFT_THRESHOLDS = (10, 20, 50)


@dataclass
class Consistency:
    samples: list[str]
    conditions: list[str]
    groups: dict  # condition → sample indices
    regions: np.ndarray  # indices of testable regions
    shape_distance: np.ndarray  # samples × samples
    meta_distance: np.ndarray
    amplitude_distance: np.ndarray
    sample_rows: list = field(default_factory=list)
    partition_rows: list = field(default_factory=list)
    region_rows: list = field(default_factory=list)
    group_rows: list = field(default_factory=list)
    flags: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def core_profiles(profiles: Profiles, config: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Background-corrected core-window signal of every sample (as in ``analysis.analyse``) and the
    regions without missing bins."""
    v = profiles.values
    complete = np.all(np.isfinite(v), axis=(0, 2))
    bg = M.background(v, profiles.edges, config.background_width, config.background_estimator)
    core = M.core_mask(profiles.edges, config.core_flank)
    e = profiles.edges
    core_edges = np.concatenate([e[:-1][core], [e[1:][core][-1]]])
    return (v - bg[..., None])[..., core], core_edges, complete


def shape_features(y: np.ndarray, centers: np.ndarray, sigma_bp: float = KERNEL_SIGMA_BP):
    """Mass-normalised profiles in the kernel feature space of ``regional.kernel_shape_test``.

    Same normalisation as the test: signed background-corrected profile divided by its mass inside
    the condition-blind support of the pooled profile; ‖φ_i − φ_j‖² is the kernel distance. Returns
    φ (samples × regions × rank), the support mass, the regions whose mass is at least
    ``MIN_RELATIVE_MASS`` of the sample's median in every sample, and the feature map L (bins × rank),
    so that φ = p L.
    """
    y = np.nan_to_num(y)
    pooled = y.mean(0)
    support = pooled > 0.2 * pooled.max(axis=1, keepdims=True)
    mass = (y * support[None]).sum(-1)
    pos = np.where(mass > 0, mass, np.nan)
    floor = MIN_RELATIVE_MASS * np.nanmedian(pos, axis=1, keepdims=True)
    ok = np.all(mass > np.nan_to_num(floor, nan=np.inf), axis=0)
    p = y[:, ok] / mass[:, ok, None]
    K = np.exp(-0.5 * ((centers[:, None] - centers[None, :]) / sigma_bp) ** 2)
    lam, V = np.linalg.eigh(K)
    keep = lam > 1e-10 * lam.max()
    L = V[:, keep] * np.sqrt(lam[keep])
    return p @ L, mass[:, ok], ok, L


def _linear_features(y: np.ndarray, mass: np.ndarray, L: np.ndarray, centers: np.ndarray):
    """Profiles scaled linearly (see :func:`linear_offsets`) in feature space, with the unit amplitude
    direction u, the aggregate norm and the translation direction t (−d/dx of the aggregate profile,
    orthogonal to u; an offset δ·t is a shift by δ bp towards higher coordinates)."""
    pooled = mass.mean(0)
    x = y / (pooled + np.median(pooled))[None, :, None]
    x = x / np.median(mass / pooled[None], axis=1)[:, None, None]
    agg_bins = x.mean((0, 1))
    agg = agg_bins @ L
    u = agg / max(np.linalg.norm(agg), 1e-300)
    t = -np.gradient(agg_bins, centers) @ L
    t = t - (t @ u) * u
    return x @ L, u, float(np.linalg.norm(agg)), t


def _remove(delta: np.ndarray, u: np.ndarray) -> np.ndarray:
    return delta - (delta @ u)[..., None] * u


def linear_offsets(y: np.ndarray, mass: np.ndarray, L: np.ndarray, centers: np.ndarray | None = None):
    """Per-region shape offset of each sample from the mean of the others, linear in the signal.

    Dividing each sample's profile by its own noisy mass (as the distances do) biases the normalised
    profile by an amount that depends on depth, the same in every region, so a shallower replicate would
    look systematically offset. Here regions are scaled by a sample-independent pooled mass (shrunk by
    its median, as in the aggregate analysis) and each sample by one global factor; the component along
    the group's aggregate profile (amplitude) is removed. Returns the offsets (samples × regions × rank),
    the norm of the aggregate profile on the same scale and the translation direction.
    """
    n = len(y)
    psi, u, norm, t = _linear_features(y, mass, L, np.arange(y.shape[-1], dtype=float) if centers is None else centers)
    tot = psi.sum(0)
    delta = np.stack([psi[k] - (tot - psi[k]) / (n - 1) for k in range(n)])
    return _remove(delta, u), norm, t


def region_blocks(regions: list[dict], ridx: np.ndarray, size: int = BLOCK_REGIONS) -> np.ndarray:
    """Block label per used region: consecutive regions along each chromosome, ``size`` per block."""
    order = sorted(range(len(ridx)), key=lambda k: (str(regions[ridx[k]].get("chrom", "")),
                                                    int(regions[ridx[k]].get("start", 0))))
    blocks = np.empty(len(ridx), int)
    b, count, chrom = -1, size, None
    for k in order:
        c = str(regions[ridx[k]].get("chrom", ""))
        if c != chrom or count == size:
            b, count, chrom = b + 1, 0, c
        blocks[k] = b
        count += 1
    return blocks


def _block_cov(x: np.ndarray, blocks: np.ndarray | None) -> np.ndarray:
    """R × Var(mean over regions) of x (regions × features), from block sums when ``blocks`` is given
    (dependence within blocks allowed), otherwise the ordinary covariance."""
    x = x.reshape(len(x), -1)
    if blocks is None:
        return np.atleast_2d(np.cov(x, rowvar=False))
    m = x.mean(0)
    ids = np.unique(blocks, return_inverse=True)[1]
    B = ids.max() + 1
    sums = np.zeros((B, x.shape[1]))
    np.add.at(sums, ids, x)
    dev = sums - np.bincount(ids, minlength=B)[:, None] * m
    return dev.T @ dev / len(x) * B / max(B - 1, 1)


def translation_part(delta: np.ndarray, t: np.ndarray, confidence: float = 0.95,
                     blocks: np.ndarray | None = None) -> dict:
    """Translation-equivalent size (bp) of per-region offsets (regions × features): per region the
    least-squares δ in delta ≈ δ·t, then the mean over regions with a normal interval from their spread
    (block variance when ``blocks`` is given). ``share`` is the part of the mean offset that this common
    translation explains (1: a pure shift; near 0: broadening, redistribution or region-specific changes).
    The projection is linear, so it is attenuated for shifts that are large relative to the peak.
    ``per_region`` keeps the per-region values (noisy: see :func:`shift_shares`)."""
    tt = float(t @ t)
    if tt <= 0:
        return dict(bp=float("nan"), low=float("nan"), high=float("nan"), share=float("nan"),
                    per_region=np.full(len(delta), np.nan))
    d = delta @ t / tt
    bp, se = float(d.mean()), float(np.sqrt(_block_cov(d, blocks)[0, 0] / len(d)))
    zq = stats.norm.ppf(0.5 + confidence / 2)
    m = delta.mean(0)
    mm = float(m @ m)
    return dict(bp=bp, low=bp - zq * se, high=bp + zq * se, share=bp * bp * tt / mm if mm > 0 else float("nan"),
                per_region=d)


def shift_shares(d: np.ndarray, thresholds=SHIFT_THRESHOLDS) -> dict:
    """Descriptive: median |per-region shift| and the share of regions whose estimated shift exceeds each
    threshold. Includes noise: in simulation with no difference at all ~40 % of regions exceed 20 bp, so
    these numbers are read against the group-mates' own values (``*_mates``), not on their own."""
    a = np.abs(d[np.isfinite(d)])
    out = {"median_abs_shift_bp": float(np.median(a)) if a.size else float("nan")}
    out.update({f"F{x}": float(np.mean(a > x)) if a.size else float("nan") for x in thresholds})
    return out


def changed_share(p: np.ndarray, lam: float = 0.5, boot: int = 200, seed: int = 0, fdr: float = 0.05,
                  blocks: np.ndarray | None = None) -> dict:
    """Experimental. Share of regions that differ between two halves of a group, π̂₁ = 1 − π̂₀ (Storey
    2002, λ = ``lam``), with a percentile interval from resampling blocks of regions, and the number of
    regions called at BH ``fdr`` (a lower bound). π̂₁ is biased downwards where per-region power is low
    (weak or broad differences): it estimates the share of *detectably* differing regions. Used only for
    splits into equal halves, where unequal depth cancels; a one-sample-versus-rest version was depth-confounded
    (simulation: π̂₁ ≈ 0.38 with no difference and one replicate at 30 % depth) and is not reported."""
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    blocks = np.arange(len(p)) if blocks is None else np.asarray(blocks)
    p, blocks = p[ok], blocks[ok]
    if p.size < 20:
        return dict(share=float("nan"), low=float("nan"), high=float("nan"), called=0, tested=int(p.size))
    est = lambda q: float(np.clip(1 - np.mean(q > lam) / (1 - lam), 0, 1))  # noqa: E731
    rng = np.random.default_rng(seed)
    ub, inv = np.unique(blocks, return_inverse=True)
    members = [np.flatnonzero(inv == b) for b in range(len(ub))]
    bs = [est(p[np.concatenate([members[b] for b in rng.integers(0, len(ub), len(ub))])]) for _ in range(boot)]
    from .regional import bh
    return dict(share=est(p), low=float(np.quantile(bs, 0.025)), high=float(np.quantile(bs, 0.975)),
                called=int(np.sum(bh(p) < fdr)), tested=int(p.size))


def region_test_p(y: np.ndarray, first: list, second: list, centers: np.ndarray, oriented: bool) -> np.ndarray:
    """Per-region p-values for ``first`` vs ``second`` (sample indices into y) from the differential shape
    tests (shape_kernel, and for oriented anchors shift, combined per region as in shape_any), with the
    same null pooled over regions and label assignments."""
    from .regional import (KERNEL_S0_QUANTILE, TEMPLATE_SMOOTH, cauchy_combination, kernel_shape_test,
                           shift_score_test)
    yo = y[list(first) + list(second)]
    ref, alt = np.arange(len(first)), np.arange(len(first), len(yo))
    pk = kernel_shape_test(yo, np.ones(yo.shape[1:], bool), ref, alt, KERNEL_SIGMA_BP, KERNEL_S0_QUANTILE, centers)
    if not oriented:
        return pk
    ps = shift_score_test(yo, ref, alt, centers, smooth=TEMPLATE_SMOOTH)["p"]
    return cauchy_combination(pk, ps)


def _median_sqdist(phi: np.ndarray) -> np.ndarray:
    """samples × samples: median over regions of the squared kernel distance (pairwise, to bound memory)."""
    n = len(phi)
    d = np.zeros((n, n))
    for i, j in itertools.combinations(range(n), 2):
        d[i, j] = d[j, i] = np.median(((phi[i] - phi[j]) ** 2).sum(-1))
    return d


def _robust_sd(x: np.ndarray, axis: int = -1) -> np.ndarray:
    return 1.4826 * np.median(np.abs(x - np.median(x, axis=axis, keepdims=True)), axis=axis)


def _strata(x: np.ndarray, n: int) -> np.ndarray:
    if n <= 1:
        return np.zeros(len(x), int)
    cuts = np.quantile(x, np.linspace(0, 1, n + 1)[1:-1])
    return np.searchsorted(cuts, x, side="right")


def _robust_z(x: np.ndarray, strata: np.ndarray) -> np.ndarray:
    z = np.full(len(x), np.nan)
    for s in np.unique(strata):
        m = strata == s
        med = np.median(x[m])
        mad = 1.4826 * np.median(np.abs(x[m] - med))
        z[m] = (x[m] - med) / mad if mad > 0 else 0.0
    return z


def group_consistency(profiles: Profiles, config: Config, groups: dict | None = None) -> Consistency:
    """Within-group consistency of every condition (or of ``groups``: name → sample indices)."""
    names = [s["sample_id"] for s in profiles.samples]
    cond = profiles.conditions
    if groups is None:
        groups = {c: np.flatnonzero(cond == c) for c in dict.fromkeys(cond.tolist())}
    groups = {g: np.asarray(i) for g, i in groups.items()}
    y, core_edges, complete = core_profiles(profiles, config)
    phi_all, mass_all, ok, L = shape_features(y[:, complete], M.centers(core_edges))
    centers = M.centers(core_edges)
    ridx = np.flatnonzero(complete)[ok]
    notes = [f"{len(ridx)} of {len(complete)} regions have no missing bins and background-corrected support mass ≥ "
             f"{MIN_RELATIVE_MASS:g} × the sample's median in every sample; all statistics use these regions"]
    flags: list[str] = []
    n_s = len(names)
    nan_matrix = np.full((n_s, n_s), np.nan)
    if len(ridx) < 20:
        notes.append("fewer than 20 usable regions: consistency not assessed")
        return Consistency(names, cond.tolist(), groups, ridx, nan_matrix, nan_matrix.copy(), nan_matrix.copy(),
                           flags=flags, notes=notes)
    shape_d = _median_sqdist(phi_all)
    meta = phi_all.mean(1)
    meta_d = np.sqrt(((meta[:, None] - meta[None]) ** 2).sum(-1))
    # Amplitude: background-corrected core mass on the log scale; the robust spread of a pair's
    # per-region log ratios ignores a global scale difference and a minority of changed regions.
    width = M.widths(core_edges)
    sig = (np.nan_to_num(y[:, ridx]) * width).sum(-1)
    pos = sig[sig > 0]
    log_sig = np.log2(np.maximum(sig, 0) + (0.1 * float(np.median(pos)) if pos.size else 1.0))
    amp_d = _robust_sd(log_sig[:, None] - log_sig[None])
    frags = getattr(profiles, "fragments", None)
    depth = (np.nanmedian(np.asarray(frags)[:, ridx], axis=1) if frags is not None else None)
    depth_label = ("median fragments per region" if frags is not None
                   else "median window signal per region (signal units, before background correction)")
    raw = (np.nan_to_num(profiles.values[:, ridx]) * M.widths(profiles.edges)).sum(-1)
    depth_v = depth if depth is not None else np.median(raw, axis=1)
    blocks = region_blocks(profiles.regions, ridx)
    oriented = bool(profiles.oriented)

    sample_rows, partition_rows, region_rows, group_rows = [], [], [], []
    for g, idx in groups.items():
        n = len(idx)
        gname = str(g)
        if n < 2:
            notes.append(f"{gname}: one sample; within-group consistency cannot be assessed")
            continue
        iu = np.triu_indices(n, 1)
        grow = dict(group=gname, samples=n, median_shape_distance=float(np.median(shape_d[np.ix_(idx, idx)][iu])),
                    median_meta_distance=float(np.median(meta_d[np.ix_(idx, idx)][iu])),
                    median_amplitude_distance=float(np.median(amp_d[np.ix_(idx, idx)][iu])))
        others = np.setdiff1d(np.arange(n_s), idx)
        if others.size:
            grow["median_shape_distance_to_other_groups"] = float(np.median(shape_d[np.ix_(idx, others)]))
            grow["median_meta_distance_to_other_groups"] = float(np.median(meta_d[np.ix_(idx, others)]))
        phi = phi_all[idx]
        # Per-region within-group spread, scored against regions of similar signal mass (sampling noise in
        # the normalised profile falls with mass, so raw spreads would rank low-signal regions first).
        spread = ((phi - phi.mean(0)) ** 2).sum((0, 2)) / (n - 1)
        gmass = np.log(mass_all[idx].mean(0))
        n_str = int(np.clip(len(ridx) // MIN_REGIONS_PER_STRATUM, 1, 10))
        z = _robust_z(np.log(np.maximum(spread, 1e-300)), _strata(gmass, n_str))
        divergent = z > REGION_Z_FLAG
        grow["divergent_regions"] = int(divergent.sum())
        grow["divergent_share"] = float(divergent.mean())
        driver = np.full(len(ridx), -1)
        driver_ratio = np.full(len(ridx), np.nan)
        tot = phi.sum(0)
        # Offset of each sample from the mean of the others (its own profile not in the reference); for two
        # replicates both offsets are the same difference with opposite sign.
        delta = np.stack([phi[k] - (tot - phi[k]) / (n - 1) for k in range(n)])
        if n >= 3:
            loo = (delta ** 2).sum(-1)
            driver = np.argmax(loo, axis=0)
            srt = np.sort(loo, axis=0)
            driver_ratio = srt[-1] / np.maximum(np.median(srt[:-1], axis=0), 1e-300)
        for k in np.flatnonzero(divergent):
            region_rows.append(dict(group=gname, region_index=int(ridx[k]), within_spread=float(spread[k]),
                                    robust_z=float(z[k]), log_mass=float(gmass[k]),
                                    driver=names[idx[driver[k]]] if driver[k] >= 0 else "",
                                    driver_ratio=float(driver_ratio[k])))
        # Per-sample systematic offset (shared across regions): mean over regions of the linear offset
        # vector, tested against its own across-region variability.
        yg = np.nan_to_num(y[idx][:, ridx])
        lin, agg_norm, t_dir = linear_offsets(yg, mass_all[idx], L, centers)
        offset = [offset_test(lin[k], blocks) for k in range(n)]
        sizes = np.array([o["size"] for o in offset])
        tr = [translation_part(lin[k], t_dir, blocks=blocks) if oriented else None for k in range(n)]
        shares = [shift_shares(tr[k]["per_region"]) if oriented else {} for k in range(n)]
        mate_shares = [dict() for _ in range(n)]
        # Typical offset of the group-mates, measured without the candidate (with it, each mate's reference
        # contains the candidate's own deviation). The factor puts a reference of n − 2 samples on the scale
        # of one of n − 1, so exchangeable replicates give a ratio near 1.
        mates = np.full(n, np.nan)
        if n >= 3:
            f = np.sqrt((1 + 1 / (n - 1)) / (1 + 1 / max(n - 2, 1)))
            for k in range(n):
                keep = np.delete(np.arange(n), k)
                d_k, _, t_k = linear_offsets(yg[keep], mass_all[idx][keep], L, centers)
                mates[k] = f * np.median(np.linalg.norm(d_k.mean(1), axis=-1))
                if oriented and float(t_k @ t_k) > 0:
                    # Noise reference for F_δ: the group-mates' own per-region shifts, on the same scale.
                    ms = [shift_shares(f * (d_k[j] @ t_k) / float(t_k @ t_k)) for j in range(n - 1)]
                    mate_shares[k] = {f"{key}_mates": float(np.median([m[key] for m in ms])) for key in ms[0]}
        ls = log_sig[idx]
        amp_loo = _robust_sd(ls - (ls.sum(0) - ls) / (n - 1))
        med = np.median((delta ** 2).sum(-1), axis=1)
        n_div = int(divergent.sum())
        for k, i in enumerate(idx):
            tk = tr[k] or {}
            row = dict(sample=names[i], group=gname, offset_size=offset[k]["size"],
                       offset_relative=offset[k]["size"] / agg_norm if agg_norm > 0 else float("nan"),
                       offset_vs_mates=float(sizes[k] / mates[k]) if n >= 3 and mates[k] > 0 else float("nan"),
                       offset_p=offset[k]["p"],
                       shape_deviation=float(med[k]), shape_deviation_ratio=float(med[k] / np.median(med)),
                       amplitude_deviation=float(amp_loo[k]), depth=float(depth_v[i]),
                       shift_bp=tk.get("bp", float("nan")), shift_bp_low=tk.get("low", float("nan")),
                       shift_bp_high=tk.get("high", float("nan")), shift_share_of_offset=tk.get("share", float("nan")),
                       **shares[k], **mate_shares[k])
            if n >= 3:
                row.update(most_deviant_share=float(np.mean(driver == k)), expected_share=1.0 / n,
                           divergent_regions_driven=int(np.sum(divergent & (driver == k))),
                           divergent_driven_share=float(np.sum(divergent & (driver == k)) / n_div) if n_div else float("nan"))
            sample_rows.append(row)
        if n == 2:
            # No group-mates to compare with: report the pair's offset, do not flag on significance alone.
            notes.append(f"{gname}: two replicates; their systematic offset is {offset[0]['size'] / agg_norm:.1%} of the "
                         f"aggregate profile (p = {offset[0]['p']:.2g}), but which one deviates cannot be told, and "
                         "with no third replicate there is no reference for how large a replicate difference is usual")
        else:
            ratio = sizes / np.maximum(mates, 1e-300)
            for k in np.flatnonzero((np.array([o["p"] for o in offset]) < FLAG_ALPHA / n) & (ratio >= OFFSET_FLAG_RATIO)):
                shallow = depth_v[idx[k]] < 0.7 * np.median(depth_v[idx])
                flags.append(f"{gname}: {names[idx[k]]} stands out from its group-mates: its sample-wide shape offset is "
                             f"{ratio[k]:.1f}× theirs ({sizes[k] / agg_norm:.1%} of the aggregate profile; p = "
                             f"{offset[k]['p']:.2g}, Bonferroni over {n} samples). Check fragment length, alignment and IP quality"
                             + ("; it is also much shallower than its group-mates" if shallow else ""))
        # Replicate structure: splits into two halves of ≥ 2 samples.
        if n >= 4:
            md = meta_d[np.ix_(idx, idx)]
            pair = {(i, j): ((phi[i] - phi[j]) ** 2).sum(-1) for i, j in itertools.combinations(range(n), 2)}
            seen, rows_g = set(), []
            psi_g, u_g, _, t_g = _linear_features(yg, mass_all[idx], L, centers)
            for size in range(2, n // 2 + 1):
                for a in itertools.combinations(range(n), size):
                    b = tuple(sorted(set(range(n)) - set(a)))
                    key = tuple(sorted((a, b)))
                    if key in seen:
                        continue
                    seen.add(key)
                    within = np.mean([np.mean([md[i, j] for i, j in itertools.combinations(h, 2)]) for h in (a, b)])
                    across = np.mean([md[i, j] for i in a for j in b])
                    # Per-region contrast of pair distances: across-half mean minus the average of the two
                    # within-half means. Each sample's own noise enters both terms with the same weight, so
                    # unequal depth cancels; under exchangeable replicates its expectation is zero.
                    c = (np.mean([pair[min(i, j), max(i, j)] for i in a for j in b], axis=0)
                         - 0.5 * sum(np.mean([pair[i, j] for i, j in itertools.combinations(h, 2)], axis=0) for h in (a, b)))
                    var_c = _block_cov(c, blocks)[0, 0]
                    zc = float(c.mean() / np.sqrt(var_c / len(c))) if var_c > 0 else 0.0
                    d_ab = _remove(psi_g[list(b)].mean(0) - psi_g[list(a)].mean(0), u_g)
                    ts = translation_part(d_ab, t_g, blocks=blocks) if oriented else {}
                    # π̂₁ only for halves of equal size: there each sample's noise enters every label assignment
                    # alike; with unequal halves a shallow replicate inflated it (evaluation v1: 0.044 with no
                    # difference, 6 replicates, one at 30 % depth).
                    ext = (changed_share(region_test_p(yg, list(a), list(b), centers, oriented), seed=len(rows_g),
                                         blocks=blocks) if len(a) == len(b) else
                           dict(share=float("nan"), low=float("nan"), high=float("nan"), called=0))
                    rows_g.append(dict(group=gname, half_a="+".join(names[idx[i]] for i in a),
                                       half_b="+".join(names[idx[i]] for i in b),
                                       meta_ratio=float(across / within) if within > 0 else float("nan"),
                                       contrast_z=zc, contrast_p=float(stats.norm.sf(zc)),
                                       shift_bp=ts.get("bp", float("nan")), shift_bp_low=ts.get("low", float("nan")),
                                       shift_bp_high=ts.get("high", float("nan")),
                                       regions_differing_share=ext["share"], regions_differing_low=ext["low"],
                                       regions_differing_high=ext["high"], regions_called=ext["called"]))
            partition_rows.extend(rows_g)
            best = min(rows_g, key=lambda r: r["contrast_p"])
            grow.update(partitions=len(rows_g), best_partition=f"{best['half_a']} | {best['half_b']}",
                        best_partition_meta_ratio=best["meta_ratio"], best_partition_p=best["contrast_p"])
            if best["contrast_p"] < FLAG_ALPHA / len(rows_g):
                flags.append(f"{gname}: replicates form two sub-groups ({best['half_a']} | {best['half_b']}; per-region "
                             f"distance contrast p = {best['contrast_p']:.2g}, Bonferroni over {len(rows_g)} splits; "
                             f"meta-profile across/within ratio {best['meta_ratio']:.2f}). Contrasts or null splits along "
                             "this partition are not between exchangeable replicates")
        group_rows.append(grow)
    return Consistency(names, cond.tolist(), groups, ridx, shape_d, meta_d, amp_d, sample_rows, partition_rows,
                       region_rows, group_rows, flags, notes + [f"depth column: {depth_label}"])


def offset_test(delta: np.ndarray, blocks: np.ndarray | None = None) -> dict:
    """Is the mean over regions of per-region offset vectors (regions × features) zero?

    R‖mean‖² has expectation tr(Σ) under a zero mean, with Σ = R·Var(mean) estimated from blocks of
    neighbouring regions (or the ordinary covariance without ``blocks``); its distribution is
    approximated by a scaled χ² with matching first two moments (Satterthwaite). Unequal depth enlarges
    Σ but not the mean, so it does not by itself give small p-values.
    """
    R = len(delta)
    m = delta.mean(0)
    S = _block_cov(delta, blocks)
    tr, tr2 = float(np.trace(S)), float((S ** 2).sum())
    T = R * float(m @ m)
    if tr <= 0 or tr2 <= 0:
        return dict(size=float(np.sqrt(m @ m)), p=float("nan"))
    g, h = tr2 / tr, tr ** 2 / tr2
    return dict(size=float(np.sqrt(m @ m)), p=float(stats.chi2.sf(T / g, h)))
