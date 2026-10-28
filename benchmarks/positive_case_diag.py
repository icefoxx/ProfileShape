"""Diagnostic supplement to the GSE218269 positive case (amendment B5; exploratory throughout).

    PYTHONPATH=benchmarks python benchmarks/positive_case_diag.py BAM_DIR BIGWIG_CACHE.npz BAM_CACHE.npz OUT_DIR

D1 A2 ablation · D2 per-region width, replicate support, detectability · D3 processing and replicate sensitivity.
Every variant is reported; nothing here updates the confirmatory verdict.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

import insertion as I
import positive_case as PC
from profileshape import metrics as M
from profileshape.analysis import analyse
from profileshape.data import Config, Profiles, extract
from profileshape.regional import TEMPLATE_SMOOTH, per_region_tests, translation_estimate

ROOT = Path(__file__).resolve().parents[1]
SEED = 9_700_000
RANGES = ((20, 70), (20, 120), (70, 120))
MODES = ("midpoint", "coverage")
FACTORS = (0.95, 0.90, 0.80)
CALLED = ("Changed", "Higher", "Lower")


def _load(cache):
    d = np.load(cache, allow_pickle=True)
    prof = Profiles(d["values"], d["edges"], json.loads(str(d["samples"])), json.loads(str(d["regions"])), True,
                    provenance=json.loads(str(d["provenance"])), fragments=d["fragments"] if d["fragments"].size else None)
    readout = {k: np.array(v, float) for k, v in json.loads(str(d["readout"])).items()}
    return d, prof, readout


def _cfg():
    return Config(reference="DMSO", alternative="CDK7i", **PC.CFG)


def _subset(prof, keep):
    return Profiles(prof.values[keep], prof.edges, [prof.samples[i] for i in keep], prof.regions, True,
                    provenance=prof.provenance, fragments=None if prof.fragments is None else prof.fragments[keep])


def _q(x, qs=(0.05, 0.25, 0.5, 0.75, 0.95)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [float(v) for v in np.quantile(x, qs)] if x.size else []


# ---------------------------------------------------------------- D1
def ablation(base, shapes: dict, y, chrom, boot=500, seed=0) -> dict:
    """CV R² of base and base + each shape set, the same rows, folds and bootstrap draws for every model."""
    allx = np.c_[tuple(shapes.values())]
    ok = np.isfinite(y) & np.all(np.isfinite(base), 1) & np.all(np.isfinite(allx), 1)
    base, y, chrom = base[ok], y[ok], chrom[ok]
    X = {k: np.c_[base, v[ok]] for k, v in shapes.items()}
    chroms = np.unique(chrom)
    fold_of = {c: i % 10 for i, c in enumerate(chroms)}
    folds = np.array([fold_of[c] for c in chrom])
    r2 = dict(base=PC._cv_r2(base, y, folds), **{k: PC._cv_r2(x, y, folds) for k, x in X.items()})
    rng = np.random.default_rng(seed)
    draws = {k: [] for k in X}
    for _ in range(boot):
        pick = rng.choice(chroms, len(chroms), replace=True)
        idx, f = PC.bootstrap_folds(pick, chrom, fold_of)
        if len(np.unique(f)) < 2:
            continue
        b0 = PC._cv_r2(base[idx], y[idx], f)
        for k, x in X.items():
            draws[k].append(PC._cv_r2(x[idx], y[idx], f) - b0)
    out = dict(regions=int(ok.sum()), r2=r2, gain={})
    for k in X:
        g = np.array(draws[k])
        out["gain"][k] = dict(estimate=r2[k] - r2["base"], low=float(np.quantile(g, 0.025)), high=float(np.quantile(g, 0.975)))
    full = np.array(draws["both"])
    for k in ("shift", "kernel"):
        dd = full - np.array(draws[k])
        out["gain"][f"both_minus_{k}"] = dict(estimate=r2["both"] - r2[k], low=float(np.quantile(dd, 0.025)),
                                              high=float(np.quantile(dd, 0.975)))
    return out


def d1(a, res, prof, readout, length):
    D = PC.a2_design(a, res, prof, readout, length)
    shapes = dict(shift=D["signed_shift"][:, None], kernel=D["kernel"][:, None], both=np.c_[D["signed_shift"], D["kernel"]])
    return {f"{m}min": ablation(*D["readouts"][m][:1], shapes, D["readouts"][m][1], D["chrom"], seed=m) for m in (30, 15, 60)}


# ---------------------------------------------------------------- D2
def sample_widths(a, idx):
    strong = np.all(a.core_mass[:, idx] > 0, axis=0)
    p = M.normalise(a.core_signal[:, idx] * M.widths(a.core_edges))
    w = M.mean_abs_deviation(p, a.core_edges)
    w[:, ~strong] = np.nan
    return w


def d2_width(a, res):
    any_ = res.tests["shape_any"]
    idx = any_["region_index"]
    called = np.isin(any_["call"], CALLED)
    w = sample_widths(a, idx)
    ref, alt = a.ref_idx, a.alt_idx
    diff = np.nanmean(w[alt], 0) - np.nanmean(w[ref], 0)
    below = np.all(np.isfinite(w), 0) & (w[alt].max(0) < w[ref].min(0))
    above = np.all(np.isfinite(w), 0) & (w[alt].min(0) > w[ref].max(0))
    out = {}
    for name, m in (("called", called), ("not_called", ~called)):
        ok = m & np.isfinite(diff)
        out[name] = dict(regions=int(ok.sum()), diff_quantiles_bp=_q(diff[ok]), share_narrower=float(np.mean(diff[ok] < 0)),
                         all_alt_below_all_ref=float(np.mean(below[ok])), all_alt_above_all_ref=float(np.mean(above[ok])),
                         exchangeable_expectation=1 / 20)
    mw = res.tests["mad_width"]
    out["mad_width_test"] = dict(tested=int(np.sum(mw["call"] != "not_tested")), fdr_calls=int(np.sum(np.asarray(mw["fdr"]) < 0.05)),
                                 p_below_001=float(np.nanmean(np.asarray(mw["p"], float) < 0.01)))
    return out


def d2_shift_subset(a, res):
    sh = res.tests["shift"]
    est = np.asarray(sh["estimate"], float)
    flags = {k: int(v) for k, v in zip(*np.unique(sh["translation_flag"], return_counts=True))}
    sub = np.isfinite(est)
    out = dict(condition="translation_consistent: cross-correlation and projection agree in sign, translation R² ≥ 0.5, "
                         "ambiguity < 0.9, not at the search edge (regional.py:848–858)",
               regions=int(sub.sum()), of_tested=int(len(est)), estimate_quantiles_bp=_q(est[sub]),
               flags=flags, shift_fdr_calls_in_subset=int(np.sum(sub & (np.asarray(sh["fdr"]) < 0.05))))
    # all regions, every 3 vs 3 label assignment (the observed one included); rank of the observed share
    sidx = sh["region_index"]
    y = np.nan_to_num(a.core_signal[:, sidx])
    ids = list(a.ref_idx) + list(a.alt_idx)
    rows = []
    for A in itertools.combinations(ids, 3):
        if ids[0] not in A:
            continue
        B = [i for i in ids if i not in A]
        s = np.asarray(translation_estimate(y, np.array(A), np.array(B), M.centers(a.core_edges), smooth=TEMPLATE_SMOOTH)["shift"], float)
        rows.append(dict(ref=[a.profiles.samples[i]["sample_id"] for i in A], median_bp=float(np.nanmedian(s)),
                         share_negative=float(np.nanmean(s < 0)), observed=sorted(A) == sorted(a.ref_idx)))
    obs = next(r for r in rows if r["observed"])
    rank = sum(r["share_negative"] >= obs["share_negative"] for r in rows)
    out["all_regions_xcorr"] = dict(assignments=rows, observed_share_negative=obs["share_negative"],
                                    one_sided_p=rank / len(rows), p_floor=1 / len(rows),
                                    rule="share of regions with a negative (upstream) cross-correlation shift, CDK7i vs DMSO; "
                                         "one-sided p = rank of the observed share among the 10 assignments (ties counted against)")
    return out


def _aggregate_iqr(a):
    p = M.normalise(np.clip(np.nanmean(np.nan_to_num(a.core_signal[a.ref_idx]), axis=(0, 1)), 0, None) * M.widths(a.core_edges))
    cdf = np.cumsum(p)
    c = M.centers(a.core_edges)
    return float(np.interp(0.75, cdf, c) - np.interp(0.25, cdf, c))


def d2_insert_job(job):
    cache, f, draw, chosen = job
    d, prof, readout = _load(cache)
    cfg = _cfg()
    alt = [i for i, s in enumerate(prof.samples) if s["condition"] == "CDK7i"]
    rng = np.random.default_rng(SEED + 100 * FACTORS.index(f) + draw)
    v = I.insert(prof.values, prof.edges, alt, np.array(chosen), f"broaden{f}", rng, cfg.background_width)
    p2 = Profiles(v, prof.edges, prof.samples, prof.regions, True, provenance=prof.provenance, fragments=prof.fragments)
    res = per_region_tests(analyse(p2, cfg, bootstrap=0))
    return f, draw, {m: _called_share(res, m, chosen) for m in ("mad_width", "shape_kernel")}


def _called_share(res, m, regions):
    t = res.tests[m]
    pos = {r: i for i, r in enumerate(t["region_index"])}
    hit = [np.asarray(t["fdr"])[pos[r]] < 0.05 for r in regions if r in pos]
    return dict(share=float(np.mean(hit)) if hit else float("nan"), regions=len(hit))


# ---------------------------------------------------------------- D3
def _summary(prof, readout, length, fragments_core=None):
    A1, A2, res = PC.case_results(prof, readout, length, _cfg(), bootstrap=0, minutes=(30,))
    out = dict(A1=A1, A2_30min=A2["30min"])
    if fragments_core is not None:
        med = np.nanmedian(fragments_core, axis=1)
        out["depth_median_over_samples"] = float(np.median(med))
        out["depth_per_sample"] = med.round(1).tolist()
    return out


def d3_variant_job(job):
    bam_dir, bigwig_cache, lo, hi, mode = job
    d = np.load(bigwig_cache, allow_pickle=True)
    regions = json.loads(str(d["regions"]))
    smp = [dict(sample_id=k, condition=k.split("_")[0], signal_file=str(Path(bam_dir) / f"{k}.bam")) for k in PC.POL2]
    opt = dict(PC.BAM_OPTIONS, bam_mode=mode, bam_min_fragment=lo, bam_max_fragment=hi)
    prof = extract(Config(**PC.CFG, **opt), smp, regions)
    core = M.core_mask(prof.edges, PC.CFG["core_flank"])
    fc = np.asarray(prof.counts)[..., core].sum(-1) if mode == "midpoint" else None
    readout = {k: np.array(v, float) for k, v in json.loads(str(d["readout"])).items()}
    return f"{lo}-{hi} {mode}", _summary(prof, readout, d["length"], fc)


def d3_loo_job(job):
    cache, drop = job
    d, prof, readout = _load(cache)
    names = [s["sample_id"] for s in prof.samples]
    keep = [i for i, n in enumerate(names) if n != drop]
    return drop, _summary(_subset(prof, keep), readout, d["length"])


def main(argv):
    bam_dir, bigwig_cache, cache, out = argv[:4]
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *x: subprocess.run(["git", *x], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]  # noqa: E731
    R = dict(provenance=dict(commit=git("rev-parse", "HEAD"), uncommitted=len(git("status", "--porcelain", "--", "src", "benchmarks").splitlines()),
                             regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"), protocol_sha256=sha(PC.PROTOCOL),
                             script_sha256=sha(__file__), positive_case_sha256=sha(PC.__file__), seeds=f"{SEED}+",
                             started=time.strftime("%Y-%m-%dT%H:%M:%S"), exploratory=True))
    d, prof, readout = _load(cache)
    a = analyse(prof, _cfg(), bootstrap=0)
    res = per_region_tests(a)
    # reproduction check against the B4 record (bootstrap only affects curves)
    A1 = PC.a1_results(a, res)
    b4 = json.loads((ROOT / "benchmarks/results/positive_case_gse218269_bam_exploratory/positive_case.json").read_text())
    R["check_B4_reproduced"] = all(np.isclose(A1[k], b4["A1"][k], equal_nan=True) for k in ("shape_any_share", "shift_significant", "centroid_bp", "centroid_low"))
    print("B4 reproduced:", R["check_B4_reproduced"], flush=True)
    R["D1"] = d1(a, res, prof, readout, d["length"])
    print("D1 done", flush=True)
    R["D2"] = dict(width=d2_width(a, res), shift_subset=d2_shift_subset(a, res))
    iqr_ref = _aggregate_iqr(a)
    iqr_diff = next(e for e in a.effects if e["metric"] == "iqr_width")["estimate"]
    R["D2"]["aggregate_iqr"] = dict(dmso_bp=iqr_ref, diff_bp=iqr_diff, equivalent_f=1 + iqr_diff / iqr_ref)
    tested = res.tests["shape_kernel"]["region_index"]
    jobs = []
    for f in FACTORS:
        for draw in range(3):
            rng = np.random.default_rng(SEED + 1000 + 100 * FACTORS.index(f) + draw)
            jobs.append((cache, f, draw, sorted(rng.choice(tested, 500, replace=False).tolist())))
    R["D2"]["insertion"] = dict(without={}, with_insertion={})
    for f, draw, chosen in [(j[1], j[2], j[3]) for j in jobs]:
        R["D2"]["insertion"]["without"][f"{f}|{draw}"] = {m: _called_share(res, m, chosen) for m in ("mad_width", "shape_kernel")}
    with ProcessPoolExecutor(6) as ex:
        for f, draw, r in ex.map(d2_insert_job, jobs):
            R["D2"]["insertion"]["with_insertion"][f"{f}|{draw}"] = r
        print("D2 done", flush=True)
        R["D3"] = dict(processing={}, leave_one_out={})
        for k, r in ex.map(d3_variant_job, [(bam_dir, bigwig_cache, lo, hi, m) for lo, hi in RANGES for m in MODES]):
            R["D3"]["processing"][k] = r
        for k, r in ex.map(d3_loo_job, [(cache, n) for n in PC.POL2]):
            R["D3"]["leave_one_out"][k] = r
    R["provenance"]["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "diag.json").write_text(json.dumps(R, indent=1, default=lambda v: v.item() if isinstance(v, np.generic) else str(v)))
    (out / "DIAG.md").write_text(report(R))
    print(report(R))


def _g(x):
    return f"{x['estimate']:+.4f} [{x['low']:+.4f}, {x['high']:+.4f}]"


def report(R):
    p = R["provenance"]
    L = ["# GSE218269 · diagnostic supplement (amendment B5) — EXPLORATORY", "",
         "S3 was not met; the confirmatory case is \"not evaluable\". Nothing here updates that verdict. "
         f"Commit `{p['commit'][:10]}` (uncommitted files {p['uncommitted']}); regional.py `{p['regional_py_sha256']}`; protocol "
         f"`{p['protocol_sha256']}`; script `{p['script_sha256']}`; {p['started']} → {p['finished']}. "
         f"B4 A1 reproduced: {R['check_B4_reproduced']}.", "",
         "## D1 · A2 ablation (same rows, folds, bootstrap draws)", "",
         "| readout | regions | R² base | + shift | + kernel | + both | gain shift | gain kernel | gain both | both − shift | both − kernel |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for k, r in R["D1"].items():
        g = r["gain"]
        L.append(f"| {k} | {r['regions']} | {r['r2']['base']:.4f} | {r['r2']['shift']:.4f} | {r['r2']['kernel']:.4f} | {r['r2']['both']:.4f} | "
                 f"{_g(g['shift'])} | {_g(g['kernel'])} | {_g(g['both'])} | {_g(g['both_minus_shift'])} | {_g(g['both_minus_kernel'])} |")
    w = R["D2"]["width"]
    L += ["", "## D2 · Per-region width (MAD of the core profile, CDK7i − DMSO)", "",
          "| regions | n | diff quantiles 5/25/50/75/95 % (bp) | share narrower | all 3 CDK7i < all 3 DMSO | all 3 > | expected under exchangeability |",
          "|---|---|---|---|---|---|---|"]
    for k in ("called", "not_called"):
        r = w[k]
        L.append(f"| `shape_any` {k.replace('_', ' ')} | {r['regions']} | {' / '.join(f'{v:+.1f}' for v in r['diff_quantiles_bp'])} | "
                 f"{r['share_narrower']:.3f} | {r['all_alt_below_all_ref']:.3f} | {r['all_alt_above_all_ref']:.3f} | 0.050 each |")
    mw = w["mad_width_test"]
    ag = R["D2"]["aggregate_iqr"]
    L += ["", f"`mad_width` test: {mw['tested']} tested, {mw['fdr_calls']} FDR calls, share p < 0.01 {mw['p_below_001']:.3f}. "
          f"Aggregate IQR: DMSO {ag['dmso_bp']:.0f} bp, difference {ag['diff_bp']:+.1f} bp → equivalent f ≈ {ag['equivalent_f']:.2f}.", "",
          "Detectability (500 random shape-tested regions narrowed by f in the three CDK7i samples; share called at FDR 0.05):", "",
          "| f | draw | `mad_width` without → with | `shape_kernel` without → with |", "|---|---|---|---|"]
    for key, r in R["D2"]["insertion"]["with_insertion"].items():
        r0 = R["D2"]["insertion"]["without"][key]
        f, dr = key.split("|")
        L.append(f"| {f} | {dr} | {r0['mad_width']['share']:.3f} → {r['mad_width']['share']:.3f} | "
                 f"{r0['shape_kernel']['share']:.3f} → {r['shape_kernel']['share']:.3f} |")
    s = R["D2"]["shift_subset"]
    x = s["all_regions_xcorr"]
    L += ["", "Reported `shift` estimate: " + s["condition"] + f"; {s['regions']} of {s['of_tested']} regions; quantiles 5/25/50/75/95 % "
          + " / ".join(f"{v:+.1f}" for v in s["estimate_quantiles_bp"]) + f" bp; flags {s['flags']}.", "",
          f"All regions, cross-correlation shift: observed share negative {x['observed_share_negative']:.3f}; one-sided p = "
          f"{x['one_sided_p']:.1f} (floor {x['p_floor']:.1f}, 10 assignments). Rule: {x['rule']}.", "",
          "## D3 · Processing and replicate sensitivity", "",
          "| variant | depth (median) | `shape_any` share | `shift` sig. | upstream share | centroid bp [95 %] | A2 30 min gain [95 %] |",
          "|---|---|---|---|---|---|---|"]
    for grp, label in (("processing", "{}"), ("leave_one_out", "drop {}")):
        for k, r in R["D3"][grp].items():
            a1, a2 = r["A1"], r["A2_30min"]
            dep = f"{r['depth_median_over_samples']:.0f}" if "depth_median_over_samples" in r else "–"
            L.append(f"| {label.format(k)} | {dep} | {a1['shape_any_share']:.3f} | {a1['shift_significant']} | {a1['upstream_share']:.2f} | "
                     f"{a1['centroid_bp']:+.1f} [{a1['centroid_low']:+.1f}, {a1['centroid_high']:+.1f}] | "
                     f"{a2['gain']:+.4f} [{a2['gain_low']:+.4f}, {a2['gain_high']:+.4f}] |")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main(sys.argv[1:])
