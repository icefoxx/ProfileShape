"""Step 3 of docs/PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md: power of each candidate's primary A1 criterion.

    PYTHONPATH=benchmarks python benchmarks/power_sim.py PILOT_ESTIMATE.json OUT_DIR [--runs 100]

Templates: condition-blind pooled pilot profiles. Depth: pilot estimate of the full run. Noise: from the PRO-seq
group of consistency_v2 (written to the result before any candidate runs). Effects and criteria: the step 3
specification of the protocol. Representations: midpoint counts, and coverage = counts convolved with the pilot
fragment lengths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from profileshape import metrics as M
from profileshape.analysis import analyse
from profileshape.data import ASSAY_PRESETS, Config, Profiles
from profileshape.regional import per_region_tests
from profileshape.simulate import SimSettings, _ar_noise, base_shapes

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md"
SEED = 9_800_000
N_REGIONS, N_REP = 3000, 3
CFG = dict(ASSAY_PRESETS["proseq"], orientation="oriented")
# candidate → (effect, size, share of regions changed, criterion); from the protocol's step 3 specification
EFFECTS = {"GSE112379 heat shock RPB1 ChIP": ("body_loss", 0.5, 0.30, "upstream"),
           "GSE205504 NELFB-AID Pol II ChIP": ("shift", 60.0, 0.50, "downstream"),
           "GSE199086 KL-2 FP-TT-seq": ("body_loss", 0.3, 0.50, "upstream")}
MODES = ("midpoint", "coverage")


def noise_params() -> dict:
    """Replicate amplitude sd and bin dispersion from the PRO-seq group of consistency_v2 (see protocol)."""
    d = np.load(ROOT / "inputs/caches/cons_proseq.npz", allow_pickle=True)
    v = np.clip(np.nan_to_num(d["values"]), 0, None)
    ids = [str(x) for x in d["sample_ids"]]
    groups = json.loads(str(d["groups"]))  # group → sample ids (2 replicates each)
    sds, cv2 = [], []
    for members in groups.values():
        x = v[[ids.index(m) for m in members]]
        tot = x.sum(-1)
        x = x / np.maximum(tot.sum(1, keepdims=True), 1e-12)[..., None] * tot.sum(1).mean()  # library-normalised
        tot = x.sum(-1)
        top = tot.mean(0) >= np.median(tot.mean(0))
        sds.append(np.median(np.std(np.log(np.maximum(tot[:, top], 1e-12)), axis=0, ddof=1)))
        m, s = x.mean(0), x.std(0, ddof=1)
        hi = m >= np.quantile(m, 0.9)
        cv2.append(np.median((s[hi] / m[hi]) ** 2))
    return dict(replicate_sd=float(np.median(sds)), dispersion=float(np.median(cv2)),
                source="consistency_v2 PRO-seq cache (cons_proseq.npz)")


def _templates(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = np.load(path)
    pooled = t["pooled"].astype(float)
    k = np.array([1.0, 1.0, 1.0]) / 3
    pooled = np.apply_along_axis(lambda r: np.convolve(r, k, mode="same"), 1, pooled)
    return pooled[pooled.sum(1) > 0], t["edges"], t["fragment_lengths"]


def _effect(peak, c, kind, size):
    out = peak.copy()
    if kind == "body_loss":
        out[:, c > 300] *= 1 - size
    elif kind == "shift":
        for r in range(len(out)):
            out[r] = np.interp(c - size, c, out[r], left=0, right=0)
    return out


def generate(job):
    """One simulated 3 vs 3 experiment with the simulator's noise model and the candidate's effect."""
    name, rep, mode, tpl_path, target, noise, seed = job
    templates, edges, lengths = _templates(tpl_path)
    s = SimSettings(n_ref=N_REP, n_alt=N_REP, n_regions=N_REGIONS, flank=1000, binsize=10,
                    replicate_sd=noise["replicate_sd"], dispersion=noise["dispersion"])
    rng = np.random.default_rng(seed)
    _, c, peak, bg, strength = base_shapes(s, rng, templates)
    kind, size, share, direction = EFFECTS[name.split(" | ")[0]]
    rows = rng.choice(N_REGIONS, int(share * N_REGIONS), replace=False)
    peak_alt = peak.copy()
    peak_alt[rows] = _effect(peak[rows], c, kind, size)
    core = M.core_mask(edges, CFG["core_flank"])
    depth = target / np.median(strength * (peak + bg[:, None])[:, core].sum(1))
    values = np.empty((2 * N_REP, N_REGIONS, len(c)))
    for i in range(2 * N_REP):
        pk = peak_alt if i >= N_REP else peak
        lib = rng.lognormal(-(s.library_sd ** 2) / 2, s.library_sd)
        rp = rng.lognormal(-(s.replicate_sd ** 2) / 2, s.replicate_sd, N_REGIONS)[:, None]
        mean = depth * lib * strength[:, None] * (rp * pk + bg[:, None])
        mean = mean * _ar_noise(rng, mean.shape, s.ar_rho, s.ar_sd)
        r = 1.0 / s.dispersion
        values[i] = rng.negative_binomial(r, r / (r + mean))
    frags = values[..., core].sum(-1)
    if mode == "coverage":  # midpoints spread over fragment spans (mean of three boxes: q25, median, q75)
        cov = np.zeros_like(values)
        for L in lengths:
            w = max(1, int(round(L / 10)))
            cov += np.apply_along_axis(lambda x: np.convolve(x, np.ones(w), mode="same"), -1, values) / len(lengths)
        values = cov
    smp = [dict(sample_id=f"{'B' if i >= N_REP else 'A'}{i % N_REP + 1}", condition="B" if i >= N_REP else "A") for i in range(2 * N_REP)]
    regions = [dict(chrom="chrSim", start=10_000 * (k + 1), end=10_000 * (k + 1) + 1, name=f"r{k}", strand="+") for k in range(N_REGIONS)]
    prof = Profiles(values, edges, smp, regions, True, provenance={"source": "bam"}, fragments=frags)
    res = per_region_tests(analyse(prof, Config(reference="A", alternative="B", **CFG), bootstrap=0))
    any_, sh = res.tests["shape_any"], res.tests["shift"]
    tested = any_["call"] != "not_tested"
    share_any = float(np.isin(any_["call"], ("Changed",)).sum() / max(tested.sum(), 1))
    sig = np.asarray(sh["fdr"]) < 0.05
    dirs = np.asarray(sh["direction"])[sig]
    n_sig = int(sig.sum())
    toward = float(np.mean(dirs < 0) if direction == "upstream" else np.mean(dirs > 0)) if n_sig else float("nan")
    met = (n_sig >= 50 and toward >= 0.80) and (share_any >= 0.10 if direction == "upstream" else True)
    return name, mode, rep, dict(share_any=share_any, shift_sig=n_sig, toward=toward, met=bool(met),
                                 median_core=float(np.median(frags)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("estimate")
    ap.add_argument("out")
    ap.add_argument("--runs", type=int, default=100)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *x: subprocess.run(["git", *x], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]  # noqa: E731
    noise = noise_params()
    sens = dict(noise, dispersion=0.15, source="sensitivity (amendment R2): simulator default dispersion, reported only")
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted=len(git("status", "--porcelain", "--", "src", "benchmarks").splitlines()),
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"), protocol_sha256=sha(PROTOCOL), script_sha256=sha(__file__),
                seeds=f"{SEED}+", runs=a.runs, started=time.strftime("%Y-%m-%dT%H:%M:%S"), noise=noise)
    (out / "noise_params.json").write_text(json.dumps(prov, indent=1))  # written before any candidate runs
    print("noise", noise, flush=True)
    est = json.loads(Path(a.estimate).read_text())
    jobs = []
    for ci, (name, e) in enumerate(est.items()):
        for mode in MODES:
            for vi, (tag, nz) in enumerate((("", noise), (" | dispersion 0.15", sens))):
                for rep in range(a.runs):
                    jobs.append((name + tag, rep, mode, e["templates"], e["estimated_full_median_over_pilots"], nz,
                                 SEED + 10_000 * ci + 1000 * vi + rep))
    rows = []
    with ProcessPoolExecutor(a.workers) as ex:
        rows = list(ex.map(generate, jobs, chunksize=2))
    prov["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    summ = {}
    for name, mode, rep, r in rows:
        summ.setdefault(name, {}).setdefault(mode, []).append(r)
    table = {n: {m: dict(power=float(np.mean([x["met"] for x in v])), runs=len(v),
                         share_any=float(np.mean([x["share_any"] for x in v])), shift_sig=float(np.mean([x["shift_sig"] for x in v])),
                         toward=float(np.nanmean([x["toward"] for x in v])) if any(np.isfinite(x["toward"]) for x in v) else float("nan"),
                         median_core=float(np.mean([x["median_core"] for x in v])))
                 for m, v in d.items()} for n, d in summ.items()}
    (out / "power_sim.json").write_text(json.dumps(dict(provenance=prov, estimate=est, table=table), indent=1, default=float))
    L = ["# Step 3 · simulated power of the primary A1 criterion", "",
         f"Protocol sha256 `{prov['protocol_sha256']}`; commit `{prov['commit'][:10]}` (uncommitted {prov['uncommitted']}); script "
         f"`{prov['script_sha256']}`; noise {noise}; {a.runs} runs per cell; {prov['started']} → {prov['finished']}.", "",
         "| candidate | representation | depth (core median) | power | mean shape_any share | mean shift-significant | mean share in predicted direction |",
         "|---|---|---|---|---|---|---|"]
    for n, d in table.items():
        for m, r in d.items():
            L.append(f"| {n} | {m} | {r['median_core']:.0f} | **{r['power']:.2f}** | {r['share_any']:.3f} | {r['shift_sig']:.0f} | {r['toward']:.2f} |")
    (out / "POWER_SIM.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
