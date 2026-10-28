"""Scenario C of docs/PROFILESHAPE_APPLICATIONS_PLAN.md: changes that window counts and aggregate profiles miss
(GSE218269, CDK7 inhibition, Pol II MNase-ChIP 3 vs 3; B2 BAM cache, midpoints). Exploratory: S3 failed.

    PYTHONPATH=benchmarks python benchmarks/app_scenario_c.py BAM_CACHE.npz OUT_DIR
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from scipy import stats

import independent_validation as IV
from app_scenario_a import diffsplice
from app_scenario_b import EDGER
from positive_case import CFG
from positive_case_diag import _load, _subset
from profileshape import metrics as M
from profileshape.analysis import analyse
from profileshape.data import Config
from profileshape.regional import bh, per_region_tests

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "docs" / "PROFILESHAPE_APPLICATIONS_PLAN.md"


def shape_measures(prof, edges, core_flank):
    """Per sample and region: width (MAD) and centroid of the mass-normalised core profile."""
    core = M.core_mask(edges, core_flank)
    e = edges[np.r_[np.flatnonzero(core), np.flatnonzero(core)[-1] + 1]]
    v = np.clip(np.nan_to_num(prof.values[..., core]), 0, None)
    p = M.normalise(v * M.widths(e))
    return M.mean_abs_deviation(p, e), M.centroid(p, e)


def calls(prof):
    res = per_region_tests(analyse(prof, Config(reference="DMSO", alternative="CDK7i", **CFG), bootstrap=0))
    t = res.tests["shape_any"]
    called = np.zeros(prof.values.shape[1], bool)
    called[np.asarray(t["region_index"])[np.isin(t["call"], ("Changed",))]] = True
    tested = np.zeros(prof.values.shape[1], bool)
    tested[np.asarray(t["region_index"])[t["call"] != "not_tested"]] = True
    return called, tested


def agreement(x, ids, called):
    """Held-out sign agreement by direction subgroup (direction from replicates 1–2, check on replicate 3)."""
    d12 = np.nanmean(x[[ids.index("CDK7i_1"), ids.index("CDK7i_2")]], 0) - np.nanmean(x[[ids.index("DMSO_1"), ids.index("DMSO_2")]], 0)
    d3 = x[ids.index("CDK7i_3")] - x[ids.index("DMSO_3")]
    ok = called & np.isfinite(d12) & np.isfinite(d3) & (d12 != 0) & (d3 != 0)
    out = {}
    for name, m in (("negative", ok & (d12 < 0)), ("positive", ok & (d12 > 0))):
        n, k = int(m.sum()), int((np.sign(d3[m]) == np.sign(d12[m])).sum())
        ci = stats.binomtest(k, n).proportion_ci(0.95, method="exact") if n else None
        out[name] = dict(regions=n, agree=k, share=k / n if n else float("nan"),
                         low=float(ci.low) if ci else float("nan"), high=float(ci.high) if ci else float("nan"))
    return out


def edger_amount(counts_core, workdir):
    fd, path = tempfile.mkstemp(suffix=".tsv", dir=workdir)
    with os.fdopen(fd, "w") as fh:
        fh.write("\t".join(f"s{i}" for i in range(counts_core.shape[0])) + "\n")
        for r in range(counts_core.shape[1]):
            fh.write(f"r{r}\t" + "\t".join(str(int(v)) for v in counts_core[:, r]) + "\n")
    subprocess.run(["Rscript", str(Path(workdir) / "edger.R"), path, "A,A,A,B,B,B"], capture_output=True, text=True, check=True)
    res = json.loads(Path(path.replace(".tsv", ".json")).read_text())
    p = np.full(counts_core.shape[1], np.nan)
    p[[int(x[1:]) for x in res["region"]]] = np.array(res["p"], float).ravel()
    os.unlink(path)
    os.unlink(path.replace(".tsv", ".json"))
    return p


def main(argv):
    cache, out = argv[0], Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *x: subprocess.run(["git", *x], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted=len(git("status", "--porcelain", "--", "src", "benchmarks/*.py").splitlines()),
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"), plan_sha256=sha(PLAN), script_sha256=sha(__file__),
                started=time.strftime("%Y-%m-%dT%H:%M:%S"), exploratory="S3 not met; confirmatory case not evaluable")
    d, prof, _ = _load(cache)
    ids = [s["sample_id"] for s in prof.samples]
    qc = json.loads(str(d["qc"]))
    lib = np.array([qc[i]["library_size"] for i in ids])
    counts = np.rint(prof.values * lib[:, None, None] / 1e6)  # exact: values are CPM per bin
    edges = prof.edges
    R = dict(provenance=prov)
    width, cen = shape_measures(prof, edges, CFG["core_flank"])
    # (i) aggregate per replicate
    core = M.core_mask(edges, CFG["core_flank"])
    e = edges[np.r_[np.flatnonzero(core), np.flatnonzero(core)[-1] + 1]]
    agg = {}
    for k, i in enumerate(ids):
        m = np.clip(np.nan_to_num(prof.values[k][:, core]), 0, None).sum(0)
        p = M.normalise(m * M.widths(e))
        cdf = np.cumsum(p)
        c = M.centers(e)
        agg[i] = dict(centroid_bp=float(M.centroid(p, e)), mad_bp=float(M.mean_abs_deviation(p, e)),
                      iqr_bp=float(np.interp(0.75, cdf, c) - np.interp(0.25, cdf, c)))
    R["aggregate"] = agg
    # (ii) 3 vs 3 shape calls against amount calls
    called, tested = calls(prof)
    work = tempfile.mkdtemp(prefix="appC_", dir=out)
    (Path(work) / "dsp.R").write_text(IV.R_SCRIPT)
    (Path(work) / "edger.R").write_text(EDGER)
    amt = bh(np.nan_to_num(edger_amount(counts[..., core].sum(-1), work), nan=1.0)) < 0.05
    dsp = bh(np.nan_to_num(diffsplice(counts, ["A"] * 3 + ["B"] * 3, edges, CFG["core_flank"], work), nan=1.0)) < 0.05
    for f in Path(work).iterdir():
        f.unlink()
    os.rmdir(work)
    R["calls_3v3"] = dict(tested=int(tested.sum()), shape_any=int(called.sum()), amount=int(amt.sum()),
                          shape_without_amount=int((called & ~amt).sum()), share_shape_without_amount=float((called & ~amt).sum() / max(called.sum(), 1)),
                          amount_without_shape=int((amt & ~called).sum()), diffsplice=int(dsp.sum()),
                          shape_and_diffsplice=int((called & dsp).sum()), diffsplice_without_amount=int((dsp & ~amt).sum()))
    # (iii) held-out sign agreement
    keep = [ids.index(i) for i in ("DMSO_1", "DMSO_2", "CDK7i_1", "CDK7i_2")]
    c12, t12 = calls(_subset(prof, keep))
    R["heldout"] = dict(called_on_reps12=int(c12.sum()), tested_on_reps12=int(t12.sum()),
                        width=agreement(width, ids, c12), centroid=agreement(cen, ids, c12),
                        width_not_called=agreement(width, ids, t12 & ~c12), centroid_not_called=agreement(cen, ids, t12 & ~c12))
    ok_amount = R["calls_3v3"]["share_shape_without_amount"] >= 0.25
    ok_agree = all(R["heldout"][m][s]["low"] > 0.5 for m in ("width", "centroid") for s in ("negative", "positive"))
    R["grade"] = "consistent" if ok_amount and ok_agree else ("partly consistent" if ok_amount or ok_agree else "not supported")
    prov["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "scenario_c.json").write_text(json.dumps(R, indent=1, default=float))
    c3 = R["calls_3v3"]
    fa = lambda a: "; ".join(f"{s} {a[s]['agree']}/{a[s]['regions']} = {a[s]['share']:.3f} [{a[s]['low']:.3f}, {a[s]['high']:.3f}]" for s in ("negative", "positive"))  # noqa: E731
    L = ["# Scenario C · changes that window counts and aggregate profiles miss (GSE218269) — EXPLORATORY", "",
         "S3 (depth) not met; the confirmatory case is not evaluable and stays so. "
         f"Plan sha256 `{prov['plan_sha256']}`; commit `{prov['commit'][:10]}` (uncommitted code files {prov['uncommitted']}); regional.py "
         f"`{prov['regional_py_sha256']}`; {prov['started']} → {prov['finished']}.", "",
         "(i) Aggregate core profile per replicate: " + "; ".join(f"{k} centroid {v['centroid_bp']:+.1f} bp, MAD {v['mad_bp']:.1f}, IQR {v['iqr_bp']:.0f}" for k, v in agg.items()), "",
         f"(ii) 3 vs 3: `shape_any` {c3['shape_any']} / {c3['tested']}; edgeR amount {c3['amount']}; shape calls without an amount call "
         f"{c3['shape_without_amount']} ({c3['share_shape_without_amount']:.2f}); amount without shape {c3['amount_without_shape']}; diffSpliceDGE "
         f"{c3['diffsplice']} ({c3['shape_and_diffsplice']} shared with `shape_any`; {c3['diffsplice_without_amount']} without an amount call).", "",
         f"(iii) Held-out sign agreement (called on replicates 1–2: {R['heldout']['called_on_reps12']} / {R['heldout']['tested_on_reps12']}; checked on replicate 3), share [95 % exact]:", "",
         f"- width, called: {fa(R['heldout']['width'])}",
         f"- centroid, called: {fa(R['heldout']['centroid'])}",
         f"- width, not called (reference): {fa(R['heldout']['width_not_called'])}",
         f"- centroid, not called (reference): {fa(R['heldout']['centroid_not_called'])}", "",
         f"**Grade: {R['grade']}** (criteria in the plan, clarification for C). Exploratory; below the planned depth; no mechanism claimed."]
    (out / "SCENARIO_C.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1:])
