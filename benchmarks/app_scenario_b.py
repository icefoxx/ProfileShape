"""Scenario B of docs/PROFILESHAPE_APPLICATIONS_PLAN.md: complex shape change at nucleosome arrays (yeast MNase).

    PYTHONPATH=benchmarks python benchmarks/app_scenario_b.py PROFILES.npz MANIFEST.tsv OUT_DIR

Development dataset (GSE280671), WT (4) vs isw1Δ chd1Δ (4); an application demonstration. Independent measure:
per-gene phasing (autocorrelation at the nucleosome repeat, lags 150–190 bp, gene-body side), not used by
ProfileShape. Held-out design: genes called on replicates 1–2, phasing change measured on replicates 3–4.
"""

from __future__ import annotations

import csv
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

import independent_validation as IV
from app_scenario_a import diffsplice
from positive_cases import YEAST_CFG
from profileshape import metrics as M
from profileshape.analysis import analyse
from profileshape.data import Config, Profiles, bigwig_count_unit
from profileshape.regional import bh, per_region_tests

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "docs" / "PROFILESHAPE_APPLICATIONS_PLAN.md"
SEED = 9_910_000
WT = [f"WT_rep{i}" for i in range(1, 5)]
DKO = [f"dKO_rep{i}" for i in range(1, 5)]
EDGER = r"""
suppressMessages(library(edgeR))
args <- commandArgs(TRUE); x <- as.matrix(read.table(args[1], header = TRUE, row.names = 1, check.names = FALSE))
g <- factor(strsplit(args[2], ",")[[1]], levels = c("A", "B")); d <- model.matrix(~g)
y <- normLibSizes(DGEList(x)); y <- estimateDisp(y, d); f <- glmQLFit(y, d); t <- glmQLFTest(f, coef = 2)
writeLines(jsonlite::toJSON(list(region = rownames(x), p = as.numeric(t$table$PValue)), digits = NA), sub("\\.tsv$", ".json", args[1]))
"""


def phasing(values, edges):
    """Per gene and sample: max autocorrelation of the mass-normalised gene-body profile at lags 150–190 bp."""
    c = M.centers(edges)
    v = np.clip(np.nan_to_num(values[..., (c > 0) & (c < 1200)]), 0, None)
    p = v / np.where(v.sum(-1, keepdims=True) > 0, v.sum(-1, keepdims=True), np.nan)
    m = p - np.nanmean(p, -1, keepdims=True)
    den = (m * m).sum(-1)
    ac = np.stack([(m[..., :-k] * m[..., k:]).sum(-1) / den for k in range(15, 20)], -1)
    return np.nanmax(ac, -1)


def boot_diff(x, called, rng, n=2000):
    a, b = x[called & np.isfinite(x)], x[~called & np.isfinite(x)]
    if not len(a) or not len(b):
        return dict(called_median=float("nan"), not_called_median=float("nan"), diff=float("nan"), low=float("nan"), high=float("nan"),
                    n_called=int(len(a)), n_not=int(len(b)))
    d = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return dict(called_median=float(np.median(a)), not_called_median=float(np.median(b)), diff=float(np.median(a) - np.median(b)),
                low=float(np.quantile(d, 0.025)), high=float(np.quantile(d, 0.975)), n_called=int(len(a)), n_not=int(len(b)))


def run_tests(values, edges, ids, ref, alt):
    use = [ids.index(i) for i in ref + alt]
    smp = [dict(sample_id=i, condition="A" if i in ref else "B") for i in ref + alt]
    regions = [dict(chrom="c", start=k, end=k + 1, name=f"r{k}", strand="+") for k in range(values.shape[1])]
    res = per_region_tests(analyse(Profiles(values[use], edges, smp, regions, True), Config(reference="A", alternative="B", **YEAST_CFG), bootstrap=0))
    called = np.zeros(values.shape[1], bool)
    t = res.tests["shape_any"]
    called[np.asarray(t["region_index"])[np.isin(t["call"], ("Changed",))]] = True
    tested = np.zeros(values.shape[1], bool)
    tested[np.asarray(t["region_index"])[t["call"] != "not_tested"]] = True
    return res, called, tested


def edger_amount(cnt, groups, edges, workdir):
    core = M.core_mask(edges, YEAST_CFG["core_flank"])
    tot = cnt[..., core].sum(-1)
    fd, path = tempfile.mkstemp(suffix=".tsv", dir=workdir)
    with os.fdopen(fd, "w") as fh:
        fh.write("\t".join(f"s{i}" for i in range(tot.shape[0])) + "\n")
        for r in range(tot.shape[1]):
            fh.write(f"r{r}\t" + "\t".join(str(int(v)) for v in tot[:, r]) + "\n")
    subprocess.run(["Rscript", str(Path(workdir) / "edger.R"), path, ",".join(groups)], capture_output=True, text=True, check=True)
    res = json.loads(Path(path.replace(".tsv", ".json")).read_text())
    p = np.full(tot.shape[1], np.nan)
    p[[int(x[1:]) for x in res["region"]]] = np.array(res["p"], float).ravel()
    os.unlink(path)
    os.unlink(path.replace(".tsv", ".json"))
    return p


def main(argv):
    cache, manifest, out = argv[0], argv[1], Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *x: subprocess.run(["git", *x], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted=len(git("status", "--porcelain", "--", "src", "benchmarks").splitlines()),
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"), plan_sha256=sha(PLAN), script_sha256=sha(__file__),
                seed=SEED, started=time.strftime("%Y-%m-%dT%H:%M:%S"))
    d = np.load(cache)
    values, edges = d["values"], d["edges"]
    rows = list(csv.DictReader(open(manifest), delimiter="\t"))
    ids = [r["sample_id"] for r in rows]
    rng = np.random.default_rng(SEED)
    ph = phasing(values, edges)  # samples × genes
    wt, dk = [ids.index(i) for i in WT], [ids.index(i) for i in DKO]
    dphase = np.nanmean(ph[dk], 0) - np.nanmean(ph[wt], 0)
    R = dict(provenance=prov)
    # (i) main 4 vs 4
    res, called, tested = run_tests(values, edges, ids, WT, DKO)
    sh = res.tests["shift"]
    flags = [str(f) for f in np.asarray(sh["translation_flag"])[sh["call"] == "Changed"]]
    R["main"] = dict(tested=int(tested.sum()), shape_any=int(called.sum()), shift_called=int(np.isin(sh["call"], ("Higher", "Lower", "Changed")).sum()),
                     shift_directional=int(np.isin(sh["call"], ("Higher", "Lower")).sum()), flags={f: flags.count(f) for f in sorted(set(flags))})
    # (ii) phasing change in called vs not-called (among tested)
    R["phasing_main"] = boot_diff(np.where(tested, dphase, np.nan), called, rng)
    R["phasing_per_sample_median"] = {i: float(np.nanmedian(ph[ids.index(i)])) for i in WT + DKO}
    # (iii) held-out: call on replicates 1–2, phasing change on 3–4
    _, c12, t12 = run_tests(values, edges, ids, WT[:2], DKO[:2])
    d34 = np.nanmean(ph[[ids.index(i) for i in DKO[2:]]], 0) - np.nanmean(ph[[ids.index(i) for i in WT[2:]]], 0)
    R["heldout"] = dict(called_on_reps12=int(c12.sum()), tested_on_reps12=int(t12.sum()), **boot_diff(np.where(t12, d34, np.nan), c12, rng))
    # comparators (fragment equivalents from the CPM coverage tracks)
    unit = {r["sample_id"]: bigwig_count_unit(r["signal_file"]) for r in rows if r["sample_id"] in WT + DKO}
    cnt = np.stack([np.rint(np.nan_to_num(values[ids.index(i)]) * np.diff(edges) / (unit[i][0] * unit[i][1])) for i in WT + DKO])
    work = tempfile.mkdtemp(prefix="appB_", dir=out)
    (Path(work) / "dsp.R").write_text(IV.R_SCRIPT)
    (Path(work) / "edger.R").write_text(EDGER)
    groups = ["A"] * 4 + ["B"] * 4
    dsp = bh(np.nan_to_num(diffsplice(cnt, groups, edges, YEAST_CFG["core_flank"], work), nan=1.0)) < 0.05
    amt = bh(np.nan_to_num(edger_amount(cnt, groups, edges, work), nan=1.0)) < 0.05
    R["comparators"] = dict(units={k: list(v) for k, v in unit.items()},
                            diffsplice=dict(called=int(dsp.sum()), overlap_with_shape_any=int((dsp & called).sum()),
                                            phasing=boot_diff(np.where(tested, dphase, np.nan), dsp, rng)),
                            edger_amount=dict(called=int(amt.sum()), overlap_with_shape_any=int((amt & called).sum()),
                                              shape_any_without_amount=int((called & ~amt).sum()),
                                              phasing=boot_diff(np.where(tested, dphase, np.nan), amt, rng)))
    for f in Path(work).iterdir():
        f.unlink()
    os.rmdir(work)
    pm, ho = R["phasing_main"], R["heldout"]
    ok2 = pm["high"] < 0
    ok3 = np.sign(ho["diff"]) == np.sign(pm["diff"])
    R["grade"] = "consistent" if ok2 and ok3 else ("partly consistent" if ok2 or ok3 else "not supported")
    prov["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "scenario_b.json").write_text(json.dumps(R, indent=1, default=float))
    f = lambda r: f"{r['diff']:+.3f} [{r['low']:+.3f}, {r['high']:+.3f}] (called {r['called_median']:+.3f}, n = {r['n_called']}; not {r['not_called_median']:+.3f}, n = {r['n_not']})"  # noqa: E731
    m = R["main"]
    L = ["# Scenario B · complex shape change at nucleosome arrays (yeast MNase, development data)", "",
         f"Plan sha256 `{prov['plan_sha256']}`; commit `{prov['commit'][:10]}` (uncommitted {prov['uncommitted']}); regional.py "
         f"`{prov['regional_py_sha256']}`; seed {SEED}; {prov['started']} → {prov['finished']}.", "",
         f"(i) WT (4) vs isw1Δ chd1Δ (4): {m['shape_any']} / {m['tested']} genes `shape_any`; `shift` {m['shift_called']} called, "
         f"{m['shift_directional']} with a direction; reasons for the rest {m['flags']}.", "",
         "Per-gene phasing change (dKO − WT; autocorrelation at 150–190 bp), difference of medians called − not called [95 % bootstrap]:", "",
         f"- (ii) main, 4 vs 4: {f(pm)}",
         f"- (iii) held-out: called on replicates 1–2 ({ho['called_on_reps12']} / {ho['tested_on_reps12']}), phasing on 3–4: {f(ho)}",
         f"- diffSpliceDGE ({R['comparators']['diffsplice']['called']} calls, {R['comparators']['diffsplice']['overlap_with_shape_any']} shared with `shape_any`): {f(R['comparators']['diffsplice']['phasing'])}",
         f"- edgeR window counts, amount ({R['comparators']['edger_amount']['called']} calls; `shape_any` calls without an amount call: "
         f"{R['comparators']['edger_amount']['shape_any_without_amount']}): {f(R['comparators']['edger_amount']['phasing'])}", "",
         "Median per-sample phasing: " + ", ".join(f"{k} {v:.2f}" for k, v in R["phasing_per_sample_median"].items()) + ".", "",
         f"**Grade: {R['grade']}** (criteria in the plan). Development data; counts for the comparators are fragment equivalents "
         "from CPM coverage tracks (approximation); held-out calls use 2 vs 2, where the yeast null is known to exceed nominal."]
    (out / "SCENARIO_B.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1:])
