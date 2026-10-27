"""Real positive cases re-run with the current code (v3.1): yeast remodeler MNase-seq and NELFB PRO-seq.

Both datasets were used during development, so these are biological-consistency cases, not independent
validation. The script replaces the v2.4-era numbers in docs/PROFILESHAPE_COMPARISON.md §6.2–6.3 with results
of the current code (covariate-stratified null, shape_any, cross-correlation translation estimate with
consistency flags), next to descriptive checks that do not use ProfileShape tests.

Inputs (profile caches built from GEO bigWigs with `profileshape.data.extract`):
  yeast   GSE280671 mononucleosome CPM bigWigs; 3,592 ORFs ≥ 1 kb, ATG-anchored, oriented, ±1.5 kb, 10 bp bins;
          cache `values` samples × regions × bins in the order of the manifest (WT 1–4, dKO 1–4, chd1 1–2, isw1 1–2)
  PRO-seq GSE233930 NELFB-FKBP12F36V HAP1, sense-strand bigWigs; RefSeq Select TSSs (genes ≥ 2 kb), ±1 kb,
          10 bp bins; cache order DMSO 1–2, VHL 2 h 1–2, VHL 48 h 1–2 (plus-strand genes then minus-strand genes)

    python benchmarks/positive_cases.py --yeast-cache Y.npz --yeast-manifest M.tsv \
        --proseq-cache P.npz --out DIR
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from scipy import stats

from profileshape.analysis import analyse
from profileshape.data import Config, Profiles
from profileshape.regional import per_region_tests

YEAST_CFG = dict(anchor="start", orientation="oriented", flank=1500, binsize=10, background_width=300, core_flank=1000,
                 min_mass_quantile=0.1)
PROSEQ_CFG = dict(anchor="start", orientation="oriented", flank=1000, binsize=10, background_width=300, core_flank=500)
PROSEQ_ORDER = [("DMSO", 1), ("DMSO", 2), ("VHL2h", 1), ("VHL2h", 2), ("VHL48h", 1), ("VHL48h", 2)]


def summarise(prof: Profiles, ref: str, alt: str, cfg: dict) -> dict:
    a = analyse(prof, Config(reference=ref, alternative=alt, **cfg), bootstrap=0)
    r = per_region_tests(a)
    out = {}
    for m in ("shape_any", "shape_kernel", "shift"):
        t = r.tests[m]
        out[m] = dict(called=int(np.isin(t["call"], ("Higher", "Lower", "Changed")).sum()), tested=int(len(t["call"])))
    t = r.tests["shift"]
    hi, lo, ch = (t["call"] == "Higher"), (t["call"] == "Lower"), (t["call"] == "Changed")
    est = np.asarray(t["estimate"], float)
    flags = [str(f) for f in np.asarray(t["translation_flag"])[ch]]
    out["shift"].update(higher=int(hi.sum()), lower=int(lo.sum()), changed_no_direction=int(ch.sum()),
                        downstream_share=float(hi.sum() / max(hi.sum() + lo.sum(), 1)),
                        median_bp_consistent=float(np.nanmedian(est[hi | lo])) if (hi | lo).any() else float("nan"),
                        iqr_bp_consistent=[float(np.nanpercentile(est[hi | lo], 25)), float(np.nanpercentile(est[hi | lo], 75))]
                        if (hi | lo).any() else [float("nan")] * 2,
                        flags={f: flags.count(f) for f in sorted(set(flags))},
                        high_r2=int((np.asarray(t["translation_fit"]) == "high_r2").sum()))
    out["notes"] = [n for n in r.notes if "strata" in n or "Changed" in n or "depth" in n.lower()]
    out["influence_flags"] = [str(f) for f in r.influence.get("flags", [])]
    return out


def yeast_job(job):
    label, ref_ids, alt_ids, cache, manifest = job
    d = np.load(cache)
    ids = [row["sample_id"] for row in csv.DictReader(open(manifest), delimiter="\t")]
    use = [ids.index(i) for i in ref_ids + alt_ids]
    smp = [dict(sample_id=i, condition="A" if i in ref_ids else "B") for i in ref_ids + alt_ids]
    regions = [dict(chrom="c", start=k, end=k + 1, name=f"r{k}", strand="+") for k in range(d["values"].shape[1])]
    return label, summarise(Profiles(d["values"][use], d["edges"], smp, regions, True), "A", "B", YEAST_CFG)


def phasing(values: np.ndarray, edges: np.ndarray) -> list[float]:
    """Autocorrelation of each sample's mean mass-normalised profile at the nucleosome repeat (150–190 bp lag),
    a phasing measure that ProfileShape does not use."""
    c = 0.5 * (edges[:-1] + edges[1:])
    out = []
    for v in values:
        p = np.clip(np.nan_to_num(v), 0, None)
        p = p / np.maximum(p.sum(-1, keepdims=True), 1e-12)
        m = p.mean(0)[(c > 0) & (c < 1200)]
        m = m - m.mean()
        ac = [float(np.dot(m[:-k], m[k:]) / np.dot(m, m)) for k in range(15, 20)]
        out.append(max(ac))
    return out


def proseq_job(job):
    alt, cache = job
    d = np.load(cache)
    use = [k for k, (cnd, _) in enumerate(PROSEQ_ORDER) if cnd in ("DMSO", alt)]
    smp = [dict(sample_id=f"{PROSEQ_ORDER[k][0]}_rep{PROSEQ_ORDER[k][1]}", condition=PROSEQ_ORDER[k][0]) for k in use]
    regions = [dict(chrom="c", start=k, end=k + 1, name=f"g{k}", strand="+") for k in range(d["values"].shape[1])]
    return alt, summarise(Profiles(d["values"][use], d["edges"], smp, regions, True), "DMSO", alt, PROSEQ_CFG)


def proseq_meta(cache: str) -> dict:
    """Pause-peak position and 0–300 bp centroid of each sample's meta-profile over the 2,000 promoters with the
    highest pooled signal (descriptive; no ProfileShape test)."""
    d = np.load(cache)
    v, e = d["values"], d["edges"]
    c = 0.5 * (e[:-1] + e[1:])
    prom = np.nansum(v[:, :, (c > 0) & (c < 300)], axis=-1).mean(0)
    top = np.argsort(prom)[-2000:]
    out = {}
    for k, (cnd, i) in enumerate(PROSEQ_ORDER):
        p = np.nan_to_num(v[k, top]); p = p / np.maximum(p.sum(-1, keepdims=True), 1e-9)
        m = p.mean(0); sel = (c > -50) & (c < 400); w = (c > 0) & (c < 300)
        out[f"{cnd}_rep{i}"] = dict(peak_bp=float(c[sel][np.argmax(m[sel])]), centroid_0_300=float((m * c)[w].sum() / m[w].sum()))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--yeast-cache", required=True)
    ap.add_argument("--yeast-manifest", required=True)
    ap.add_argument("--proseq-cache", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    git = lambda *a: subprocess.run(["git", *a], cwd=root, capture_output=True, text=True).stdout.strip()  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                regional_py_sha256=hashlib.sha256((root / "src/profileshape/regional.py").read_bytes()).hexdigest()[:16],
                started=datetime.datetime.now().isoformat(timespec="seconds"))
    wt = [f"WT_rep{i}" for i in range(1, 5)]
    dk = [f"dKO_rep{i}" for i in range(1, 5)]
    yj = [("WT(4) vs isw1Δ chd1Δ(4)", wt, dk), ("WT(4) vs chd1Δ(2)", wt, ["chd1_rep1", "chd1_rep2"]),
          ("WT(4) vs isw1Δ(2)", wt, ["isw1_rep1", "isw1_rep2"])]
    yj += [(f"null WT {a}+{b} vs rest", [a, b], [x for x in wt if x not in (a, b)]) for a, b in (("WT_rep1", "WT_rep2"), ("WT_rep1", "WT_rep3"), ("WT_rep1", "WT_rep4"))]
    jobs_y = [(lab, r, a, args.yeast_cache, args.yeast_manifest) for lab, r, a in yj]
    with ProcessPoolExecutor(4) as ex:
        yres = dict(ex.map(yeast_job, jobs_y))
        pres = dict(ex.map(proseq_job, [("VHL2h", args.proseq_cache), ("VHL48h", args.proseq_cache)]))
    d = np.load(args.yeast_cache)
    ids = [row["sample_id"] for row in csv.DictReader(open(args.yeast_manifest), delimiter="\t")]
    ph = dict(zip(ids, phasing(d["values"], d["edges"]), strict=True))
    meta = proseq_meta(args.proseq_cache)
    prov["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    (out / "positive_cases.json").write_text(json.dumps(dict(provenance=prov, yeast=yres, yeast_phasing=ph, proseq=pres,
                                                             proseq_meta=meta), indent=1))
    L = ["# Real positive cases with ProfileShape v3.1 (development datasets; biological consistency, not validation)", "",
         f"Commit `{prov['commit'][:10]}`; uncommitted src files {prov['uncommitted_src']}; regional.py `{prov['regional_py_sha256']}`; "
         f"{prov['started']} → {prov['finished']}.", "",
         "## Yeast MNase-seq (GSE280671), ATG-anchored ORFs", "",
         "| comparison | shape_any | shape_kernel | shift (called) | consistent: downstream / upstream | Changed (no direction) | median consistent shift, bp [IQR] |",
         "|---|---|---|---|---|---|---|"]
    for lab, r in yres.items():
        s = r["shift"]
        L.append(f"| {lab} | {r['shape_any']['called']} / {r['shape_any']['tested']} | {r['shape_kernel']['called']} | {s['called']} | "
                 f"{s['higher']} / {s['lower']} | {s['changed_no_direction']} | {s['median_bp_consistent']:.1f} [{s['iqr_bp_consistent'][0]:.1f}, {s['iqr_bp_consistent'][1]:.1f}] |")
    groups = {}
    for k, v in ph.items():
        groups.setdefault(k.rsplit("_rep", 1)[0], []).append(v)
    L += ["", "Phasing (meta-profile autocorrelation at the nucleosome repeat; not used by ProfileShape): "
          + "; ".join(f"{g} {np.round(v, 2).tolist()}" for g, v in groups.items()), "",
          "## PRO-seq (GSE233930), NELFB degradation vs DMSO, 2 vs 2, sense strand at TSSs", "",
          "| comparison | shape_any | shape_kernel | shift (called) | consistent: downstream / upstream (share downstream) | Changed (no direction) | median consistent shift, bp [IQR] |",
          "|---|---|---|---|---|---|---|"]
    for alt, r in pres.items():
        s = r["shift"]
        L.append(f"| DMSO vs {alt} | {r['shape_any']['called']} / {r['shape_any']['tested']} | {r['shape_kernel']['called']} | {s['called']} | "
                 f"{s['higher']} / {s['lower']} ({s['downstream_share']:.2f}) | {s['changed_no_direction']} | "
                 f"{s['median_bp_consistent']:.1f} [{s['iqr_bp_consistent'][0]:.1f}, {s['iqr_bp_consistent'][1]:.1f}] |")
    L += ["", "Meta-profile over the top 2,000 promoters (descriptive; not used by ProfileShape): "
          + "; ".join(f"{k} peak {v['peak_bp']:.0f} bp, centroid {v['centroid_0_300']:.1f} bp" for k, v in meta.items()), ""]
    L += ["Notes: " + " | ".join(f"{k}: " + "; ".join(v["notes"]) for k, v in {**yres, **pres}.items() if v["notes"]), ""]
    (out / "POSITIVE_CASES.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
