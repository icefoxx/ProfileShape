"""Prerequisite R of docs/PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md: midpoint vs fragment coverage under real
nulls, on the same reads (amendment R1).

    PYTHONPATH=benchmarks python benchmarks/representation_null.py --cutrun-bams DIR --cutrun-anchors cutrun_profiles.npz \
        --yeast-bams DIR --yeast-bed orf_atg.bed --atac-cache psv_cache.npz --out DIR

CUT&RUN A673: ten 3 vs 3 splits; yeast MNase WT: three 2 vs 2 splits; ATAC Sth1 0 h / 1 h: three 2 vs 2 splits each,
three anchor sets (as in its validation). Per split and test: whole-list calls (FDR 0.05) and share of p < 0.01.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

import independent_validation as IV
from profileshape.analysis import analyse
from profileshape.data import Config, Profiles, extract, read_regions
from profileshape.regional import per_region_tests

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md"
TESTS = ("shape_any", "shape_kernel", "shift")
MODES = ("midpoint", "coverage")
CUTRUN_CFG = dict(anchor="center", orientation="oriented", flank=1000, binsize=10, background_width=300, core_flank=600)
YEAST_CFG = dict(anchor="start", orientation="oriented", flank=1500, binsize=10, background_width=300, core_flank=1000)
P01_MARGIN = 0.005  # coverage "higher" on CUT&RUN / ATAC if its mean p < 0.01 share exceeds midpoints' by more


def evaluate(prof: Profiles, cfg: dict) -> dict:
    res = per_region_tests(analyse(prof, Config(reference="A", alternative="B", **cfg), bootstrap=0))
    out = {}
    for m in TESTS:
        t = res.tests[m]
        p = np.asarray(t["p"], float)
        out[m] = dict(calls=int(np.isin(t["call"], ("Higher", "Lower", "Changed")).sum()), tested=int(np.isfinite(p).sum()),
                      p01=float(np.mean(p[np.isfinite(p)] < 0.01)) if np.isfinite(p).any() else float("nan"))
    return out


def _relabel(prof: Profiles, a_ids, b_ids) -> Profiles:
    ids = [s["sample_id"] for s in prof.samples]
    use = [ids.index(i) for i in list(a_ids) + list(b_ids)]
    smp = [dict(prof.samples[i], condition="A" if prof.samples[i]["sample_id"] in a_ids else "B") for i in use]
    return Profiles(prof.values[use], prof.edges, smp, prof.regions, True, provenance=prof.provenance,
                    fragments=None if prof.fragments is None else prof.fragments[use])


def _splits(ids, k):
    return [(list(a), [x for x in ids if x not in a]) for a in itertools.combinations(ids, k) if ids[0] in a]


# ------------------------------------------------------------------ BAM datasets (CUT&RUN, yeast)
def bam_profiles(job):
    name, mode, samples, regions, cfg = job
    prof = extract(Config(**cfg, bam_mode=mode), samples, regions)
    return name, mode, prof


def bam_split_job(job):
    name, mode, prof, cfg, a_ids, b_ids = job
    return name, mode, "+".join(a_ids), evaluate(_relabel(prof, a_ids, b_ids), cfg)


# ------------------------------------------------------------------ ATAC (fragments from the validation cache)
class CovData(IV.Data):
    """Validation fragments binned as coverage: each fragment's span (midpoint ± length / 2), bp overlap per bin."""

    def counts(self, names, anchor, size, modify=None):
        ok = np.flatnonzero(np.isfinite(anchor))
        idx = np.full(len(self.genes), -1)
        idx[ok] = np.arange(len(ok))
        F, B = IV.FLANK, IV.BIN
        nb = 2 * F // B
        out = np.zeros((len(names), len(ok), nb))
        for i, n in enumerate(names):
            reg, pos, ln = self.frag[n]
            m = (ln >= size[0]) & (ln <= size[1]) & (idx[reg] >= 0)
            reg, pos, ln = reg[m], pos[m].astype(float), ln[m].astype(float)
            rel = pos - anchor[reg]
            lo = np.clip(np.floor(rel - ln / 2).astype(int) + F, 0, 2 * F)
            hi = np.clip(np.floor(rel + ln / 2).astype(int) + F, 0, 2 * F)
            diff = np.zeros((len(ok), 2 * F + 1))
            np.add.at(diff, (idx[reg], lo), 1.0)
            np.add.at(diff, (idx[reg], hi), -1.0)
            cov = np.cumsum(diff, axis=1)[:, : 2 * F]
            out[i] = cov.reshape(len(ok), nb, B).mean(-1)
        return out, ok.tolist()


def atac_job(job):
    cache, mode, anchor_set, cond, a_ids, b_ids = job
    D = IV.Data(cache) if mode == "midpoint" else CovData(cache)
    anc = IV.Data.anchors(D)[anchor_set]  # anchors always from midpoints (condition-blind, as in the validation)
    size = IV.NUC if anchor_set != "ndr" else IV.SUB
    prof, _ = D.profiles(a_ids + b_ids, ["A"] * 2 + ["B"] * 2, anc, size)
    return "atac", mode, f"{anchor_set}|{cond}|{'+'.join(a_ids)}", evaluate(prof, IV.CFG)


# ------------------------------------------------------------------ verdict and report
def verdict(rows: list[dict]) -> dict:
    out = {}
    for ds in sorted({r["dataset"] for r in rows}):
        v = {}
        for m in TESTS:
            mid = {r["split"]: r["res"][m] for r in rows if r["dataset"] == ds and r["mode"] == "midpoint"}
            cov = {r["split"]: r["res"][m] for r in rows if r["dataset"] == ds and r["mode"] == "coverage"}
            keys = sorted(mid)
            more_calls = sum(cov[k]["calls"] > mid[k]["calls"] for k in keys)
            dp01 = float(np.nanmean([cov[k]["p01"] for k in keys]) - np.nanmean([mid[k]["p01"] for k in keys]))
            if ds == "yeast":
                higher = more_calls >= 2
            else:
                higher = dp01 > P01_MARGIN
            if ds == "cutrun":
                zero = sum(cov[k]["calls"] == 0 for k in keys)
                higher = higher or zero < 9
                v[m] = dict(splits=len(keys), coverage_more_calls=more_calls, delta_p01=dp01, coverage_zero_call_splits=zero, higher=higher)
            else:
                v[m] = dict(splits=len(keys), coverage_more_calls=more_calls, delta_p01=dp01, higher=higher)
        out[ds] = dict(tests=v, coverage_higher=any(x["higher"] for x in v.values()))
    out["coverage_eligible"] = not any(out[ds]["coverage_higher"] for ds in out if ds != "coverage_eligible")
    return out


def report(R) -> str:
    p = R["provenance"]
    L = ["# Prerequisite R · midpoint vs coverage under real nulls (same reads)", "",
         f"Protocol `docs/PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md` sha256 `{p['protocol_sha256']}` (amendment R1); commit "
         f"`{p['commit'][:10]}` (uncommitted {p['uncommitted']}); regional.py `{p['regional_py_sha256']}`; script `{p['script_sha256']}`; "
         f"{p['started']} → {p['finished']}. Regions: {p['regions']}.", "",
         "Per split: whole-list calls at FDR 0.05 / share p < 0.01, for midpoint | coverage.", "",
         "| dataset | split | " + " | ".join(f"{m} mid" for m in TESTS) + " | " + " | ".join(f"{m} cov" for m in TESTS) + " |",
         "|---|---|" + "---|" * (2 * len(TESTS))]
    by = {}
    for r in R["rows"]:
        by.setdefault((r["dataset"], r["split"]), {})[r["mode"]] = r["res"]
    for (ds, sp), d in sorted(by.items()):
        cells = [f"{d[mo][m]['calls']} / {d[mo][m]['p01']:.4f}" for mo in MODES for m in TESTS]
        L.append(f"| {ds} | {sp} | " + " | ".join(cells) + " |")
    L += ["", "## Verdict", "", "| dataset | test | coverage more calls (splits) | Δ p<0.01 share (cov − mid) | other | coverage higher |",
          "|---|---|---|---|---|---|"]
    for ds, d in R["verdict"].items():
        if ds == "coverage_eligible":
            continue
        for m, v in d["tests"].items():
            other = f"zero-call splits {v['coverage_zero_call_splits']}/10" if "coverage_zero_call_splits" in v else ""
            L.append(f"| {ds} | {m} | {v['coverage_more_calls']} / {v['splits']} | {v['delta_p01']:+.4f} | {other} | {'**yes**' if v['higher'] else 'no'} |")
    L += ["", f"**Coverage eligible for stage 2: {'yes' if R['verdict']['coverage_eligible'] else 'no'}** "
          "(not higher than midpoints on any dataset and test). Yeast midpoints are known to exceed nominal at 2 vs 2; "
          "this comparison is relative. CUT&RUN uses a 4 M-pair, chr1 + chr2 subsample (R1)."]
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    for a in ("--cutrun-bams", "--cutrun-anchors", "--yeast-bams", "--yeast-bed", "--atac-cache", "--out"):
        ap.add_argument(a, required=True)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *x: subprocess.run(["git", *x], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted=len(git("status", "--porcelain", "--", "src", "benchmarks").splitlines()),
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"), protocol_sha256=sha(PROTOCOL),
                script_sha256=sha(__file__), started=time.strftime("%Y-%m-%dT%H:%M:%S"))
    names = [str(x) for x in np.load(a.cutrun_anchors, allow_pickle=True)["anchors"]]
    cut_regions = [dict(chrom=n.split(":")[0], start=int(n.split(":")[1]), end=int(n.split(":")[1]) + 1, name=n, strand="+")
                   for n in names if n.split(":")[0] in ("chr1", "chr2")]
    cut_samples = [dict(sample_id=f"A673_{k}", condition="A", signal_file=str(Path(a.cutrun_bams) / f"A673_{k}.bam")) for k in range(1, 7)]
    y_regions = read_regions(Path(a.yeast_bed))
    y_samples = [dict(sample_id=f"WT_rep{k}", condition="A", signal_file=str(Path(a.yeast_bams) / f"WT_rep{k}.bam")) for k in range(1, 5)]
    build = [("cutrun", m, cut_samples, cut_regions, CUTRUN_CFG) for m in MODES] + [("yeast", m, y_samples, y_regions, YEAST_CFG) for m in MODES]
    rows = []
    with ProcessPoolExecutor(6) as ex:
        profs = {(n, m): p for n, m, p in ex.map(bam_profiles, build)}
        prov["regions"] = {n: int(p.values.shape[1]) for (n, m), p in profs.items() if m == "midpoint"}
        prov["depth_median_fragments"] = {n: float(np.median(np.nanmedian(p.fragments, axis=1))) for (n, m), p in profs.items()
                                          if m == "midpoint" and p.fragments is not None}
        print("profiles built", prov["regions"], prov["depth_median_fragments"], flush=True)
        jobs = []
        for (n, m), p in profs.items():
            cfg = CUTRUN_CFG if n == "cutrun" else YEAST_CFG
            ids = [s["sample_id"] for s in p.samples]
            jobs += [(n, m, p, cfg, x, y) for x, y in _splits(ids, 3 if n == "cutrun" else 2)]
        for n, m, sp, res in ex.map(bam_split_job, jobs):
            rows.append(dict(dataset=n, mode=m, split=sp, res=res))
        print("BAM datasets done", flush=True)
        ajobs = [(a.atac_cache, m, s, cond, x, y) for m in MODES for s in ("plus1", "minus1", "ndr") for cond in ("0h", "1h")
                 for x, y in _splits(IV.STH1[cond], 2)]
        for n, m, sp, res in ex.map(atac_job, ajobs):
            rows.append(dict(dataset=n, mode=m, split=sp, res=res))
    prov["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    R = dict(provenance=prov, rows=rows, verdict=verdict(rows))
    (out / "representation_null.json").write_text(json.dumps(R, indent=1, default=float))
    (out / "REPRESENTATION_NULL.md").write_text(report(R))
    print(report(R))


if __name__ == "__main__":
    main()
