"""Stress tests of the pooled null beyond the calibration study (protocol: docs/PROFILESHAPE_STRESS_PROTOCOL.md).

    PYTHONPATH=benchmarks python benchmarks/stress_null.py --out DIR --templates REAL_TEMPLATES.npz [--reps 100]
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

import calibration_boundary as CB
from profileshape.analysis import analyse
from profileshape.data import Config, Profiles
from profileshape.regional import per_region_tests
from profileshape.simulate import SimSettings, simulate

ROOT = Path(__file__).resolve().parents[1]
SEED = 9_600_000
CELLS = ([("proportion_high", p) for p in (0.7, 0.9, 1.0)] + [("same_direction", p) for p in (0.5, 0.9)]
         + [("sample_effect", t) for t in (5, 10, 20)] + [("batch_partial", k) for k in ("undeclared", "declared")]
         + [("real_shapes", "GSE112560 templates"), ("tiling", "3 offsets")])
N = 600


def _shift_rows(v: np.ndarray, edges: np.ndarray, d: float) -> np.ndarray:
    c = 0.5 * (edges[:-1] + edges[1:])
    out = v.copy()
    for r in range(v.shape[0]):
        x = np.nan_to_num(v[r])
        out[r] = np.interp(c - d, c, x, left=x[0], right=x[-1])
    return out


def _templ_part(name, n, k, seed, templates):
    if not k:
        return None
    prof, truth = simulate(name, SimSettings(n_regions=k, n_ref=n, n_alt=n), seed=seed, templates=templates)
    lib = np.asarray(truth["library"], float)
    prof.fragments = np.rint(np.nan_to_num(prof.values) * lib[:, None, None]).sum(-1)
    return prof


def build(factor, level, n, seed, templates_path):
    frac = {"proportion_high": None, "same_direction": None}.get(factor, 0.2)
    if factor in ("proportion_high", "same_direction"):
        frac = float(level)
    m = int(round(frac * N))
    if factor == "same_direction":
        plan = (("shift", m), ("null", N - m))
    else:
        plan = (("shift", m // 2), ("broaden", m - m // 2), ("null", N - m))
    templates = np.load(templates_path)["templates"] if factor == "real_shapes" else None
    parts, changed = [], []
    for j, (name, cnt) in enumerate(plan):
        if cnt:
            parts.append(_templ_part(name, n, cnt, seed + 101 * j, templates) if templates is not None
                         else CB._part(name, n, cnt, seed + 101 * j))
            changed += [name != "null"] * cnt
    prof = CB._concat(parts)
    changed = np.array(changed)
    batches = None
    if factor == "tiling":
        step = float(np.median(np.diff(prof.edges)))
        copies = [prof]
        for off in (150, 300):
            sh = int(round(off / step))
            dup = np.concatenate([prof.values[..., sh:], np.repeat(prof.values[..., -1:], sh, axis=-1)], axis=-1)
            copies.append(Profiles(dup, prof.edges, prof.samples, prof.regions, True, fragments=prof.fragments))
        prof = CB._concat(copies)
        changed = np.r_[changed, changed, changed]
    if factor == "sample_effect":
        rng = np.random.default_rng(seed + 7)
        v = prof.values.copy()
        for s in range(v.shape[0]):
            v[s] = _shift_rows(v[s], prof.edges, rng.normal(0, float(level)))
        prof = Profiles(v, prof.edges, prof.samples, prof.regions, True, fragments=prof.fragments)
    if factor == "batch_partial":
        batches = ["A", "A", "B", "A", "B", "B"] if n == 3 else ["A", "A", "A", "B"]
        v = prof.values.copy()
        for s, b in enumerate(batches):
            if b == "B":
                v[s] = _shift_rows(v[s], prof.edges, 15.0)
        samples = [dict(s, batch=b) for s, b in zip(prof.samples, batches, strict=True)] if level == "declared" else prof.samples
        prof = Profiles(v, prof.edges, samples, prof.regions, True, fragments=prof.fragments)
    return prof, changed


def run(job):
    factor, level, n, rep, templates = job
    seed = SEED + 10_000 * [c[0] for c in CELLS].index(factor) + 1000 * rep + 100 * n + zlib.crc32(str(level).encode()) % 97
    prof, changed = build(factor, level, n, seed, templates)
    cfg = Config(reference="WT", alternative="KO",
                 batch_permutation="within" if (factor == "batch_partial" and level == "declared") else "free")
    a = analyse(prof, cfg, bootstrap=0)
    reg = per_region_tests(a)
    out = {}
    for m in CB.TESTS:
        t = reg.tests[m]
        idx = t["region_index"]
        called = np.isin(t["call"], ("Higher", "Lower", "Changed"))
        tr = changed[idx]
        nc = int(called.sum())
        pv = np.asarray(t["p"], float)[~tr]
        out[m] = dict(fdp=float((called & ~tr).sum() / nc) if nc else 0.0,
                      p01=float(np.mean(pv[np.isfinite(pv)] < 0.01)) if np.isfinite(pv).any() else float("nan"),
                      power=float((called & tr).sum() / tr.sum()) if tr.any() else float("nan"),
                      untested=float(np.mean(t["call"] == "not_tested")))
    return factor, str(level), n, out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--templates", required=True)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *x: subprocess.run(["git", *x], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"),
                protocol_sha256=sha(ROOT / "docs/PROFILESHAPE_STRESS_PROTOCOL.md"), script_sha256=sha(__file__),
                seeds=f"{SEED}+", reps=a.reps, started=datetime.datetime.now().isoformat(timespec="seconds"))
    jobs = [(f, lv, n, r, a.templates) for f, lv in CELLS for n in (2, 3) for r in range(a.reps)]
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(run, jobs, chunksize=4))
    agg: dict = {}
    for f, lv, n, o in res:
        agg.setdefault((f, lv, n), []).append(o)
    prov["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    cells = {}
    L = ["# Stress tests of the pooled null (nominal FDR 0.05)", "",
         f"Protocol `docs/PROFILESHAPE_STRESS_PROTOCOL.md` (sha256 `{prov['protocol_sha256']}`); commit `{prov['commit'][:10]}`, "
         f"uncommitted src files {prov['uncommitted_src']} (diff `{prov['src_diff_sha256']}`); regional.py "
         f"`{prov['regional_py_sha256']}`; script `{prov['script_sha256']}`; seeds {prov['seeds']}; {a.reps} runs per cell; "
         f"{prov['started']} → {prov['finished']}.", "",
         "FDP mean [95 % Monte Carlo CI] / power / p01 (unchanged regions with p < 0.01). Verdict as in the calibration "
         "study: EXCEEDS = CI above 0.05; WATCH = estimate above 0.075; ok = no clear exceedance.", "",
         "| factor | level | design | " + " | ".join(f"{t}" for t in CB.TESTS) + " | verdict |", "|---|---|---|" + "---|" * (len(CB.TESTS) + 1)]
    for (f, lv, n), rows in agg.items():
        parts, bad, cell = [], [], {}
        for t in CB.TESTS:
            fd = [r[t]["fdp"] for r in rows]
            m_, lo, hi = CB.ci_mean(fd)
            v = CB.verdict(m_, lo, hi) if float(lv if f in ("proportion_high",) else 0) < 1.0 else "n/a"
            pw = float(np.nanmean([r[t]["power"] for r in rows]))
            p01 = float(np.nanmean([r[t]["p01"] for r in rows])) if any(np.isfinite(r[t]["p01"]) for r in rows) else float("nan")
            cell[t] = dict(fdp=m_, low=lo, high=hi, verdict=v, power=pw, p01=p01)
            parts.append(f"{m_:.3f} [{max(lo, 0):.3f}, {hi:.3f}] / {pw:.2f} / {p01:.4f}" if v != "n/a" else f"– / {pw:.2f} / –")
            if v not in ("ok", "n/a"):
                bad.append(f"{t}:{v}")
        cells[f"{f}|{lv}|{n}v{n}"] = cell
        L.append(f"| {f} | {lv} | {n}v{n} | " + " | ".join(parts) + f" | {'; '.join(bad) or 'ok'} |")
    (out / "stress_null.json").write_text(json.dumps(dict(provenance=prov, cells=cells), indent=1))
    (out / "STRESS_NULL.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
