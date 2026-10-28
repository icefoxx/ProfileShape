"""Evaluation of consistency verdicts against user tolerance margins (protocol: docs/PROFILESHAPE_MARGIN_EVALUATION.md).

    python benchmarks/margin_eval.py run INPUTS_DIR OUT_DIR [--workers 10]
"""

from __future__ import annotations

import argparse
import hashlib
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

sys.path.insert(0, str(Path(__file__).parent))
import consistency_eval as CE  # noqa: E402

from profileshape.consistency import group_consistency, region_blocks  # noqa: E402
from profileshape.data import Config  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "PROFILESHAPE_MARGIN_EVALUATION.md"
SIM_SEED, SEMI_SEED = 9_990_000, 9_995_000
SHIFTS = (0, 5, 10, 15, 20)
STRETCHES = (2.5, 5, 7.5, 10)
KEYS = ("shift_bp_est", "shift_bp_low90", "shift_bp_high90", "shift_verdict", "stretch_pct", "stretch_pct_low90",
        "stretch_pct_high90", "stretch_verdict", "offset_pct", "offset_pct_lower", "offset_pct_upper", "offset_verdict")


def _kind(alter):
    kind, size = alter
    return f"shift{size}" if kind == "shift" else f"broaden{1 + size / 100:.4f}"


def _row(c):
    r = c.sample_rows[0]
    return {k: r.get(k) for k in KEYS}


def job_sim(args):
    from profileshape.simulate import SimSettings, simulate
    cell, n, depth, alter, run = args
    seed = SIM_SEED + 1000 * cell + run
    cfg = Config(consistency_margin_shift_bp=10, consistency_margin_stretch_pct=5, consistency_margin_offset_pct=5)
    s = SimSettings(n_ref=n, n_alt=1, n_regions=2000, depth=depth)
    out = {}
    for label, noiseless in (("obs", False), ("truth", True)):
        p, _ = simulate("null", s, seed=seed, noiseless=noiseless)
        v = p.values[:n]
        if alter[1]:
            v = CE.insert(v, p.edges, [0], _kind(alter), np.random.default_rng(seed + 7), cfg.background_width)
        c = group_consistency(CE._profiles(v, p.edges, [f"r{i + 1}" for i in range(n)], p.regions), cfg)
        out[label] = _row(c)
    return dict(arm="A", cell=cell, n=n, depth=depth, kind=alter[0], size=alter[1], run=run, obs=out["obs"],
                truth_offset_pct=out["truth"]["offset_pct"])


def job_semi(args):
    inputs, ds, gi, group, alter, run = args
    v, edges, members, regions, cfg, _, _ = CE._group_data(inputs, ds, group)
    M = max(10.0, float(np.median(np.diff(edges))))
    size = alter[1] * M if alter[0] == "shift" else alter[1]
    cfg2 = Config(**{**{k: getattr(cfg, k) for k in Config.__dataclass_fields__ if hasattr(cfg, k)},
                     "consistency_margin_shift_bp": M, "consistency_margin_stretch_pct": 5.0})
    seed = SEMI_SEED + 1000 * list(CE.DATASETS).index(ds) + 100 * gi + run
    rng = np.random.default_rng(seed)
    (v,) = CE.shuffle_blocks(v, region_blocks(regions, np.arange(len(regions))), rng)
    if size:
        v = CE.insert(v, edges, [0], _kind((alter[0], size)), np.random.default_rng(seed + 7), cfg.background_width)
    c = group_consistency(CE._profiles(v, edges, [f"{group}_p{i + 1}" for i in range(len(members))], regions), cfg2)
    return dict(arm="B", dataset=ds, group=group, n=len(members), kind=alter[0], multiple=alter[1], size=size,
                margin_shift=M, run=run, obs=_row(c))


def job_real(args):
    inputs, ds, group = args
    v, edges, members, regions, cfg, _, _ = CE._group_data(inputs, ds, group)
    c = group_consistency(CE._profiles(v, edges, members, regions), cfg)
    return dict(arm="C", dataset=ds, group=group, n=len(members),
                samples=[{"sample": r["sample"], **{k: r.get(k) for k in KEYS}} for r in c.sample_rows])


def jobs_for(inputs):
    sim, cell = [], 0
    for n in (2, 3, 4):
        for depth in (12, 40, 160):
            for alter in [("shift", d) for d in SHIFTS] + [("stretch", f) for f in STRETCHES]:
                sim += [(cell, n, depth, alter, r) for r in range(100)]
                cell += 1
    semi, real = [], []
    for ds, (name, _, _) in CE.DATASETS.items():
        path = Path(inputs) / "caches" / name
        if not path.exists():
            continue
        groups = json.loads(str(np.load(path, allow_pickle=True)["groups"]))
        for gi, g in enumerate(groups):
            real.append((inputs, ds, g))
            for alter in [("shift", m) for m in (0, 0.5, 1, 1.5, 2)] + [("stretch", f) for f in STRETCHES]:
                semi += [(inputs, ds, gi, g, alter, r) for r in range(20)]
    return sim, semi, real


def _wrong_within(rows, kind, truth, margin):
    key = {"shift": "shift_verdict", "stretch": "stretch_verdict", "offset": "offset_verdict"}[kind]
    return float(np.mean([r["obs"][key] == "within" for r in rows])) if rows else float("nan")


def report(prov, R) -> str:
    sim, semi, real = R["sim"], R["semireal"], R["real"]
    L = ["# Consistency against tolerance margins — evaluation results", "",
         f"Protocol sha256 `{prov['protocol_sha256']}`; commit `{prov['commit'][:10]}`, uncommitted src files "
         f"{prov['uncommitted_src']} (diff `{prov['src_diff_sha256']}`); consistency.py `{prov['consistency_py_sha256']}`; "
         f"script `{prov['script_sha256']}`; seeds simulation {SIM_SEED}+, semi-real {SEMI_SEED}+; "
         f"{prov['started']} → {prov['finished']}.", ""]
    V = {}
    # ---- A
    L += ["## A · Simulation (100 runs per cell; margins: shift 10 bp, stretch 5 %, offset 5 %)", "",
          "| n | depth | alteration | estimate (mean) | 90 % CI covers truth | within | inconclusive | beyond | offset % truth / est / upper | offset within |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    m1, cover_s, cover_w, cover_o = [], [], [], []
    for key in sorted({(r["n"], r["depth"], r["kind"], r["size"]) for r in sim}):
        rows = [r for r in sim if (r["n"], r["depth"], r["kind"], r["size"]) == key]
        n, depth, kind, size = key
        q = "shift_bp" if kind == "shift" else "stretch_pct"
        est = np.mean([r["obs"][q + ("_est" if kind == "shift" else "")] for r in rows])
        lo, hi = (q + "_low90", q + "_high90")
        cov = np.mean([r["obs"][lo] <= size <= r["obs"][hi] for r in rows])
        (cover_s if kind == "shift" else cover_w).append(cov)
        vkey = "shift_verdict" if kind == "shift" else "stretch_verdict"
        frac = {v: np.mean([r["obs"][vkey] == v for r in rows]) for v in ("within", "inconclusive", "beyond")}
        margin = 10 if kind == "shift" else 5
        if size >= margin:
            m1.append(frac["within"] <= (0.10 if size == margin else 0.05))
        to = np.mean([r["truth_offset_pct"] for r in rows])
        oe = np.mean([r["obs"]["offset_pct"] for r in rows])
        ou = np.mean([r["obs"]["offset_pct_upper"] for r in rows])
        cover_o += [r["obs"]["offset_pct_upper"] >= r["truth_offset_pct"] for r in rows]
        ow = np.mean([r["obs"]["offset_verdict"] == "within" for r in rows])
        if to >= 5:
            m1.append(ow <= (0.10 if to < 7.5 else 0.05))
        L.append(f"| {n} | {depth} | {kind} {size:g} | {est:+.2f} | {cov:.2f} | {frac['within']:.2f} | {frac['inconclusive']:.2f} | "
                 f"{frac['beyond']:.2f} | {to:.2f} / {oe:.2f} / {ou:.2f} | {ow:.2f} |")
    V["M1"] = all(m1)
    V["M3"] = (min(cover_s + cover_w) >= 0.85 and np.mean(cover_s) >= 0.88 and np.mean(cover_w) >= 0.88
               and np.mean(cover_o) >= 0.93)
    L += ["", f"Pooled coverage: shift {np.mean(cover_s):.3f}, stretch {np.mean(cover_w):.3f}; offset upper bound ≥ truth "
          f"{np.mean(cover_o):.3f}."]
    # ---- B
    L += ["", "## B · Semi-real (20 shuffles; margins: shift M = max(10 bp, bin), stretch 5 %)", "",
          "| dataset | group | n | alteration | estimate (mean) | 90 % CI covers truth | within | inconclusive | beyond |",
          "|---|---|---|---|---|---|---|---|---|"]
    m2, cov_b = [], []
    for key in sorted({(r["dataset"], r["group"], r["kind"], r["multiple"]) for r in semi}):
        rows = [r for r in semi if (r["dataset"], r["group"], r["kind"], r["multiple"]) == key]
        ds, g, kind, mult = key
        size, M = rows[0]["size"], rows[0]["margin_shift"]
        q = "shift_bp" if kind == "shift" else "stretch_pct"
        est = np.mean([r["obs"][q + ("_est" if kind == "shift" else "")] for r in rows])
        cov = [r["obs"][q + "_low90"] <= size <= r["obs"][q + "_high90"] for r in rows if r["obs"][q + "_low90"] is not None
               and np.isfinite(r["obs"][q + "_low90"])]
        cov_b += cov
        vkey = "shift_verdict" if kind == "shift" else "stretch_verdict"
        nw = sum(r["obs"][vkey] == "within" for r in rows)
        ni = sum(r["obs"][vkey] == "inconclusive" for r in rows)
        nb = sum(r["obs"][vkey] == "beyond" for r in rows)
        margin = M if kind == "shift" else 5
        if size >= margin and (nw + ni + nb):
            m2.append(nw <= (3 if np.isclose(size, margin) else 1))
        L.append(f"| {ds} | {g} | {rows[0]['n']} | {kind} {size:g}{' bp' if kind == 'shift' else ' %'} | {est:+.2f} | "
                 f"{np.mean(cov) if cov else float('nan'):.2f} | {nw} | {ni} | {nb} |")
    V["M2"] = all(m2)
    V["M4"] = bool(cov_b) and np.mean(cov_b) >= 0.85
    L += ["", f"Pooled semi-real coverage {np.mean(cov_b):.3f} (m = {len(cov_b)})."]
    # ---- C
    L += ["", "## C · Real replicates as they are (reference ranges, no truth)", "",
          "| dataset | group | n | offset % (upper bound) range | shift bp range | stretch % range |", "|---|---|---|---|---|---|"]
    for r in real:
        s = r["samples"]
        rng = lambda k: f"{min(x[k] for x in s):.2f} – {max(x[k] for x in s):.2f}" if all(  # noqa: E731
            x[k] is not None and np.isfinite(x[k]) for x in s) else "—"
        L.append(f"| {r['dataset']} | {r['group']} | {r['n']} | {rng('offset_pct')} ({rng('offset_pct_upper')}) | "
                 f"{rng('shift_bp_est')} | {rng('stretch_pct')} |")
    L += ["", "## Criteria", "", "| id | result |", "|---|---|"]
    L += [f"| {k} | {'pass' if V[k] else '**FAIL**'} |" for k in ("M1", "M2", "M3", "M4")]
    return "\n".join(L) + "\n"


def run(inputs, out, workers):
    import gzip
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()[:16]  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                consistency_py_sha256=sha(ROOT / "src/profileshape/consistency.py"), protocol_sha256=sha(PROTOCOL),
                script_sha256=sha(Path(__file__)), started=time.strftime("%Y-%m-%dT%H:%M:%S"))
    sim, semi, real = jobs_for(inputs)
    print(f"jobs: sim {len(sim)}, semi-real {len(semi)}, real {len(real)}", flush=True)
    R = {}
    with ProcessPoolExecutor(workers) as ex:
        for arm, fn, jobs in (("real", job_real, real), ("semireal", job_semi, semi), ("sim", job_sim, sim)):
            t0 = time.time()
            R[arm] = list(ex.map(fn, jobs, chunksize=4))
            print(f"  {arm}: {len(jobs)} jobs in {time.time() - t0:.0f} s", flush=True)
            with gzip.open(out / f"{arm}.json.gz", "wt") as f:
                json.dump(R[arm], f, default=CE._json)
    prov["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "provenance.json").write_text(json.dumps(prov, indent=2))
    (out / "MARGIN_EVALUATION.md").write_text(report(prov, R))
    print((out / "MARGIN_EVALUATION.md").read_text())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("run",))
    ap.add_argument("inputs")
    ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=10)
    a = ap.parse_args()
    run(a.inputs, a.out, a.workers)
