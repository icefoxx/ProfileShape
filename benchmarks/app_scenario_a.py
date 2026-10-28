"""Scenario A of docs/PROFILESHAPE_APPLICATIONS_PLAN.md: directional redistribution at the pause (PRO-seq NELFB).

    PYTHONPATH=benchmarks python benchmarks/app_scenario_a.py PROSEQ_DIR OLD_CACHE.npz OUT_DIR

Development dataset (GSE233930); an application demonstration, not independent confirmation. The named cache is
rebuilt from the bigWigs as the original cache was (RefSeq Select genes ≥ 2 kb, plus strand then minus strand)
and checked against the original before use.
"""

from __future__ import annotations

import gzip
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
from positive_cases import PROSEQ_CFG
from profileshape import metrics as M
from profileshape.analysis import analyse
from profileshape.data import Config, Profiles, bigwig_count_unit, extract
from profileshape.regional import bh, per_region_tests, support_mask

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "docs" / "PROFILESHAPE_APPLICATIONS_PLAN.md"
FILES = {"DMSO_rep1": "GSM7472402_NELFB-DMSO-{}-dedup-Rep1.bw", "DMSO_rep2": "GSM7472403_NELFB-DMSO-{}-dedup-Rep2.bw",
         "VHL2h_rep1": "GSM7472407_NELFB-VHL-2h-{}-dedup-Rep1.bw", "VHL2h_rep2": "GSM7472408_NELFB-VHL-2h-{}-dedup-Rep2.bw",
         "VHL48h_rep1": "GSM7472398_NELFB-VHL-48h-{}-dedup-Rep1.bw", "VHL48h_rep2": "GSM7472399_NELFB-VHL-48h-{}-dedup-Rep2.bw"}


def build(directory: Path, old: Path) -> dict:
    rows = [line.split("\t") for line in gzip.open(directory / "hg38_ncbiRefSeqSelect.txt.gz", "rt")]
    main = {f"chr{i}" for i in list(range(1, 23)) + ["X"]}
    regions, seen = {"+": [], "-": []}, set()
    for r in rows:
        if not r[1].startswith("NM_") or r[2] not in main or int(r[5]) - int(r[4]) < 2000 or r[12] in seen:
            continue
        seen.add(r[12])
        tss = int(r[4]) if r[3] == "+" else int(r[5]) - 1
        regions[r[3]].append(dict(chrom=r[2], start=tss, end=tss + 1, name=r[12], strand=r[3]))
    blocks = []
    for strand, kind in (("+", "FW"), ("-", "RV")):
        smp = [dict(sample_id=k, condition=k.split("_")[0], signal_file=str(directory / f.format(kind))) for k, f in FILES.items()]
        blocks.append(extract(Config(**PROSEQ_CFG), smp, regions[strand]).values)
    values = np.concatenate(blocks, axis=1)
    regs = regions["+"] + regions["-"]
    ref = np.load(old)
    same = values.shape == ref["values"].shape and bool(np.allclose(np.nan_to_num(values), np.nan_to_num(ref["values"])))
    quanta = {k: [bigwig_count_unit(str(directory / f.format(s)))[0] for s in ("FW", "RV")] for k, f in FILES.items()}
    return dict(values=values, edges=ref["edges"], regions=regs, identical_to_original=same, quanta=quanta)


def counts_of(values, regions, quanta, edges):
    """Integer counts per bin: track area / one read's value (the strand's quantum)."""
    w = np.diff(edges)
    q = np.array([[quanta[k][0] if r["strand"] == "+" else quanta[k][1] for r in regions] for k in FILES])
    return np.rint(np.nan_to_num(values) * w / q[..., None])


def diffsplice(cnt, groups, edges, core_flank, workdir):
    """edgeR QL diffSpliceDGE on core-window counts, support bins merged to 60 bp (round-6 comparator)."""
    core = M.core_mask(edges, core_flank)
    c = cnt[..., core]
    pooled = c.sum(0).astype(float)
    nb = max(1, pooled.shape[1] // 10)
    bgc = np.median(np.concatenate([pooled[:, :nb], pooled[:, -nb:]], axis=1), axis=1, keepdims=True)
    sup = support_mask((pooled - bgc)[None], 0.2)
    fd, path = tempfile.mkstemp(suffix=".tsv", dir=workdir)
    with os.fdopen(fd, "w") as fh:
        fh.write("\t".join(f"s{i}" for i in range(c.shape[0])) + "\n")
        for r in range(c.shape[1]):
            bins = np.flatnonzero(sup[r])
            for g, s in enumerate(range(0, len(bins), 6)):
                col = c[:, r][:, bins[s:s + 6]].sum(1)
                if col.sum() > 0:
                    fh.write(f"r{r}_b{g}\t" + "\t".join(str(int(v)) for v in col) + "\n")
    subprocess.run(["Rscript", str(Path(workdir) / "dsp.R"), path, ",".join(groups)], capture_output=True, text=True, check=True)
    res = json.loads(Path(path.replace(".tsv", ".json")).read_text())
    p = np.full(c.shape[1], np.nan)
    p[[int(x[1:]) for x in res["region"]]] = np.array(res["p"], float).ravel()
    os.unlink(path)
    os.unlink(path.replace(".tsv", ".json"))
    return p


def meta(values, edges):
    """Per-sample pause peak and 0–300 bp centroid over the 2000 promoters with the highest pooled signal."""
    c = M.centers(edges)
    prom = np.nansum(values[:, :, (c > 0) & (c < 300)], axis=-1).mean(0)
    top = np.argsort(prom)[-2000:]
    out = {}
    for k, name in enumerate(FILES):
        p = np.nan_to_num(values[k, top])
        p = p / np.maximum(p.sum(-1, keepdims=True), 1e-9)
        m = p.mean(0)
        sel, w = (c > -50) & (c < 400), (c > 0) & (c < 300)
        out[name] = dict(peak_bp=float(c[sel][np.argmax(m[sel])]), centroid_0_300=float((m * c)[w].sum() / m[w].sum()))
    return out


def contrast(values, edges, regions, cnt, alt, workdir):
    use = [k for k, n in enumerate(FILES) if n.startswith(("DMSO", alt))]
    names = [list(FILES)[k] for k in use]
    smp = [dict(sample_id=n, condition=n.split("_")[0]) for n in names]
    prof = Profiles(values[use], edges, smp, regions, True)
    res = per_region_tests(analyse(prof, Config(reference="DMSO", alternative=alt, **PROSEQ_CFG), bootstrap=0))
    sh, any_ = res.tests["shift"], res.tests["shape_any"]
    idx = np.asarray(sh["region_index"])
    hi, lo, ch = sh["call"] == "Higher", sh["call"] == "Lower", sh["call"] == "Changed"
    est = np.asarray(sh["estimate"], float)
    flags = [str(f) for f in np.asarray(sh["translation_flag"])[ch]]
    pd = diffsplice(cnt[use], ["A" if n.startswith("DMSO") else "B" for n in names], edges, PROSEQ_CFG["core_flank"], workdir)
    dsp = bh(np.nan_to_num(pd, nan=1.0)) < 0.05
    any_called = np.zeros(len(regions), bool)
    any_called[np.asarray(any_["region_index"])[np.isin(any_["call"], ("Changed",))]] = True
    down = np.zeros(len(regions), bool)
    down[idx[hi]] = True
    cons = hi | lo
    q = np.nanpercentile(est[cons], [25, 50, 75]) if cons.any() else [np.nan] * 3
    reps = [regions[idx[cons][np.nanargmin(np.abs(est[cons] - v))]]["name"] for v in q] if cons.any() else []
    return dict(tested=int(len(idx)), shape_any=int(any_called.sum()), shift_called=int((hi | lo | ch).sum()),
                downstream=int(hi.sum()), upstream=int(lo.sum()), downstream_share=float(hi.sum() / max(cons.sum(), 1)),
                changed_no_direction=int(ch.sum()), flags={f: flags.count(f) for f in sorted(set(flags))},
                consistent_share_of_shift=float(cons.sum() / max((hi | lo | ch).sum(), 1)),
                estimate_quantiles_bp=[float(np.nanpercentile(est[cons], x)) for x in (5, 25, 50, 75, 95)] if cons.any() else [],
                representatives=reps, diffsplice_called=int(dsp.sum()),
                diffsplice_of_shape_any=float(np.mean(dsp[any_called])) if any_called.any() else float("nan"),
                diffsplice_of_downstream=float(np.mean(dsp[down])) if down.any() else float("nan"),
                shape_any_of_diffsplice=float(np.mean(any_called[dsp])) if dsp.any() else float("nan"))


def figure(values, edges, regions, names, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    c = M.centers(edges)
    by = {r["name"]: i for i, r in enumerate(regions)}
    fig, axes = plt.subplots(1, len(names), figsize=(4 * len(names), 3), sharey=False)
    for ax, g in zip(np.atleast_1d(axes), names, strict=True):
        for k, n in enumerate(FILES):
            if n.startswith("VHL48h"):
                continue
            p = np.nan_to_num(values[k, by[g]])
            ax.plot(c, p / max(p.sum(), 1e-9), color="#4a6fa5" if n.startswith("DMSO") else "#c0504d",
                    ls="-" if n.endswith("1") else "--", lw=1, label=n)
        ax.set_xlim(-300, 800)
        ax.set_title(g)
        ax.set_xlabel("bp from TSS")
    np.atleast_1d(axes)[0].set_ylabel("share of window signal")
    np.atleast_1d(axes)[0].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def grade(r2, r48, m) -> str:
    d_ok = r2["downstream_share"] >= 0.90 and r48["downstream_share"] >= 0.90
    dm = [m[k]["centroid_0_300"] for k in m if k.startswith("DMSO")]
    m_ok = all(m[k]["centroid_0_300"] > max(dm) for k in m if k.startswith("VHL"))
    return "consistent" if d_ok and m_ok else ("partly consistent" if d_ok or m_ok else "not supported")


def main(argv):
    directory, old, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *x: subprocess.run(["git", *x], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted=len(git("status", "--porcelain", "--", "src", "benchmarks").splitlines()),
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"), plan_sha256=sha(PLAN), script_sha256=sha(__file__),
                started=time.strftime("%Y-%m-%dT%H:%M:%S"))
    b = build(directory, old)
    print("named cache identical to the original:", b["identical_to_original"], flush=True)
    cnt = counts_of(b["values"], b["regions"], b["quanta"], b["edges"])
    work = tempfile.mkdtemp(prefix="appA_", dir=out)
    (Path(work) / "dsp.R").write_text(IV.R_SCRIPT)
    R = dict(provenance=prov, identical_to_original=b["identical_to_original"], quanta=b["quanta"])
    R["VHL2h"] = contrast(b["values"], b["edges"], b["regions"], cnt, "VHL2h", work)
    R["VHL48h"] = contrast(b["values"], b["edges"], b["regions"], cnt, "VHL48h", work)
    R["meta"] = meta(b["values"], b["edges"])
    R["grade"] = grade(R["VHL2h"], R["VHL48h"], R["meta"])
    figure(b["values"], b["edges"], b["regions"], R["VHL2h"]["representatives"], out / "representatives_2h.png")
    prov["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "scenario_a.json").write_text(json.dumps(R, indent=1, default=float))
    L = ["# Scenario A · directional redistribution at the pause (PRO-seq NELFB, development data)", "",
         f"Plan sha256 `{prov['plan_sha256']}`; commit `{prov['commit'][:10]}` (uncommitted {prov['uncommitted']}); regional.py "
         f"`{prov['regional_py_sha256']}`; named cache identical to the original: {b['identical_to_original']}; {prov['started']} → {prov['finished']}.", "",
         "| contrast | tested | shape_any | shift called | consistent down / up | downstream share | Changed (no direction) | consistent shift bp 5/25/50/75/95 % | diffSplice calls | diffSplice among shape_any | diffSplice among downstream | shape_any among diffSplice |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k in ("VHL2h", "VHL48h"):
        r = R[k]
        L.append(f"| DMSO vs {k} | {r['tested']} | {r['shape_any']} | {r['shift_called']} | {r['downstream']} / {r['upstream']} | {r['downstream_share']:.3f} | "
                 f"{r['changed_no_direction']} {r['flags']} | {' / '.join(f'{v:+.0f}' for v in r['estimate_quantiles_bp'])} | {r['diffsplice_called']} | "
                 f"{r['diffsplice_of_shape_any']:.2f} | {r['diffsplice_of_downstream']:.2f} | {r['shape_any_of_diffsplice']:.2f} |")
    L += ["", "Per-replicate meta-profile (2000 promoters with the highest pooled signal): " +
          "; ".join(f"{k} peak {v['peak_bp']:.0f} bp, centroid {v['centroid_0_300']:.1f} bp" for k, v in R["meta"].items()), "",
          f"Representative genes (25 / 50 / 75 % of the consistent 2 h shift): {', '.join(R['VHL2h']['representatives'])} (`representatives_2h.png`).", "",
          f"**Grade: {R['grade']}** (criteria in the plan). Development data, seen before the plan was written: the grade was "
          "largely foreseeable; this is a demonstration, not a test."]
    (out / "SCENARIO_A.md").write_text("\n".join(L) + "\n")
    (Path(work) / "dsp.R").unlink()
    os.rmdir(work)
    print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1:])
