"""Independent positive case GSE218269 (CDK7 inhibition, Pol II MNase-ChIP-seq; TT-seq readout).

Protocol: docs/PROFILESHAPE_POSITIVE_CASE_PROTOCOL.md (stage 2). Implements it as written:

    python benchmarks/positive_case.py build INPUTS_DIR CACHE.npz
    python benchmarks/positive_case.py build-bam BAM_DIR BIGWIG_CACHE.npz BAM_CACHE.npz   # amendment B2
    python benchmarks/positive_case.py run   CACHE.npz OUT_DIR
"""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from profileshape import metrics as M
from profileshape.analysis import analyse
from profileshape.consistency import group_consistency
from profileshape.data import ASSAY_PRESETS, Config, Profiles, extract
from profileshape.regional import per_region_tests

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "PROFILESHAPE_POSITIVE_CASE_PROTOCOL.md"
MAIN = [f"chr{i}" for i in list(range(1, 23)) + ["X"]]
CFG = dict(ASSAY_PRESETS["proseq"], orientation="oriented")  # flank 1000, bin 10, background 300, core 500
POL2 = {f"{c}_{r}": f"GSM67381{97 + k}_HEK293CDK7as_MNaseChIPseq_PolII_RPB1_{c2}_30min_Rep{r}.fragment.smaller70bp.norm.coverage.bw"
        for k, (c, c2, r) in enumerate([("DMSO", "DMSO", 1), ("DMSO", "DMSO", 2), ("DMSO", "DMSO", 3)])}
POL2.update({f"CDK7i_{r}": f"GSM67382{r - 1:02d}_HEK293CDK7as_MNaseChIPseq_PolII_RPB1_1-NM-PP1_30min_Rep{r}.fragment.smaller70bp.norm.coverage.bw"
             for r in (1, 2, 3)})
READ_WINDOW = (1000, 4000)  # TT-seq readout window relative to the TSS, sense strand


def _tt_files(folder: Path) -> dict:
    out = {}
    for p in sorted(folder.glob("*TTseq_LabeledRNA_*plus.bw")):
        cond = "DMSO" if "_DMSO_" in p.name else "CDK7i"
        minutes = p.name.split("min_")[0].rsplit("_", 1)[-1]  # …_DMSO_15min_Rep1… → 15
        rep = p.name.split("_Rep")[1][0]
        out[(cond, int(minutes), int(rep))] = (p, Path(str(p).replace("plus.bw", "minus.bw")))
    return out


def _mean(h, chrom, a, b):
    size = h.chroms().get(chrom)
    a, b = max(int(a), 0), min(int(b), size or 0)
    if not size or b <= a:
        return np.nan
    v = h.stats(chrom, a, b, type="mean", exact=True)[0]
    return float(v) if v is not None else 0.0


def build(inputs: Path, cache: Path) -> None:
    import pyBigWig
    folder = inputs / "cdk7_gse218269"
    rows = [line.split("\t") for line in gzip.open(inputs / "proseq_gse233930" / "hg38_ncbiRefSeqSelect.txt.gz", "rt")]
    genes, seen = [], set()
    for r in rows:
        if not r[1].startswith("NM_") or r[2] not in MAIN or int(r[5]) - int(r[4]) < 5000 or r[12] in seen:
            continue
        seen.add(r[12])
        tss = int(r[4]) if r[3] == "+" else int(r[5]) - 1
        genes.append(dict(chrom=r[2], start=tss, end=tss + 1, name=r[12], strand=r[3], length=int(r[5]) - int(r[4])))
    print(f"genes ≥ 5 kb: {len(genes)}", flush=True)
    # condition-blind selection 1: pooled Pol II (all six samples, each scaled to its genome-wide mean) in TSS ± 300
    hs = {k: pyBigWig.open(str(folder / f)) for k, f in POL2.items()}
    pooled = np.zeros(len(genes))
    for h in hs.values():
        s = np.array([_mean(h, g["chrom"], g["start"] - 300, g["start"] + 301) for g in genes])
        pooled += np.nan_to_num(s) / max(np.nanmean(s), 1e-12)
    keep = pooled >= np.median(pooled)
    # TT-seq readout (sense strand), all time points; selection 2 uses the pooled 30-min samples
    tt = _tt_files(folder)
    lo, hi = READ_WINDOW
    readout = {}
    for key, (plus, minus) in tt.items():
        hp, hm = pyBigWig.open(str(plus)), pyBigWig.open(str(minus))
        vals = []
        for g in genes:
            if g["strand"] == "+":
                vals.append(_mean(hp, g["chrom"], g["start"] + lo, g["start"] + hi))
            else:
                vals.append(_mean(hm, g["chrom"], g["start"] - hi + 1, g["start"] - lo + 1))
        readout["{}_{}min_{}".format(*key)] = np.abs(np.array(vals))  # minus-strand tracks may be stored negative
    pooled_tt = sum(v for k, v in readout.items() if "_30min_" in k)
    keep &= np.nan_to_num(pooled_tt) > 0
    sel = [g for g, k in zip(genes, keep, strict=True) if k]
    print(f"selected regions: {len(sel)}", flush=True)
    smp = [dict(sample_id=k, condition=k.split("_")[0], signal_file=str(folder / f)) for k, f in POL2.items()]
    prof = extract(Config(**CFG), smp, sel)
    np.savez_compressed(cache, values=prof.values, edges=prof.edges, samples=json.dumps(prof.samples),
                        regions=json.dumps(prof.regions), fragments=prof.fragments if prof.fragments is not None else np.array([]),
                        provenance=json.dumps(prof.provenance), readout=json.dumps({k: v[keep].tolist() for k, v in readout.items()}),
                        length=np.array([g["length"] for g in sel]))


BAM_OPTIONS = dict(bam_mode="midpoint", bam_min_fragment=20, bam_max_fragment=120, bam_min_mapq=10,
                   bam_proper_pairs=True, bam_remove_duplicates=False)  # amendment B2


def build_bam(bam_dir: Path, bigwig_cache: Path, cache: Path) -> None:
    """Amendment B2: the same regions and readout as the bigWig cache, Pol II from re-aligned reads."""
    d = np.load(bigwig_cache, allow_pickle=True)
    regions = json.loads(str(d["regions"]))
    smp = [dict(sample_id=k, condition=k.split("_")[0], signal_file=str(bam_dir / f"{k}.bam")) for k in POL2]
    prof = extract(Config(**CFG, **BAM_OPTIONS), smp, regions)
    core = M.core_mask(prof.edges, CFG["core_flank"])
    core_counts = np.asarray(prof.counts)[..., core].sum(-1)  # exact fragments with a midpoint in TSS ± 500 bp
    np.savez_compressed(cache, values=prof.values, edges=prof.edges, samples=json.dumps(prof.samples),
                        regions=json.dumps(prof.regions), fragments=prof.fragments, fragments_core=core_counts,
                        provenance=json.dumps(prof.provenance), readout=d["readout"], length=d["length"],
                        qc=json.dumps(prof.qc, default=str))


def _cv_r2(X, y, folds):
    pred = np.empty_like(y)
    for f in np.unique(folds):
        tr, te = folds != f, folds == f
        A = np.c_[np.ones(tr.sum()), X[tr]]
        beta, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
        pred[te] = np.c_[np.ones(te.sum()), X[te]] @ beta
    return 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)


def bootstrap_folds(pick: np.ndarray, chrom: np.ndarray, fold_of: dict) -> tuple[np.ndarray, np.ndarray]:
    """Rows and folds of one chromosome bootstrap sample: every copy of a drawn chromosome keeps that
    chromosome's original fold, so the same genes never sit in a training and a test fold at once
    (amendment B3; folds by draw position had split repeated chromosomes across folds)."""
    idx = np.concatenate([np.flatnonzero(chrom == c) for c in pick])
    return idx, np.array([fold_of[c] for c in chrom[idx]])


def added_value(base: np.ndarray, shape: np.ndarray, y: np.ndarray, chrom: np.ndarray, boot: int = 500, seed: int = 0) -> dict:
    """A2: cross-validated R² gain of base + shape over base, folds by chromosome, bootstrap over chromosomes."""
    ok = np.isfinite(y) & np.all(np.isfinite(base), 1) & np.all(np.isfinite(shape), 1)
    base, shape, y, chrom = base[ok], shape[ok], y[ok], chrom[ok]
    chroms = np.unique(chrom)
    fold_of = {c: i % 10 for i, c in enumerate(chroms)}
    folds = np.array([fold_of[c] for c in chrom])
    full = np.c_[base, shape]
    gain = _cv_r2(full, y, folds) - _cv_r2(base, y, folds)
    rng = np.random.default_rng(seed)
    gains = []
    for _ in range(boot):
        pick = rng.choice(chroms, len(chroms), replace=True)
        idx, f = bootstrap_folds(pick, chrom, fold_of)
        if len(np.unique(f)) < 2:  # all draws in one fold: no held-out data
            continue
        gains.append(_cv_r2(full[idx], y[idx], f) - _cv_r2(base[idx], y[idx], f))
    A = np.c_[np.ones(len(y)), (full - full.mean(0)) / full.std(0)]
    beta = np.linalg.lstsq(A, y, rcond=None)[0]
    return dict(regions=int(ok.sum()), r2_base=float(_cv_r2(base, y, folds)), r2_full=float(_cv_r2(full, y, folds)),
                gain=float(gain), gain_low=float(np.quantile(gains, 0.025)), gain_high=float(np.quantile(gains, 0.975)),
                std_coef_signed_shift=float(beta[1 + base.shape[1]]), std_coef_kernel=float(beta[2 + base.shape[1]]))


def run(cache: Path, out: Path, exploratory: bool = False) -> None:
    out.mkdir(parents=True, exist_ok=True)
    git = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()[:16]  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"), protocol_sha256=sha(PROTOCOL),
                script_sha256=sha(Path(__file__)), started=time.strftime("%Y-%m-%dT%H:%M:%S"),
                exploratory=exploratory)
    d = np.load(cache, allow_pickle=True)
    frags = d["fragments"] if d["fragments"].size else None
    prof = Profiles(d["values"], d["edges"], json.loads(str(d["samples"])), json.loads(str(d["regions"])), True,
                    provenance=json.loads(str(d["provenance"])), fragments=frags)
    readout = {k: np.array(v, float) for k, v in json.loads(str(d["readout"])).items()}
    cfg = Config(reference="DMSO", alternative="CDK7i", **CFG)
    R = dict(provenance=prov)
    # ---- A0: depth (S3) and consistency
    core = M.core_mask(prof.edges, cfg.core_flank)
    exact = "fragments_core" in d.files
    if exact or frags is not None:
        med = np.nanmedian(d["fragments_core"] if exact else frags, axis=1)
        R["depth"] = dict(median_fragments_per_sample=dict(zip([s["sample_id"] for s in prof.samples], med.round(1).tolist(), strict=True)),
                          median_over_samples=float(np.median(med)), S3_met=bool(np.median(med) >= 100),
                          source="BAM, exact fragments with a midpoint in the core window (amendment B2)" if exact
                          else "bigWig count quantum (estimate, whole window)")
    else:
        R["depth"] = dict(S3_met=None, note="count quantum not recoverable from the tracks: depth not estimable")
    cons = group_consistency(prof, cfg)
    R["consistency"] = dict(flags=cons.flags, samples=cons.sample_rows, groups=cons.group_rows)
    if R["depth"].get("S3_met") is False and not exploratory:
        R["verdict"] = "not evaluable: S3 (depth) not met"
        _write(out, R)
        return
    # ---- A1
    a = analyse(prof, cfg, bootstrap=200)
    res = per_region_tests(a)
    any_ = res.tests["shape_any"]
    sh = res.tests["shift"]
    called = np.isin(any_["call"], ("Changed",))
    tested = any_["call"] != "not_tested"
    sig_shift = sh["fdr"] < 0.05
    upstream = sh["direction"][sig_shift] < 0
    cen = next(e for e in a.effects if e["metric"] == "centroid")
    R["A1"] = dict(tested=int(tested.sum()), shape_any_share=float(called.sum() / max(tested.sum(), 1)),
                   P1=bool(called.sum() / max(tested.sum(), 1) >= 0.20),
                   shift_significant=int(sig_shift.sum()), upstream_share=float(upstream.mean()) if upstream.size else float("nan"),
                   P2=bool(upstream.size and upstream.mean() >= 0.80),
                   centroid_bp=cen["estimate"], centroid_low=cen["welch_low"], centroid_high=cen["welch_high"],
                   P3=bool(cen["welch_high"] < 0),
                   translation_consistent=int(np.sum(sh["translation_consistent"] & sig_shift)),
                   median_consistent_shift_bp=float(np.nanmedian(sh["estimate"][sh["translation_consistent"] & sig_shift]))
                   if np.any(sh["translation_consistent"] & sig_shift) else float("nan"))
    # ---- A2: per-region covariates on the shape-tested regions (same order as the shape tests)
    idx = any_["region_index"]
    ref, alt = a.ref_idx, a.alt_idx
    mass = a.core_mass[:, idx]
    base_pol = np.log(np.maximum(mass[: len(ref)].mean(0), 1e-12))
    amount = np.log2(np.maximum(mass[len(ref):].mean(0), 1e-12)) - np.log2(np.maximum(mass[: len(ref)].mean(0), 1e-12))
    raw = (np.nan_to_num(prof.values[:, idx]) * M.widths(prof.edges)).sum(-1)
    depth = np.log(np.maximum(raw.mean(0), 1e-12))
    pooled = np.clip(np.nan_to_num(a.core_signal[:, idx]).mean(0), 0, None)
    width = M.mean_abs_deviation(M.normalise(pooled * M.widths(a.core_edges)), a.core_edges)
    length = np.log(d["length"][idx])
    chrom = np.array([prof.regions[i]["chrom"] for i in idx])
    signed_shift = sh["direction"] * -np.log10(np.maximum(sh["p"], 1e-300))
    kern = -np.log10(np.maximum(res.tests["shape_kernel"]["p"], 1e-300))
    R["A2"] = {}
    for minutes in (30, 15, 60):
        dm = np.mean([readout[k][idx] for k in readout if k.startswith(f"DMSO_{minutes}min")], 0)
        cd = np.mean([readout[k][idx] for k in readout if k.startswith(f"CDK7i_{minutes}min")], 0)
        pc = 0.01 * np.nanmedian(np.r_[dm, cd][np.r_[dm, cd] > 0])
        y = np.log2(cd + pc) - np.log2(dm + pc)
        base = np.c_[base_pol, amount, depth, np.log(dm + pc), length, width]
        R["A2"][f"{minutes}min"] = added_value(base, np.c_[signed_shift, kern], y, chrom, seed=minutes)
    R["A2"]["criterion_met_30min"] = bool(R["A2"]["30min"]["gain_low"] > 0)
    R["notes"] = res.notes
    _write(out, R)


def _write(out: Path, R: dict) -> None:
    R["provenance"]["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out / "positive_case.json").write_text(json.dumps(R, indent=2, default=lambda v: v.item() if isinstance(v, np.generic) else str(v)))
    L = ["# Independent positive case · GSE218269 (CDK7 inhibition, Pol II MNase-ChIP-seq)", ""]
    if R["provenance"].get("exploratory"):
        L += ["**EXPLORATORY (amendment B4): S3 (depth) was not met; the confirmatory case is \"not evaluable\". "
              "This run is not a confirmation and is not counted as one.**", ""]
    L += [
         "Protocol `docs/PROFILESHAPE_POSITIVE_CASE_PROTOCOL.md` (stage 2) sha256 `{protocol_sha256}`; commit `{commit:.10}`, "
         "uncommitted src files {uncommitted_src} (diff `{src_diff_sha256}`); regional.py `{regional_py_sha256}`; script "
         "`{script_sha256}`; {started} → {finished}.".format(**R["provenance"]), "",
         "## A0 · Depth and consistency", "", f"```\n{json.dumps(R['depth'], indent=1)}\n```", "",
         "Consistency flags: " + ("; ".join(R["consistency"]["flags"]) or "none"), ""]
    if "A1" in R:
        a1 = R["A1"]
        L += ["## A1 · Expected direction", "",
              "| prediction | result | met |", "|---|---|---|",
              f"| P1 `shape_any` share ≥ 0.20 | {a1['shape_any_share']:.3f} of {a1['tested']} | {'yes' if a1['P1'] else '**no**'} |",
              f"| P2 upstream share among `shift`-significant ≥ 0.80 | {a1['upstream_share']:.3f} of {a1['shift_significant']} | {'yes' if a1['P2'] else '**no**'} |",
              f"| P3 aggregate centroid upstream (Welch 95 % < 0) | {a1['centroid_bp']:+.1f} bp [{a1['centroid_low']:+.1f}, {a1['centroid_high']:+.1f}] | {'yes' if a1['P3'] else '**no**'} |",
              "", f"Consistent translations among `shift`-significant regions: {a1['translation_consistent']} "
              f"(median {a1['median_consistent_shift_bp']:+.1f} bp).", ""]
        L += ["## A2 · Added explanatory value for TT-seq (primary at 30 min)", "",
              "| readout | regions | CV R² amount-only | CV R² + shape | gain [95 % CI] | std. coef. signed shift | std. coef. kernel |",
              "|---|---|---|---|---|---|---|"]
        for k in ("30min", "15min", "60min"):
            r = R["A2"][k]
            L.append(f"| TT-seq {k} | {r['regions']} | {r['r2_base']:.3f} | {r['r2_full']:.3f} | {r['gain']:+.4f} "
                     f"[{r['gain_low']:+.4f}, {r['gain_high']:+.4f}] | {r['std_coef_signed_shift']:+.3f} | {r['std_coef_kernel']:+.3f} |")
        L += ["", f"Primary criterion (30 min, lower bound of gain > 0): **{'met' if R['A2']['criterion_met_30min'] else 'not met'}**."]
    else:
        L += [f"**{R.get('verdict', '')}**"]
    (out / "POSITIVE_CASE.md").write_text("\n".join(L) + "\n")
    print((out / "POSITIVE_CASE.md").read_text())


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "build":
        build(Path(sys.argv[2]), Path(sys.argv[3]))
    elif cmd == "build-bam":
        build_bam(Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
    else:
        run(Path(sys.argv[2]), Path(sys.argv[3]), exploratory="--exploratory" in sys.argv[4:])
