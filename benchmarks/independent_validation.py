"""Independent real-data validation of ProfileShape v3.1 on GSE346229 (yeast ATAC-seq, Sth1 / Snf2 anchor-away).

Protocol, criteria and amendments: docs/PROFILESHAPE_INDEPENDENT_VALIDATION.md (written before the data were seen).

    python benchmarks/independent_validation.py extract BAM_DIR ANNOTATION.txt.gz CACHE.npz
    python benchmarks/independent_validation.py run CACHE.npz OUT_DIR
    python benchmarks/independent_validation.py aggregate CACHE.npz OUT.md   # exploratory

`extract` collects, per sample, the midpoint (oriented, relative to the gene start) and length of every fragment
within ±1300 bp of each gene start, with the ProfileShape BAM-reader filters (MAPQ ≥ 10, proper pairs, duplicates
removed upstream, each fragment once from its leftmost mate, TLEN ≤ 1000). `run` builds condition-blind anchors,
then performs analyses A (real-null splits), B (fragment-level insertion) and C (biology), with ProfileShape and
edgeR diffSpliceDGE on the same counts.
"""

from __future__ import annotations

import datetime
import gzip
import hashlib
import itertools
import json
import os
import subprocess
import sys
import tempfile
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
from scipy import stats

from profileshape.analysis import analyse
from profileshape.data import Config, Profiles
from profileshape.regional import bh, per_region_tests, support_mask

EXTRACT = 1300
FLANK, BIN, CORE, BGW = 600, 10, 300, 200
NUC = (140, 200)
SUB = (0, 120)
SEED = 9_900_000
STH1 = {"0h": [f"sth1_0h_rep{i}" for i in range(1, 5)], "1h": [f"sth1_1h_rep{i}" for i in range(1, 5)]}
SNF2 = {"0h": [f"snf2_0h_rep{i}" for i in range(1, 4)], "1h": [f"snf2_1h_rep{i}" for i in range(1, 4)]}
CFG = dict(anchor="start", orientation="oriented", flank=FLANK, binsize=BIN, background_width=BGW, core_flank=CORE)
TESTS = ("shape_any", "shape_kernel", "shift")
R_SCRIPT = r"""
suppressMessages(library(edgeR))
args <- commandArgs(TRUE); f <- args[1]; groups <- strsplit(args[2], ",")[[1]]
x <- as.matrix(read.table(f, header = TRUE, row.names = 1, check.names = FALSE))
region <- sub("_b[0-9]+$", "", rownames(x))
cond <- factor(groups, levels = c("A", "B"))
design <- model.matrix(~cond)
y <- normLibSizes(DGEList(x, genes = data.frame(GeneID = region)))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
sp <- diffSpliceDGE(fit, coef = 2, geneid = "GeneID", verbose = FALSE)
writeLines(jsonlite::toJSON(list(region = as.character(sp$gene.genes$GeneID), p = as.numeric(sp$gene.p.value)), digits = NA),
           sub("\\.tsv$", ".json", f))
"""


# ----------------------------------------------------------------------------------------------- extraction
def gene_starts(annotation: str) -> list[dict]:
    rows = {}
    with gzip.open(annotation, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            name, chrom, strand = f[1], f[2], f[3]
            cds_s, cds_e = int(f[6]), int(f[7])
            if chrom in ("chrM", "chrmt") or cds_e <= cds_s:
                continue
            pos = cds_s if strand == "+" else cds_e - 1
            rows.setdefault(name.split(".")[0], dict(name=name, chrom=chrom, strand=strand, pos=pos))
    genes = sorted(rows.values(), key=lambda r: (r["chrom"], r["pos"]))
    keep = []
    for g in genes:  # ≥ 500 bp from any other gene start on the same strand
        near = [h for h in genes if h is not g and h["chrom"] == g["chrom"] and h["strand"] == g["strand"] and abs(h["pos"] - g["pos"]) < 500]
        if not near:
            keep.append(g)
    return keep


def extract_sample(job):
    bam_path, genes = job
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    reg, pos, length = [], [], []
    for k, g in enumerate(genes):
        lo, hi = max(0, g["pos"] - EXTRACT - 1000), g["pos"] + EXTRACT + 1000
        for r in bam.fetch(g["chrom"], lo, hi):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or not r.is_proper_pair:
                continue
            if r.mapping_quality < 10 or r.template_length <= 0 or r.template_length > 1000:
                continue
            mid = r.reference_start + r.template_length // 2
            rel = (mid - g["pos"]) if g["strand"] == "+" else (g["pos"] - mid)
            if -EXTRACT <= rel < EXTRACT:
                reg.append(k); pos.append(rel); length.append(r.template_length)
    lib = sum(s.mapped for s in bam.get_index_statistics()) / 2.0
    return np.array(reg, np.int32), np.array(pos, np.int32), np.array(length, np.int16), lib


def extract(bam_dir: str, annotation: str, cache: str) -> None:
    genes = gene_starts(annotation)
    names = STH1["0h"] + STH1["1h"] + SNF2["0h"] + SNF2["1h"]
    jobs = [(str(Path(bam_dir) / f"{n}.bam"), genes) for n in names]
    out = dict(genes=json.dumps(genes), samples=json.dumps(names))
    with ProcessPoolExecutor(7) as ex:
        for n, (reg, pos, ln, lib) in zip(names, ex.map(extract_sample, jobs), strict=True):
            out[f"{n}__reg"], out[f"{n}__pos"], out[f"{n}__len"], out[f"{n}__lib"] = reg, pos, ln, np.array(lib)
            print(n, len(reg), "fragments near gene starts; library", int(lib), flush=True)
    np.savez_compressed(cache, **out)


# ----------------------------------------------------------------------------------------------- anchors and binning
class Data:
    def __init__(self, cache: str):
        d = np.load(cache, allow_pickle=True)
        self.genes = json.loads(str(d["genes"]))
        self.samples = json.loads(str(d["samples"]))
        self.frag = {n: (d[f"{n}__reg"], d[f"{n}__pos"], d[f"{n}__len"]) for n in self.samples}
        self.lib = {n: float(d[f"{n}__lib"]) for n in self.samples}

    def anchors(self) -> dict:
        """Condition-blind +1 / −1 dyads from pooled nucleosomal midpoints of the eight Sth1 samples."""
        G = len(self.genes)
        grid = np.arange(-EXTRACT, EXTRACT + 1, BIN)
        dens = np.zeros((G, len(grid) - 1))
        for n in STH1["0h"] + STH1["1h"]:
            reg, pos, ln = self.frag[n]
            m = (ln >= NUC[0]) & (ln <= NUC[1])
            h = np.zeros_like(dens)
            np.add.at(h, (reg[m], np.clip(np.searchsorted(grid, pos[m], side="right") - 1, 0, len(grid) - 2)), 1.0)
            dens += h / max(m.sum(), 1)
        kern = np.exp(-0.5 * (np.arange(-6, 7) * BIN / 20.0) ** 2); kern /= kern.sum()
        dens = np.apply_along_axis(lambda r: np.convolve(r, kern, mode="same"), 1, dens)
        c = 0.5 * (grid[:-1] + grid[1:])
        p1 = np.full(G, np.nan); m1 = np.full(G, np.nan)
        for g in range(G):
            w = (c >= -200) & (c <= 100)
            if dens[g, w].max() <= 0:
                continue
            p1[g] = c[w][np.argmax(dens[g, w])]
            w2 = (c >= p1[g] - 450) & (c <= p1[g] - 200)
            if w2.any() and dens[g, w2].max() > 0:
                m1[g] = c[w2][np.argmax(dens[g, w2])]
        return {"plus1": p1, "minus1": m1, "ndr": 0.5 * (p1 + m1)}

    def counts(self, names: list[str], anchor: np.ndarray, size: tuple, modify=None) -> tuple[np.ndarray, list[int]]:
        """samples × regions × bins counts over anchor ± FLANK (oriented); `modify(name, reg, pos, len, rng)` may
        alter fragments of a sample before binning."""
        ok = np.flatnonzero(np.isfinite(anchor))
        idx = {g: j for j, g in enumerate(ok)}
        edges = np.arange(-FLANK, FLANK + BIN, BIN)
        out = np.zeros((len(names), len(ok), len(edges) - 1))
        for i, n in enumerate(names):
            reg, pos, ln = self.frag[n]
            m = (ln >= size[0]) & (ln <= size[1])
            reg, pos, ln = reg[m], pos[m].astype(float), ln[m]
            if modify is not None:
                reg, pos, ln = modify(n, reg, pos, ln)
            keep = np.isin(reg, ok)
            reg, pos = reg[keep], pos[keep]
            rel = pos - anchor[reg]
            inside = (rel >= -FLANK) & (rel < FLANK)
            j = np.array([idx[g] for g in reg[inside]], int)
            b = np.searchsorted(edges, rel[inside], side="right") - 1
            np.add.at(out[i], (j, b), 1.0)
        return out, ok.tolist()

    def profiles(self, names, groups, anchor, size, modify=None) -> tuple[Profiles, np.ndarray]:
        cnt, ok = self.counts(names, anchor, size, modify)
        lib = np.array([self.lib[n] for n in names])
        values = cnt / lib[:, None, None] * 1e6
        regions = [dict(chrom="c", start=10_000 * (k + 1), end=10_000 * (k + 1) + 1, name=self.genes[g]["name"], strand="+")
                   for k, g in enumerate(ok)]
        smp = [dict(sample_id=n, condition=c) for n, c in zip(names, groups, strict=True)]
        prof = Profiles(values, np.arange(-FLANK, FLANK + BIN, BIN, dtype=float), smp, regions, True,
                        provenance={"source": "bam"},  # exact fragment counts (the depth note reads this)
                        fragments=cnt[..., (np.arange(-FLANK, FLANK, BIN) >= -CORE) & (np.arange(-FLANK, FLANK, BIN) < CORE)].sum(-1))
        return prof, cnt


# ----------------------------------------------------------------------------------------------- tests
def diffsplice(cnt: np.ndarray, groups: list[str], workdir: str) -> np.ndarray:
    """edgeR QL diffSpliceDGE on core-window counts, support bins of the background-corrected pooled profile
    merged to 60 bp (round-6 comparator); returns a p-value per region (NaN where not returned)."""
    edges = np.arange(-FLANK, FLANK + BIN, BIN)
    core = (edges[:-1] >= -CORE) & (edges[:-1] < CORE)
    c = cnt[..., core]
    pooled = c.sum(0).astype(float)
    nb = max(1, pooled.shape[1] // 10)
    bgc = np.median(np.concatenate([pooled[:, :nb], pooled[:, -nb:]], axis=1), axis=1, keepdims=True)
    sup = support_mask((pooled - bgc)[None], 0.2)
    merge = 6
    fd, path = tempfile.mkstemp(suffix=".tsv", dir=workdir)
    with os.fdopen(fd, "w") as fh:
        fh.write("\t".join(f"s{i}" for i in range(c.shape[0])) + "\n")
        for r in range(c.shape[1]):
            bins = np.flatnonzero(sup[r])
            for g, s in enumerate(range(0, len(bins), merge)):
                col = c[:, r][:, bins[s:s + merge]].sum(1)
                if col.sum() > 0:
                    fh.write(f"r{r}_b{g}\t" + "\t".join(str(int(v)) for v in col) + "\n")
    subprocess.run(["Rscript", str(Path(workdir) / "dsp.R"), path, ",".join("A" if g == groups[0] else "B" for g in groups)],
                   capture_output=True, text=True, check=True)
    res = json.loads(Path(path.replace(".tsv", ".json")).read_text())
    p = np.full(c.shape[1], np.nan)
    p[[int(x[1:]) for x in res["region"]]] = np.array(res["p"], float).ravel()
    os.unlink(path); os.unlink(path.replace(".tsv", ".json"))
    return p


def profileshape(prof: Profiles, ref: str, alt: str) -> dict:
    a = analyse(prof, Config(reference=ref, alternative=alt, **CFG), bootstrap=0)
    reg = per_region_tests(a)
    out = {}
    for m in TESTS:
        t = reg.tests[m]
        called = np.zeros(prof.values.shape[1], bool)
        called[np.asarray(t["region_index"])[np.isin(t["call"], ("Higher", "Lower", "Changed"))]] = True
        out[m] = called
    t = reg.tests["shift"]
    idx = np.asarray(t["region_index"])
    est = np.full(prof.values.shape[1], np.nan); est[idx] = t["estimate"]
    lab = np.full(prof.values.shape[1], "", object); lab[idx] = t["call"]
    out["shift_estimate"], out["shift_call"] = est, lab
    out["notes"] = reg.notes
    return out


# ----------------------------------------------------------------------------------------------- jobs
def _data():
    return Data(os.environ["PSV_CACHE"])


def job_null(args):
    anchor_set, size_name, cond, a_ids, b_ids, workdir = args
    D = _data(); anc = D.anchors()[anchor_set]
    size = NUC if size_name == "nuc" else SUB
    names = a_ids + b_ids; groups = ["A"] * 2 + ["B"] * 2
    prof, cnt = D.profiles(names, groups, anc, size)
    ps = profileshape(prof, "A", "B")
    pd = diffsplice(cnt, groups, workdir)
    res = {m: float(ps[m].mean()) for m in TESTS}
    res["diffSplice"] = float(np.mean(bh(np.nan_to_num(pd, nan=1.0)) < 0.05))
    res["regions"] = int(prof.values.shape[1])
    return dict(kind="null", anchor=anchor_set, size=size_name, cond=cond, split=f"{'+'.join(a_ids)} vs {'+'.join(b_ids)}", **res)


def job_insert(args):
    anchor_set, size_name, cond, a_ids, b_ids, insertion, rep, workdir = args
    D = _data(); anc = D.anchors()[anchor_set]
    size = NUC if size_name == "nuc" else SUB
    names = a_ids + b_ids; groups = ["A"] * 2 + ["B"] * 2
    ok = np.flatnonzero(np.isfinite(anc))
    rng = np.random.default_rng(SEED + zlib.crc32(f"{anchor_set}|{cond}|{a_ids}|{insertion}|{rep}".encode()) % 100_000)
    chosen = rng.choice(ok, int(0.2 * len(ok)), replace=False)
    sign = dict(zip(chosen.tolist(), rng.choice([-1.0, 1.0], len(chosen)).tolist()))
    sel = np.zeros(len(D.genes), bool); sel[chosen] = True

    def modify(name, reg, pos, ln):
        if name not in b_ids:
            return reg, pos, ln
        r = np.random.default_rng(zlib.crc32(f"{name}|{insertion}|{rep}|{anchor_set}".encode()))
        hit = sel[reg]
        if insertion.startswith("shift"):
            d = float(insertion[5:])
            pos = pos.copy(); pos[hit] += d * np.array([sign[g] for g in reg[hit]])
        elif insertion == "broaden":
            pos = pos.copy(); pos[hit] += r.normal(0, 30, hit.sum())
        elif insertion == "amplitude":
            keep = ~hit | (r.uniform(size=len(reg)) < 0.5)
            reg, pos, ln = reg[keep], pos[keep], ln[keep]
        return reg, pos, ln

    prof, cnt = D.profiles(names, groups, anc, size, modify)
    truth = np.isin(np.array(ok), chosen) if insertion != "amplitude" else np.zeros(len(ok), bool)
    ps = profileshape(prof, "A", "B")
    pd = diffsplice(cnt, groups, workdir)
    calls = {m: ps[m] for m in TESTS}
    calls["diffSplice"] = bh(np.nan_to_num(pd, nan=1.0)) < 0.05
    res = {}
    for m, c in calls.items():
        nc = int(c.sum())
        res[m] = dict(called=float(c.mean()), power=float((c & truth).sum() / truth.sum()) if truth.any() else float("nan"),
                      fdp=float((c & ~truth).sum() / nc) if nc else 0.0)
    if insertion.startswith("shift"):
        d = float(insertion[5:])
        okr = np.array(ok)
        sgn = np.array([sign.get(g, 0.0) for g in okr])
        dirl = np.isin(ps["shift_call"], ("Higher", "Lower")) & truth
        est = ps["shift_estimate"]
        res["translation"] = dict(coverage=float(dirl.sum() / truth.sum()),
                                  direction_ok=float(np.mean(np.sign(est[dirl]) == sgn[dirl])) if dirl.any() else float("nan"),
                                  ratio=float(np.median(est[dirl] * sgn[dirl] / d)) if dirl.any() else float("nan"))
    return dict(kind="insert", anchor=anchor_set, size=size_name, cond=cond, insertion=insertion, rep=rep,
                split=f"{'+'.join(a_ids)} vs {'+'.join(b_ids)}", **res)


def job_bio(args):
    anchor_set, size_name, ref_ids, alt_ids, label, workdir = args
    D = _data(); anc = D.anchors()[anchor_set]
    size = NUC if size_name == "nuc" else SUB
    names = ref_ids + alt_ids; groups = ["ref"] * len(ref_ids) + ["alt"] * len(alt_ids)
    prof, cnt = D.profiles(names, groups, anc, size)
    ps = profileshape(prof, "ref", "alt")
    pd = diffsplice(cnt, ["A"] * len(ref_ids) + ["B"] * len(alt_ids), workdir)
    dsp = bh(np.nan_to_num(pd, nan=1.0)) < 0.05
    lab = ps["shift_call"]
    up, down = int((lab == "Lower").sum()), int((lab == "Higher").sum())
    est = ps["shift_estimate"]
    return dict(kind="bio", anchor=anchor_set, size=size_name, label=label, regions=int(prof.values.shape[1]),
                **{f"{m}_share": float(ps[m].mean()) for m in TESTS}, diffSplice_share=float(dsp.mean()),
                overlap_any_diffSplice=float((ps["shape_any"] & dsp).sum() / max((ps["shape_any"] | dsp).sum(), 1)),
                consistent_upstream=up, consistent_downstream=down,
                median_consistent_shift_bp=float(np.nanmedian(est[np.isin(lab, ("Higher", "Lower"))])) if (up + down) else float("nan"),
                notes=[n for n in ps["notes"] if "strata" in n or "depth" in n.lower() or "Changed" in n])


def splits(ids: list[str]) -> list[tuple[list[str], list[str]]]:
    out = []
    for a in itertools.combinations(ids, 2):
        if ids[0] not in a:
            continue
        out.append((list(a), [x for x in ids if x not in a]))
    return out


def run(cache: str, outdir: str) -> None:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    os.environ["PSV_CACHE"] = str(Path(cache).resolve())
    root = Path(__file__).resolve().parents[1]
    git = lambda *a: subprocess.run(["git", *a], cwd=root, capture_output=True, text=True).stdout.strip()  # noqa: E731
    prov = dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                regional_py_sha256=hashlib.sha256((root / "src/profileshape/regional.py").read_bytes()).hexdigest()[:16],
                protocol_sha256=hashlib.sha256((root / "docs/PROFILESHAPE_INDEPENDENT_VALIDATION.md").read_bytes()).hexdigest()[:16],
                started=datetime.datetime.now().isoformat(timespec="seconds"), seed=SEED)
    work = tempfile.mkdtemp(prefix="psv_", dir=out)
    (Path(work) / "dsp.R").write_text(R_SCRIPT)
    D = Data(cache)
    anc = D.anchors()
    prov["anchors"] = {k: int(np.isfinite(v).sum()) for k, v in anc.items()}
    sets = [("plus1", "nuc"), ("minus1", "nuc"), ("ndr", "sub")]
    jobs_null = [(s, z, cond, a, b, work) for s, z in sets for cond in ("0h", "1h") for a, b in splits(STH1[cond])]
    jobs_ins = [(s, "nuc", cond, a, b, ins, rep, work) for s in ("plus1", "minus1") for cond in ("0h", "1h") for a, b in splits(STH1[cond])
                for ins in ("shift20", "shift40", "broaden", "amplitude") for rep in range(2)]
    jobs_bio = [(s, z, STH1["0h"], STH1["1h"], "Sth1 4v4", work) for s, z in sets]
    jobs_bio += [(s, z, SNF2["0h"], SNF2["1h"], "Snf2 3v3", work) for s, z in sets]
    jobs_bio += [(s, z, list(r0), list(r1), f"Sth1 3v3 {'+'.join(x[-1] for x in r0)}|{'+'.join(x[-1] for x in r1)}", work)
                 for s, z in sets[:2] for r0 in itertools.combinations(STH1["0h"], 3) for r1 in itertools.combinations(STH1["1h"], 3)]
    results = []
    with ProcessPoolExecutor(6) as ex:
        for fn, jobs in ((job_null, jobs_null), (job_bio, jobs_bio), (job_insert, jobs_ins)):
            results += list(ex.map(fn, jobs))
            print(fn.__name__, "done", datetime.datetime.now().isoformat(timespec="seconds"), flush=True)
    prov["finished"] = datetime.datetime.now().isoformat(timespec="seconds")
    (Path(work) / "dsp.R").unlink(); os.rmdir(work)
    (out / "independent_validation.json").write_text(json.dumps(dict(provenance=prov, results=results), indent=1, default=float))
    (out / "INDEPENDENT_VALIDATION.md").write_text(report(prov, results))
    print((out / "INDEPENDENT_VALIDATION.md").read_text())


def report(prov, results) -> str:
    L = ["# Independent validation on GSE346229 (results)", "",
         f"Commit `{prov['commit'][:10]}`; uncommitted src files {prov['uncommitted_src']}; regional.py `{prov['regional_py_sha256']}`; "
         f"protocol sha256 `{prov['protocol_sha256']}`; anchors {prov['anchors']}; {prov['started']} → {prov['finished']}.", ""]
    nulls = [r for r in results if r["kind"] == "null"]
    L += ["## A. Real-null splits (share of regions called, BH 5 %)", "",
          "| anchor | fragments | condition | split | regions | shape_any | shape_kernel | shift | diffSplice |", "|---|---|---|---|---|---|---|---|---|"]
    for r in nulls:
        L.append(f"| {r['anchor']} | {r['size']} | {r['cond']} | {r['split']} | {r['regions']} | {r['shape_any']:.4f} | {r['shape_kernel']:.4f} | {r['shift']:.4f} | {r['diffSplice']:.4f} |")
    L += ["", "| anchor | A1: mean shape_any share ≤ 0.005 | A2: max split ≤ 0.01 |", "|---|---|---|"]
    for s in ("plus1", "minus1", "ndr"):
        v = [r["shape_any"] for r in nulls if r["anchor"] == s]
        L.append(f"| {s} | {np.mean(v):.4f} → {'pass' if np.mean(v) <= 0.005 else 'FAIL'} | {max(v):.4f} → {'pass' if max(v) <= 0.01 else 'FAIL'} |")
    ins = [r for r in results if r["kind"] == "insert"]
    L += ["", "## B. Fragment-level insertion (20 % of regions; 6 splits × 2 sets)", "",
          "| anchor | insertion | shape_any power / FDP | shape_kernel power / FDP | shift power / FDP | diffSplice power / FDP | translation coverage · direction · estimate ÷ true | criterion |",
          "|---|---|---|---|---|---|---|---|"]
    for s in ("plus1", "minus1"):
        for ins_name in ("shift20", "shift40", "broaden", "amplitude"):
            rr = [r for r in ins if r["anchor"] == s and r["insertion"] == ins_name]
            m = lambda t, k: float(np.nanmean([r[t][k] for r in rr]))  # noqa: E731
            if ins_name == "amplitude":
                cells = [f"called {m(t, 'called'):.4f}" for t in ("shape_any", "shape_kernel", "shift", "diffSplice")]
                crit = f"B2 {'pass' if m('shape_any', 'called') <= 0.01 else 'FAIL'}"
                tr = "–"
            else:
                cells = [f"{m(t, 'power'):.2f} / {m(t, 'fdp'):.3f}" for t in ("shape_any", "shape_kernel", "shift", "diffSplice")]
                crit = f"B1 {'pass' if m('shape_any', 'fdp') <= 0.075 else 'FAIL'}"
                tr = "–"
                if ins_name.startswith("shift"):
                    tr = " · ".join(f"{np.nanmean([r['translation'][k] for r in rr]):.2f}" for k in ("coverage", "direction_ok", "ratio"))
            L.append(f"| {s} | {ins_name} | " + " | ".join(cells) + f" | {tr} | {crit} |")
    bio = [r for r in results if r["kind"] == "bio"]
    L += ["", "## C. Biology", "",
          "| comparison | anchor | fragments | regions | shape_any | shape_kernel | shift | diffSplice | overlap (Jaccard) | consistent upstream / downstream | median consistent shift (bp) |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in bio:
        if r["label"].startswith("Sth1 3v3"):
            continue
        L.append(f"| {r['label']} | {r['anchor']} | {r['size']} | {r['regions']} | {r['shape_any_share']:.3f} | {r['shape_kernel_share']:.3f} | {r['shift_share']:.3f} | "
                 f"{r['diffSplice_share']:.3f} | {r['overlap_any_diffSplice']:.2f} | {r['consistent_upstream']} / {r['consistent_downstream']} | {r['median_consistent_shift_bp']:.1f} |")
    def binom(up, down, toward):
        n = up + down
        k = up if toward == "up" else down
        return k / n if n else float("nan"), stats.binomtest(k, n, 0.5, alternative="greater").pvalue if n else float("nan")
    s4 = {r["anchor"]: r for r in bio if r["label"] == "Sth1 4v4"}
    sh1, p1 = binom(s4["plus1"]["consistent_upstream"], s4["plus1"]["consistent_downstream"], "up")
    sh2, p2 = binom(s4["minus1"]["consistent_upstream"], s4["minus1"]["consistent_downstream"], "down")
    sn = {r["anchor"]: r for r in bio if r["label"] == "Snf2 3v3"}
    sub = {s: np.mean([r["shape_any_share"] for r in bio if r["label"].startswith("Sth1 3v3") and r["anchor"] == s]) for s in ("plus1", "minus1")}
    L += ["", "| criterion | result |", "|---|---|",
          f"| C1: +1 nucleosomes, consistent calls mostly upstream (share > 0.5, binomial p < 0.05) | share {sh1:.2f}, p = {p1:.2g} → {'pass' if sh1 > 0.5 and p1 < 0.05 else 'FAIL'} |",
          f"| C2: −1 nucleosomes, consistent calls mostly downstream | share {sh2:.2f}, p = {p2:.2g} → {'pass' if sh2 > 0.5 and p2 < 0.05 else 'FAIL'} |",
          f"| C3: Sth1 (mean over 3 vs 3 subsets) has more shape_any calls than Snf2 3 vs 3 | +1: {sub['plus1']:.3f} vs {sn['plus1']['shape_any_share']:.3f}; −1: {sub['minus1']:.3f} vs {sn['minus1']['shape_any_share']:.3f} → "
          f"{'pass' if sub['plus1'] > sn['plus1']['shape_any_share'] and sub['minus1'] > sn['minus1']['shape_any_share'] else 'FAIL'} |", ""]
    L += ["Notes from the Sth1 4 vs 4 runs: " + " / ".join(f"{a}: " + "; ".join(r["notes"]) for a, r in s4.items()), ""]
    return "\n".join(L) + "\n"


def aggregate(cache: str, out_md: str) -> None:
    """Exploratory (not in the protocol; added after the per-region results): replicate-level test of the
    biological expectation on per-sample meta-profiles over all regions (4 vs 4 Sth1, 3 vs 3 Snf2; Welch)."""
    D = Data(cache); A = D.anchors()
    edges = np.arange(-FLANK, FLANK + BIN, BIN); c = 0.5 * (edges[:-1] + edges[1:])

    def meta(names, anchor, size):
        cnt, _ = D.counts(names, anchor, size)
        lib = np.array([D.lib[n] for n in names])
        return cnt.sum(1) / lib[:, None] * 1e6 / cnt.shape[1]

    def centroid(m, lo, hi):
        w = (c >= lo) & (c <= hi)
        x = m[:, w] - m[:, w].min(1, keepdims=True)
        return (x * c[w]).sum(1) / x.sum(1)

    L = ["# Exploratory aggregate analysis (GSE346229; not part of the pre-specified protocol)", "",
         "Per-sample meta-profiles over all anchored regions; replicate-level Welch test 1 h vs 0 h.", "",
         "| comparison | quantity | 0 h | 1 h | difference | Welch p |", "|---|---|---|---|---|---|"]
    for lab, grp in (("Sth1 4v4", STH1), ("Snf2 3v3", SNF2)):
        rows = []
        for q, anc, size, f in (("+1 dyad centroid, bp (negative = towards NDR)", A["plus1"], NUC, lambda m: centroid(m, -80, 80)),
                                ("−1 dyad centroid, bp (positive = towards NDR)", A["minus1"], NUC, lambda m: centroid(m, -80, 80)),
                                ("NDR sub-nucleosomal signal, log2 CPM (|x| ≤ 100 bp)", A["ndr"], SUB, lambda m: np.log2(m[:, np.abs(c) <= 100].sum(1))),
                                ("−1 to +1 dyad distance, bp", A["ndr"], NUC, lambda m: centroid(m, 30, 250) - centroid(m, -250, -30))):
            a, b = f(meta(grp["0h"], anc, size)), f(meta(grp["1h"], anc, size))
            t = stats.ttest_ind(b, a, equal_var=False)
            rows.append(f"| {lab} | {q} | {np.round(a, 1).tolist()} | {np.round(b, 1).tolist()} | {b.mean() - a.mean():+.2f} | {t.pvalue:.3g} |")
        L += rows
    Path(out_md).write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    if sys.argv[1] == "extract":
        extract(sys.argv[2], sys.argv[3], sys.argv[4])
    elif sys.argv[1] == "aggregate":
        aggregate(sys.argv[2], sys.argv[3])
    else:
        run(sys.argv[2], sys.argv[3])
