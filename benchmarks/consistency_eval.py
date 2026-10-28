"""Evaluation of the within-group consistency module (protocol: docs/PROFILESHAPE_CONSISTENCY_EVALUATION.md).

    python benchmarks/consistency_eval.py build  INPUTS_DIR          # caches with region coordinates
    python benchmarks/consistency_eval.py run    INPUTS_DIR OUT_DIR  # arms A–E, results JSON + report

Arms: A simulation; B semi-real null (replicate labels permuted within blocks of 10 neighbouring regions);
C insertion on B; D real data as is; E functional consistency (real-null splits). Seeds and criteria are
fixed in the protocol; this script only implements them.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import insertion as INS  # noqa: E402

from profileshape.consistency import group_consistency  # noqa: E402
from profileshape.data import Config, Profiles, extract, read_regions  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "PROFILESHAPE_CONSISTENCY_EVALUATION.md"
# v1 (9 950 000 / 9 960 000) found a defect (near-empty windows); v2 repeats the protocol on new seeds after the
# amendment recorded in the protocol.
SEEDS = {"v1": (9_950_000, 9_960_000), "v2": (9_970_000, 9_980_000)}
SIM_SEED, SEMI_SEED = SEEDS[os.environ.get("CONSISTENCY_EVAL_SEEDS", "v2")]
MOUSE = [f"chr{i}" for i in list(range(1, 20)) + ["X"]]
SUMMIT = dict(anchor="center", orientation="oriented")
K27 = dict(SUMMIT, flank=3000, binsize=50, background_width=500, core_flank=2500)
K4 = dict(flank=2000, binsize=25, background_width=400, core_flank=1500)
# id → cache name, config, development data (excluded from the functional criterion)
DATASETS = {
    "mnase": ("cons_mnase.npz", INS.YEAST["config"], True),
    "atac": ("cons_atac.npz", dict(anchor="center", orientation="oriented", flank=600, binsize=10, background_width=200,
                                   core_flank=300), False),
    "tf_chip": ("cons_tf_chip.npz", dict(SUMMIT, flank=2000, binsize=20, background_width=500, core_flank=1000), False),
    "tf_cutrun": ("cons_tf_cutrun.npz", dict(SUMMIT, flank=1000, binsize=10, background_width=300, core_flank=600), False),
    "proseq": ("cons_proseq.npz", INS.PROSEQ["config"], False),
    "k27ac_chip": ("cons_k27ac_chip.npz", K27, False),
    "k4me3_chip_summit": ("cons_k4me3_chip_summit.npz", dict(SUMMIT, **K4), False),
    "k4me3_chip_tss": ("cons_k4me3_chip_tss.npz", dict(anchor="start", orientation="oriented", **K4), False),
    "k27ac_cutrun": ("cons_k27ac_cutrun.npz", K27, False),
    "k4me3_cutrun_summit": ("cons_k4me3_cutrun_summit.npz", dict(SUMMIT, **K4), False),
    "k4me3_cutrun_tss": ("cons_k4me3_cutrun_tss.npz", dict(anchor="start", orientation="oriented", **K4), False),
}


# ------------------------------------------------------------------------------------------------ caches
def _save(path: Path, prof: Profiles, groups: dict, counts=None, lib=None) -> None:
    extra = {} if counts is None else dict(counts=counts, lib=lib)
    np.savez_compressed(path, values=prof.values, edges=prof.edges, sample_ids=np.array([s["sample_id"] for s in prof.samples]),
                        groups=json.dumps(groups), regions=json.dumps(prof.regions), **extra)
    print(f"  {path.name}: {prof.values.shape}", flush=True)


def _bw_profiles(cfg: dict, files: dict, regions: list[dict]) -> Profiles:
    smp = [dict(sample_id=k, condition=k.rsplit("_", 1)[0], signal_file=str(f)) for k, f in files.items()]
    return extract(Config(**cfg), smp, regions)


def _groups(ids) -> dict:
    out: dict = {}
    for s in ids:
        out.setdefault(s.rsplit("_", 1)[0], []).append(s)
    return out


def _summits(paths: list[Path], sep: int) -> list[dict]:
    import cutrun_null as CR
    CR.CHROMS = MOUSE  # the caller is genome-specific; mouse main chromosomes here
    return CR.call_anchors(paths, n=5000, step=200, sep=sep)


def _tss(annotation: Path, paths: list[Path], cfg: dict) -> list[dict]:
    """RefSeq TSS (NM_, main chromosomes, one per gene), the 5000 with the highest pooled signal in ±1 kb."""
    import pyBigWig
    rows = [line.split("\t") for line in gzip.open(annotation, "rt")]
    regions, seen = [], set()
    for r in rows:
        if not r[1].startswith("NM_") or r[2] not in MOUSE or r[12] in seen:
            continue
        seen.add(r[12])
        tss = int(r[4]) if r[3] == "+" else int(r[5]) - 1
        regions.append(dict(chrom=r[2], start=tss, end=tss + 1, name=r[12], strand=r[3]))
    hs = [pyBigWig.open(str(p)) for p in paths]
    score = np.zeros(len(regions))
    for h in hs:
        size = h.chroms()
        s = np.array([np.nan_to_num(h.stats(g["chrom"], max(g["start"] - 1000, 0), min(g["start"] + 1000, size[g["chrom"]]),
                                            type="mean", exact=False)[0] or 0.0) for g in regions])
        score += s / max(s.mean(), 1e-12)
    keep = np.argsort(score)[::-1][:5000]
    return [regions[i] for i in sorted(keep)]


def build(inputs: Path, only: set | None = None) -> None:
    import pyBigWig
    cache = inputs / "caches"
    want = lambda *k: only is None or bool(set(k) & only)  # noqa: E731
    # MNase (development data): reuse the insertion cache, coordinates from the BED it was built from.
    if want("mnase"):
        _build_mnase(inputs, cache)
    if want("atac"):
        _build_atac(cache)
    if want("tf_chip"):
        _build_tf_chip(inputs, cache)
    if want("tf_cutrun"):
        _build_tf_cutrun(cache)
    if want("proseq"):
        _build_proseq(inputs, cache)
    if want("chip", "cutrun"):
        _build_histone(inputs, cache, only)


def _build_mnase(inputs, cache):
    d = np.load(cache / "ins_yeast.npz", allow_pickle=True)
    regs = read_regions(inputs / "yeast_mnase_gse280671/workspace/orf_atg.bed")
    assert len(regs) == d["values"].shape[1]
    ids = [str(s) for s in d["sample_ids"]]
    _save(cache / "cons_mnase.npz", Profiles(d["values"], d["edges"], [dict(sample_id=s) for s in ids], regs, True),
          {g: list(m) for g, m in json.loads(str(d["groups"])).items()})


def _build_atac(cache):
    # +1 dyads, nucleosomal fragments; counts kept for thinning.
    import independent_validation as IV
    D = IV.Data(str(cache / "psv_cache.npz"))
    anc = D.anchors()["plus1"]
    names = IV.STH1["0h"] + IV.STH1["1h"] + IV.SNF2["0h"] + IV.SNF2["1h"]
    prof, cnt = D.profiles(names, [n.rsplit("_", 1)[0] for n in names], anc, IV.NUC)
    ok = np.flatnonzero(np.isfinite(anc))
    prof.regions = [dict(chrom=D.genes[g]["chrom"], start=int(D.genes[g]["pos"]), end=int(D.genes[g]["pos"]) + 1,
                         name=D.genes[g]["name"], strand="+") for g in ok]
    lib = np.array([D.lib[n] for n in names])
    _save(cache / "cons_atac.npz", prof, _groups(names), counts=cnt, lib=lib)


def _build_tf_chip(inputs, cache):
    # consensus peaks
    chip = sorted((inputs / "chip_gse112560").glob("*.bw"))
    files = {}
    for f in chip:
        g = "ZNF384WT" if "ZNF384WT" in f.name else "TAF15" if "TAF15" in f.name else "TCF3"
        files[f"{g}_{f.name.split('Rep')[-1][0]}"] = f
    regs = [dict(r, strand="+") for r in read_regions(inputs / "annotation_tss/consensus_top5000.bed")]  # along the genome
    _save(cache / "cons_tf_chip.npz", _bw_profiles(DATASETS["tf_chip"][1], files, regs), _groups(files))


def _build_tf_cutrun(cache):
    # the existing profile cache (anchors 'chr:pos')
    d = np.load(cache / "cutrun_profiles.npz", allow_pickle=True)
    ids = [str(s) for s in d["sample_ids"]]
    regs = [dict(chrom=a.split(":")[0], start=int(a.split(":")[1]), end=int(a.split(":")[1]) + 1, name=a, strand="+")
            for a in d["anchors"]]
    groups = {"A673": ids[:6], "A2.2": ids[6:]}
    _save(cache / "cons_tf_cutrun.npz", Profiles(d["values"], d["edges"], [dict(sample_id=s) for s in ids], regs, True), groups)


def _build_proseq(inputs, cache):
    # as insertion.build_proseq, keeping the coordinates
    pro = inputs / "proseq_gse233930"
    files = {"NELFB_1": "GSM7472402_NELFB-DMSO-{}-dedup-Rep1.bw", "NELFB_2": "GSM7472403_NELFB-DMSO-{}-dedup-Rep2.bw",
             "PAF1_1": "GSM7472394_PAF1-DMSO-{}-dedup-Rep1.bw", "PAF1_2": "GSM7472395_PAF1-DMSO-{}-dedup-Rep2.bw",
             "RTF1_1": "GSM7472405_RTF1-DMSO-{}-dedup-Rep1.bw", "RTF1_2": "GSM7472406_RTF1-DMSO-{}-dedup-Rep2.bw"}
    rows = [line.split("\t") for line in gzip.open(pro / "hg38_ncbiRefSeqSelect.txt.gz", "rt")]
    main = {f"chr{i}" for i in list(range(1, 23)) + ["X"]}
    by, seen = {"+": [], "-": []}, set()
    for r in rows:
        if not r[1].startswith("NM_") or r[2] not in main or int(r[5]) - int(r[4]) < 2000 or r[12] in seen:
            continue
        seen.add(r[12])
        tss = int(r[4]) if r[3] == "+" else int(r[5]) - 1
        by[r[3]].append(dict(chrom=r[2], start=tss, end=tss + 1, name=r[12], strand=r[3]))
    parts = [_bw_profiles(DATASETS["proseq"][1], {k: pro / f.format(kind) for k, f in files.items()}, by[strand])
             for strand, kind in (("+", "FW"), ("-", "RV"))]
    values = np.concatenate([p.values for p in parts], axis=1)
    c = 0.5 * (parts[0].edges[:-1] + parts[0].edges[1:])
    prom = np.nansum(values[..., (c >= 0) & (c < 300)], -1).mean(0)
    keep = prom > np.quantile(prom, 0.5)
    regs = [r for r, k in zip(by["+"] + by["-"], keep, strict=True) if k]
    _save(cache / "cons_proseq.npz", Profiles(values[:, keep], parts[0].edges, parts[0].samples, regs, True), _groups(files))


def _build_histone(inputs, cache, only):
    # GSE260577 ChIP-seq, GSE285245 CUT&RUN
    import pyBigWig
    for assay, folder in (("chip", "histone_chip_gse260577"), ("cutrun", "histone_cutrun_gse285245")):
        if only is not None and assay not in only:
            continue
        for mark in ("H3K27ac", "H3K4me3"):
            paths = sorted(p for p in (inputs / folder).glob("*.bw") if mark in p.name)
            if len(paths) != 8:
                print(f"  {folder} {mark}: {len(paths)} files, expected 8 — dataset dropped", flush=True)
                continue
            files = {}
            for p in paths:
                if assay == "chip":  # GSM…_Sham1-H3K27ac-1.bw
                    cond = "Sham" if "Sham" in p.name else "CLP"
                    files[f"{cond}_{p.stem.split('-')[-1]}"] = p
                else:  # GSM…_Naive_P14_Day30_H3K27ac_Rep1_S1.scale.bw
                    cond = "Naive" if "_Naive_" in p.name else "Arm"
                    files[f"{cond}_{p.name.split('_Rep')[1][0]}"] = p
            key = "k27ac" if mark == "H3K27ac" else "k4me3"
            if key == "k27ac":
                _save(cache / DATASETS[f"k27ac_{assay}"][0],
                      _bw_profiles(DATASETS[f"k27ac_{assay}"][1], files, _summits(paths, 5000)), _groups(files))
                continue
            _save(cache / DATASETS[f"k4me3_{assay}_summit"][0],
                  _bw_profiles(DATASETS[f"k4me3_{assay}_summit"][1], files, _summits(paths, 3000)), _groups(files))
            chr1 = pyBigWig.open(str(paths[0])).chroms().get("chr1")
            build_id = {195471971: "mm10", 195154279: "mm39"}.get(chr1)
            ann = [p for p in (inputs / "annotation_mouse").glob(f"{build_id}_ncbiRefSeq*.txt.gz")] if build_id else []
            if not ann:
                print(f"  {folder} {mark}: assembly not matched (chr1 = {chr1}) — TSS set dropped", flush=True)
                continue
            _save(cache / DATASETS[f"k4me3_{assay}_tss"][0],
                  _bw_profiles(DATASETS[f"k4me3_{assay}_tss"][1], files, _tss(ann[0], paths, K4)), _groups(files))
            print(f"  {folder} {mark}: assembly {build_id}", flush=True)


# ------------------------------------------------------------------------------------------------ perturbations
def shuffle_blocks(v: np.ndarray, blocks: np.ndarray, rng, *extra):
    """Permute replicate labels independently within each block of regions (same permutation for extras)."""
    out = [v.copy()] + [e.copy() for e in extra]
    n = len(v)
    for b in np.unique(blocks):
        m = blocks == b
        perm = rng.permutation(n)
        for o, src in zip(out, (v,) + extra, strict=True):
            o[:, m] = src[perm][:, m]
    return out


def insert(v: np.ndarray, edges: np.ndarray, samples, kind: str, rng, bgw: int, frac: float = 1.0,
           same_sign: bool = True) -> np.ndarray:
    """Shift (kind 'shift<bp>') or broaden ('broaden<factor>') the profiles of ``samples`` in a random share
    ``frac`` of regions (one draw for all listed samples), as insertion.insert; ``same_sign``: +d everywhere,
    else a random sign per region."""
    v = v.copy()
    R = v.shape[1]
    regions = np.sort(rng.choice(R, int(round(frac * R)), replace=False)) if frac < 1 else np.arange(R)
    c = 0.5 * (edges[:-1] + edges[1:])
    nb = max(1, int(bgw / float(np.median(np.diff(edges)))))
    signs = np.ones(len(regions)) if same_sign else rng.choice([-1.0, 1.0], len(regions))
    for r, sg in zip(regions, signs, strict=True):
        for s in samples:
            x = np.nan_to_num(v[s, r])
            if kind.startswith("shift"):
                v[s, r] = np.interp(c - float(kind[5:]) * sg, c, x, left=x[0], right=x[-1])
            else:
                f = float(kind[7:])
                bg = np.r_[x[:nb], x[-nb:]].mean()
                y = bg + np.interp(c / f, c, x - bg)
                v[s, r] = bg + (y - bg) * max((x - bg).sum(), 1e-9) / max((y - bg).sum(), 1e-9)
    return v


INSERTIONS = [  # (name, kind, samples: 'one' | 'half', share of regions, same sign)
    ("one_shift20_all", "shift20", "one", 1.0, True),
    ("one_shift50_all", "shift50", "one", 1.0, True),
    ("one_shift50_10pc_rand", "shift50", "one", 0.1, False),
    ("one_shift50_30pc_rand", "shift50", "one", 0.3, False),
    ("one_shift50_30pc_same", "shift50", "one", 0.3, True),
    ("one_broaden1.3_all", "broaden1.3", "one", 1.0, True),
    ("half_shift20_all", "shift20", "half", 1.0, True),
    ("half_shift50_20pc_rand", "shift50", "half", 0.2, False),
]


# ------------------------------------------------------------------------------------------------ one evaluation
def summarise(prof: Profiles, cfg: Config) -> dict:
    n = len(prof.samples)
    c = group_consistency(prof, cfg, groups={"g": np.arange(n)})
    half = "+".join(s["sample_id"] for s in prof.samples[: n // 2])
    keys = ("offset_p", "offset_vs_mates", "shift_bp", "shift_bp_low", "shift_bp_high", "F10", "F20", "F50",
            "F10_mates", "F20_mates", "F50_mates", "median_abs_shift_bp")
    return dict(
        sample_flag=any("stands out" in f for f in c.flags), split_flag=any("sub-groups" in f for f in c.flags),
        flagged_samples=[f.split(": ")[1].split(" stands")[0] for f in c.flags if "stands out" in f],
        samples=[{k: r.get(k, float("nan")) for k in keys} | {"sample": r["sample"]} for r in c.sample_rows],
        splits=[dict(half_a=r["half_a"], contrast_p=r["contrast_p"], contrast_z=r["contrast_z"], pi1=r["regions_differing_share"],
                     called=r["regions_called"], shift_bp=r["shift_bp"], target=r["half_a"] == half or r["half_b"] == half)
                for r in c.partition_rows],
        regions=len(c.regions))


def _profiles(values, edges, ids, regions, oriented=True) -> Profiles:
    return Profiles(values, edges, [dict(sample_id=s, condition="g") for s in ids], regions, oriented)


def job_sim(args):
    from profileshape.simulate import SimSettings, simulate
    cell, n, depth, shallow, ins, run = args
    seed = SIM_SEED + 1000 * cell + run
    s = SimSettings(n_ref=n, n_alt=1, n_regions=2000, depth=depth)
    p, _ = simulate("null", s, seed=seed)
    if shallow:
        q, _ = simulate("null", SimSettings(n_ref=n, n_alt=1, n_regions=2000, depth=depth * 0.3), seed=seed)
        p.values[0] = q.values[0]
    v = p.values[:n]
    if ins:
        name, kind, who, frac, same = next(x for x in INSERTIONS if x[0] == ins)
        v = insert(v, p.edges, [0] if who == "one" else list(range(n // 2, n)), kind, np.random.default_rng(seed + 7),
                   Config().background_width, frac, same)
    ids = [f"r{i + 1}" for i in range(n)]
    out = summarise(_profiles(v, p.edges, ids, p.regions), Config())
    return dict(arm="A", cell=cell, n=n, depth=depth, shallow=shallow, insertion=ins or "none", run=run, **out)


@lru_cache(maxsize=4)
def _cache(path: str):
    d = np.load(path, allow_pickle=True)
    return {k: d[k] for k in d.files}


def _group_data(inputs: str, ds: str, group: str):
    name, cfgd, _ = DATASETS[ds]
    d = _cache(str(Path(inputs) / "caches" / name))
    ids = [str(s) for s in d["sample_ids"]]
    members = json.loads(str(d["groups"]))[group]
    use = [ids.index(m) for m in members]
    regions = json.loads(str(d["regions"]))
    counts = d["counts"][use] if "counts" in d else None
    lib = d["lib"][use] if "lib" in d else None
    return d["values"][use], d["edges"], members, regions, Config(**cfgd), counts, lib


def job_semi(args):
    from profileshape.consistency import region_blocks
    inputs, ds, gi, group, ins, thin, run = args
    v, edges, members, regions, cfg, counts, lib = _group_data(inputs, ds, group)
    seed = SEMI_SEED + 1000 * list(DATASETS).index(ds) + 100 * gi + run
    rng = np.random.default_rng(seed)
    blocks = region_blocks(regions, np.arange(len(regions)))
    if thin:
        libm = np.broadcast_to(lib[:, None], v.shape[:2]).copy()
        v, counts, libm = shuffle_blocks(v, blocks, rng, counts, libm)
        kept = rng.binomial(counts[0].astype(np.int64), 0.3)
        v[0] = kept / (0.3 * libm[0][:, None]) * 1e6
    else:
        (v,) = shuffle_blocks(v, blocks, rng)
    n = len(members)
    if ins:
        name, kind, who, frac, same = next(x for x in INSERTIONS if x[0] == ins)
        v = insert(v, edges, [0] if who == "one" else list(range(n // 2, n)), kind, np.random.default_rng(seed + 7),
                   cfg.background_width, frac, same)
    out = summarise(_profiles(v, edges, [f"{group}_p{i + 1}" for i in range(n)], regions), cfg)
    return dict(arm="C" if ins else "B", dataset=ds, group=group, n=n, insertion=ins or "none", thinned=thin, run=run, **out)


def job_real(args):
    inputs, ds, group = args
    v, edges, members, regions, cfg, _, _ = _group_data(inputs, ds, group)
    c = group_consistency(_profiles(v, edges, members, regions), cfg, groups={"g": np.arange(len(members))})
    return dict(arm="D", dataset=ds, group=group, n=len(members), flags=c.flags, samples=c.sample_rows,
                splits=c.partition_rows, groups=c.group_rows, regions=len(c.regions))


def job_func(args):
    from profileshape.analysis import analyse
    from profileshape.regional import per_region_tests
    inputs, ds, group, a_idx = args
    v, edges, members, regions, cfg, _, _ = _group_data(inputs, ds, group)
    n = len(members)
    b_idx = [i for i in range(n) if i not in a_idx]
    order = list(a_idx) + b_idx
    smp = [dict(sample_id=members[i], condition="A" if k < len(a_idx) else "B") for k, i in enumerate(order)]
    cfgd = dict(vars(cfg))
    cfgd.update(reference="A", alternative="B")
    cfg2 = Config(**{k: val for k, val in cfgd.items() if k in Config.__dataclass_fields__})
    a = analyse(Profiles(v[order], edges, smp, [regions[i] for i in range(len(regions))], True), cfg2, bootstrap=0)
    res = per_region_tests(a)
    share = float(np.isin(res.tests["shape_any"]["call"], ("Changed",)).mean()) if "shape_any" in res.tests else float("nan")
    return dict(arm="E", dataset=ds, group=group, half_a="+".join(members[i] for i in a_idx),
                half_b="+".join(members[i] for i in b_idx), shape_any_null_share=share)


# ------------------------------------------------------------------------------------------------ driver
def provenance() -> dict:
    git = lambda *a: subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()  # noqa: E731
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()[:16]  # noqa: E731
    return dict(commit=git("rev-parse", "HEAD"), uncommitted_src=len(git("status", "--porcelain", "--", "src").splitlines()),
                src_diff_sha256=hashlib.sha256(git("diff", "HEAD", "--", "src").encode()).hexdigest()[:16],
                regional_py_sha256=sha(ROOT / "src/profileshape/regional.py"),
                consistency_py_sha256=sha(ROOT / "src/profileshape/consistency.py"),
                protocol_sha256=sha(PROTOCOL), script_sha256=sha(Path(__file__)),
                seeds=f"simulation {SIM_SEED}+, semi-real {SEMI_SEED}+", started=time.strftime("%Y-%m-%dT%H:%M:%S"))


def jobs_for(inputs: str):
    sim, cell = [], 0
    for n in (3, 4, 6):
        for depth in (12, 40, 160):
            for shallow in (False, True):
                sim += [(cell, n, depth, shallow, None, r) for r in range(100)]
                cell += 1
        for name, _, who, _, _ in INSERTIONS:
            if who == "half" and n < 4:
                continue
            sim += [(cell, n, 40, False, name, r) for r in range(100)]
            cell += 1
    semi, real, func = [], [], []
    for ds, (name, _, dev) in DATASETS.items():
        path = Path(inputs) / "caches" / name
        if not path.exists():
            continue
        groups = json.loads(str(np.load(path, allow_pickle=True)["groups"]))
        for gi, (g, members) in enumerate(groups.items()):
            n = len(members)
            real.append((inputs, ds, g))
            semi += [(inputs, ds, gi, g, None, False, r) for r in range(20)]
            if ds == "atac" and n == 4:
                semi += [(inputs, ds, gi, g, None, True, r) for r in range(20)]
            if n >= 3:
                for ins_name, _, who, _, _ in INSERTIONS:
                    if who == "half" and n < 4:
                        continue
                    semi += [(inputs, ds, gi, g, ins_name, False, r) for r in range(20)]
            if n in (4, 6):
                for a in itertools.combinations(range(n), n // 2):
                    if 0 in a:
                        func.append((inputs, ds, g, a))
    return sim, semi, real, func


def run(inputs: str, out: str, workers: int = 10) -> None:
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    prov = provenance()
    sim, semi, real, func = jobs_for(inputs)
    print(f"jobs: sim {len(sim)}, semi-real {len(semi)}, real {len(real)}, functional {len(func)}", flush=True)
    results = {}
    with ProcessPoolExecutor(workers) as ex:
        for arm, fn, jobs in (("real", job_real, real), ("functional", job_func, func), ("semireal", job_semi, semi),
                              ("sim", job_sim, sim)):
            t0 = time.time()
            results[arm] = list(ex.map(fn, jobs, chunksize=4))
            print(f"  {arm}: {len(jobs)} jobs in {time.time() - t0:.0f} s", flush=True)
            (out_dir / f"{arm}.json").write_text(json.dumps(results[arm], default=_json))
    prov["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    (out_dir / "provenance.json").write_text(json.dumps(prov, indent=2))
    (out_dir / "CONSISTENCY_EVALUATION.md").write_text(report(prov, results))
    print((out_dir / "CONSISTENCY_EVALUATION.md").read_text())


def _json(v):
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist()
    raise TypeError(type(v))


# ------------------------------------------------------------------------------------------------ report
def _rate(rows, key):
    return float(np.mean([r[key] for r in rows])) if rows else float("nan")


def report(prov: dict, R: dict) -> str:
    from scipy import stats
    sim, semi, real, func = R["sim"], R["semireal"], R["real"], R["functional"]
    L = ["# Within-group consistency — evaluation results", "",
         f"Protocol `docs/PROFILESHAPE_CONSISTENCY_EVALUATION.md` (sha256 `{prov['protocol_sha256']}`); commit "
         f"`{prov['commit'][:10]}`, uncommitted src files {prov['uncommitted_src']} (diff `{prov['src_diff_sha256']}`); "
         f"consistency.py `{prov['consistency_py_sha256']}`; regional.py `{prov['regional_py_sha256']}`; script "
         f"`{prov['script_sha256']}`; seeds {prov.get('seeds', '')}; {prov['started']} → {prov['finished']}.", ""]
    verdict = {}

    # ---- A: null cells
    L += ["## A · Simulation, null cells (100 runs each, 2000 regions)", "",
          "| n | depth | replicate 1 at 30 % depth | any flag | sample flag | split flag | offset p < 0.01 | split p < 0.01 | bp CI covers 0 | excess F20 | split π̂₁ |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    null = [r for r in sim if r["insertion"] == "none"]
    worst = dict(flag=0.0, cover=1.0, f20=0.0, pi=0.0)
    pools = {False: dict(off=[], spl=[]), True: dict(off=[], spl=[])}
    for key in sorted({(r["n"], r["depth"], r["shallow"]) for r in null}):
        rows = [r for r in null if (r["n"], r["depth"], r["shallow"]) == key]
        anyf = np.mean([r["sample_flag"] or r["split_flag"] for r in rows])
        offp = [s["offset_p"] for r in rows for s in r["samples"]]
        splp = [s["contrast_p"] for r in rows for s in r["splits"]]
        pools[key[2]]["off"] += offp
        pools[key[2]]["spl"] += splp
        cover = np.mean([s["shift_bp_low"] <= 0 <= s["shift_bp_high"] for r in rows for s in r["samples"]])
        f20 = np.mean([s["F20"] - s["F20_mates"] for r in rows for s in r["samples"] if np.isfinite(s["F20_mates"])])
        pi = np.nanmean([s["pi1"] for r in rows for s in r["splits"]]) if splp else float("nan")
        worst["flag"] = max(worst["flag"], anyf)
        worst["cover"] = min(worst["cover"], cover)
        if not key[2]:
            worst["f20"] = max(worst["f20"], abs(f20))
        if np.isfinite(pi):
            worst["pi"] = max(worst["pi"], pi)
        L.append(f"| {key[0]} | {key[1]} | {'yes' if key[2] else 'no'} | {anyf:.2f} | {_rate(rows, 'sample_flag'):.2f} | "
                 f"{_rate(rows, 'split_flag'):.2f} | {np.mean(np.array(offp) < 0.01):.3f} | "
                 f"{np.mean(np.array(splp) < 0.01) if splp else float('nan'):.3f} | {cover:.2f} | {f20:+.3f} | {pi:.3f} |")
    verdict["K1"] = worst["flag"] <= 0.04
    k3 = []
    for sh, pool in pools.items():
        for kind, vals in pool.items():
            x = float(np.mean(np.array(vals) < 0.01))
            k3.append(0.005 <= x <= 0.02)
            L.append(f"\nPooled {'unequal' if sh else 'equal'}-depth {('offset' if kind == 'off' else 'split contrast')} "
                     f"P(p < 0.01) = {x:.4f} (m = {len(vals)})")
    verdict["K3"] = all(k3)

    # ---- B: semi-real null
    L += ["", "## B · Semi-real null (labels permuted within blocks of 10 regions; 20 shuffles)", "",
          "| dataset | group | n | regions | runs with any flag | offset p < 0.01 | split p < 0.01 | bp CI covers 0 | excess F20 | split π̂₁ |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    bnull = [r for r in semi if r["arm"] == "B"]
    k2, k4b, k6b, k7b = [], [], [], []
    for key in sorted({(r["dataset"], r["group"], r["thinned"]) for r in bnull}):
        rows = [r for r in bnull if (r["dataset"], r["group"], r["thinned"]) == key]
        nf = int(sum(r["sample_flag"] or r["split_flag"] for r in rows))
        k2.append(nf <= 2)
        offp = np.array([s["offset_p"] for r in rows for s in r["samples"]])
        splp = np.array([s["contrast_p"] for r in rows for s in r["splits"]])
        cov = [s["shift_bp_low"] <= 0 <= s["shift_bp_high"] for r in rows for s in r["samples"] if np.isfinite(s["shift_bp"])]
        f20 = [s["F20"] - s["F20_mates"] for r in rows for s in r["samples"] if np.isfinite(s["F20_mates"])]
        pi = [s["pi1"] for r in rows for s in r["splits"] if np.isfinite(s["pi1"])]
        if cov:
            k4b.append(np.mean(cov) >= 0.90)
        if f20 and not key[2]:
            k6b.append(abs(np.mean(f20)) <= 0.03)
        if pi:
            k7b.append(np.mean(pi) <= 0.02)
        L.append(f"| {key[0]}{' (thinned)' if key[2] else ''} | {key[1]} | {rows[0]['n']} | {rows[0]['regions']} | {nf} / {len(rows)} | "
                 f"{np.mean(offp < 0.01):.3f} | {np.mean(splp < 0.01) if splp.size else float('nan'):.3f} | "
                 f"{np.mean(cov) if cov else float('nan'):.2f} | {np.mean(f20) if f20 else float('nan'):+.3f} | "
                 f"{np.mean(pi) if pi else float('nan'):.3f} |")
    verdict["K2"] = all(k2)
    verdict["K4"] = worst["cover"] >= 0.90 and all(k4b)
    verdict["K6"] = worst["f20"] <= 0.03 and all(k6b)
    verdict["K7"] = worst["pi"] <= 0.02 and all(k7b)

    # ---- A and C: insertions
    def ins_table(rows_all, by):
        out = []
        for key in sorted({by(r) for r in rows_all}):
            rows = [r for r in rows_all if by(r) == key]
            s0 = [r["samples"][0] for r in rows]
            tgt = [s for r in rows for s in r["splits"] if s["target"]]
            out.append((key, rows, s0, tgt))
        return out
    L += ["", "## A · Simulation, insertions (depth 40; replicate 1 or the second half altered)", "",
          "| n | insertion | sample flag (replicate 1 named) | split flag | shift bp (replicate 1) | F20 − mates | F50 − mates | target split π̂₁ | target split called |",
          "|---|---|---|---|---|---|---|---|---|"]
    k5, k8, k9 = [], [], []
    ins_sim = [r for r in sim if r["insertion"] != "none"]
    for (n, ins), rows, s0, tgt in ins_table(ins_sim, lambda r: (r["n"], r["insertion"])):
        named = np.mean(["r1" in r["flagged_samples"] for r in rows])
        bp = np.mean([s["shift_bp"] for s in s0])
        if ins in ("one_shift20_all", "one_shift50_all"):
            k5.append(0.8 <= bp / float(ins.split("shift")[1].split("_")[0]) <= 1.1)
        if ins == "one_shift20_all":
            k9.append(named >= 0.8)
        if ins == "half_shift20_all":
            k9.append(_rate(rows, "split_flag") >= 0.8)
        if ins == "half_shift50_20pc_rand":
            k8.append(np.nanmean([s["pi1"] for s in tgt]) <= 0.25)
        L.append(f"| {n} | {ins} | {named:.2f} | {_rate(rows, 'split_flag'):.2f} | {bp:+.1f} | "
                 f"{np.mean([s['F20'] - s['F20_mates'] for s in s0]):+.3f} | {np.mean([s['F50'] - s['F50_mates'] for s in s0]):+.3f} | "
                 f"{np.mean([s['pi1'] for s in tgt]) if tgt else float('nan'):.3f} | "
                 f"{np.mean([s['called'] for s in tgt]) if tgt else float('nan'):.0f} |")
    verdict["K5"], verdict["K8"], verdict["K9"] = all(k5), all(k8), all(k9)
    L += ["", "## C · Semi-real insertions (20 shuffles each; no criterion)", "",
          "| dataset | group | n | insertion | sample flag (pseudo-replicate 1 named) | split flag | shift bp | F20 − mates | target split π̂₁ |",
          "|---|---|---|---|---|---|---|---|---|"]
    for (ds, g, ins), rows, s0, tgt in ins_table([r for r in semi if r["arm"] == "C"], lambda r: (r["dataset"], r["group"], r["insertion"])):
        named = np.mean([f"{g}_p1" in r["flagged_samples"] for r in rows])
        L.append(f"| {ds} | {g} | {rows[0]['n']} | {ins} | {named:.2f} | {_rate(rows, 'split_flag'):.2f} | "
                 f"{np.mean([s['shift_bp'] for s in s0]):+.1f} | {np.mean([s['F20'] - s['F20_mates'] for s in s0]):+.3f} | "
                 f"{np.mean([s['pi1'] for s in tgt]) if tgt else float('nan'):.3f} |")

    # ---- D: real data
    L += ["", "## D · Real data as is (descriptive; no truth)", ""]
    for r in real:
        L.append(f"**{r['dataset']} · {r['group']}** (n = {r['n']}, {r['regions']} regions): "
                 + ("; ".join(r["flags"]) if r["flags"] else "no flag"))
        for s in r["samples"]:
            L.append(f"- {s['sample']}: offset {100 * s['offset_relative']:.1f} % of profile ({s.get('offset_vs_mates', float('nan')):.2f}× mates), "
                     f"shift {s['shift_bp']:+.1f} bp [{s['shift_bp_low']:+.1f}, {s['shift_bp_high']:+.1f}], "
                     f"F20 {s['F20']:.2f} (mates {s.get('F20_mates', float('nan')):.2f})")
        for s in sorted(r["splits"], key=lambda x: x["contrast_p"])[:3]:
            L.append(f"- split {s['half_a']} | {s['half_b']}: z {s['contrast_z']:+.2f}, p {s['contrast_p']:.2g}, "
                     f"π̂₁ {s['regions_differing_share']:.3f}, called {s['regions_called']}, shift {s['shift_bp']:+.1f} bp")
        L.append("")

    # ---- E: functional
    L += ["## E · Functional consistency (real-null splits)", "",
          "| dataset | group | split | contrast z | split flagged | shape_any null-call share |", "|---|---|---|---|---|---|"]
    z_of = {}
    flagged_split = {}
    for r in real:
        best = min(r["splits"], key=lambda s: s["contrast_p"]) if r["splits"] else None
        for s in r["splits"]:
            z_of[(r["dataset"], r["group"], s["half_a"], s["half_b"])] = s["contrast_z"]
            flagged_split[(r["dataset"], r["group"], s["half_a"], s["half_b"])] = (
                any("sub-groups" in f for f in r["flags"]) and s is best)
    pairs = []
    for f in func:
        key = (f["dataset"], f["group"], f["half_a"], f["half_b"])
        if key not in z_of:  # the split listed with halves swapped
            key = (f["dataset"], f["group"], f["half_b"], f["half_a"])
        z, fl = z_of.get(key, float("nan")), flagged_split.get(key, False)
        dev = DATASETS[f["dataset"]][2]
        if not dev:
            pairs.append((z, f["shape_any_null_share"], fl))
        L.append(f"| {f['dataset']}{' (development)' if dev else ''} | {f['group']} | {f['half_a']} vs {f['half_b']} | {z:+.2f} | "
                 f"{'yes' if fl else ''} | {f['shape_any_null_share']:.4f} |")
    if len(pairs) >= 5:
        zs, sh, fl = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs]), np.array([p[2] for p in pairs])
        rho, p_two = stats.spearmanr(zs, sh)
        p_one = p_two / 2 if rho > 0 else 1 - p_two / 2
        part2 = (sh[fl].mean() > sh[~fl].mean()) if fl.any() and (~fl).any() else None
        # No variation in the outcome (no null call in any split): the association cannot be assessed.
        verdict["K10"] = None if np.ptp(sh) == 0 else bool(rho > 0 and p_one < 0.05 and (part2 is None or part2))
        L += ["", f"New datasets: {len(pairs)} splits; Spearman(contrast z, null-call share) = {rho:+.2f}, one-sided p = {p_one:.3g}; "
              f"flagged splits {int(fl.sum())} (mean share {sh[fl].mean() if fl.any() else float('nan'):.4f}) vs unflagged "
              f"{int((~fl).sum())} (mean {sh[~fl].mean():.4f})."]
    else:
        verdict["K10"] = None
    L += ["", "## Criteria", "", "| id | result |", "|---|---|"]
    for k in ("K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K10"):
        v = verdict.get(k)
        L.append(f"| {k} | {'pass' if v else 'not evaluable' if v is None else '**FAIL**'} |")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("build", "run", "report"))
    ap.add_argument("inputs")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--only", nargs="*", help="build: dataset ids or 'chip' / 'cutrun'")
    a = ap.parse_args(argv)
    if a.command == "build":
        build(Path(a.inputs), set(a.only) if a.only else None)
    elif a.command == "run":
        run(a.inputs, a.out, a.workers)
    else:
        d = Path(a.out)
        R = {k: json.loads(gzip.open(d / f"{k}.json.gz", "rt").read() if (d / f"{k}.json.gz").exists()
                           else (d / f"{k}.json").read_text()) for k in ("sim", "semireal", "real", "functional")}
        (d / "CONSISTENCY_EVALUATION.md").write_text(report(json.loads((d / "provenance.json").read_text()), R))
    return 0


if __name__ == "__main__":
    sys.exit(main())
