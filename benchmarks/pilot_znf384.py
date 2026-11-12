"""Step 2 of docs/PROFILESHAPE_ZNF384_SCENARIO.md: pilot depth for the two ZNF384 designs.

    python benchmarks/pilot_znf384.py PILOT_DIR PEAK_DIR OUT.json

(ii) GSE181617 +e8 / −e8: pilot BAMs of `pilot_depth.sh`, regions = summits of the two samples' GEO narrowPeak files
merged within 200 bp, `tf_chip` preset, condition-blind top half by pooled core fragments; full-run depth by
region-level complexity (R5), with library-wide complexity and linear scaling alongside (as `pilot_estimate.py`).
(i) GSE112560 WT / TCF3-ZNF384: processed bigWigs at `consensus_top5000.bed`, fragments = area / (quantum × footprint).
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np

from pilot_estimate import PILOT_N, complexity
from profileshape import metrics as M
from profileshape.data import ASSAY_PRESETS, Config, bigwig_count_unit, extract

ROOT = Path(__file__).resolve().parents[1]
CFG = dict(ASSAY_PRESETS["tf_chip"])  # flank 2000, bin 20, core 1000
RUNS = {"SRR15365234": ("plus_e8", "GSM5508181", 46_769_491), "SRR15365232": ("minus_e8", "GSM5508179", 39_746_088)}
MAIN = {f"chr{i}" for i in list(range(1, 20)) + ["X"]}


def summits(peak_dir: Path) -> list[dict]:
    pts = []
    for _, gsm, _ in RUNS.values():
        for r in (x.split("\t") for x in gzip.open(peak_dir / f"{gsm}.narrowPeak.gz", "rt")):
            if r[0] in MAIN:
                pts.append((r[0], int(r[1]) + int(r[9])))
    pts.sort()
    out, cur = [], []
    for c, p in pts + [("", 0)]:  # merge summits within 200 bp; the merged centre is their mean
        if cur and (c != cur[-1][0] or p - cur[-1][1] > 200):
            m = int(np.mean([q for _, q in cur]))
            out.append(dict(chrom=cur[0][0], start=m, end=m + 1, name=f"{cur[0][0]}_{m}", strand="."))
            cur = []
        cur.append((c, p))
    return out


def gse181617(pilot: Path, peak_dir: Path) -> dict:
    regs = summits(peak_dir)
    smp = [dict(sample_id=r, condition=c, signal_file=str(pilot / f"{r}.bam")) for r, (c, _, _) in RUNS.items()]
    prof = extract(Config(**CFG), smp, regs)
    c = np.asarray(prof.counts)
    core = c[..., M.core_mask(prof.edges, CFG["core_flank"])].sum(-1)
    pooled = core.sum(0)
    keep = pooled >= np.median(pooled)
    out = dict(regions=len(regs), selected=int(keep.sum()), samples={})
    win = [(g["chrom"], g["start"] - CFG["core_flank"], g["start"] + CFG["core_flank"]) for g, k in zip(regs, keep, strict=True) if k]
    for i, (r, (cond, gsm, total)) in enumerate(RUNS.items()):
        med = float(np.median(core[i, keep]))
        cx = complexity(str(pilot / f"{r}.bam"))
        cr = complexity(str(pilot / f"{r}.bam"), windows=win)
        n_p, u_p = cx["points"][-1]
        u_full = cx["C"] * (1 - np.exp(-(n_p * total / PILOT_N) / cx["C"]))
        nr, ur = cr["points"][-1]
        ur_full = cr["C"] * (1 - np.exp(-(nr * total / PILOT_N) / cr["C"]))
        q = prof.qc.get(r, {})
        out["samples"][r] = dict(condition=cond, gsm=gsm, run_pairs=total, pilot_median=med,
                                 linear_full_median=med * total / PILOT_N, estimated_full_median=med * u_full / u_p,
                                 region_full_median=med * ur_full / ur, complexity=cx, region_complexity=cr,
                                 pilot_core_by_condition_selected=float(core[i, keep].sum()),
                                 fragment_median=q.get("median"), fragment_q25=q.get("q25"), fragment_q75=q.get("q75"))
    s = out["samples"].values()
    for k in ("region_full_median", "estimated_full_median", "linear_full_median"):
        out[k + "_over_pilots"] = float(np.median([x[k] for x in s]))
    return out


def gse112560() -> dict:
    import pyBigWig
    regs = [x.split("\t") for x in open(ROOT / "inputs/annotation_tss/consensus_top5000.bed")]
    out = {}
    for f in sorted((ROOT / "inputs/chip_gse112560").glob("*.bw")):
        unit = bigwig_count_unit(str(f))
        h = pyBigWig.open(str(f))
        fr = []
        for r in regs:
            mid = int(r[1])
            lo, hi = max(mid - CFG["core_flank"], 0), min(mid + CFG["core_flank"], h.chroms().get(r[0], 0))
            v = np.nan_to_num(np.asarray(h.values(r[0], lo, hi), float)) if hi > lo else np.zeros(0)
            fr.append(v.sum() / (unit[0] * unit[1]) if unit else np.nan)
        h.close()
        out[f.name.split("_")[0] + "_" + f.name.split("_")[2].rsplit("-", 1)[0]] = dict(file=f.name, quantum=unit[0] if unit else None, footprint_bp=unit[1] if unit else None,
                                         core_median=float(np.nanmedian(fr)), core_q25=float(np.nanpercentile(fr, 25)))
    return out


def main(argv):
    pilot, peaks, out = Path(argv[0]), Path(argv[1]), Path(argv[2])
    res = dict(GSE181617=gse181617(pilot, peaks), GSE112560=gse112560())
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=str))
    a = res["GSE181617"]
    print("GSE181617 regions", a["regions"], "selected", a["selected"], "| region-level", round(a["region_full_median_over_pilots"]),
          "| library-wide", round(a["estimated_full_median_over_pilots"]), "| linear", round(a["linear_full_median_over_pilots"]),
          {r: (s["condition"], round(s["region_full_median"]), s["fragment_median"]) for r, s in a["samples"].items()})
    print("GSE112560", {k: (round(v["core_median"]), v["footprint_bp"]) for k, v in res["GSE112560"].items()})


if __name__ == "__main__":
    main(sys.argv[1:])
