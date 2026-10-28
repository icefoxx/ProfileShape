"""Step 2 of docs/PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md: full-run depth estimated from pilot alignments.

    python benchmarks/pilot_estimate.py PILOT_DIR OUT.json

For each candidate: genes ≥ 5 kb (RefSeq Select, one TSS per gene); the condition-blind top half by pooled pilot
fragments in TSS ± 300 bp (as in the first attempt); fragments with a midpoint in the core window (TSS ± 500 bp)
per region and pilot sample, scaled by the run's read pairs (reads) / the pilot's. Same `extract` path and filters
as the full analysis would use.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import numpy as np

from profileshape import metrics as M
from profileshape.data import ASSAY_PRESETS, Config, extract

ROOT = Path(__file__).resolve().parents[1]
PILOT_N = 2_000_000
CFG = dict(ASSAY_PRESETS["proseq"], orientation="oriented")  # flank 1000, bin 10, core 500
ANNOT = {"hg38": ROOT / "inputs/proseq_gse233930/hg38_ncbiRefSeqSelect.txt.gz",
         "mm10": ROOT / "inputs/annotation_mouse/mm10_ncbiRefSeqSelect.txt.gz"}
CANDIDATES = {
    "GSE112379 heat shock RPB1 ChIP": ("hg38", {"SRR6904142": ("NHS", 47_710_698), "SRR6904139": ("HS", 59_369_910)}),
    "GSE199086 KL-2 FP-TT-seq": ("hg38", {"SRR18401340": ("DMSO", 7_909_263), "SRR18401334": ("KL2", 7_529_382)}),
    "GSE205504 NELFB-AID Pol II ChIP": ("mm10", {"SRR19544320": ("EV_IAA", 47_115_975), "SRR19544316": ("NELFB_IAA", 30_118_664)}),
}


def genes(path: Path) -> list[dict]:
    main = {f"chr{i}" for i in list(range(1, 23)) + ["X"]}
    out, seen = [], set()
    for r in (line.split("\t") for line in gzip.open(path, "rt")):
        if not r[1].startswith("NM_") or r[2] not in main or int(r[5]) - int(r[4]) < 5000 or r[12] in seen:
            continue
        seen.add(r[12])
        tss = int(r[4]) if r[3] == "+" else int(r[5]) - 1
        out.append(dict(chrom=r[2], start=tss, end=tss + 1, name=r[12], strand=r[3]))
    return out


def complexity(bam: str, seed: int = 1) -> dict:
    """Distinct fragments at 25 / 50 / 100 % of the pilot and a least-squares Poisson complexity C (amendment R3):
    distinct(n) = C (1 − exp(−n / C))."""
    import pysam
    from scipy.optimize import minimize_scalar
    keys = []
    for a in pysam.AlignmentFile(bam):
        if a.is_secondary or a.is_supplementary or a.is_unmapped:
            continue
        if a.is_paired:
            if a.is_read1:
                keys.append(hash((a.reference_id, a.reference_start, a.next_reference_start, a.template_length)))
        else:
            keys.append(hash((a.reference_id, a.reference_end if a.is_reverse else a.reference_start, a.is_reverse)))
    keys = np.array(keys)
    np.random.default_rng(seed).shuffle(keys)
    pts = [(int(f * len(keys)), len(np.unique(keys[: int(f * len(keys))]))) for f in (0.25, 0.5, 1.0)]
    n, u = np.array(pts, float).T
    loss = lambda logc: float(np.sum((np.exp(logc) * (1 - np.exp(-n / np.exp(logc))) - u) ** 2))  # noqa: E731
    C = float(np.exp(minimize_scalar(loss, bounds=(np.log(u[-1]), np.log(1e11)), method="bounded").x))
    fit = C * (1 - np.exp(-n / C))
    return dict(points=pts, C=C, max_rel_residual=float(np.max(np.abs(fit - u) / u)))


def estimate(pilot: Path, genome: str, runs: dict) -> dict:
    regs = genes(ANNOT[genome])
    smp = [dict(sample_id=r, condition=c, signal_file=str(pilot / f"{r}.bam")) for r, (c, _) in runs.items()]
    prof = extract(Config(**CFG), smp, regs)
    c = np.asarray(prof.counts)
    mids = M.centers(prof.edges)
    prom = c[..., np.abs(mids) <= 300].sum(-1)
    pooled = (prom / np.maximum(prom.sum(1, keepdims=True), 1)).sum(0)
    keep = pooled >= np.median(pooled)
    core = c[..., M.core_mask(prof.edges, CFG["core_flank"])].sum(-1)[:, keep]
    out = dict(genome=genome, genes=len(regs), selected=int(keep.sum()), samples={})
    for i, (r, (cond, total)) in enumerate(runs.items()):
        med = float(np.median(core[i]))
        q = prof.qc.get(r, {})
        cx = complexity(str(pilot / f"{r}.bam"))
        n_pilot, u_pilot = cx["points"][-1]
        n_full = n_pilot * total / PILOT_N  # pilot reads after alignment and filters, scaled to the run
        u_full = cx["C"] * (1 - np.exp(-n_full / cx["C"]))
        out["samples"][r] = dict(condition=cond, pilot_median=med, run_reads=total, linear_scale=total / PILOT_N,
                                 linear_full_median=med * total / PILOT_N, complexity=cx, complexity_scale=u_full / u_pilot,
                                 estimated_full_median=med * u_full / u_pilot, paired=q.get("paired"),
                                 fragment_length=q.get("fragment_length"), fragment_median=q.get("median"))
    out["estimated_full_median_over_pilots"] = float(np.median([s["estimated_full_median"] for s in out["samples"].values()]))
    out["linear_full_median_over_pilots"] = float(np.median([s["linear_full_median"] for s in out["samples"].values()]))
    # condition-blind templates for step 3: pooled pilot counts of the selected regions; fragment lengths
    qc = [prof.qc[r] for r in runs]
    if all(x.get("paired") for x in qc):  # paired: quartiles and median of the observed fragment lengths
        lengths = [float(np.median([x[k] for x in qc])) for k in ("q25", "median", "q75")]
    else:  # single-end: the estimated fragment length
        lengths = [float(np.median([x["fragment_length"] for x in qc]))] * 3
    np.savez_compressed(pilot / f"templates_{genome}_{list(runs)[0]}.npz", pooled=c[:, keep].sum(0), edges=prof.edges,
                        fragment_lengths=np.array(lengths))
    out["templates"] = str(pilot / f"templates_{genome}_{list(runs)[0]}.npz")
    return out


def main(argv):
    pilot, out = Path(argv[0]), Path(argv[1])
    res = {k: estimate(pilot, g, runs) for k, (g, runs) in CANDIDATES.items() if all((pilot / f"{r}.bam").exists() for r in runs)}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=str))
    for k, v in res.items():
        print(k, "| selected", v["selected"], "| complexity-adjusted", round(v["estimated_full_median_over_pilots"]),
              {r: (round(s["estimated_full_median"]), round(s["complexity"]["C"] / 1e6, 1), round(s["complexity"]["max_rel_residual"], 4))
               for r, s in v["samples"].items()}, "| linear", round(v["linear_full_median_over_pilots"]))


if __name__ == "__main__":
    main(sys.argv[1:])
