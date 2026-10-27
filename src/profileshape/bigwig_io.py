"""bigWig / bedGraph window reader (from GenomeCanvas `pipeline/genome_bins.py`, MIT licence).

`window_matrix` returns binned signal for many windows with the semantics of per-window
``pyBigWig.stats(nBins=…)``: ``mean`` over covered bases, 0 where a bin has no data, NaN beyond chromosome ends.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

DENSE_BYTES_PER_BASE = 0.1  # 1-bp BigWigs are ~5 B/base; 10 bp ~0.5; sparse bins ~0.005


class SignalReader:
    """Interval access for BigWig or bedGraph files."""

    def __init__(self, path: str | Path, read_mode: str = "auto"):
        self.path = str(path)
        self.bigwig = self.path.lower().endswith((".bw", ".bigwig"))
        self.dense = False
        if self.bigwig:
            import pyBigWig

            self.handle = pyBigWig.open(self.path)
            if self.handle is None:
                raise ValueError(f"Cannot open BigWig {self.path}")
            self.sizes = dict(self.handle.chroms())
            density = Path(self.path).stat().st_size / max(1, sum(self.sizes.values()))
            self.dense = read_mode == "values" or (
                read_mode == "auto" and density > DENSE_BYTES_PER_BASE
            )
            self.dense = self.dense and bool(pyBigWig.numpy)
        else:
            self._load_bedgraph()

    def _load_bedgraph(self):
        rows: dict[str, list] = {}
        with open(self.path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                chrom, start, end, value = line.split()[:4]
                rows.setdefault(chrom, []).append((int(start), int(end), float(value)))
        self.data = {}
        for chrom, items in rows.items():
            array = np.array(sorted(items), dtype=float)
            self.data[chrom] = (array[:, 0], array[:, 1], array[:, 2])
        self.sizes = {chrom: int(arrays[1].max()) for chrom, arrays in self.data.items()}

    def intervals(self, chrom: str, start: int, end: int):
        """Sorted (starts, ends, values) clipped to [start, end); NaN values dropped."""
        if self.bigwig:
            raw = (
                self.handle.intervals(chrom, int(start), int(end)) if chrom in self.sizes else None
            )
            if not raw:
                return _empty()
            array = np.asarray(raw, dtype=float)
            s, e, v = array[:, 0], array[:, 1], array[:, 2]
        else:
            if chrom not in self.data:
                return _empty()
            s, e, v = self.data[chrom]
            lo = np.searchsorted(e, start, side="right")
            hi = np.searchsorted(s, end, side="left")
            s, e, v = s[lo:hi], e[lo:hi], v[lo:hi]
        keep = np.isfinite(v)
        s, e, v = np.clip(s[keep], start, end), np.clip(e[keep], start, end), v[keep]
        nonempty = e > s
        return s[nonempty], e[nonempty], v[nonempty]

    def values(self, chrom: str, start: int, end: int):
        """Per-base values (NaN = no data) for [start, end) of a BigWig."""
        start, end = int(start), int(end)
        out = np.full(end - start, np.nan)
        length = self.sizes.get(chrom)
        if length is None:
            return out
        lo, hi = max(0, start), min(end, length)  # windows may run past a chromosome edge
        if hi > lo:
            out[lo - start : hi - start] = self.handle.values(chrom, lo, hi, numpy=True)
        return out

    def close(self):
        if self.bigwig:
            self.handle.close()


def _empty():
    z = np.zeros(0)
    return z, z, z


def _integrals(s, e, v, positions):
    """∫signal and ∫coverage from the first interval to each position (piecewise linear)."""
    if not s.size:
        zero = np.zeros(positions.shape)
        return zero, zero
    width = e - s
    cum_sum = np.concatenate([[0.0], np.cumsum(v * width)])
    cum_cov = np.concatenate([[0.0], np.cumsum(width)])
    k = np.searchsorted(s, positions, side="right") - 1
    inside = k >= 0
    kk = np.where(inside, k, 0)
    part = np.clip(positions - s[kk], 0, width[kk])
    total = np.where(inside, cum_sum[kk] + v[kk] * part, 0.0)
    cover = np.where(inside, cum_cov[kk] + part, 0.0)
    return total, cover


def _alias(reader, chrom):
    if chrom in reader.sizes:
        return chrom
    key = chrom.lower().removeprefix("chr")
    return next((c for c in reader.sizes if c.lower().removeprefix("chr") == key), chrom)


def window_matrix(path: str | Path, windows, n_bins: int, stat: str = "mean") -> np.ndarray:
    """Rows of binned signal for many windows (heatmap / profile matrices).

    Same semantics as per-window ``pyBigWig.stats(nBins=…)``: ``mean`` is over covered bases and
    0 where a bin has no data; bins beyond a chromosome edge are NaN. Windows may differ in
    length (scale-regions); bin edges are then fractional, exactly like pyBigWig.
    """
    if stat not in ("mean", "max", "sum"):
        raise ValueError("window_matrix supports mean, max or sum")
    reader = SignalReader(path)
    out = np.full((len(windows), n_bins), np.nan)
    try:
        by_chrom: dict[str, list[int]] = {}
        for index, (chrom, _start, _end) in enumerate(windows):
            by_chrom.setdefault(_alias(reader, chrom), []).append(index)
        frac = np.arange(n_bins + 1) / n_bins
        for chrom, rows in by_chrom.items():
            length = reader.sizes.get(chrom)
            if length is None:
                continue
            starts = np.array([windows[i][1] for i in rows], dtype=float)
            ends = np.array([windows[i][2] for i in rows], dtype=float)
            edges = starts[:, None] + (ends - starts)[:, None] * frac[None, :]
            inside = (edges[:, :-1] >= 0) & (edges[:, 1:] <= length)
            clipped = np.clip(edges, 0, length)
            if reader.dense or stat == "max":
                for k, row in enumerate(rows):
                    lo, hi = int(clipped[k, 0]), int(np.ceil(clipped[k, -1]))
                    if hi <= lo:
                        continue
                    block = reader.values(chrom, lo, hi)
                    finite = np.isfinite(block)
                    cut = np.clip(np.round(clipped[k] - lo).astype(np.int64), 0, hi - lo)
                    empty = cut[1:] <= cut[:-1]  # zero-width bins, e.g. beyond a chromosome end
                    starts_ = np.minimum(cut[:-1], hi - lo - 1)  # reduceat needs in-range indices
                    total = np.add.reduceat(np.where(finite, block, 0.0), starts_)
                    cover = np.add.reduceat(finite.astype(float), starts_)
                    total[empty], cover[empty] = 0.0, 0.0
                    if stat == "max":
                        peak = np.maximum.reduceat(np.where(finite, block, -np.inf), starts_)
                        peak[empty] = -np.inf
                        values = np.where(np.isfinite(peak), peak, 0.0)
                    elif stat == "sum":
                        values = total
                    else:
                        values = np.where(cover > 0, total / np.maximum(cover, 1e-12), 0.0)
                    out[row] = np.where(inside[k], values, np.nan)
                continue
            lo, hi = int(max(0, starts.min())), int(min(length, ends.max()))
            s, e, v = reader.intervals(chrom, lo, hi)
            f_sum, f_cov = _integrals(s, e, v, clipped.ravel())
            total = np.diff(f_sum.reshape(clipped.shape), axis=1)
            cover = np.diff(f_cov.reshape(clipped.shape), axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                values = (
                    total
                    if stat == "sum"
                    else np.where(cover > 0, total / np.maximum(cover, 1e-12), 0.0)
                )
            out[rows] = np.where(inside, values, np.nan)
    finally:
        reader.close()
    return out
