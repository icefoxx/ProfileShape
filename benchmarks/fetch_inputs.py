"""Manifest of benchmark inputs (source URL, size, sha256) and re-download of the GEO / UCSC files.

    python benchmarks/fetch_inputs.py manifest inputs docs/INPUT_MANIFEST.tsv   # record what is on disk
    python benchmarks/fetch_inputs.py fetch    inputs docs/INPUT_MANIFEST.tsv   # download missing files, check sha256

GEO supplementary files (bigWigs) and UCSC annotations are re-downloaded from their URLs; BAMs and caches are
derived (alignment and build commands in docs/DATA.md) and carry a checksum only. Downloads use the system
curl (see CLAUDE.md: TLS-intercepting gateway).
"""

from __future__ import annotations

import hashlib
import re
import subprocess
import sys
from pathlib import Path

UCSC = {"hg38_ncbiRefSeqSelect.txt.gz": "hg38", "mm10_ncbiRefSeqSelect.txt.gz": "mm10",
        "mm39_ncbiRefSeqSelect.txt.gz": "mm39", "sacCer3_ncbiRefSeq.txt.gz": "sacCer3"}
SKIP = re.compile(r"\.(bai|log)$|/caches/|/workspace/out_|\.DS_Store")


def source(path: Path) -> str:
    name = path.name
    m = re.match(r"(GSM\d+)_", name)
    if m and name.lower().endswith((".bw", ".bigwig")):
        gsm = m.group(1)
        return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"
    if name in UCSC:
        table = name.split("_", 1)[1]
        return f"https://hgdownload.soe.ucsc.edu/goldenPath/{UCSC[name]}/database/{table}"
    return "derived"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest(inputs: Path, out: Path) -> None:
    rows = ["path\tbytes\tsha256\tsource"]
    for p in sorted(x for x in inputs.rglob("*") if x.is_file() and not SKIP.search(str(x))):
        rows.append(f"{p.relative_to(inputs)}\t{p.stat().st_size}\t{sha256(p)}\t{source(p)}")
    out.write_text("\n".join(rows) + "\n")
    print(f"{len(rows) - 1} files → {out}")


def fetch(inputs: Path, table: Path) -> int:
    bad = 0
    for line in table.read_text().splitlines()[1:]:
        rel, size, digest, url = line.split("\t")
        p = inputs / rel
        if not p.exists():
            if url == "derived":
                print(f"missing (derived; rebuild per docs/DATA.md): {rel}")
                bad += 1
                continue
            p.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(["/usr/bin/curl", "-sfL", "--retry", "3", "-o", str(p), url], check=True)
        if sha256(p) != digest:
            print(f"checksum mismatch: {rel}")
            bad += 1
    print("all files present and verified" if not bad else f"{bad} problem(s)")
    return int(bad > 0)


if __name__ == "__main__":
    cmd, inputs, table = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    sys.exit(manifest(inputs, table) if cmd == "manifest" else fetch(inputs, table))
