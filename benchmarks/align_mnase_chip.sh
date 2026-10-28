#!/usr/bin/env bash
# Align the GSE218269 Pol II MNase-ChIP-seq runs (protocol amendment B2): bowtie2 --very-sensitive, concordant
# pairs only, MAPQ >= 10 and proper pairs kept, no duplicate removal (MNase cuts give identical fragments).
#   bash benchmarks/align_mnase_chip.sh FASTQ_DIR INDEX_PREFIX OUT_DIR [THREADS]
set -euo pipefail
FQ=$1; IDX=$2; OUT=$3; T=${4:-10}
BIN=$HOME/anaconda3/envs/psbench/bin
mkdir -p "$OUT"
# run:sample pairs (plain list: macOS ships bash 3.2, without associative arrays)
for pair in SRR22327380:DMSO_1 SRR22327379:DMSO_2 SRR22327378:DMSO_3 \
            SRR22327377:CDK7i_1 SRR22327376:CDK7i_2 SRR22327375:CDK7i_3; do
  run=${pair%%:*}; s=${pair##*:}
  [ -s "$OUT/$s.bam.bai" ] && continue
  for r in 1 2; do  # wait for a download still in progress: both files must pass the gzip integrity check
    until gzip -t "$FQ/${run}_$r.fastq.gz" 2>/dev/null; do sleep 60; done
  done
  "$BIN/bowtie2" --very-sensitive -X 1000 --no-mixed --no-discordant -p "$T" -x "$IDX" \
      -1 "$FQ/${run}_1.fastq.gz" -2 "$FQ/${run}_2.fastq.gz" 2> "$OUT/$s.bowtie2.log" \
    | "$BIN/samtools" view -u -q 10 -f 2 - \
    | "$BIN/samtools" sort -@ 4 -m 1G -o "$OUT/$s.bam" -
  "$BIN/samtools" index "$OUT/$s.bam"
  "$BIN/samtools" markdup -s -n "$OUT/$s.bam" /dev/null 2> "$OUT/$s.dupstats.log" || true  # duplicate share only
  echo "$s done: $(tail -1 "$OUT/$s.bowtie2.log")"
done
