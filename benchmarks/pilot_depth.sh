#!/bin/bash
# Step 2 of docs/PROFILESHAPE_POSITIVE_CASE_SELECTION_V2.md: pilot depth from the first N read pairs of one run.
#   bash benchmarks/pilot_depth.sh RUN INDEX OUTDIR [N=2000000] [THREADS=8]
# Streams the first N reads (pairs) from ENA, aligns as the full analysis would (bowtie2 --very-sensitive,
# -X 1000 --no-mixed --no-discordant for pairs; MAPQ >= 10, proper pairs), sorts and indexes.
set -euo pipefail
RUN=$1; INDEX=$2; OUT=$3; N=${4:-2000000}; T=${5:-8}
BIN=~/anaconda3/envs/psbench/bin
mkdir -p "$OUT"
urls=$(/usr/bin/curl -sg "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$RUN&result=read_run&fields=fastq_ftp" | tail -1 | cut -f2 | tr ';' ' ')
i=0
for u in $urls; do
  i=$((i + 1))
  for try in 1 2 3 4; do  # streamed downloads occasionally return nothing; retry until the expected lines arrive
    /usr/bin/curl -sg --retry 3 "https://$u" | gzip -dc 2>/dev/null | head -n $((4 * N)) > "$OUT/$RUN.$i.fq" || true
    lines=$(wc -l < "$OUT/$RUN.$i.fq")
    [ "$lines" -gt 0 ] && [ $((lines % 4)) -eq 0 ] && break
    sleep 5
  done
done
if [ -s "$OUT/$RUN.2.fq" ] && [ "$(wc -l < "$OUT/$RUN.1.fq")" -ne "$(wc -l < "$OUT/$RUN.2.fq")" ]; then
  echo "$RUN: mates differ in length after download" >&2; exit 1
fi
if [ -s "$OUT/$RUN.2.fq" ]; then
  $BIN/bowtie2 -p $T --very-sensitive -X 1000 --no-mixed --no-discordant -x "$INDEX" \
    -1 "$OUT/$RUN.1.fq" -2 "$OUT/$RUN.2.fq" 2> "$OUT/$RUN.bowtie2.log" \
    | $BIN/samtools view -b -q 10 -f 2 - | $BIN/samtools sort -@ 4 -o "$OUT/$RUN.bam" -
else
  $BIN/bowtie2 -p $T --very-sensitive -x "$INDEX" -U "$OUT/$RUN.1.fq" 2> "$OUT/$RUN.bowtie2.log" \
    | $BIN/samtools view -b -q 10 - | $BIN/samtools sort -@ 4 -o "$OUT/$RUN.bam" -
fi
$BIN/samtools index "$OUT/$RUN.bam"
rm -f "$OUT/$RUN".*.fq
echo "$RUN done: $(tail -1 "$OUT/$RUN.bowtie2.log")"
