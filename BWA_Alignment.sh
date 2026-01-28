#!/usr/bin/env bash
set -euo pipefail

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025 Brian Behnken

############################################
# Whole-genome alignment to Labor Ovalle
# - Uses full reference FASTA (all chromosomes/contigs)
# - Indexes (samtools faidx + bwa index)
# - Aligns samples, produces:
#     1) dedup BAM (all alignments kept)
#     2) more stringent BAM (primary, properly paired, MAPQ>=20)
############################################

THREADS="${THREADS:-32}"

# Full reference FASTA
REF_FA="PvulgarisLaborOvalle_670_v1.0.fa"

# Semi-stringent filter settings
MIN_MAPQ="${MIN_MAPQ:-20}"        # default
REQ_FLAGS="${REQ_FLAGS:-2}"       # require properly paired
EXCL_FLAGS="${EXCL_FLAGS:-2304}"  # 256 (secondary) + 2048 (supplementary)

# Pv250 = PI 311785
# Pv63h = W6 13807

SAMPLES=(
  Pv250
  Pv63H
  Pv4B
  Pv4C
  Pv7A
  Pv7D
  Pv8A
  Pv8C
)

# ---- checks ----
for exe in samtools bwa; do
  command -v "$exe" >/dev/null || { echo "ERROR: $exe not found"; exit 1; }
done

[[ -s "$REF_FA" ]] || { echo "ERROR: missing reference: $REF_FA"; exit 1; }

echo "[*] Indexing reference (samtools faidx)..."
[[ -s "${REF_FA}.fai" ]] || samtools faidx "$REF_FA"

echo "[*] Building BWA index for full reference (one-time)..."
if [[ -s "${REF_FA}.bwt" ]]; then
  echo "    BWA index exists; skipping."
else
  bwa index "$REF_FA"
fi

echo "[*] Aligning samples to full reference..."
for S in "${SAMPLES[@]}"; do
  echo "→ Processing $S"

  R1="${S}_R1_001.fastq.gz"
  R2="${S}_R2_001.fastq.gz"
  [[ -s "$R1" ]] || { echo "ERROR: missing $R1"; exit 1; }
  [[ -s "$R2" ]] || { echo "ERROR: missing $R2"; exit 1; }

  # Outputs
  DEDUP_BAM="${S}.LaborOvalle.dedup.bam"
  SEMI_BAM="${S}.LaborOvalle.dedup.primary.pp.MQ${MIN_MAPQ}.bam"

  # 1) Align -> BAM
  bwa mem -t "$THREADS" "$REF_FA" "$R1" "$R2" \
    | samtools view -b -o "${S}.unsorted.bam" -

  # 2) Coordinate sort
  samtools sort -@ "$THREADS" -o "${S}.sorted.bam" "${S}.unsorted.bam"

  # 3) Name sort for fixmate
  samtools sort -n -@ "$THREADS" -o "${S}.namesorted.bam" "${S}.sorted.bam"

  # 4) Fix mates
  samtools fixmate -m -@ "$THREADS" "${S}.namesorted.bam" "${S}.fixmate.bam"

  # 5) Coordinate sort
  samtools sort -@ "$THREADS" -o "${S}.fixmate.sorted.bam" "${S}.fixmate.bam"

  # 6) Mark/remove duplicates
  samtools markdup -r -@ "$THREADS" "${S}.fixmate.sorted.bam" "$DEDUP_BAM"

  # 7) Index dedup BAM
  samtools index -@ "$THREADS" "$DEDUP_BAM"

  # 8) more stringent filtered BAM:
  #    -f 2    require properly paired
  #    -F 2304 drop secondary+supplementary
  #    -q MIN  keep MAPQ >= MIN_MAPQ
  samtools view -@ "$THREADS" -b \
    -f "$REQ_FLAGS" -F "$EXCL_FLAGS" -q "$MIN_MAPQ" \
    "$DEDUP_BAM" > "$SEMI_BAM"
  samtools index -@ "$THREADS" "$SEMI_BAM"

  # Quick sanity metrics
  mapped_dedup=$(samtools view -c -F 4 "$DEDUP_BAM")
  mapped_semi=$(samtools view -c -F 4 "$SEMI_BAM")
  echo "    mapped reads (dedup): $mapped_dedup"
  echo "    mapped reads (semi):  $mapped_semi"

  # Optional cleanup
rm -f "${S}.unsorted.bam" "${S}.sorted.bam" "${S}.namesorted.bam" "${S}.fixmate.bam" "${S}.fixmate.sorted.bam"
done

echo "[*] Done."
echo "    Dedup BAMs: *.LaborOvalle.dedup.bam"
echo "    BAMs: *.LaborOvalle.dedup.primary.pp.MQ${MIN_MAPQ}.bam"
