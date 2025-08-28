#!/usr/bin/env bash
set -euo pipefail

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025 Brian Behnken

# Reference to use
REF=Pvulgaris_442_v2.0.fa

# Sample prefixes (omit _R1/_R2 suffixes)
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

for S in "${SAMPLES[@]}"; do
  echo "→ Processing $S"

  # 1) Align → unsorted BAM using Burrows-Wheeler Alignment. The R1/R2 file names were default from the service. 
  bwa mem -t 30 "$REF" \
    "${S}_R1_001.fastq.gz" "${S}_R2_001.fastq.gz" \
  | samtools view -bS - \
    > "${S}.unsorted.bam"

  # 2) Coordinate‐sort
  samtools sort -@30 -o "${S}.sorted.bam" "${S}.unsorted.bam"

  # 3) Name‐sort for fixmate
  samtools sort -n -@30 -o "${S}.namesorted.bam" "${S}.sorted.bam"

  # 4) Fix mate fields
  samtools fixmate -m -@30 "${S}.namesorted.bam" "${S}.fixmate.bam"

  # 5) Coordinate‐sort the fixmated BAM
  samtools sort -@30 -o "${S}.fixmate.sorted.bam" "${S}.fixmate.bam"

  # 6) Mark & remove duplicates
  samtools markdup -r -@30 "${S}.fixmate.sorted.bam" "${S}.dedup.bam"

  # 7) Index deduped BAM
  samtools index "${S}.dedup.bam"

  echo "✔ $S complete"
done
