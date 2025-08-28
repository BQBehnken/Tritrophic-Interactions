#!/usr/bin/env bash
set -euo pipefail

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025 Brian Behnken


##############################
# User settings
##############################
REF="Pvulgaris_442_v2.0.fa"            # path to reference FASTA; This is P. vulgaris G19833
REGION="Chr07"                         # whole chromosome 7
THREADS="${THREADS:-8}"                # allow override: THREADS=4 ./script.sh

# Your 8 BAMs (deduplicated) in the current directory
BAMS="Pv250.dedup.bam Pv63H.dedup.bam Pv4B.dedup.bam Pv4C.dedup.bam Pv7A.dedup.bam Pv7D.dedup.bam Pv8A.dedup.bam Pv8C.dedup.bam"

# Outputs
RAWBCF="NILs.chr07.raw.bcf"
SNPNORM="NILs.chr07.snps.norm.bcf"
HCBCF="NILs.chr07.hc.bcf"
HCVGZ="NILs_chr07.hc.vcf.gz"
ZOOMRGN="Chr07:7407744-7413510"
ZOOMVCF="NILs_focus.vcf.gz"

##############################
# Checks
##############################
echo "[*] Tool checks..."
command -v samtools >/dev/null || { echo "ERROR: samtools not found in PATH"; exit 1; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found in PATH"; exit 1; }

echo "[*] Reference checks..."
[ -s "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
[ -s "${REF}.fai" ] || { echo "[*] Indexing reference..."; samtools faidx "$REF"; }

echo "[*] BAM index checks (create if missing)..."
for b in $BAMS; do
  [ -s "$b" ] || { echo "ERROR: missing BAM: $b"; exit 1; }
  samtools idxstats "$b" >/dev/null 2>&1 || { echo "    Indexing $b ..."; samtools index -@ "$THREADS" "$b"; }
done

##############################
# 1) mpileup + call (biallelic SNPs)
##############################
if [ ! -s "$RAWBCF" ] || [ ! -s "${RAWBCF}.csi" ]; then
  echo "[*] Running mpileup + call on ${REGION} ..."
  # Read-level filters: -q 30 (mapQ), -Q 30 (baseQ), cap depth -d 800
  bcftools mpileup -Ou -r "$REGION" -f "$REF" -q 30 -Q 30 -d 800 \
    -a DP,AD $BAMS \
  | bcftools call -mv -Ob --threads "$THREADS" -o "$RAWBCF"
  bcftools index "$RAWBCF"
else
  echo "[*] Skipping mpileup: $RAWBCF exists."
fi

##############################
# 2) Normalize and keep clean bi-allelic SNPs
##############################
if [ ! -s "$SNPNORM" ] || [ ! -s "${SNPNORM}.csi" ]; then
  echo "[*] Normalizing and keeping bi-allelic SNPs..."
  bcftools view -v snps -m2 -M2 -Ob --threads "$THREADS" "$RAWBCF" \
  | bcftools norm -f "$REF" -Ob --threads "$THREADS" -o "$SNPNORM"
  bcftools index "$SNPNORM"
else
  echo "[*] Skipping normalize: $SNPNORM exists."
fi

##############################
# 3) Per-genotype DP masking + site-level filters
#    (No GQ, no AD gating to avoid version quirks)
##############################
# Determine sample count to set missingness threshold: NS >= ceil(N/2)
N_SAMPLES=$(bcftools query -l "$SNPNORM" | wc -l | awk '{print $1}')
NS_MIN=$(( (N_SAMPLES + 1) / 2 ))
echo "[*] Samples detected: $N_SAMPLES (requiring NS>=${NS_MIN})"

if [ ! -s "$HCBCF" ] || [ ! -s "${HCBCF}.csi" ]; then
  echo "[*] Applying DP-based genotype masking and site filters..."
  bcftools +setGT "$SNPNORM" -- -t q -n . -i 'FMT/DP<15 || FMT/GQ<30' \
  | bcftools +fill-tags -- -t AC,AN,AF,NS,DP \
  | bcftools view -i "QUAL>30 && AC>0 && NS>=${NS_MIN}" \
  -Ob --threads "$THREADS" -o "$HCBCF"
  bcftools index "$HCBCF"
else
  echo "[*] Skipping filtering: $HCBCF exists."
fi

##############################
# 4) Exports: VCF for R + zoom window
##############################
if [ ! -s "$HCVGZ.tbi" ]; then
  echo "[*] Writing final VCF for R..."
  bcftools view "$HCBCF" -Oz --threads "$THREADS" -o "$HCVGZ"
  tabix -p vcf "$HCVGZ"
else
  echo "[*] Skipping VCF export: $HCVGZ exists."
fi

if [ ! -s "$ZOOMVCF.tbi" ]; then
  echo "[*] Extracting zoom window $ZOOMRGN ..."
  bcftools view -r "$ZOOMRGN" -Oz --threads "$THREADS" -o "$ZOOMVCF" "$HCBCF"
  tabix -p vcf "$ZOOMVCF"
else
  echo "[*] Skipping zoom: $ZOOMVCF exists."
fi

##############################
# 5) Summaries
##############################
echo "[*] Done. Outputs:"
ls -lh "$RAWBCF" "${RAWBCF}.csi" \
       "$SNPNORM" "${SNPNORM}.csi" \
       "$HCBCF" "${HCBCF}.csi" \
       "$HCVGZ" "${HCVGZ}.tbi" \
       "$ZOOMVCF" "${ZOOMVCF}.tbi"

echo "[*] Sanity checks:"
echo "  - Chromosomes in final VCF:"
bcftools query -f '%CHROM\n' "$HCVGZ" | sort -u
echo "  - Samples in final VCF:"
bcftools query -l "$HCVGZ"
