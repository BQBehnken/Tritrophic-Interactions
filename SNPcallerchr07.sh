#!/usr/bin/env bash
set -euo pipefail

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025 Brian Behnken


REF="PvulgarisLaborOvalle_670_v1.0.fa"     
THREADS="${THREADS:-32}"                

OUTDIR="Labor Ovalle VCFs"
mkdir -p "$OUTDIR"

# 8 BAMs (deduplicated) in the current directory
BAMS="Pv250.LaborOvalle.dedup.primary.pp.MQ20.bam Pv63H.LaborOvalle.dedup.primary.pp.MQ20.bam Pv4B.LaborOvalle.dedup.primary.pp.MQ20.bam Pv4C.LaborOvalle.dedup.primary.pp.MQ20.bam Pv7D.LaborOvalle.dedup.primary.pp.MQ20.bam Pv7A.LaborOvalle.dedup.primary.pp.MQ20.bam Pv8A.LaborOvalle.dedup.primary.pp.MQ20.bam Pv8C.LaborOvalle.dedup.primary.pp.MQ20.bam"


# Checks
##############################
echo "[*] Tool checks..."
command -v samtools  >/dev/null || { echo "ERROR: samtools not found in PATH"; exit 1; }
command -v bcftools  >/dev/null || { echo "ERROR: bcftools not found in PATH"; exit 1; }
command -v tabix     >/dev/null || { echo "ERROR: tabix not found in PATH"; exit 1; }

echo "[*] Reference checks..."
[ -s "$REF" ] || { echo "ERROR: reference not found: $REF"; exit 1; }
[ -s "${REF}.fai" ] || { echo "[*] Indexing reference..."; samtools faidx "$REF"; }

echo "[*] BAM index checks (create if missing)..."
for b in $BAMS; do
  [ -s "$b" ] || { echo "ERROR: missing BAM: $b"; exit 1; }
  samtools idxstats "$b" >/dev/null 2>&1 || { echo "    Indexing $b ..."; samtools index -@ "$THREADS" "$b"; }
done


# Main loop: Chr01–Chr11
##############################
for i in $(seq 1 11); do
  CHR=$(printf "Chr%02d" "$i")

  # Outputs (per chromosome)
  RAWBCF="$OUTDIR/F4.${CHR}.raw.bcf"
  SNPNORM="$OUTDIR/F4.${CHR}.snps.norm.bcf"
  HCBCF="$OUTDIR/F4.${CHR}.hc.bcf"
  HCVGZ="$OUTDIR/F4_${CHR}.hc.vcf.gz"

  echo "============================================================"
  echo "[*] Processing $CHR"
  echo "============================================================"


  # 1) mpileup + call (biallelic SNPs)
  ##############################
  if [ ! -s "$RAWBCF" ] || [ ! -s "${RAWBCF}.csi" ]; then
    echo "[*] Running mpileup + call on ${CHR} ..."
	bcftools mpileup -Ou -r "$CHR" -f "$REF" -q 30 -Q 30 -d 800 -a DP,AD $BAMS \
	| bcftools call -mv -Ob --threads "$THREADS" -o "$RAWBCF"
	bcftools index "$RAWBCF"
  else
    echo "[*] Skipping mpileup: $RAWBCF exists."
  fi


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


  # 3) Per-genotype DP masking + site-level filters
  ##############################
  N_SAMPLES=$(bcftools query -l "$SNPNORM" | wc -l | awk '{print $1}')
  NS_MIN=$(( (N_SAMPLES + 1) / 2 ))
  echo "[*] Samples detected: $N_SAMPLES (requiring NS>=${NS_MIN})"

  if [ ! -s "$HCBCF" ] || [ ! -s "${HCBCF}.csi" ]; then
    echo "[*] Applying DP-based genotype masking and site filters..."
    bcftools +setGT "$SNPNORM" -- -t q -n . -i 'FMT/DP<15' \
    | bcftools +fill-tags -- -t AC,AN,AF,NS \
    | bcftools view -i "QUAL>30 && AC>0 && NS>=${NS_MIN}" \
      -Ob --threads "$THREADS" -o "$HCBCF"
    bcftools index "$HCBCF"
  else
    echo "[*] Skipping filtering: $HCBCF exists."
  fi


  # 4) Export final VCF.gz (per chromosome)
  ##############################
  if [ ! -s "$HCVGZ.tbi" ]; then
    echo "[*] Writing final VCF for R..."
    bcftools view "$HCBCF" -Oz --threads "$THREADS" -o "$HCVGZ"
    tabix -p vcf "$HCVGZ"
  else
    echo "[*] Skipping VCF export: $HCVGZ exists."
  fi


  # 5) Summaries
  ##############################
  echo "[*] Done for $CHR. Outputs:"
  ls -lh "$RAWBCF" "${RAWBCF}.csi" \
         "$SNPNORM" "${SNPNORM}.csi" \
         "$HCBCF" "${HCBCF}.csi" \
         "$HCVGZ" "${HCVGZ}.tbi"

  echo "[*] Sanity checks for $CHR:"
  echo "  - Chromosomes in final VCF:"
  bcftools query -f '%CHROM\n' "$HCVGZ" | sort -u
  echo "  - Samples in final VCF:"
  bcftools query -l "$HCVGZ"

done

echo "[*] All requested chromosomes completed."
