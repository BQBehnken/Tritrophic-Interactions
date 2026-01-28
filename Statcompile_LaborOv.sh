#!/usr/bin/env bash
set -euo pipefail


# --- inputs, can edit ---
REF=${REF:-PvulgarisLaborOvalle_670_v1.0.fa}
THREADS=${THREADS:-32}
BAMS=(${BAMS:-Pv250.LaborOvalle.dedup.bam Pv63H.LaborOvalle.dedup.bam Pv4B.LaborOvalle.dedup.bam Pv4C.LaborOvalle.dedup.bam Pv7A.LaborOvalle.dedup.bam Pv7D.LaborOvalle.dedup.bam Pv8A.LaborOvalle.dedup.bam Pv8C.LaborOvalle.dedup.bam})
VCF=${VCF:-"Labor Ovalle VCFs"/F4_genome.hc.vcf.gz}   

# --- checks ---
command -v samtools >/dev/null || { echo "samtools missing"; exit 1; }
[ -x ./mosdepth ] || { echo "mosdepth missing - download from https://github.com/brentp/mosdepth/releases"; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools missing"; exit 1; }
[ -s "${REF}.fai" ] || samtools faidx "$REF"
# total reference length
REFLEN=$(awk '{s+=$2} END{print s}' "${REF}.fai")

# header
OUT=qc_per_sample.tsv
echo -e "sample\tref_len_bp\tmean_read_len_bp\tread_N50_bp\treads_total\treads_mapped\tmap_rate_pct\tmean_depth\tcov_>=1x_pct\tcov_>=4x_pct\tcov_>=10x_pct\tcov_>=20x_pct\tmean_baseQ_phred\tmean_mapQ\tmismatch_rate_pct" > "$OUT"
# helper: compute mean MAPQ from histogram in samtools stats
mean_mapq() { awk '$1=="MAPQ"{sum+=$2*$3; n+=$3} END{if(n>0) printf "%.2f",sum/n; else print "NA"}'; }
# helper: compute read N50 from RL histogram
read_n50() { awk '
  $1=="RL"{len=$2; n=$3; tot+=len*n; c[len]+=n}
  END{
    target=tot/2; run=0
    n50="NA"
    PROCINFO["sorted_in"]="@ind_num_desc"
    for (L in c){
      run += L*c[L]
      if (run>=target){ n50=L; break }
    }
    print n50
  }'
}

# loop samples
for BAM in "${BAMS[@]}"; do
  [ -s "$BAM" ] || { echo "WARN: missing $BAM – skipping"; continue; }
  SAMPLE=$(basename "$BAM" .bam)
  # samtools stats (one pass)
  samtools stats -@ "$THREADS" "$BAM" > "${SAMPLE}.samstats.txt"
  MEAN_RLEN=$(awk -F'\t' '$2=="average length:"{print $3}' "${SAMPLE}.samstats.txt")
  TOTAL_READS=$(awk -F'\t' '$2=="raw total sequences:"{print $3}' "${SAMPLE}.samstats.txt")
  MAPPED_READS=$(awk -F'\t' '$2=="reads mapped:"{print $3}' "${SAMPLE}.samstats.txt")
  MAP_RATE=$(awk -F'\t' '$2=="reads mapped:"{printf "%.2f", ($3/'"$TOTAL_READS"')*100}' "${SAMPLE}.samstats.txt")
  MEAN_BASEQ=$(awk -F'\t' '$2=="average quality:"{print $3}' "${SAMPLE}.samstats.txt")
  ERR_RATE=$(awk -F'\t' '$2=="error rate:"{print $3*100}' "${SAMPLE}.samstats.txt" 2>/dev/null || echo NA)
  MEAN_MAPQ=$(mean_mapq < "${SAMPLE}.samstats.txt")
  READ_N50=$(read_n50 < "${SAMPLE}.samstats.txt")
  # mosdepth (fast) – thresholds give breadth ≥1/4/10/20×
  ./mosdepth -t "$THREADS" --fast-mode --by 1000 --thresholds 1,4,10,20 "$SAMPLE" "$BAM" 1>/dev/null
  # mean depth from summary:
  MEAN_DP=$(awk 'NR>1{sum+=$4; n++} END{if(n>0) printf "%.2f",sum/n; else print "NA"}' "${SAMPLE}.mosdepth.summary.txt")
  # breadth from thresholds bed.gz (columns: chrom start end >=1 >=4 >=10 >=20)
  read GE1 GE4 GE10 GE20 < <(zcat "${SAMPLE}.thresholds.bed.gz" | \
    awk -v L="$REFLEN" '{ge1+=$5; ge4+=$6; ge10+=$7; ge20+=$8} END{printf "%.2f %.2f %.2f %.2f\n", (ge1/L)*100,(ge4/L)*100,(ge10/L)*100,(ge20/L)*100}')
  echo -e "${SAMPLE}\t${REFLEN}\t${MEAN_RLEN}\t${READ_N50}\t${TOTAL_READS}\t${MAPPED_READS}\t${MAP_RATE}\t${MEAN_DP}\t${GE1}\t${GE4}\t${GE10}\t${GE20}\t${MEAN_BASEQ}\t${MEAN_MAPQ}\t${ERR_RATE}" >> "$OUT"
done

# VCF-level variant stats (SNP totals, Ti/Tv, per-chrom counts)
if [ -s "$VCF" ]; then
  bcftools stats -s - "$VCF" > variant.stats.txt
  echo "[*] Wrote: $OUT  and variant.stats.txt"
else
  echo "[*] Wrote: $OUT  (VCF not found, skipping variant stats)"
fi
