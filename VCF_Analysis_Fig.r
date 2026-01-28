# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025 Brian Behnken

cat("\014"); rm(list=ls())
gc() # working with these files uses a lot of scratch disk
pacman::p_load(vcfR, tidyverse)

capabilities("cairo")
# -----------------------
# Inputs
# -----------------------

# per-chromosome VCFs
vcf_files <- list.files( pattern = "^F4_Chr\\d{2}\\.hc\\.vcf\\.gz$", full.names = TRUE)

# Parents
parentA <- "Pv250.LaborOvalle.dedup.primary.pp.MQ20.bam"
parentB <- "Pv63H.LaborOvalle.dedup.primary.pp.MQ20.bam"

# Samples to show in parents, then family 4, 7, and 8, showing inr-1/inr-1 first
sample_order <- c(parentA, parentB, "Pv4B.LaborOvalle.dedup.primary.pp.MQ20.bam", "Pv4C.LaborOvalle.dedup.primary.pp.MQ20.bam", "Pv7D.LaborOvalle.dedup.primary.pp.MQ20.bam", "Pv7A.LaborOvalle.dedup.primary.pp.MQ20.bam", "Pv8A.LaborOvalle.dedup.primary.pp.MQ20.bam", "Pv8C.LaborOvalle.dedup.primary.pp.MQ20.bam")

# -----------------------
# Helper functions
# -----------------------
allele_from_gt <- function(x){
  ifelse(is.na(x) | x=="./.", NA,
         ifelse(x %in% c("0/0","0|0"), "REF",
                ifelse(x %in% c("1/1","1|1"), "ALT",
                       ifelse(x %in% c("0/1","1/0","0|1","1|0"), "HET", NA)))) # sorting informative SNPs
}

read_and_classify <- function(vcf_path, parentA, parentB, samples_keep){
  vcf <- read.vcfR(vcf_path, verbose = FALSE)

  fix <- as_tibble(getFIX(vcf)) |>
    mutate(POS = as.integer(POS),
           CHROM = as.character(CHROM),
           Mb = POS/1e6)

  gt <- extract.gt(vcf, element = "GT")
  samples <- colnames(gt)

  if(!(parentA %in% samples && parentB %in% samples)){
    stop("Parents not found in VCF samples for: ", basename(vcf_path),
         "\nSamples are: ", paste(samples, collapse=", "))
  }

  # Keep only requested samples that actually exist in this VCF
  keep <- intersect(samples_keep, samples)

  # Parents’ alleles at each site
  PA <- allele_from_gt(gt[, parentA])
  PB <- allele_from_gt(gt[, parentB])

  # Informative sites = parents homozygous and different
  informative <- PA %in% c("REF","ALT") & PB %in% c("REF","ALT") & (PA != PB)

  # classify relative to parent-of-origin
  classify_vec <- function(s){
    S <- allele_from_gt(gt[, s])
    dplyr::case_when(
      informative & S == PA ~ "A",      # Pv250-like
      informative & S == PB ~ "B",      # Pv63H-like
      informative & S == "HET" ~ "HET", # heterozygous
      TRUE ~ NA_character_             # not informative or missing
    )
  }

  # Build long table across retained samples
  long <- map_dfr(keep, \(s){
    tibble(sample = s, cls = classify_vec(s)) |>
      bind_cols(fix)
  }) |>
    filter(!is.na(cls)) |>
    mutate(
      sample = factor(sample, levels = samples_keep),
      CHROM  = factor(CHROM, levels = sprintf("Chr%02d", 1:11)),
      cls    = factor(cls, levels = c("A","B","HET"))
    )

  list(long = long, fix = fix, gt = gt, informative = informative)
}

# -----------------------
# Read all chromosomes + combine
# -----------------------
res_list <- lapply(vcf_files, read_and_classify,
                   parentA = parentA, parentB = parentB, samples_keep = sample_order)

all_long <- bind_rows(lapply(res_list, `[[`, "long"))

n_distinct(all_long[, c("CHROM", "POS")]) # find out how many informative sites were used across all 8 NILs for reporting in caption

### Now Onto Plotting: 4 Plots
# Plot A: Large overview of all chromosomes from all parents and NILs to demonstrate introgression was only on chromosome 7
# Plot B: Large overview of only chromosome 7 to demonstrate introgression boundaries
# Plot C: Zoom overview of the 7.8 Mb region to demonstrate that INR is in the correct respective introgression
# Plot D: Read depth plot to demonstrate the 103-bp presence or absence in respective parent or NIL

# -----------------------
# (1) Genome-wide visual comparison; all chromosomes, all parents and NILs
# -----------------------

track_df <- all_long |>
  mutate(x_pos = 0)

# Labor Ovalle INR ID: PvLabOv.07G078300
# Labor Ovalle INR: 7,831,835 - 7,834,672 from Phytozome (https://phytozome-next.jgi.doe.gov/report/transcript/PvulgarisLaborOvalle_v1_1/PvLabOv.07G078300.1)
# Labor Ovalle INR Locus: Mb 7.830999, 7.834826

roi <- c(7.830999, 7.834826)  

# Create annotation data for Chr07 only to mark black line of INR across all chromosome 7s. 

gene_annotation <- tibble(
  CHROM = factor("Chr07", levels = levels(track_df$CHROM)),
  ymin = roi[1],
  ymax = roi[2]
)

############# PLOTTING

p_all <- ggplot(track_df, aes(x = x_pos, y = Mb, color = cls)) +
  ggrastr::rasterise(geom_jitter(shape = 15, size = 0.3, alpha = 0.9, 
                                  width = 0.3, height = 0), dpi = 600) +
  geom_rect(data = gene_annotation,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "yellow", alpha = 0.3, color = "black", linewidth = 0.3,
            inherit.aes = FALSE) +  # Don't inherit color aesthetic
  facet_grid(sample ~ CHROM, scales = "free_y", space = "free_y") +
  scale_color_manual(values = c(A="#595a5b", B="#bf7049", HET="#9fb1d6"),
                     labels = c("Pv250-like","Pv63H-like","HET"),
                     name = NULL) +
  scale_x_continuous(NULL, breaks = NULL, limits = c(-0.6, 0.6)) +
  labs(y = "Position (Mb)", title = "Genome-wide parent-of-origin SNP calls (informative sites only)") +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA), # Ensure the panel itself is clear
    plot.background = element_rect(fill = "transparent", color = NA),  # Ensure the whole plot is clear
    strip.text = element_text(size = 9),
    panel.grid = element_blank(),
    panel.spacing.x = unit(3, "mm")
  )

p_all
ggsave("Introgression_AllChromosomes_F4.pdf", p_all, width = 16, height = 18)

# -----------------------
# (2) Chr07 Full Chromosome View with 7.8 Mb ROI highlight
# -----------------------
roi_chr <- "Chr07"
roi_7p8_mb <- c(7.831835, 7.834672)

# Labor Ovalle INR 7.831835, 7.834672
# INR Locus: 7.830999, 7.834826

df_chr07 <- all_long |>
  filter(CHROM == roi_chr, 
         !is.na(cls),
         !sample %in% c(parentA, parentB)) |>
  mutate(sample = droplevels(sample))  # Remove unused parent levels

############# PLOTTING

p_chr07 <- ggplot(df_chr07, aes(x = sample, y = Mb, color = cls)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = roi_7p8_mb[1], ymax = roi_7p8_mb[2],
           fill = "grey50", alpha = 0.5, color = NA) +
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 1, shape = 15,
                                position = position_jitter(width = 0.35, height = 0)), 
                     dpi = 1200) +
  scale_color_manual(values = c(A="#595a5b", B="#bf7049", HET="#9fb1d6"),
                     labels = c("Pv250-like","Pv63H-like","HET"),
                     name = NULL) +
  labs(x = NULL, y = "Position (Mb)", 
       title = sprintf("%s SNPs by parent-of-origin", roi_chr)) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA), # Ensure the panel itself is clear
    plot.background = element_rect(fill = "transparent", color = NA),  # Ensure the whole plot is clear
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_chr07
ggsave("Introgression_Chr07_Full_F4.pdf", p_chr07, width = 7, height = 8)

# -----------------------
# (3) Zoom to 7.8 Mb region (broad scan window)
# -----------------------
roiscan_mb <- c(6.700000, 8.500000)
inr_roi_mb <- c(7.831835, 7.834672)

# Locus 7830999-7834826

# Define progeny samples (excluding parents)
progeny_samples <- c("Pv4B.LaborOvalle.dedup.primary.pp.MQ20.bam",
                     "Pv4C.LaborOvalle.dedup.primary.pp.MQ20.bam",
                     "Pv7D.LaborOvalle.dedup.primary.pp.MQ20.bam",
                     "Pv7A.LaborOvalle.dedup.primary.pp.MQ20.bam",
                     "Pv8A.LaborOvalle.dedup.primary.pp.MQ20.bam",
                     "Pv8C.LaborOvalle.dedup.primary.pp.MQ20.bam")

df_zoom_broad <- all_long |>
  filter(CHROM == roi_chr, 
         Mb >= roiscan_mb[1], 
         Mb <= roiscan_mb[2],
         !is.na(cls),
         sample %in% progeny_samples) |>
  mutate(sample = factor(sample, levels = progeny_samples))

############# PLOTTING

p_zoom_broad <- ggplot(df_zoom_broad, aes(x = sample, y = Mb, color = cls)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = inr_roi_mb[1], ymax = inr_roi_mb[2],
           fill = "grey50", alpha = 1, color = NA) +
  geom_point(size = 0.6, alpha = 1, shape = 15,
             position = position_jitter(width = 0.35, height = 0)) +
  scale_color_manual(values = c(A="#595a5b", B="#bf7049", HET="#9fb1d6"),
                     labels = c("Pv250-like","Pv63H-like","HET"),
                     name = NULL) +
  labs(x = NULL, y = "Position (Mb)",
       title = sprintf("%s ROI scan: %.2f–%.2f Mb", roi_chr, roiscan_mb[1], roiscan_mb[2])) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA), # Ensure the panel itself is clear
    plot.background = element_rect(fill = "transparent", color = NA),  # Ensure the whole plot is clear
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_zoom_broad
ggsave("Introgression_Chr07_7p4Mb_Broad_F4.pdf", p_zoom_broad, width = 7, height = 8)

# -----------------------
# (4) TSV Depth Plot for INR 103-bp deletion
# -----------------------

# Bash Script for retrieving TSV files for the Depth Plot
# on Labor Ovalle (PhaseolusLaborOvalle 670) PvLabOv.07G078300

#mkdir -p depth_roi
#REGION="Chr07:7830999-7834826"
#for s in Pv4B Pv4C Pv7A Pv7D Pv8A Pv8C Pv250 Pv63H; do
#samtools depth -a -r "$REGION" "${s}.LaborOvalle.dedup.primary.pp.MQ20.bam" > depth_roi/${s}.tsv
#done

pacman::p_load(readr, purrr)

# 1) TSVs are in working directory. 

# 2) Find the files (accept .tsv or .tsv.gz)
depth_files <- list.files(depth_dir, pattern = "\\.tsv(\\.gz)?$", full.names = TRUE)
if (length(depth_files) == 0) stop("No depth files found. Check depth_dir or pattern.")

# 3) Reader that forces 3 columns (no header!) and adds sample + Mb
read_depth <- function(f){
  read_tsv(
    f,
    col_names = c("CHROM","POS","DP"),   # <-- critical: NO HEADER in  files
    col_types = cols(CHROM = col_character(),
                     POS   = col_integer(),
                     DP    = col_double()),
    progress = FALSE
  ) |>
    mutate(
      sample = str_remove(basename(f), "\\.tsv(\\.gz)?$"),
      Mb = POS / 1e6
    )
}

cov <- map_dfr(depth_files, read_depth)

# Quick sanity checks
print(depth_files)
print(head(cov))
cat("Chrom names in files:", paste(unique(cov$CHROM), collapse = ", "), "\n")
cov %>% count(sample) %>% print(n = Inf)

# Phytozome ROI for Labor Ovalle: 7830999-7834826; INR: 7.831835, 7.834672
roi <- c(7.830999, 7.834826)
inrroi <- c(7.831835, 7.834672) # 7.408696, 7.411351 from the G19833 reference file

cov_roi <- cov %>%
  filter(CHROM == "Chr07", Mb >= roi[1], Mb <= roi[2])

# If this is empty, check spelling & case of CHROM and roi units (Mb)
if (nrow(cov_roi) == 0) {
  message("cov_roi is empty. Debug info:")
  print(unique(cov$CHROM))
  print(range(cov$Mb))
  stop("CHROM name or ROI range mismatch. Make sure CHROM=='Chr07' and roi is in Mb.")
}

# Order tracks to match SNP panels from the previous 3
cov_roi$sample <- factor(
  cov_roi$sample,
  levels = c("Pv250","Pv63H","Pv4B","Pv4C","Pv7D","Pv7A","Pv8A","Pv8C")
)

#### INR only

order_samples <- c("Pv250","Pv63H","Pv4B","Pv4C","Pv7D","Pv7A","Pv8A","Pv8C")
# If cov sample names have suffixes (e.g. .dedup.bam), adapt:

cov_roi <- cov %>%                     # cov = full depth table  created earlier
  filter(CHROM == "Chr07", Mb >= inrroi[1], Mb <= inrroi[2]) %>%
  mutate(sample = factor(sample, levels = order_samples)) %>%
  # drop any samples not present (defensive)
  filter(!is.na(sample))

# quick sanity check
if(nrow(cov_roi)==0) stop("cov_roi is empty — check CHROM name, Mb range, or cov data")

# Coverage panel with read depth INR
depthplot <- ggplot(cov_roi, aes(x = Mb, y = DP)) +
  geom_line(size = 0.4) +                     # line per sample
  facet_grid(sample ~ ., scales = "free_y", switch = "y") +
  labs(x = "Position (Mb)", y = "Depth",
       title = sprintf("Chr07 coverage: %.6f–%.6f Mb (INR)", inrroi[1], inrroi[2])) +
  theme_linedraw() +
  theme(
    strip.text.y = element_text(angle = 0),   # nicer sample labels
    panel.grid.minor = element_blank()
  )

# Print it:
print(depthplot)

ggsave("coverage_inr.pdf", width = 12, height = 6, depthplot)


# -----------------------
# Lastly, find introgression boundaries for reporting
# -----------------------

# Function to find genotype blocks for a given sample and chromosome for reporting
find_introgression_blocks <- function(data, sample_name, chrom) {
  
  df <- data |>
    filter(sample == sample_name, CHROM == chrom, !is.na(cls)) |>
    arrange(POS)
  
  if(nrow(df) == 0) return(NULL)
  
  # Run-length encoding to find contiguous blocks
  rle_result <- rle(as.character(df$cls))
  
  # Calculate block boundaries
  block_ends <- cumsum(rle_result$lengths)
  block_starts <- c(1, block_ends[-length(block_ends)] + 1)
  
  blocks <- tibble(
    sample = sample_name,
    chrom = chrom,
    genotype = rle_result$values,
    n_snps = rle_result$lengths,
    start_idx = block_starts,
    end_idx = block_ends,
    start_pos = df$POS[block_starts],
    end_pos = df$POS[block_ends],
    start_Mb = df$Mb[block_starts],
    end_Mb = df$Mb[block_ends],
    length_bp = df$POS[block_ends] - df$POS[block_starts],
    length_Mb = round((df$POS[block_ends] - df$POS[block_starts]) / 1e6, 3)
  )
  
  return(blocks)
}

# Samples with introgressions (columns 3 and 5 = Pv7D and Pv8A)
introgression_samples <- c("Pv7D.LaborOvalle.dedup.primary.pp.MQ20.bam",
                           "Pv8A.LaborOvalle.dedup.primary.pp.MQ20.bam")

# Find blocks for all chromosomes for these samples
all_blocks <- map_dfr(introgression_samples, function(s) {

  map_dfr(sprintf("Chr%02d", 1:11), function(chr) {
    find_introgression_blocks(all_long, s, chr)
  })
})

# Filter to show only B (introgressed) blocks with decent size
# (filtering out tiny blocks that might be noise)
introgression_blocks <- all_blocks |>
  filter(genotype == "B", n_snps >= 10) |>  # At least 10 SNPs to be considered real

  arrange(sample, chrom, start_pos)

# View the results
print(introgression_blocks, n = 50)

# Export to CSV for paper
write_csv(introgression_blocks, "introgression_boundaries.csv")

# Summary: just the key columns for the paper
introgression_summary <- introgression_blocks |>
  select(sample, chrom, start_Mb, end_Mb, length_Mb, n_snps) |>
  mutate(
    start_Mb = round(start_Mb, 3),
    end_Mb = round(end_Mb, 3)
  )

print(introgression_summary)
write_csv(introgression_summary, "introgression_summary.csv")
