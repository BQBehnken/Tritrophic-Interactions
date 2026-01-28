

pacman::p_load(tidyverse, openxlsx)

# --- Configuration ---
# Set WD
stats_dir <- "."

file_pattern <- "\\.LaborOvalle\\.dedup\\.samstats\\.txt$"
# Output file
output_file <- "sequencing_summary.xlsx"

# --- Functions ---
parse_samstats <- function(filepath) {
  lines <- readLines(filepath)
  
  # Extract SN (Summary Numbers) lines
  sn_lines <- lines[grepl("^SN\t", lines)]
  
  # Helper to extract value by key
  get_sn <- function(key) {
  line <- sn_lines[grepl(paste0("^SN\t", key, ":"), sn_lines)]
  if (length(line) == 0) return(NA_real_)
  as.numeric(str_extract(line, "\\d+\\.?\\d*([eE][+-]?\\d+)?"))
}
  
  tibble(
    raw_total_sequences = get_sn("raw total sequences"),
    reads_mapped = get_sn("reads mapped"),
    reads_properly_paired = get_sn("reads properly paired"),
    average_length = get_sn("average length"),
    average_quality = get_sn("average quality"),
    error_rate = get_sn("error rate"),
    insert_size_average = get_sn("insert size average")
  )
}

# --- Main ---
# Find all samstats files
stats_files <- list.files(stats_dir, pattern = file_pattern, full.names = TRUE)

if (length(stats_files) == 0) {
  stop("No samstats files found matching pattern: ", file_pattern)
}

cat("Found", length(stats_files), "samstats files\n")

# Parse all files
results <- map_dfr(stats_files, function(f) {
  sample <- basename(f) |> 
    str_remove(file_pattern) |>
    str_remove("_LaborOvalle_dedup$|_dedup$|\\.LaborOvalle\\.dedup$")
  
  parse_samstats(f) |>
    mutate(sample = sample, .before = 1)
})

# Calculate derived metrics for reporting in M&M; figure caption
summary_df <- results |>
  mutate(
    paired_end_fragments = raw_total_sequences / 2,
    total_gb = (raw_total_sequences * average_length) / 1e9,
    mapping_rate_pct = (reads_mapped / raw_total_sequences) * 100,
    proper_pair_rate_pct = (reads_properly_paired / raw_total_sequences) * 100
  ) |>
  select(
    sample,
    raw_total_sequences,
    reads_mapped,
    paired_end_fragments,
    total_gb,
    mapping_rate_pct,
    proper_pair_rate_pct,
    average_length,
    average_quality,
    insert_size_average,
    error_rate
  )

# Calculate totals
totals <- summary_df |>
  summarise(
    sample = "TOTAL",
    raw_total_sequences = sum(raw_total_sequences),
    reads_mapped = sum(reads_mapped),
    paired_end_fragments = sum(paired_end_fragments),
    total_gb = sum(total_gb),
    mapping_rate_pct = (sum(reads_mapped) / sum(raw_total_sequences)) * 100,
    proper_pair_rate_pct = NA_real_,
    average_length = mean(average_length),
    average_quality = mean(average_quality),
    insert_size_average = mean(insert_size_average),
    error_rate = mean(error_rate)
  )

# Combine
final_df <- bind_rows(summary_df, totals)

# --- Export to xlsx ---
wb <- createWorkbook()
addWorksheet(wb, "Sequencing Summary")

# Write data
writeData(wb, 1, final_df, startRow = 1, startCol = 1, headerStyle = createStyle(textDecoration = "bold"))

# Format columns
# Number formats
addStyle(wb, 1, style = createStyle(numFmt = "#,##0"), rows = 2:(nrow(final_df)+1), cols = 2:4, gridExpand = TRUE)
addStyle(wb, 1, style = createStyle(numFmt = "0.00"), rows = 2:(nrow(final_df)+1), cols = 5:7, gridExpand = TRUE)
addStyle(wb, 1, style = createStyle(numFmt = "0.0"), rows = 2:(nrow(final_df)+1), cols = 8:10, gridExpand = TRUE)
addStyle(wb, 1, style = createStyle(numFmt = "0.00E+00"), rows = 2:(nrow(final_df)+1), cols = 11, gridExpand = TRUE)

# Bold the totals row
addStyle(wb, 1, style = createStyle(textDecoration = "bold"), rows = nrow(final_df)+1, cols = 1:11, gridExpand = TRUE)

# Auto-width columns
setColWidths(wb, 1, cols = 1:11, widths = "auto")

# Save
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n=== Summary ===\n")
cat("Total paired-end fragments:", format(totals$paired_end_fragments, big.mark = ","), "\n")
cat("Total data:", round(totals$total_gb, 2), "Gb\n")
cat("Overall mapping rate:", round(totals$mapping_rate_pct, 2), "%\n")
cat("\nExported to:", output_file, "\n")

# In excel format in case I need to upload directly for supplemental file
