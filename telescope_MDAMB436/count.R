library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# List your Telescope files
files <- c("GSK_rep1-telescope_report.tsv",
           "GSK_rep2-telescope_report.tsv",
           "LLY_rep1-telescope_report.tsv",
           "LLY_rep2-telescope_report.tsv")

# Function to extract transcript and count columns, skipping first row and using second as header
extract_counts <- function(file) {
  df <- read_tsv(file, skip = 1, show_col_types = FALSE)
  df %>% select(transcript = 1, count = 3)
}

# Extract counts from each file and name the count column uniquely
counts_list <- map2(files, files, ~{
  df <- extract_counts(.x)
  colnames(df)[2] <- gsub("-telescope_report.tsv", "", basename(.y)) # sample name
  return(df)
})

# Merge all counts by transcript using full_join
merged_counts <- reduce(counts_list, full_join, by = "transcript")

# Replace NA counts with zero
merged_counts[is.na(merged_counts)] <- 0

# Save merged counts to a single tab-delimited text file
write_tsv(merged_counts, "merged_transcript_counts_all_replicates.txt")

cat("All replicates' counts merged and saved to 'merged_transcript_counts_all_replicates.txt'.\n")
