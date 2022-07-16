# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#IN_FILES = list("/hps/nobackup/birney/users/ian/somites/dpABs/F2/hdrr/F1_het_min_DP/628.txt",
#                "/hps/nobackup/birney/users/ian/somites/dpABs/F2/hdrr/F1_het_min_DP/639.txt")
#BIN_LENGTH = 5000
#MAX_READS = "None"

## True
IN_FILES = snakemake@input
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
MAX_READS = snakemake@params[["max_reads"]]
OUT_FILE = snakemake@output[[1]]

# Convert MAX_READS to numeric if not "None"

if (MAX_READS != "None"){
  MAX_READS = as.numeric(MAX_READS)
}

# Name IN_FILES by sample name
names(IN_FILES) = IN_FILES %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".txt")

# Read in files and process

final = purrr::map(IN_FILES, function(FILE) {
  # Read in files
  out = FILE %>%
    readr::read_tsv(.,
                    col_names = c("CHROM", "POS", "F0_1_ALLELE", "F0_1_COUNT", "F0_2_ALLELE", "F0_2_COUNT"),
                    col_types = c("iicici"))
  
  # Filter for rows where both F0_1 and F0_2 allele counts are less than MAX_READS
  if (MAX_READS != "None"){
    out = out %>% 
      dplyr::filter(F0_1_COUNT < MAX_READS & F0_2_COUNT < MAX_READS)
  }
  
  # Continue to process
  out = out %>% 
    # bin
    dplyr::mutate(BIN = floor(POS / BIN_LENGTH)) %>% 
    # group by CHROM and BIN
    dplyr::group_by(CHROM, BIN) %>% 
    # get proportion of F0_2 (i.e. Kaga) so that proportion 0 -> state 0 and proportion 1 -> state 2
    dplyr::summarise(CAB_COUNT = sum(F0_1_COUNT),
                     KAGA_COUNT = sum(F0_2_COUNT),
                     TOTAL_COUNT = CAB_COUNT + KAGA_COUNT,
                     PROP_KAGA = KAGA_COUNT / TOTAL_COUNT) %>% 
    dplyr::ungroup()
  
  return(out)
}) %>% 
  # bind into single DF
  dplyr::bind_rows(.id = "SAMPLE") %>% 
  # replace NAs in PROP_KAGA (from counts of 0 for both Cab and Kaga) with 0.5
  dplyr::mutate(PROP_KAGA = tidyr::replace_na(PROP_KAGA, 0.5))

# Write to file

readr::write_csv(final, OUT_FILE)
