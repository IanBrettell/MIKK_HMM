# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/F2/hdrr/hmmlearn_true/5000/0.8.csv"

## True
IN = snakemake@input[["hmm"]]
OUT = snakemake@output[[1]]

# Read in data

df = readr::read_csv(IN) %>% 
  tidyr::separate(col = "SAMPLE",
                  into = c("SAMPLE", "PAT_LINE", "MAT_LINE"),
                  sep = "_")

# Split by sample

dat_list = df %>% 
  split(., f = .$SAMPLE)

# Save each df

OUT_DIR = dirname(OUT)

lapply(dat_list, function(DF){
  # Set up path to write to
  SAMPLE = DF$SAMPLE[1]
  PAT = DF$PAT_LINE[1]
  MAT = DF$MAT_LINE[1]
  OUT_PATH = file.path(OUT_DIR,
                       paste(SAMPLE, "_", PAT, "_", MAT, ".csv", sep = ""))
  # Write df to file
  readr::write_csv(DF, OUT_PATH)
})
