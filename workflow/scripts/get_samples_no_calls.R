# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8", 
                pattern = "^7.csv",
                recursive = T, 
                full.names = T)
TARGET_CHROM = "7"

## True
IN = snakemake@input
OUT = snakemake@output[["csv"]]
TARGET_CHROM = snakemake@params[["contig"]]

# Read in files
dat_list = purrr::map(IN, function(PATH){
  df = readr::read_csv(PATH)
  
  # if the tibble is empty (because the sample had no SNPs), then create an empty one with a sample name
  ## get sample name
  if (nrow(df) == 0){
    sample = PATH %>% 
      stringr::str_split(pattern = "/",simplify = T) %>% 
      .[length(.) - 1] %>% 
      stringr::str_split('_', simplify = T) %>% 
      .[1]
    df = tibble::tibble(SAMPLE = sample,
                        CHROM = NA_real_,
                        POS = NA_real_,
                        GENO_AB_REF_ALT = NA_character_)    
  }
  
  # Get sample
  ID = as.character(df %>% 
                      dplyr::distinct(SAMPLE) %>% 
                      dplyr::filter(!is.na(SAMPLE)) %>% 
                      dplyr::pull(SAMPLE))
  # Remove sample column and rename GENO_NT column
  df = df %>% 
    dplyr::rename(GT = GENO_AB_REF_ALT) %>% 
    #    # convert AB to 012
    #    dplyr::mutate(GT = dplyr::case_when(GENO_AB_REF_ALT == "AA" ~ 0,
    #                                        GENO_AB_REF_ALT == "AB" ~ 1,
    #                                        GENO_AB_REF_ALT == "BB" ~ 2)) %>% 
    dplyr::select(CHROM, POS, {{ID}} := GT)
})

# If there are no samples with no calls for this chromosome, write an empty file
if (!any(purrr::map_lgl(dat_list, ~nrow(.x) == 1))){
  file.create(OUT)
} else {
  # Find samples that have only 1 row (indicating that they have no SNP calls for this chromosome)
  MISS_SAMPLES = purrr::keep(dat_list, ~nrow(.x) == 1) %>% 
    purrr::reduce(dplyr::full_join, by = c("CHROM", "POS")) %>% 
    dplyr::select(-c(CHROM, POS)) %>% 
    colnames(.)
  # Write to file
  readr::write_lines(MISS_SAMPLES, OUT)
}


