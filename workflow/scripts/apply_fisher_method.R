# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.05/15/5000/0.8/dge/invnorm/novel_object/1/time-quadrant.mlma",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.05/15/5000/0.8/dge/invnorm/novel_object/2/time-quadrant.mlma",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.05/15/5000/0.8/dge/invnorm/novel_object/3/time-quadrant.mlma")

## True
IN = snakemake@input
ADJ_METHOD = snakemake@params[["adj_method"]]
OUT = snakemake@output[["csv"]]

# Read in files

## Name files as states
names(IN) = IN %>% 
  unlist() %>% 
  stringr::str_split(., "/", simplify = T) %>% 
  .[,20]

## Read into single DF
df = purrr::map_dfr(IN, ~readr::read_tsv(.x,
                                         col_types = c("iciccdddd")),
                    .id = "STATE")

# Create correlation matrix for the different states between each set of SNPs
cor_mat = df %>% 
  #dplyr::filter(SNP %in% c("1:8664", "1:8720", "1:13479", "1:13523", "1:19602")) %>% 
  dplyr::select(STATE, SNP, b) %>% 
  tidyr::pivot_wider(id_cols = SNP,
                     names_from = STATE,
                     values_from = b) %>% 
  tibble::column_to_rownames(var = "SNP") %>% 
  cor(.)

## Group by SNP and apply Fisher method

out = df %>% 
  #dplyr::filter(SNP %in% c("1:8664", "1:8720", "1:13479", "1:13523", "1:19602")) %>% 
  dplyr::group_by(SNP) %>% 
  dplyr::mutate(FISHER_P = poolr::fisher(p,
                                         adjust = ADJ_METHOD,
                                         R = cor_mat)$p) %>% 
  dplyr::ungroup() %>% 
  # get distinct SNPs (each SNP should have the same FISHER_P)
  dplyr::distinct(Chr, SNP, bp, FISHER_P)


# Write to file

readr::write_csv(out, OUT)