# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8", pattern = "14.csv",recursive = T, full.names = T)
IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/2_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/1_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/3_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/18_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/17_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/20_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/194_50-2_18-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/193_50-2_18-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/196_50-2_18-2/10.csv")


## True
IN = snakemake@input
OUT_PED = snakemake@output[["ped"]]
OUT_MAP = snakemake@output[["map"]]



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

# Join

df = dat_list %>%
  purrr::reduce(full_join, by=c("CHROM", "POS")) %>% 
  # remove NAs in CHROM (and POS)
  dplyr::filter(!is.na(CHROM)) %>% 
  # order by CHROM and POS
  dplyr::arrange(CHROM, POS) %>% 
  # add SNP column
  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
  # send SNP to rowname
  tibble::column_to_rownames(var = "SNP") %>% 
  # replace NA with 00 as required by Plink https://www.cog-genomics.org/plink/1.9/input#plink_irreg
  dplyr::mutate(dplyr::across(-c(CHROM, POS),
                              ~tidyr::replace_na(., "00"))) %>% 
  # remove CHROM and POS columns
  dplyr::select(-c(CHROM, POS)) %>% 
  # remove all rows where all samples are HOM
  dplyr::filter(!dplyr::if_all(dplyr::everything(),
                               ~ . == "AA"),
                !dplyr::if_all(dplyr::everything(),
                               ~ . == "BB"))
  
# Create .ped file
ped = df %>% 
#  # remove non-sample-genotype columns
#  dplyr::select(-c(CHROM, POS)) %>% 
  # convert to matrix
  as.matrix() %>% 
  # transpose to put SNPs as columns
  t()


# Make map
map = df %>% 
  tibble::rownames_to_column(var = "SNP") %>% 
  tidyr::separate(SNP,
                  into = c("CHROM", "POS"),
                  sep = ":",
                  remove = F) %>% 
  dplyr::select(CHROM, SNP, POS)



## Create .ped file
#
#ped = df %>% 
#  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
#  #dplyr::select(-c(CHROM, POS)) %>% 
#  # pivot into 3 column
#  tidyr::pivot_longer(cols = -c(SNP, CHROM, POS), 
#                      names_to = "SAMPLE", 
#                      values_to = "GT") %>% 
#  # convert SAMPLE to numeric
#  dplyr::mutate(SAMPLE = as.numeric(SAMPLE)) %>% 
#  # order by SAMPLE, CHROM, and POS
#  dplyr::arrange(SAMPLE, CHROM, POS) %>% 
#  # remove CHROM and POS columns
#  dplyr::select(-c(CHROM, POS)) %>% 
#  # replace NA with 0 as required by Plink https://www.cog-genomics.org/plink/1.9/input#plink_irreg
#  dplyr::mutate(GT = tidyr::replace_na(GT, "00")) %>% 
#  # pivot wide into .ped format (samples to rows, SNPs to columns)
#  tidyr::pivot_wider(id_cols = SAMPLE, 
#                     names_from = SNP, 
#                     values_from = GT)
#
## Create .map file
#
#map = tibble::tibble(SNP = colnames(ped)) %>% 
#  dplyr::filter(SNP != "SAMPLE") %>% 
#  tidyr::separate(col = "SNP",into = c("CHROM", "POS"),sep = ":",remove = F) %>% 
#  dplyr::select(CHROM, SNP, POS)

# Write to file

## Need to use write.table for matrixes
## Saves the row names (SAMPLE) to file
write.table(ped, OUT_PED, quote = F, sep = "\t", row.names = T, col.names = F)
readr::write_tsv(map, OUT_MAP, col_names = F)
