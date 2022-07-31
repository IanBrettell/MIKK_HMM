# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/2_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/1_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/3_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/18_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/17_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/20_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/194_50-2_18-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/193_50-2_18-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/196_50-2_18-2/10.csv")

# Read in files

raw = purrr::map_dfr(IN, readr::read_csv)


# Process

df_genos = raw %>% 
  # convert AB to 012
  dplyr::mutate(GT = dplyr::case_when(GENO_AB_REF_ALT == "AA" ~ 0,
                                      GENO_AB_REF_ALT == "AB" ~ 1,
                                      GENO_AB_REF_ALT == "BB" ~ 2)) %>% 
  dplyr::select(SAMPLE, CHROM, POS, GT) %>% 
  # pivot wider to put samples as separate columns
  tidyr::pivot_wider(id_cols = c("CHROM", "POS"),
                     names_from = SAMPLE,
                     values_from = GT) %>% 
  # keep only complete cases
  dplyr::filter(complete.cases(.)) %>% 
  # add SNP column
  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
  # send to rowname
  tibble::column_to_rownames(var = "SNP") %>% 
  # remove CHROM and POS columns
  dplyr::select(-c(CHROM, POS)) %>% 
  # remove all rows where all samples are HOM
  dplyr::filter(!dplyr::if_all(dplyr::everything(),
                               ~ . == 0),
                !dplyr::if_all(dplyr::everything(),
                               ~ . == 2))


x = df_genos %>% 
  # take first 100
  #dplyr::slice_head(n = 1000) %>% 
  # transpose
  t(.) %>% 
  # convert back to data frame
  as.data.frame(.)

# Compute GRM "manually"
# Following guidance here: https://zjuwhw.github.io/2021/08/20/GRM.html

n = dim(x)[1]
m = dim(x)[2]
# For each SNP, sum the ALT alleles and divide by 2n to get the ALT allele frequency
p_hat = apply(x, 2, sum)/(2*n)
# 
w = apply(rbind(x,p_hat), 2, function(x) (x-2*x[length(x)])/sqrt(2*x[length(x)]*(1-x[length(x)])))[1:n,]
w

# Remove all columns with NaN
empties = apply(w, 2, function(COL) all(is.na(COL)))

w_new = w[, !empties]

# Cacluate the GRM
A = w_new %*% t(w_new) / m

# Order
ord = hclust(dist(A, method = "euclidean"), method = "ward.D")$order
labs = hclust(dist(A, method = "euclidean"), method = "ward.D")$labels

# Order matrix
A_ord = A[ord, ord]

# Convert to Df
df_fig = A_ord %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "SAMPLE_1") %>% 
  tidyr::pivot_longer(-c(SAMPLE_1), names_to = "SAMPLE_2", values_to = "VALUE") %>% 
  dplyr::mutate(dplyr::across(c("SAMPLE_1", "SAMPLE_2"),
                              ~factor(., levels = labs)))

# Inspect
df_fig %>% 
  dplyr::arrange(desc(VALUE)) %>% 
  dplyr::filter(!SAMPLE_1 == SAMPLE_2) %>% 
  View()

fig = df_fig %>% 
  ggplot() +
  geom_tile(aes(x = SAMPLE_1, y = SAMPLE_2, fill = VALUE)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(aspect.ratio = 1) +
  theme(axis.text.x = element_text(angle = 90))

fig
