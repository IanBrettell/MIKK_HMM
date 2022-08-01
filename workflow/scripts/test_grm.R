# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
### Select samples
IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/2_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/1_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/3_38-2_21-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/18_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/17_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/20_8-2_40-1/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/194_50-2_18-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/193_50-2_18-2/10.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8/196_50-2_18-2/10.csv")
IN = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos_AB_chr/hdrr/5000/0.8",pattern = "10.csv", recursive = T, full.names = T)
# Try with all SNPs
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8"

F2_SAMPLES = "/hps/software/users/birney/ian/repos/MIKK_HMM/config/F2_samples_converted.csv"


# Read bed
fam = genio::read_fam(IN_PREF)
bim = genio::read_bim(IN_PREF)
df_genos = genio::read_bed(IN_PREF,
                           names_loci = bim$id,
                           names_ind = fam$id) %>% 
  as.data.frame()

# Read in F2 samples

f2 = readr::read_csv(F2_SAMPLES,
                     col_types = c("cccc")) %>% 
  dplyr::select(SAMPLE = finclip_id,
                PAT = pat_line,
                MAT = mat_line) %>% 
  dplyr::mutate(PAT_MAT = paste(PAT, "x", MAT, sep = ""))

## Read in files
#
#raw = purrr::map_dfr(IN, readr::read_csv)
#
## How many unique SNPs before removing NAs?
#TOT_SNPS = unique(paste(raw$CHROM, raw$POS, sep = ":"))
#
## Process
#
#df_genos = raw %>% 
#  # convert AB to 012
#  dplyr::mutate(GT = dplyr::case_when(GENO_AB_REF_ALT == "AA" ~ 0,
#                                      GENO_AB_REF_ALT == "AB" ~ 1,
#                                      GENO_AB_REF_ALT == "BB" ~ 2)) %>% 
#  dplyr::select(SAMPLE, CHROM, POS, GT) %>% 
#  # pivot wider to put samples as separate columns
#  tidyr::pivot_wider(id_cols = c("CHROM", "POS"),
#                     names_from = SAMPLE,
#                     values_from = GT) %>% 
#  # keep only complete cases
#  dplyr::filter(complete.cases(.)) %>% 
#  # add SNP column
#  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
#  # send to rowname
#  tibble::column_to_rownames(var = "SNP") %>% 
#  # remove CHROM and POS columns
#  dplyr::select(-c(CHROM, POS)) %>% 
#  # remove all rows where all samples are HOM
#  dplyr::filter(!dplyr::if_all(dplyr::everything(),
#                               ~ . == 0),
#                !dplyr::if_all(dplyr::everything(),
#                               ~ . == 2))

x = df_genos %>% 
  # take first 100
  #dplyr::slice_head(n = 1000) %>% 
  # transpose
  t(.) #%>% 
  # convert back to data frame
  #as.data.frame(.)

# Compute GRM "manually"
# Following guidance here: https://zjuwhw.github.io/2021/08/20/GRM.html

n = dim(x)[1]
m = dim(x)[2]
# For each SNP, sum the ALT alleles and divide by 2n to get the ALT allele frequency
p_hat = apply(x, 2, sum)/(2*n)
# 
w = apply(rbind(x,p_hat), 2, function(x) (x-2*x[length(x)])/sqrt(2*x[length(x)]*(1-x[length(x)])))[1:n,]
#w

# Remove all columns with NaN
#empties = apply(w, 2, function(COL) all(is.na(COL)))

#w_new = w[, !empties]

# Caclulate the GRM
A = w %*% t(w) / m

# Order
## By sample to compare with manual GRM
ord = order(match(rownames(A), f2$SAMPLE))
## By cluster
ord = hclust(dist(A, method = "euclidean"), method = "ward.D")$order
labs = rownames(A)[ord]
# Get labs with cross
labs_x = tibble(SAMPLE = labs) %>% 
  dplyr::left_join(f2 %>% 
                     dplyr::select(SAMPLE, PAT_MAT),
                   by = "SAMPLE") %>% 
  # combine
  dplyr::mutate(S_X = paste(SAMPLE, PAT_MAT, sep = "_"))

# Order matrix
A_ord = A[ord, ord]

# Convert to DF
df_fig = A_ord %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "SAMPLE_1") %>% 
  tidyr::pivot_longer(-c(SAMPLE_1), names_to = "SAMPLE_2", values_to = "VALUE") %>% 
  # bind with f2
  dplyr::left_join(f2 %>% 
                     dplyr::select(SAMPLE_1 = SAMPLE,
                                   S1_PAR = PAT_MAT),
                   by = "SAMPLE_1") %>% 
  dplyr::left_join(f2 %>% 
                     dplyr::select(SAMPLE_2 = SAMPLE,
                                   S2_PAR = PAT_MAT),
                   by = "SAMPLE_2") %>% 
  dplyr::mutate(S1_X = paste(SAMPLE_1, S1_PAR, sep = "_"),
                S2_X = paste(SAMPLE_2, S2_PAR, sep = "_")) %>% 
  dplyr::mutate(dplyr::across(c("S1_X", "S2_X"),
                              ~factor(., levels = labs_x$S_X)))

# Inspect
df_fig %>% 
  dplyr::arrange(desc(VALUE)) %>% 
  dplyr::filter(!SAMPLE_1 == SAMPLE_2) %>% 
  View()

fig = df_fig %>% 
  ggplot() +
  geom_tile(aes(x = S1_X, y = S2_X, fill = VALUE)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(aspect.ratio = 1) +
  theme(axis.text.x = element_text(angle = 90))

fig

ggsave(OUT,
       fig,
       device = "png",
       width = 30,
       height = 30,
       units = "in",
       dpi = 400)


# Check if ID files are all the same

IDS = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/peds_contigs_sample_ids/F2/hdrr/5000/0.8", full.names = T)

ids_list = purrr::map_dfc(IDS, ~read_csv(.x, col_names = "SAMPLE"))
