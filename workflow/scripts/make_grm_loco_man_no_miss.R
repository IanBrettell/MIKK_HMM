# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)
library(genio)

# Set variables

## Debug
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds_all/F2/hdrr/5000/0.8/22"
F2_SAMPLES = "/hps/software/users/birney/ian/repos/MIKK_HMM/config/F2_samples_converted.csv"

## True
IN_PREF = snakemake@params[["in_pref"]]
F2_SAMPLES = snakemake@input[["F2_samples"]]
OUT_PREF = snakemake@params[["out_pref"]]
OUT_PNG = snakemake@output[["png"]]
OUT_PDF = snakemake@output[["pdf"]]

# Read in files

## Read bed
fam = genio::read_fam(IN_PREF)
bim = genio::read_bim(IN_PREF)
df_genos = genio::read_bed(IN_PREF,
                           names_loci = bim$id,
                           names_ind = fam$id) %>% 
  as.data.frame()

## Read in F2 samples

f2 = readr::read_csv(F2_SAMPLES,
                     col_types = c("cccc")) %>% 
  dplyr::select(SAMPLE = finclip_id,
                PAT = pat_line,
                MAT = mat_line) %>% 
  dplyr::mutate(PAT_MAT = paste(PAT, "x", MAT, sep = ""))


# Convert genos to matrix

x = df_genos %>% 
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
# Standardise matrix by allele frequency
w = apply(rbind(x,p_hat), 2, function(x) (x-2*x[length(x)])/sqrt(2*x[length(x)]*(1-x[length(x)])))[1:n,]

# Cacluate the GRM
A = w %*% t(w) / m

# Write to file
genio::write_grm(OUT_PREF,
                 kinship = A,
                 fam = fam %>% 
                   dplyr::select(fam, id))

#########################
# Plot GRM
#########################

# Order
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

## Inspect
#df_fig %>% 
#  dplyr::arrange(desc(VALUE)) %>% 
#  dplyr::filter(!SAMPLE_1 == SAMPLE_2) %>% 
#  View()

fig = df_fig %>% 
  ggplot() +
  geom_tile(aes(x = S1_X, y = S2_X, fill = VALUE)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(aspect.ratio = 1) +
  theme(axis.text.x = element_text(angle = 90))

fig

ggsave(OUT_PNG,
       fig,
       device = "png",
       width = 30,
       height = 30,
       units = "in",
       dpi = 400)

ggsave(OUT_PDF,
       fig,
       device = "pdf",
       width = 30,
       height = 30,
       units = "in",
       dpi = 400)

