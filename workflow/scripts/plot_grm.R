# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
## Original
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms/hdrr/5000/0.8"
## Inbred
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms_inbred/hdrr/5000/0.8"
## Original no missing
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms_no_miss/hdrr/5000/0.8"
## Inbred no missing
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms_inbred_no_miss/hdrr/5000/0.8"
## Chr10
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms_per_chr/hdrr/5000/0.8/10"
## Inbred chr10
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms_inbred_per_chr/hdrr/5000/0.8/10"

F2_SAMPLES = "/hps/software/users/birney/ian/repos/MIKK_HMM/config/F2_samples_converted.csv"

## True
IN_PREF = snakemake@params[["grm_pref"]]
F2_SAMPLES = snakemake@input[["F2_samples"]]
OUT = snakemake@output[["png"]]

# Read in F2 samples

f2 = readr::read_csv(F2_SAMPLES,
                     col_types = c("cccc")) %>% 
  dplyr::select(SAMPLE = finclip_id,
                PAT = pat_line,
                MAT = mat_line) %>% 
  dplyr::mutate(PAT_MAT = paste(PAT, "x", MAT, sep = ""))

# Read in binary files using script provided here: https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

grm = ReadGRMBin(IN_PREF)

# Convert to symmetric matrix
# Guidance here: https://stackoverflow.com/questions/23040676/how-to-fill-in-a-matrix-given-diagonal-and-off-diagonal-elements-in-r
A_GCTA <- matrix(NA, ncol = length(grm$diag), nrow = length(grm$diag))
A_GCTA[lower.tri(A_GCTA)] <- grm$off
A_GCTA[upper.tri(A_GCTA)] <- t(A_GCTA)[upper.tri(t(A_GCTA))]
diag(A_GCTA) <- grm$diag

rownames(A_GCTA) = grm$id$V1
colnames(A_GCTA) = grm$id$V2

# Order
## By sample to compare with manual GRM
ord = order(match(rownames(A_GCTA), f2$SAMPLE))
## by cluster
ord = hclust(dist(A_GCTA, method = "euclidean"), method = "ward.D")$order
labs = rownames(A_GCTA)[ord]
# Get labs with cross
labs_x = tibble(SAMPLE = labs) %>% 
  dplyr::left_join(f2 %>% 
                     dplyr::select(SAMPLE, PAT_MAT),
                   by = "SAMPLE") %>% 
  # combine
  dplyr::mutate(S_X = paste(SAMPLE, PAT_MAT, sep = "_"))

# Order matrix
A_GCTA = A_GCTA[ord, ord]

# Convert to Df
df_fig = A_GCTA %>% 
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


ggsave(OUT,
       fig,
       device = "png",
       width = 30,
       height = 30,
       units = "in",
       dpi = 400)
