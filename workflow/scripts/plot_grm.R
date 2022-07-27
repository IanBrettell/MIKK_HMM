# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms/hdrr/5000/0.8"
IN_PREF = "/hps/nobackup/birney/users/ian/MIKK_HMM/grms_inbred/hdrr/5000/0.8"
PED = "/hps/nobackup/birney/users/ian/MIKK_HMM/peds/F2/hdrr/5000/0.8.ped"


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
m <- matrix(NA, ncol = length(grm$diag), nrow = length(grm$diag))
m[lower.tri(m)] <- grm$off
m[upper.tri(m)] <- t(m)[upper.tri(t(m))]
diag(m) <- grm$diag

rownames(m) = grm$id$V1
colnames(m) = grm$id$V2

# Get order
ord = hclust( dist(m, method = "euclidean"))$order

# Plot

test = m %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "f_id") %>% 
  tidyr::pivot_longer(-c(f_id), names_to = "samples", values_to = "values") %>% 
  dplyr::mutate(dplyr::across(c("f_id", "samples"),
                              ~factor(., levels = ord))) %>% 
  ggplot() +
  geom_tile(aes(x = samples, y = f_id, fill = values)) +
  scale_fill_viridis_c(option = "plasma") +
  theme(aspect.ratio = 1)

######################
# Compute "manually"
######################

IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos/hdrr/5000/0.8/1_38-2_21-2.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos/hdrr/5000/0.8/2_38-2_21-2.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos/hdrr/5000/0.8/99_38-2_40-1.csv")

# Read in files
dat_list = purrr::map(IN, function(DF){
  df = readr::read_csv(DF)
  # Get sample
  ID = as.character(df$SAMPLE[1])
  # Remove sample column and rename GENO_NT column
  df = df %>% 
    dplyr::select(-SAMPLE, {{ID}} := GENO_NT)
})

# Join

df = dat_list %>%
  purrr::reduce(full_join, by=c("CHROM", "POS")) %>% 
  dplyr::arrange(CHROM, POS)

# Create .ped file

ped = df %>% 
  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
  #dplyr::select(-c(CHROM, POS)) %>% 
  # pivot into 3 column
  tidyr::pivot_longer(cols = -c(SNP, CHROM, POS), 
                      names_to = "SAMPLE", 
                      values_to = "GT") %>% 
  # convert SAMPLE to numeric
  dplyr::mutate(SAMPLE = as.numeric(SAMPLE)) %>% 
  # order by SAMPLE, CHROM, and POS
  dplyr::arrange(SAMPLE, CHROM, POS) %>% 
  # remove CHROM and POS columns
  dplyr::select(-c(CHROM, POS)) %>% 
  # replace NA with 0 as required by Plink https://www.cog-genomics.org/plink/1.9/input#plink_irreg
  dplyr::mutate(GT = tidyr::replace_na(GT, "00")) %>% 
  # pivot wide into .ped format (samples to rows, SNPs to columns)
  tidyr::pivot_wider(id_cols = SAMPLE, 
                     names_from = SNP, 
                     values_from = GT)



